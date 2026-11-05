"""独立于页面代码的测试：用 Playwright 打开 index.html 取设计结果，再用 Python 另写的
PCR 闭包模拟 / Tm 去核对，并对故意改坏的引物做“变异测试”，确认验证器会拦下。

运行：  <venv>/bin/python -m pytest -q tests      （需要 playwright + chromium）
"""
import csv, math, os, pathlib, pytest
from playwright.sync_api import sync_playwright

ROOT = pathlib.Path(__file__).resolve().parent.parent
UP, DOWN = "GGAGAGGGACGCGT", "TGGCCATTAACCGC"

# ---------- Python 侧的独立实现（不与页面共享代码） ----------
def rc(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))
NN = {"AA":(-7.9,-22.2),"TT":(-7.9,-22.2),"AT":(-7.2,-20.4),"TA":(-7.2,-21.3),"CA":(-8.5,-22.7),"TG":(-8.5,-22.7),
      "GT":(-8.4,-22.4),"AC":(-8.4,-22.4),"CT":(-7.8,-21.0),"AG":(-7.8,-21.0),"GA":(-8.2,-22.2),"TC":(-8.2,-22.2),
      "CG":(-10.6,-27.2),"GC":(-9.8,-24.4),"GG":(-8.0,-19.9),"CC":(-8.0,-19.9)}
def tm(s):
    """SantaLucia 1998 最近邻 + 熵项盐校正，50 mM Na+、0.9 µM 引物（用 SnapGene 读数拟合，见 README）。"""
    dH = dS = 0.0
    for b in (s[0], s[-1]):
        dH, dS = (dH+0.1, dS-2.8) if b in "GC" else (dH+2.3, dS+4.1)
    for i in range(len(s)-1):
        h, d = NN[s[i:i+2]]; dH += h; dS += d
    dS += 0.368*(len(s)-1)*math.log(0.05)
    return dH*1000/(dS + 1.987*math.log(0.9e-6/4)) - 273.15

def py_assemble(oligos, nmin=12):
    """返回 PCR 闭包得到的全部单链集合。与页面实现算法相同，但独立编写。"""
    pool = set(oligos)
    while True:
        new = set()
        for x in pool:
            key = rc(x[-nmin:])
            for y in pool:
                start = 0
                while (j := y.find(key, start)) != -1:
                    cand = x + rc(y[:j])
                    if cand not in pool: new.add(cand)
                    start = j + 1
        if not new: return pool
        pool |= new
        if len(pool) > 400: return pool

def py_ok(F, M, R, sp, up=UP, down=DOWN):
    T = up + sp + down
    pool = py_assemble([F, M, R])
    # 角色检查：F 必须是上链前缀、R 必须是下链前缀（F/R 只是标签，互换后物理上是同一组寡核苷酸，
    # 所以闭包本身分辨不出，必须单独核对）
    return (T in pool and rc(T) in pool and all(s in T or s in rc(T) for s in pool)
            and T.startswith(F) and rc(T).startswith(R)
            and max(len(F), len(M), len(R)) <= 60)

# ---------- 固件 ----------
@pytest.fixture(scope="module")
def page():
    with sync_playwright() as p:
        b = p.chromium.launch()
        pg = b.new_page()
        pg.goto((ROOT/"index.html").as_uri())
        yield pg
        b.close()

def golden():
    with open(ROOT/"data"/"golden_sp.csv") as f:
        return list(csv.DictReader(f))

def js_design(page, name, sp, cfg=None):
    return page.evaluate("([n,s,c])=>SPDesigner.design(n,s,c)", [name, sp, cfg or {}])

def js_verify(page, **spec):
    return page.evaluate("s=>SPDesigner.verify(s)", spec)

# ---------- 测试 ----------
SNAPGENE = [  # (序列, SnapGene 显示的 Tm 区间)；读数取整，所以允许落在区间两端各 0.5 内，再放宽 0.4 容差
    ("AAATCATCTCAAAAAAATGGGT", (51, 51)), ("AGTCTTTTAAGTAAGTCTACTCTGAAC", (54, 54)),
    ("AAGCTAGACAAAACGGACAAAATAAAAATTGGC", (61, 61)), ("TTTGGAAAATGGGAGTATAATGGAAGTAAGAAATCCTTTATTG", (60, 61)),
    ("TGCAATATAAAAAAATTAATAAAAATC", (47, 47)), ("AAATAATAAAGATAATTCTTATAATCAGCACATGCATGTAT", (58, 58)),
    ("CTAACAATGTGTTTTCTGGGGCTCA", (59, 59)), ("TTTTCTGGGGCTCATACTTTTTGTGC", (59, 59)),
]
def test_tm_matches_snapgene_readings(page):
    for seq, (lo, hi) in SNAPGENE:
        for v in (page.evaluate("s=>SPDesigner.tm(s)", seq), tm(seq)):
            assert lo - 0.9 <= v <= hi + 0.9, (seq, v, lo, hi)
    assert abs(page.evaluate("s=>SPDesigner.tm(s)", SNAPGENE[3][0]) - tm(SNAPGENE[3][0])) < 1e-6

@pytest.mark.parametrize("row", golden(), ids=lambda r: r["name"])
def test_golden_design_and_two_verifiers_agree(page, row):
    d = js_design(page, row["name"], row["sp_dna"])
    assert d["ok"], d.get("error")
    F, M, R = d["F"], d["M"], d["R"]
    assert max(len(F), len(M), len(R)) <= 60
    assert F.startswith(UP) and rc(R).endswith(DOWN)              # 引物自带同源臂
    v = js_verify(page, up=UP, down=DOWN, sp=row["sp_dna"], F=F, M=M, R=R, expectedAA=row["expected_aa"])
    assert v["ok"], [c for c in v["checks"] if not c["ok"]]
    assert py_ok(F, M, R, row["sp_dna"])                           # Python 独立判定也通过
    # 重叠区 Tm 满足厂家门槛，并与 Python 重算一致
    for ov in (d["ovF"], d["ovR"]):
        assert ov["tm"] >= 55            # 硬下限
        assert abs(tm(ov["seq"]) - ov["tm"]) < 1e-6
    # M 从 F 的 5′ 侧起、到 R 的 3′ 侧止（F 与 R 允许直接重叠）
    assert d["ovF"]["a"] < d["ovR"]["a"] and d["ovF"]["b"] < d["ovR"]["b"]

def mutate_cases(F, M, R, d):
    ovf = d['ovF']['b'] - d['ovF']['a']; ovr = d['ovR']['b'] - d['ovR']['a']
    flip = lambda s, i: s[:i] + {"A":"C","C":"G","G":"T","T":"A"}[s[i]] + s[i+1:]
    return {
        "F 中间替换一个碱基": (flip(F, len(F)//2), M, R),
        "M 中间替换一个碱基": (F, flip(M, len(M)//2), R),
        "R 中间替换一个碱基": (F, M, flip(R, len(R)//2)),
        "M 中间删一个碱基": (F, M[:len(M)//2] + M[len(M)//2+1:], R),
        "F 的上游臂首碱基改错": (flip(F, 0), M, R),
        "R 的下游臂末碱基改错": (F, M, flip(R, 0)),
        "F 与 R 互换": (R, M, F),
        "M 5′ 端截短到与 F 只剩 8 nt 重叠": (F, M[ovf-8:], R),
        "M 3′ 端截短到与 R 只剩 8 nt 重叠": (F, M[:len(M)-(ovr-8)], R),
        "漏掉 M（用 F 代替）": (F, F, R),
    }

@pytest.mark.parametrize("row", golden(), ids=lambda r: r["name"])
def test_verifier_catches_mutations(page, row):
    d = js_design(page, row["name"], row["sp_dna"])
    for label, (F, M, R) in mutate_cases(d["F"], d["M"], d["R"], d).items():
        v = js_verify(page, up=UP, down=DOWN, sp=row["sp_dna"], F=F, M=M, R=R)
        assert not v["ok"], f"JS 验证器没拦住：{row['name']} {label}"
        assert not py_ok(F, M, R, row["sp_dna"]), f"Python 验证器没拦住：{row['name']} {label}"

def test_wrong_arm_is_caught(page):
    row = golden()[0]
    d = js_design(page, row["name"], row["sp_dna"])
    v = js_verify(page, up="GGAGAGGGACGCGA", down=DOWN, sp=row["sp_dna"], F=d["F"], M=d["M"], R=d["R"])
    assert not v["ok"]

def test_bad_sp_inputs_caught(page):
    row = golden()[0]; d = js_design(page, "x", row["sp_dna"])
    sp = row["sp_dna"]
    bad = {"不以 ATG 开头": "TTG" + sp[3:], "长度不是 3 的倍数": sp + "A", "内部终止密码子": sp[:9] + "TAA" + sp[12:]}
    for label, s in bad.items():
        dd = js_design(page, "x", s)
        v = js_verify(page, up=UP, down=DOWN, sp=s, F=dd["F"], M=dd["M"], R=dd["R"])
        assert not v["ok"], label

def test_infeasible_reports_clear_error(page):
    d = js_design(page, "long", row_sp := ("ATG" + "GCT"*120))
    assert d["ok"] is False and "找不到可行" in d["error"]
