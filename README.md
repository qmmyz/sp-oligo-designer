# sp-oligo-designer

**Signal-peptide DNA → three ≤60-nt oligos (F / M / R) with homology arms → in-page verification → order sheet.**
**信号肽 DNA → 带同源臂的 3 条 ≤60 nt 寡核苷酸（F / M / R）→ 页面内独立验证 → 订购表。**

A single static page (`index.html`, no build step, no server, nothing is uploaded). Open it in a browser, paste signal-peptide (SP) sequences, and get oligos you can order directly. The UI is in Chinese.
单文件网页，浏览器直接打开，序列不会上传。界面为中文。

![screenshot](docs/screenshot.png)

## Why three oligos / 为什么是 3 条

A signal peptide plus two 14-nt homology arms is ~100–115 bp. Two primers of ≤60 nt can overlap by only 5–17 nt, which is too short for a 60 °C anneal. Adding a middle oligo **M** gives two independent overlaps of 18–45 nt:

```
 upstream arm   signal peptide (SP)                  downstream arm
 |-----------|=================================|--------------|   full product (top strand)
 F  →→→→→→→→→→→→→→→→→→
                 M  →→→→→→→→→→→→→→→→→→→→→→→→→→→
                                         ←←←←←←←←←←←←←←←←←←  R
        F∩M overlap ┘                  └ M∩R overlap
```

- **F** = upstream arm + SP 5′ part. **M** = middle of the sense strand (no arms). **R** = reverse complement of the SP 3′ part, with the downstream arm at its 5′ end.
- PCR: F and R at ~0.4 µM, M at ~0.02 µM, no template, **one** annealing temperature (60 °C).

## Usage / 使用

1. Open `index.html` in a browser (it loads a 9-SP example and designs it on load).
2. Paste one SP per line: `name  DNA` (optionally a third column with the expected amino acids), or a FASTA.
3. Set the homology arms and limits if needed (defaults: 14-nt arms, ≤60 nt, overlap 18–45 nt, overlap Tm window 59–61 °C).
4. **Design and verify** → check each design's map, warnings and PCR-program card → copy the oligos.

The page also exposes `window.SPDesigner = {design, verify, tm, …}` for scripting.

## Verification / 验证

`verify()` looks only at the final three oligos: it simulates PCR as a closure (a 3′ end of ≥12 nt that pairs to another strand is extended to the template's 5′ end), then asserts that

- the pool assembles `upstream arm + SP + downstream arm` and its reverse complement,
- every strand in the pool is a piece of that product (no off-target strands),
- arms, ATG start, reading frame, no internal stop and the expected amino acids are right.

`tests/test_design.py` re-implements Tm and the closure independently in Python and runs **mutation tests** (substituted/deleted bases, swapped F/R, broken overlaps, wrong arms) to confirm the verifier rejects bad oligos.

## Tm model / Tm 算法

SantaLucia 1998 nearest-neighbour with the entropy-term salt correction `ΔS += 0.368·(N−1)·ln[Na⁺]`, 50 mM Na⁺, 0.9 µM oligo. The salt/concentration parameters were **fitted to 8 SnapGene readings** (all within ±0.9 °C); they are not SnapGene's published parameters. Details and caveats: [docs/design-notes.md](docs/design-notes.md).

## Tests

```bash
pip install -r requirements-dev.txt
playwright install chromium
python -m pytest -q tests
```

The tests open `index.html` in headless Chromium.

## Limitations / 未覆盖

- No vector sequence: the vector reading frame and primer mispriming sites on the vector are **not** checked.
- The Tm model is calibrated on only 8 SnapGene readings; the polymerase vendor's own Tm scale (Primer Premier) was not reconciled.
- No wet-lab validation yet of the M amount, cycle number, or single 60 °C annealing for 18–39 nt overlaps.

## Data

`data/golden_sp.csv` holds the 9 signal-peptide sequences used as golden test input; `data/*.fasta` are the same sequences for opening in SnapGene; provenance and hashes are in `data/MANIFEST.md`.

## Method attribution

SantaLucia, J. (1998) *PNAS* 95:1460–1465 (nearest-neighbour parameters). The 60 °C single-anneal rule follows the Vazyme Phanta UniFi master-mix manual (P516/P526 family).

## License

MIT — see [LICENSE](LICENSE).
