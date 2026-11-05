# Short-DNA oligo assembler

**Need a short piece of DNA (up to ~130 bp) and don't want to pay gene-synthesis prices? Order three ordinary oligos instead, and let one PCR join them.**
**想合成一段较短的 DNA（约 130 bp 以内），又不想付基因合成的价格？订 3 条普通寡核苷酸，一次 PCR 拼起来。**

This tool designs the three oligos for you, **checks its own answer** with an independent PCR simulation, and exports an order sheet.
It is one static web page — no install, no server, **your sequences never leave your browser** — in English and Chinese (中 / EN switch, top right).

> **[▶ Open the tool](https://qmmyz.github.io/oligo-assembly-designer/)** · or download `index.html` and open it locally.

![The tool, designing nine signal peptides](docs/screenshot.png)

## Why do it this way

Oligos are priced per base, but most vendors quote that price only up to about **60 nt**. Past that you move to a different synthesis scale or purification and the per-base price climbs — and because coupling efficiency compounds over every step, a longer oligo also carries a smaller fraction of full-length product. Ordering the piece as a synthetic gene instead usually runs into a minimum per-fragment charge, so a ~100 bp fragment can cost about as much as a far larger one.

Splitting the target into **three ordinary oligos of ≤60 nt and joining them in a single PCR** is normally much cheaper and faster. It suits short, defined sequences: signal peptides, linkers and tags, short promoters/RBS/terminators. Compare against your own vendor's price list — the numbers differ, the shape of the curve usually doesn't.

## How it works

![F, M and R overlap to build the fragment](docs/assembly.svg)

Two primers of ≤60 nt can only overlap by a few bases on a ~100 bp fragment, which is too short for a 60 °C anneal. A third oligo **M** in the middle gives two overlaps of 18–45 nt instead. So that **one PCR with a single annealing temperature (59–61 °C)** works, the tool picks overlaps whose Tm falls in that window, keeps every oligo ≤60 nt, and prefers G/C 3′ ends and no dimer-prone ends.

- **F** = (upstream arm +) 5′ part of your sequence · **M** = middle, sense strand · **R** = reverse complement of the 3′ part (+ downstream arm).
- PCR: F and R at ~0.4 µM, M at ~0.02 µM, no template, any high-fidelity master mix. The page prints a starting recipe and program.

## Quick start

1. Open the page. It loads a 9-sequence example and designs it immediately.
2. Replace the text with your sequences, one per line: `name  DNA` (or paste FASTA).
3. **Homology arms:** put yours in the two arm boxes if you will clone by recombination, or clear both if you just need the bare DNA. Keep *coding sequence* ticked for a protein-coding insert (checks ATG, reading frame, internal stops); untick it for non-coding DNA.
4. **Design & verify.** Each design shows a map, the three oligos, any warnings, and the full check list. **Download CSV / Copy** gives the order sheet.

Typical limit: product (including arms) up to about **130 bp** works reliably; at 140 bp only some sequences have a solution (random-sequence test: 30/30 designs up to 130 bp, 9/30 at 140 bp). Longer needs a different strategy.

## How the answer is checked

The checker looks only at the final three oligos, simulates PCR (any 3′ end of ≥12 nt that pairs with another strand is extended), and requires that:

- the full product and its complement are formed, and **nothing else** is (no mispriming products among the three oligos);
- the product is exactly *arm + your sequence + arm*, and F/R start the right strands;
- for coding sequences: ATG start, frame, no internal stop, expected amino acids.

`tests/` holds 32 tests: an independent Python re-implementation of the Tm and the PCR simulation, plus **mutation tests** (substituted/deleted bases, swapped F/R, broken overlaps, wrong arms) that confirm the checker rejects bad oligos — for the signal-peptide set and for random non-coding DNA of 60–130 bp.

## Please read: what this does not do

- **In-silico only.** The design has not been wet-lab tested yet. The single 60 °C anneal, the amount of M and the cycle number are principled starting points, not a validated protocol. Run a gradient the first time.
- **Sequence-verify your clone.** Oligos are synthesized with errors, and PCR assembly does not correct them. Clone the product and sequence it before you rely on it.
- **No vector checks.** Reading frame across the vector junction and primer mispriming on the vector are not checked.
- **Tm is calibrated, not authoritative.** Nearest-neighbour Tm (SantaLucia 1998), tuned to agree with SnapGene readings to within about 1 °C. If your own design tool disagrees, trust the gradient. Details: [docs/design-notes.md](docs/design-notes.md).

## For developers

`window.SPDesigner = {design, verify, tm, setLang, …}` is exposed for scripting. Tests open `index.html` in headless Chromium:

```bash
pip install -r requirements-dev.txt
playwright install chromium
python -m pytest -q tests
```

`data/` holds the nine signal-peptide sequences used as golden test input (and as the page's example); provenance and hashes are in `data/MANIFEST.md`.

## Attribution and license

Nearest-neighbour parameters: SantaLucia, J. (1998) *PNAS* 95:1460–1465. MIT licensed, see [LICENSE](LICENSE).
