# src/ — search code index

All scripts are plain Python (numpy, scipy, mpmath, python-flint). They are research code: the verified claims
are re-checked by the independent scripts in `../certify/` and `../verify.py`.

| Script | Purpose | Hypothesis |
|---|---|---|
| `digits.py` | Carry automaton (carries in {0,1}): growth rate λ of digit sums of S in base m = n+1 | H0 |
| `digits2.py` | General carry automaton (any base m ≥ 2, arbitrary carries) | H0, H6 |
| `test_digits.py` | Automaton vs brute force; λ of the Coppersmith–Davenport seeds | H0 |
| `zero_patterns.py` | Annealing over zero patterns of p² (base n+1) | H1 |
| `enum_sym.py` | Exhaustive enumeration of symmetric patterns with square systems | H2, H5 |
| `symsolve.py` | Random-start Newton solver for the symmetric seed system | H2 |
| `homotopy.py` | Total-degree homotopy: all isolated solutions of the symmetric system | H6a |
| `tower.py` | Exact (integer) tower experiments with python-flint | H4 |
| `mix.py` | Mixed-radix / different seeds per level | H3 |
| `miracle_sym.py` | Search for overdetermined ("miracle") symmetric solutions | H5 |
| `enum_overlap.py` | Exhaustive symmetric patterns × overlapping bases | H6 |
| `solve_targets.py` | Solve selected patterns (random starts) | H6a |
| `superseed.py` | High-precision refinement; super-seed analysis | H6c |
| `overlap_check.py` | Density / consecutive-ratio diagnostics for overlapping towers | H6b |
| `minds.py` | Brute-force check of the unique-minimal-digit-sum lemma | H6b |
| `anneal_overlap.py` | Annealing over patterns and bases, n ≤ 32 | H7 |
| `rational_hunt.py`, `rational_strict.py` | Rationality test of all isolated solutions (float screen / 70-digit strict) | H8 |
| `rational_big.py` | Same for n = 19, 20, 22 | H12 |
| `slack_tab.py` | Exponents of slack-1 (curve) patterns | H9 |
| `curve_scan.py` | Rational-point sweep in ρ on the n = 18 curve | H9 |
| `rational_overlap.py` | Known rational seeds with overlapping bases | H11 |
| `ratbeam.py` | Rational beam search over products A(x)(v + u x^k) | H13 |
| `rational_point_sampler.py` | Rational-point sampler on slack-1 families (long runs) | H14 |
| `sparsify.py` | Structure-agnostic sparsification optimizer | H15 |
