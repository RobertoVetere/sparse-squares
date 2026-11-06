# Hypothesis log

Every idea tried, what it gave, and why it was kept or dropped. Spanish original: [HYPOTHESES.es.md](HYPOTHESES.es.md).

| # | Hypothesis / experiment | Result | Decision |
|---|---|---|---|
| H0 | The exponent comes from towers p(x)p(x^b)…; supp(P_K²) is the set of carry digit-sums of supp(p²). | Carry automaton agrees with brute force (700 random cases). The Coppersmith–Davenport degree-12 seed gives λ = 9.87, exponent 0.893. | Basis for everything below. |
| H0′ | Correction to my first summary: I said Verdenius' seeds cannot be dense. | False. His bound is for complete polynomials; the seed is dense and works through the H0 mechanism. | Corrected. |
| H1 | Anneal over zero patterns Z of p² with base m = n+1. | \|Z\| = n reproduces exactly log 6/log 9 (n = 8) and log 8/log 13 (n = 12). | Verdenius' mechanism identified. |
| H2 | Symmetric ansatz x^n p(1/x) = C p(ρx); exhaustive search over square systems, n ≤ 20, m = n+1. | Optimum is Verdenius n = 12; the family has λ = n/2+2. Integer seed n = 8 reproduced exactly. | Family exhausted. |
| H3 | Different seeds at different levels (mixed radix). | ρ(M₁M₂) never beats the best single seed (best 0.81071). | Dropped. |
| H4 | Coefficient cancellations in P_K² beyond the support bound. | Verdenius n = 8: exact count equals the bound up to K = 5. Same for CD12 and A6. | Dropped for these seeds. |
| H5 | Symmetric systems overdetermined by one ("miracles"); potential exponents 0.747–0.827. | Complex Levenberg–Marquardt, 150–300 starts each: no non-degenerate solutions. | Dropped (not exhaustive, but no signal). |
| H6 | **Overlapping base m < n+1.** | Exhaustive n ≤ 18: n = 14, m = 12 → 0.78309. Three-block family S = A + m{0,1,2}; optimum n = 20, m = 17 → **0.77552**. | **Pursued.** |
| H6a | Do real seeds exist for these patterns? | Total-degree homotopy: n = 14 → 28 complex solutions, 8 real and dense; n = 20 → 259/262 complex, 42 real and dense. Krawczyk certifies one of each. | Confirmed. |
| H6b | Density of the overlapping tower. | Consecutive-ratio / Möbius route: inconclusive (ratios approach the forbidden values). ℓ-adic route plus the unique-minimal-digit-sum lemma: **proved for all K**. | Solved. |
| H6c | Super-seed P_K stacked without overlap. | Converges to 0.7755, but would need K ≈ 8–9 computed explicitly. | Unnecessary after H6b. |
| H7 | Wider annealing (n ≤ 32, all bases). | Always lands in the three-block family; best 0.77552. | Optimum of this framework (heuristic). |
| H8 | **Rational** seeds with exponent < 0.8155 (integer record). | All isolated solutions for n ≤ 16: only Verdenius n = 8 is rational. | Negative. |
| H9 | Curve n = 18, m = 15 (exponent 0.81137): rational points via a sweep over ρ ∈ ℚ. | ~20 values of ρ, 70–78 complex and ~20 real solutions each, none rational (the apparent hits were tolerance false positives). | Negative (limited search). |
| H10 | Can this framework reach 1/2 (the Epoch target)? | Heuristic counting argument: no; the optimum is ≈ 0.77. | Needs a different idea (README). |
| H11 | Known rational seeds (V8, A6, CD12 with α ∈ {−110, −253, 15625}) with overlapping bases. | Best is V8 at m = 9 (0.8155); overlap never helps (A6 0.847, CD12 0.850). | Negative. |
| H12 | Rationality of all isolated solutions for n = 19 and n = 20 (patterns giving 0.7925 and 0.7755). | 894 + 521 complex solutions, 0 rational. | Negative. |
| H13 | Rational beam search: P = A(x)·(v + u x^k), u/v a rational root of a cancellation condition; 4 rounds, width 40, degree ≤ 40. | ~16 000 candidates; best 0.850 (= CD12); gets worse with depth. | Negative. |
| H14 | Rational-point sampling on the five slack-1 curves (n = 18, 22, 23, 24, 26; potential exponents 0.800–0.815): fix ρ ∈ ℚ, solve the square system by homotopy, test rationality at 60 digits (`src/rational_point_sampler.py`; the control run recovers Verdenius' seed). | Local test run: ~100 s per sample per core; 5 samples (n = 22, 23) gave 2247 solutions, 252 real and dense, 0 rational. | Dropped: too slow here and low prior probability. Runnable on one's own machine. |
| H15 | Structure-agnostic search: minimise a smooth surrogate of \|supp p²\| without prescribing the zero pattern (`src/sparsify.py`). | Degree 12, ~9 min: best \|p²\| = 14 = n+2 (generic count). Does **not** rediscover the Coppersmith–Davenport examples (12). | Weaker than theory as is; would need structured seeding. |
