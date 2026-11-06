# sparse-squares

**Dense polynomials whose squares are sparse.** This repository is a computational attack on the Rényi–Erdős
problem of how few terms the square of a polynomial with many terms can have (Erdős Problem
[#485](https://www.erdosproblems.com/485); Epoch AI, *FrontierMath: Open Problems*,
["Dense Polynomials Whose Squares Are Sparse"](https://epoch.ai/frontiermath/open-problems/polynomial-square-few-coefficients)).

*Versión en español: [README.es.md](README.es.md).*

> **Status (November 2026): draft, not yet reviewed by an expert.**
> The main claim below has been machine-checked (interval-arithmetic certificate and independent
> recomputation), and the proof is written out in [`note/note.pdf`](note/note.pdf). Corrections are very welcome.

## TL;DR

| Question | Status |
|---|---|
| **Epoch target:** an integer, dense p of degree n with \|p²\| < √n (exponent < 1/2) | **Not solved.** No candidate. A heuristic argument (below) suggests the methods here stall near 0.77. |
| **Integer/rational record:** Verdenius 1949, exponent log 6/log 9 ≈ 0.8155 | **Not improved.** Reproduced exactly; all searches for better rational seeds were negative. |
| **Real record for complete polynomials:** Verdenius 1949, exponent log 8/log 13 ≈ 0.81071 | **Improved to log 9/log 17 ≈ 0.77552** (claimed; see the theorem and caveats below). |

Notation: ‖f‖ is the number of nonzero coefficients. Q̂(k) = min ‖f²‖ over real polynomials f of degree k−1 with all k coefficients nonzero ("complete" or "dense" polynomials).

## Main result (claimed)

**Theorem.** For every K ≥ 1 there is a real polynomial of degree D_K = 20(17^K − 1)/16 with all D_K + 1
coefficients nonzero whose square has at most 21·9^(K−1) nonzero coefficients. Consequently

  Q̂(k) ≤ (21/9) · k^(log 9 / log 17) ≈ 2.33 · k^0.77552,  for k = D_K + 1.

Verdenius' 1949 exponent is log 8/log 13 ≈ 0.81071. As far as we can tell from the literature (Schinzel–Zannier
2009, Mantova 2016, the #485 thread in 2026), it had not been improved.

### Caveats

1. **Real coefficients only.** The seed is irrational. For ℚ[x], which is how #485 is stated and also the setting of the Epoch problem, the best exponent is still Verdenius' log 6/log 9.
2. **Asymptotic.** The explicit bound 21·9^(K−1) drops below k^0.81071 only from K = 7, which is about 5·10⁸ terms.
3. **Along a sequence of k.** This is the same form as Verdenius' statement. Extending it to all k should follow from Erdős' 1949 interpolation argument, but that is not written out here.
4. **No expert review yet.**

## How it works

**1. Towers and carries (Verdenius' mechanism, reconstructed).** If p has degree n and S = supp(p²), then
P_K(x) = ∏_{j<K} p(x^{m^j}) satisfies supp(P_K²) ⊆ {Σ s_j m^j : s_j ∈ S}. Because of carries, this set grows like
λ^K with λ much smaller than |S|. We compute λ exactly with a carry automaton (`src/digits.py`, `src/digits2.py`,
cross-checked against brute force).

With m = n + 1, P_K is automatically dense. An exhaustive search over symmetric zero patterns for n ≤ 20
returns **exactly Verdenius' seeds**:
- n = 8, rational, λ = 6, p = [1, −2, −2, −4, −10, −28, −98, −686, 2401];
- n = 12, real, λ = 8.

That family has λ = n/2 + 2, which is optimal at n = 12. This explains why the record stood.

**2. Overlapping bases.** Take m < n + 1. The best family has three blocks:

  S = [0, r] ∪ [m, m+r] ∪ [2m, 2m+r],  n = 3r + 2,  m = n − r/2,

which gives λ = r + 3 with no further carries. The optimum is at r = 6: **n = 20, m = 17, λ = 9**. Annealing over all symmetric patterns with n ≤ 32 and all bases always lands in this family.

**3. A real seed exists.** A palindromic p of degree 20 has its zeros on [7, 16] imposed by 10 quadratic equations. Total-degree homotopy finds 259 complex solutions, 42 of them real with all coefficients nonzero. One of them is certified by the **Krawczyk test** in 110-digit interval arithmetic: there is a unique exact solution in a box of radius 10⁻⁶⁰, all coefficients are nonzero, and supp(p²) = S exactly.

**4. The tower stays dense for every K.** With overlapping bases, cancellations are possible a priori. The proof:
- Choose a prime ℓ at which all coefficients of p are units, and replace p(x) by p(ℓx).
- Each representation N = Σ i_j m^j with digits 0 ≤ i_j ≤ n then contributes a term of ℓ-adic valuation proportional to Σ i_j.
- **Lemma:** if m − 1 ≤ n < 2m, the representation with minimal digit sum is unique. The proof is a carry-chain argument.
- So exactly one term has minimal valuation, and every coefficient is nonzero.

## Verification (independent of the search code)

| Check | Script | Output |
|---|---|---|
| Exact integer check of the Epoch criterion | `verify.py` | `results/verify_*.txt` |
| Krawczyk certificate for the seeds n = 14 (m = 12) and n = 20 (m = 17) | `certify/krawczyk.py` | `results/krawczyk_*.json` (all `true`) |
| Tower check: re-derives the seed, counts digit sums by direct enumeration (21·9^(K−1) up to K = 6), and computes P₂² directly (dense, ‖P₂²‖ = 189) | `certify/tower_check.py` | `results/tower_check.txt` |
| Brute-force check of the unique-minimal-digit-sum lemma (K ≤ 4) | `src/minds.py` | `results/unique_min_digitsum_check.txt` |

### Reproduce

```bash
pip install -r requirements.txt
./reproduce.sh            # a few seconds on a laptop
```

## Negative results

All of these are logged in [HYPOTHESES.md](HYPOTHESES.md):
- mixing different seeds across levels;
- extra coefficient cancellations in P_K²;
- "miracle" solutions of overdetermined symmetric systems;
- rational points on all isolated solutions (n ≤ 16, and n = 19, 20);
- rational beam search over products A(x)·(v + u·x^k), in the style of Coppersmith–Davenport;
- rational-point sampling on the five relevant curves;
- a structure-agnostic sparsification optimizer (it does not even rediscover the degree-12 examples).

## Why not exponent 1/2? (heuristic)

A solvable pattern needs |S| ≈ n + 1 by dimension counting. Writing S = A + m·E gives λ ≈ |A + E| ≥ |A| + |E| − 1
and m ≲ 2n / max E. Balancing these leaves the exponent near 0.77. Going below 1/2 would need patterns with
|S| ≪ n that are still solvable, i.e. a new source of "free" cancellations of the kind Coppersmith–Davenport
exploited at degree 12, or a construction that is not a single-seed tower.

## Repository layout

```
README.md / README.es.md   this file (English / Spanish)
HYPOTHESES.md              log of every hypothesis tried, with outcome (Spanish original: HYPOTHESES.es.md)
note/                      technical note (LaTeX + PDF) with full proofs
verify.py                  independent exact integer verifier (Epoch JSON format)
certify/                   interval-arithmetic certificate and independent tower check
seeds/                     seeds as JSON (verdenius_integer_n8.json is in Epoch format)
results/                   outputs of all verification runs
src/                       search code; see src/README.md for an index
reproduce.sh               re-runs all verifications
```

## References

- A. Rényi, *On the minimal number of terms of the square of a polynomial*, Hungarica Acta Math. 1 (1947).
- P. Erdős, *On the number of terms of the square of a polynomial*, Nieuw Arch. Wisk. 23 (1949).
- W. Verdenius, *On the number of terms of the square and the cube of polynomials*, Indag. Math. 11 (1949).
- D. Coppersmith, J. Davenport, *Polynomials whose powers are sparse*, Acta Arith. 58 (1991).
- J. Abbott, *Sparse squares of polynomials*, Math. Comp. 71 (2002).
- A. Schinzel, Acta Arith. 49 (1987); A. Schinzel, U. Zannier, Rend. Lincei Mat. Appl. 20 (2009).
- V. Mantova, survey, Riv. Mat. Univ. Parma 7 (2016).

## Authorship and AI disclosure

Directed by Roberto Vetere. The search, the certificates and the drafts were produced with an AI assistant
(Claude, Anthropic). Every claim marked as verified is machine-checked; the proofs still await human expert review.
