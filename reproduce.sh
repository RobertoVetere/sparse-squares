#!/usr/bin/env bash
# Re-run every verification behind the claims in README.md. Outputs go to results/.
set -euo pipefail
cd "$(dirname "$0")"
echo "== 1. carry automaton vs brute force"
(cd src && python3 test_digits.py)
echo "== 2. exact integer checks (Verdenius seed, Epoch-format example)"
python3 verify.py seeds/verdenius_integer_n8.json
python3 verify.py results/epoch_format_example_verdenius_K2.json
echo "== 3. Krawczyk certificates (interval arithmetic, 110 digits)"
for s in seed_n14_m12 seed_n20_m17; do
  python3 certify/krawczyk.py seeds/$s.json > results/krawczyk_$s.json
  grep -E '"(krawczyk_inclusion|all_coefficients_nonzero|square_support_exactly_S)"' results/krawczyk_$s.json
done
echo "== 4. independent tower check (digit-sum counts, P_2 dense, |P_2^2|)"
python3 certify/tower_check.py seeds/seed_n20_m17.json 3 2 | tee results/tower_check.txt
echo "== 5. unique-minimal-digit-sum lemma, brute force"
(cd src && python3 minds.py) | tee results/unique_min_digitsum_check.txt
echo "done"
