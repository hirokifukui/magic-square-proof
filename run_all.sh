#!/usr/bin/env bash
#
# run_all.sh — execute the full 17-script verification suite for
# "Nonexistence of 3x3 Magic Squares of Squares" (Fukui 2026).
#
# Usage: bash run_all.sh
#
# Requires: Python 3.9+ with SymPy installed.
#
# Total runtime: ~15 seconds on standard hardware.

set -u

# Resolve the directory of this script so the harness can be run from anywhere.
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

SCRIPTS=(
  # Group A: Pillar A — case k=1 (split primes)
  "scripts/phase_A_pillar_A_verification/verify_P2_identity.py"
  "scripts/phase_A_pillar_A_verification/verify_P2_nu_p.py"
  "scripts/phase_A_pillar_A_verification/verify_P2_nu_p_symbolic.py"
  "scripts/phase_A_pillar_A_verification/verify_P2_row_sum_reduction.py"
  "scripts/phase_A_pillar_A_verification/verify_c_value.py"
  "scripts/phase_A_pillar_A_verification/verify_rep_count_v4.py"
  # Group B: Phase C — Gaussian shadow and dichotomy
  "scripts/phase_C_supplementary/verify_branch1_dichotomy.py"
  "scripts/phase_C_supplementary/verify_c12_vacuousness.py"
  "scripts/phase_C_supplementary/verify_k3_cardinality.py"
  # Group C: Main case-analysis scripts (root level)
  "verify_k3base.py"
  "verify_birational.py"
  "verify_cor64.py"
  "verify_prop81.py"
  "verify_caseIV.py"
  "verify_k3_via_thm75.py"
  "generate_k3_table.py"
  "verify_appendixA.py"
)

TOTAL=${#SCRIPTS[@]}
pass=0
fail=0
declare -a failures

mkdir -p .verify_logs

echo "Running $TOTAL verification scripts..."
echo

i=0
for s in "${SCRIPTS[@]}"; do
  i=$((i + 1))
  base=$(basename "$s" .py)
  log=".verify_logs/${i}_${base}.log"
  t0=$(date +%s)
  if python3 "$s" > "$log" 2>&1; then
    t1=$(date +%s)
    printf "  [%2d/%d] PASS  (%2ds)  %s\n" "$i" "$TOTAL" "$((t1 - t0))" "$s"
    pass=$((pass + 1))
  else
    t1=$(date +%s)
    printf "  [%2d/%d] FAIL  (%2ds)  %s\n" "$i" "$TOTAL" "$((t1 - t0))" "$s"
    failures+=("$s")
    fail=$((fail + 1))
  fi
done

echo
echo "=== Summary ==="
echo "PASS: $pass / $TOTAL"
echo "FAIL: $fail / $TOTAL"

if [ "$fail" -gt 0 ]; then
  echo
  echo "Failed scripts (see .verify_logs/ for details):"
  for f in "${failures[@]}"; do
    echo "  - $f"
  done
  exit 1
fi

exit 0
