#!/usr/bin/env python3
"""
Verify whether the λ-scaled parity sum Σ λ_ε σ_ε s_{jε} can vanish.

Setup (Lemma A.3):
  λ_x = (D̄/D_x)⁴  for x ∈ {a,b,c}
  D̄ = lcm(D_a, D_b, D_c)
  D_a, D_b, D_c are products of inert primes (q ≡ 3 mod 4).

The sum in question:
  S_j = λ_a σ_a s_{ja} + λ_b σ_b s_{jb} + λ_c σ_c s_{jc}
      = λ_a·(±1) + λ_b·(±1) + λ_c·(±1)     (each ±1 is σ_ε s_{jε})

With σ_a = σ_b = +1, σ_c = -1 (the convention R_a + R_b = R_c):
  σ_a s_{ja} ∈ {±1}, σ_b s_{jb} ∈ {±1}, σ_c s_{jc} ∈ {∓1}
  (the sign of σ_c s_{jc} = -s_{jc})

We need S_j ≠ 0 for the determinant D to be nonzero.
"""
from math import gcd, lcm
from itertools import product as iprod

print("=" * 60)
print("  Structure of λ_ε = (D̄/D_ε)⁴")
print("=" * 60)

# D_a, D_b, D_c are products of primes ≡ 3 (mod 4).
# Inert primes up to 50: 3, 7, 11, 19, 23, 31, 43, 47
inert_primes = [3, 7, 11, 19, 23, 31, 43, 47]

# D_ε is a product of a subset of these primes (or 1 for empty set).
# For efficiency, use small D values that are products of inert primes.
def gen_D_values(primes, max_val=200):
    """Generate all products of subsets of primes up to max_val."""
    vals = {1}
    for p in primes:
        new = set()
        for v in vals:
            if v * p <= max_val:
                new.add(v * p)
        vals |= new
    return sorted(vals)

D_values = gen_D_values(inert_primes, 200)
print(f"  D values (products of inert primes ≤ 200): {D_values}")
print(f"  Count: {len(D_values)}")

print(f"\n  λ_ε = (D̄/D_ε)⁴ where D̄ = lcm(D_a, D_b, D_c).")
print(f"  Since D̄/D_ε is a POSITIVE INTEGER, λ_ε is a perfect 4th power.")
print(f"  Possible λ values: 1⁴=1, 2⁴=16, 3⁴=81, 4⁴=256, 5⁴=625, ...")

print(f"\n{'='*60}")
print(f"  Search for S_j = 0")
print(f"{'='*60}")

# For each triple (D_a, D_b, D_c) and each sign pattern:
# S_j = ε_a · λ_a + ε_b · λ_b + ε_c · λ_c
# where ε_a, ε_b, ε_c ∈ {+1, -1} are the signs σ_x s_{jx}.

# Sign patterns: there are 8, but due to σ_a=σ_b=+1, σ_c=-1,
# the actual signs are (s_{ja}, s_{jb}, -s_{jc}).
# All 8 combinations of s_{ja}, s_{jb}, s_{jc} ∈ {+1,-1} give
# 8 sign patterns for (ε_a, ε_b, ε_c).

found_zero = False
zero_examples = []

for Da in D_values:
    for Db in D_values:
        for Dc in D_values:
            Dbar = lcm(Da, lcm(Db, Dc))
            la = (Dbar // Da) ** 4
            lb = (Dbar // Db) ** 4
            lc = (Dbar // Dc) ** 4

            for signs in iprod([1, -1], repeat=3):
                ea, eb, ec = signs
                S = ea * la + eb * lb + ec * lc
                if S == 0:
                    found_zero = True
                    zero_examples.append((Da, Db, Dc, la, lb, lc, signs))

if found_zero:
    print(f"  FOUND {len(zero_examples)} zero cases!")
    for Da, Db, Dc, la, lb, lc, signs in zero_examples[:20]:
        ea, eb, ec = signs
        print(f"    D=({Da},{Db},{Dc}), λ=({la},{lb},{lc}), "
              f"signs=({ea:+d},{eb:+d},{ec:+d}), "
              f"S={ea*la}+{eb*lb}+{ec*lc}=0")
else:
    print(f"  No S_j = 0 found for any D triple and any sign pattern.")
    print(f"  (Searched {len(D_values)}³ × 8 = {len(D_values)**3 * 8} cases)")

# Theoretical analysis
print(f"\n{'='*60}")
print(f"  Theoretical analysis")
print(f"{'='*60}")

print("""
  S_j = ε_a·n_a⁴ + ε_b·n_b⁴ + ε_c·n_c⁴
  where n_x = D̄/D_x ∈ Z_{>0} and ε_x ∈ {+1,-1}.

  For S_j = 0: need ε_a·n_a⁴ + ε_b·n_b⁴ = -ε_c·n_c⁴.
  Since n_x ≥ 1, this means a sum/difference of two 4th powers
  equals ±(a 4th power).

  Cases:
  (+,+,-): n_a⁴ + n_b⁴ = n_c⁴  → Fermat for exponent 4 → IMPOSSIBLE
  (+,-,+): n_a⁴ - n_b⁴ = n_c⁴  → n_a⁴ = n_b⁴ + n_c⁴ → Fermat → IMPOSSIBLE
  (-,+,+): same as (+,-,+) by symmetry → IMPOSSIBLE
  (+,-,-): n_a⁴ = n_b⁴ + n_c⁴  → Fermat → IMPOSSIBLE
  (-,+,-): n_b⁴ = n_a⁴ + n_c⁴  → Fermat → IMPOSSIBLE
  (-,-,+): n_c⁴ = n_a⁴ + n_b⁴  → Fermat → IMPOSSIBLE
  (+,+,+): n_a⁴ + n_b⁴ + n_c⁴ > 0 → NEVER ZERO
  (-,-,-): -(n_a⁴ + n_b⁴ + n_c⁴) < 0 → NEVER ZERO

  ALL CASES: S_j ≠ 0, by Fermat's Last Theorem for exponent 4!
  (FLT for n=4: x⁴ + y⁴ = z⁴ has no positive integer solutions,
   proved by Fermat himself using infinite descent.)
""")

print(f"  CONCLUSION:")
print(f"  S_j = Σ λ_ε σ_ε s_{{jε}} ≠ 0 for ALL valid (D_a,D_b,D_c)")
print(f"  because λ_ε = (D̄/D_ε)⁴ are FOURTH POWERS of positive")
print(f"  integers, and no sum/difference of three fourth powers")
print(f"  of positive integers can vanish (by FLT for n=4).")
print(f"")
print(f"  The reviewer's counterexample λ_c = 2 is NOT a valid λ_ε")
print(f"  because 2 is not a perfect fourth power.")
print(f"")
print(f"  Remark 4.10 should be revised to note this structure.")
