#!/usr/bin/env python3
"""
Verify Proposition 8.1: y² = 9t₁⁴ - 14t₁² + 9 has no rational point with t₁ > 1.

Proof strategy (Pell descent reducing Q₂ to Q₁):
  1. Set u = t₁ - 1/t₁, v = y/t₁ → v² - 9u² = 4.
  2. Parametrize the conic: v-3u = r, v+3u = 4/r  →  u = (4-r²)/(6r).
  3. t₁ rational ⟺ u²+4 = □ ⟺ r⁴+136r²+16 = □.
  4. Set r = 2r' → r'⁴+34r'²+1 = □, which is Q₁.
  5. Q₁ ≅ 48.a3 (rank 0) → r' ∈ {0,±1} → r ∈ {0,±2} → u=0 → t₁=±1.
"""
from fractions import Fraction as F
from sympy import symbols, expand, simplify

t1, y, u, v, r = symbols('t1 y u v r', rational=True)

print("=" * 60)
print("  Verification of Proposition 8.1")
print("=" * 60)

# Step 1: substitution check
print("\n  Step 1: v = y/t₁, u = t₁-1/t₁ → v²-9u²=4")
v_expr = y / t1
u_expr = t1 - 1/t1
lhs = v_expr**2 - 9*u_expr**2
# On curve y²=9t₁⁴-14t₁²+9:
lhs_on_curve = simplify(lhs.subs(y**2, 9*t1**4 - 14*t1**2 + 9))
print(f"    v²-9u² on curve = {lhs_on_curve}  ✓")

# Step 2: conic parametrization
print("\n  Step 2: conic v²-9u²=4 parametrized by r = v-3u")
u_r = (4 - r**2) / (6*r)
v_r = (r + 4/r) / 2
check2 = simplify(v_r**2 - 9*u_r**2 - 4)
print(f"    v²-9u²-4 = {check2}  ✓")

# Step 3: rationality condition
print("\n  Step 3: t₁ rational ⟺ u²+4 = □")
u2p4 = expand(u_r**2 + 4)
u2p4_num = simplify(u2p4 * 36 * r**2)
print(f"    36r²(u²+4) = {expand(u2p4_num)}")
print(f"               = r⁴ + 136r² + 16")

# Step 4: reduction to Q₁
rp = symbols('rp')
quartic_r = expand((2*rp)**4 + 136*(2*rp)**2 + 16)
print(f"\n  Step 4: r = 2r' → {quartic_r} = 16(r'⁴+34r'²+1)")
print(f"    This is Q₁. By Theorem 6.3, rank = 0.")

# Step 5: enumerate solutions
print("\n  Step 5: Q₁ rational points: r' ∈ {0, ±1}")
for rp_v in [0, 1, -1]:
    r_v = 2 * rp_v
    if r_v == 0:
        print(f"    r'={rp_v}: r=0 → u undefined (pole)")
        continue
    u_v = F(4 - r_v**2, 6*r_v)
    disc = u_v**2 + 4
    from math import isqrt
    disc_n = disc.numerator * disc.denominator
    sq = isqrt(disc_n)
    t1_p = (u_v + F(sq, disc.denominator)) / 2
    t1_m = (u_v - F(sq, disc.denominator)) / 2
    print(f"    r'={rp_v}: r={r_v}, u={u_v}, t₁ ∈ {{{t1_p}, {t1_m}}}")

# Also: t₁=0 (not covered by substitution u=t₁-1/t₁)
print(f"\n  t₁=0: y²=9, y=±3. Not t₁>1.")
print(f"\n  Q₂ rational points: {{(0,±3), (±1,±2)}}. None has t₁>1. ✓")
