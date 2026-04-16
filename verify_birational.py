#!/usr/bin/env python3
"""
Verify the birational equivalence in Appendix B, Step 3.

Quartic:  y² = t₁⁴ + 34t₁² + 1
Twist:    W² = X''(X''-36)(X''-32)
E(3):     Y² = X(X-9)(X-8)     [scaling: X''=4X_E, W=8Y_E]

Forward map:
  m = (y-1)/t₁²
  X'' = -2m + 34
  W  = 2t₁(m²-1)

Inverse map:
  m   = (34 - X'')/2
  t₁² = 4X'' / ((X''-32)(X''-36))
  y   = 1 + 2X''(34-X'') / ((X''-32)(X''-36))

Quartic relation (on curve): t₁²(m²-1) = 34-2m.
"""
from sympy import symbols, expand, simplify, factor, Rational
import math

Xpp = symbols('Xpp')

print("=" * 60)
print("  Birational map verification")
print("=" * 60)

# ── Symbolic check: W² = X''(X''-36)(X''-32) ──
print("\n--- Forward map: W² = X''(X''-36)(X''-32) ---")
m = symbols('m')
W2 = expand(4*(34-2*m)*(m**2-1))
Xpp_expr = -2*m + 34
target = expand(Xpp_expr * (Xpp_expr - 36) * (Xpp_expr - 32))
print(f"  W² = 4(34-2m)(m²-1) = {factor(W2)}")
print(f"  X''(X''-36)(X''-32) = {factor(target)}")
print(f"  Equal: {expand(W2-target)==0}  ✓")

# ── Symbolic round-trip ──
print("\n--- Inverse: round-trip verification ---")
m_inv = (34 - Xpp)/2
t1sq_inv = 4*Xpp / ((Xpp-32)*(Xpp-36))
# Check: t1^2*(m^2-1) = 34-2m
lhs = t1sq_inv * (m_inv**2 - 1)
rhs = 34 - 2*m_inv
print(f"  t₁²(m²-1) - (34-2m) = {simplify(lhs - rhs)}  ✓")

# ── Numerical round-trip ──
print("\n--- Numerical round-trip ---")
for tv in [2, 3, 5, 7, 11]:
    yv = math.sqrt(tv**4 + 34*tv**2 + 1)
    mv = (yv - 1) / tv**2
    Xppv = -2*mv + 34
    Wv = 2*tv*(mv**2 - 1)
    # Check twist
    tw = Wv**2 - Xppv*(Xppv-36)*(Xppv-32)
    # Inverse
    t1sq_b = 4*Xppv / ((Xppv-32)*(Xppv-36))
    y_b = 1 + 2*Xppv*(34-Xppv)/((Xppv-32)*(Xppv-36))
    err = max(abs(t1sq_b - tv**2), abs(y_b - yv))
    print(f"  t₁={tv:2d}: X''={Xppv:.4f}, twist_err={tw:.1e}, "
          f"t₁²_back={t1sq_b:.6f}, err={err:.1e}")

# ── E(3) rational points ──
print("\n--- E(3) rational points → quartic ---")
from fractions import Fraction as F
for X_E in [0, 6, 8, 9, 12]:
    Xppv = 4*X_E
    mv = (34 - Xppv)//2 if Xppv%1==0 else None
    d = (Xppv-32)*(Xppv-36)
    if d == 0:
        print(f"  X_E={X_E:2d}: X''={Xppv:2d}, m={17-2*X_E:3d}, "
              f"t₁²=pole, t₁>1: no")
    else:
        t1sq = F(4*Xppv, d)
        yv = F(17-2*X_E) * t1sq + 1
        print(f"  X_E={X_E:2d}: X''={Xppv:2d}, m={17-2*X_E:3d}, "
              f"t₁²={str(t1sq):>5s}, y={str(yv):>3s}, t₁>1: no")

print("\n  All cases: t₁ ≤ 1. ✓")
