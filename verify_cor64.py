#!/usr/bin/env python3
"""
Verify Corollary 6.4: E(3) rational points → quartic (t₁, y).

E(3):  Y² = X(X-9)(X-8),  twist: W² = X''(X''-36)(X''-32)
Map:   X'' = 4·X_E, W = 8·Y_E  (since roots scale by 4)

Correct birational map (quartic ↔ twist):
  Forward:  m = (y-1)/t₁², X'' = -2m+34, W = 2t₁(m²-1)
  Inverse:  m = (34-X'')/2 = 17-2X_E
            t₁² = X_E / ((X_E-8)(X_E-9))
            y = m·t₁² + 1
"""
from fractions import Fraction as F

E3_points = [(0,0),(6,6),(6,-6),(8,0),(9,0),(12,12),(12,-12)]

print(f"{'X_E':>4}  {'Y_E':>4}  {'m':>4}  {'t₁²':>12}  {'y':>5}  {'quartic✓':>8}  {'t₁>1':>5}")
print("-" * 55)
for X_E, Y_E in E3_points:
    m = 17 - 2*X_E
    d = (X_E - 8)*(X_E - 9)
    if d == 0:
        print(f"{X_E:4}  {Y_E:4}  {m:4}  {'pole':>12}  {'—':>5}  {'—':>8}  {'no':>5}")
        continue
    t1sq = F(X_E, d)
    y = m * t1sq + 1
    ok = (y*y == t1sq**2 + 34*t1sq + 1)
    adm = "no" if t1sq <= 1 else "YES"
    print(f"{X_E:4}  {Y_E:4}  {m:4}  {str(t1sq):>12}  {str(y):>5}  "
          f"{'✓' if ok else '✗':>8}  {adm:>5}")

print(f"{'∞':>4}  {'∞':>4}  {'—':>4}  {'t₁=0':>12}  {'±1':>5}  {'✓':>8}  {'no':>5}")
print()
print("All cases: t₁ ∈ {0, ±1, pole}. None has t₁ > 1.  ✓")
