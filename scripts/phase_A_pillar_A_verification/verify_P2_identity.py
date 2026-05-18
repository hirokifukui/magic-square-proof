"""
Pillar A Step P2 verification — Stage 1: algebraic identity Delta_l = 2 p^(2l) J_(e-l)

Setup:
  p ≡ 1 mod 4, split as p = π·\bar{π} in Z[i], π = u + vi with u^2 + v^2 = p.
  N := 2 p^(2e). Then |w|^2 = N has solutions w (up to Gaussian units and conjugation)
  parameterised by a ∈ {0, 1, ..., 2e}, with w = (1+i) · π^a · \bar{π}^(2e-a).
  
  Symmetry: a and 2e-a give conjugate W (same {U^2, V^2}, hence same R_l = min(a, 2e-a) class).
  Centre: a = e gives w = (1+i) π^e \bar{π}^e = (1+i) p^e, so W = p^e (1+i),
  i.e. (U, V) = (p^e, p^e), the degenerate "centre rep" R_e excluded by P1.

  For l = min(a, 2e-a) < e (i.e. a < e WLOG), write a = e - m with m = e - a ≥ 1.
  Then w = (1+i) π^(e-m) \bar{π}^(e+m) = (1+i) p^(e-m) \bar{π}^(2m)
  So U + V i = w / (some unit), and writing W := U^2, W' := V^2 we want to show
  
      Delta_l := U^2 - V^2 = 2 p^(2l) · J_(e-l)
  
  where l = e - m, so e - l = m. The claim becomes
      Delta = 2 p^(2(e-m)) · J_m,  J_m := |Im(\bar{π}^(4m))|.

This script:
  (1) constructs π symbolically as u + v i with u^2 + v^2 = p as a formal relation,
  (2) computes w = (1+i) π^(e-m) \bar{π}^(e+m),
  (3) extracts U = Re(w), V = Im(w) (with sign convention U ≥ V ≥ 0 left implicit),
  (4) computes U^2 - V^2 and shows it equals 2 p^(2(e-m)) · Im(\bar{π}^(4m)) (up to sign).

We verify for (e, m) ∈ {(2,1), (3,1), (3,2), (4,1), (4,2), (4,3),
  (5,1), (5,2), (5,3), (5,4), (6,1), (6,3), (6,5)}
symbolically (treating u, v as free symbols with p = u^2 + v^2 as a free parameter).
"""

import sympy as sp

u, v = sp.symbols('u v', real=True)
p_sym = u**2 + v**2  # p as polynomial in u, v
i = sp.I

pi = u + v*i           # Gaussian prime
pi_bar = u - v*i       # conjugate

def compute_Delta_and_RHS(e, m):
    """
    Returns (LHS, RHS, diff_simplified) where
      LHS = U^2 - V^2 with w = (1+i) · pi^(e-m) · pi_bar^(e+m), W = w (no unit normalisation here)
      RHS = 2 p^(2(e-m)) · Im(pi_bar^(4m))
    """
    w = (1 + i) * (pi**(e - m)) * (pi_bar**(e + m))
    w = sp.expand(w)
    U = sp.re(w)
    V = sp.im(w)
    LHS = sp.expand(U**2 - V**2)
    
    # RHS: 2 p^(2(e-m)) · Im(pi_bar^(4m))
    pi_bar_4m = sp.expand(pi_bar**(4*m))
    J_m_signed = sp.im(pi_bar_4m)   # Im(\bar{π}^(4m)), signed
    RHS = 2 * (p_sym ** (2*(e - m))) * J_m_signed
    RHS = sp.expand(RHS)
    
    # Test both signs (the absolute value in J_m means LHS = ± RHS is OK)
    diff_plus = sp.expand(LHS - RHS)
    diff_minus = sp.expand(LHS + RHS)
    
    return LHS, RHS, diff_plus, diff_minus


print("=" * 80)
print("Stage 1: Verify Delta_l = ±2 p^(2l) · Im(\\bar{π}^(4m))  symbolically")
print("=" * 80)
print()
print(f"{'(e,m)':>8} | {'LHS - RHS':>15} | {'LHS + RHS':>15} | conclusion")
print("-" * 75)

results = []
for (e, m) in [(2,1), (3,1), (3,2), (4,1), (4,2), (4,3),
               (5,1), (5,2), (5,3), (5,4),
               (6,1), (6,3), (6,5)]:
    LHS, RHS, dplus, dminus = compute_Delta_and_RHS(e, m)
    if dplus == 0:
        sign = "+"
        ok = True
    elif dminus == 0:
        sign = "-"
        ok = True
    else:
        sign = "?"
        ok = False
    print(f"  ({e},{m})  | {str(dplus)[:15]:>15} | {str(dminus)[:15]:>15} |  Delta = {sign} RHS  [{'OK' if ok else 'FAIL'}]")
    results.append((e, m, ok))

all_ok = all(r[2] for r in results)
print()
print(f"All cases pass: {all_ok}")
print()

# Numerical check with actual split primes p ≡ 1 mod 4
print("=" * 80)
print("Stage 2: Numerical confirmation with real split primes")
print("=" * 80)
print()
# p = 5: π = 2 + i (u=2, v=1)
# p = 13: π = 3 + 2i (u=3, v=2)
# p = 17: π = 4 + i (u=4, v=1)
# p = 29: π = 5 + 2i
# p = 37: π = 6 + i

primes_data = [(5, 2, 1), (13, 3, 2), (17, 4, 1), (29, 5, 2), (37, 6, 1), (41, 5, 4), (53, 7, 2)]

print(f"{'p':>4} {'(u,v)':>10} {'(e,m)':>8} {'|Delta|':>20} {'2 p^(2l) |J_m|':>20} {'match'}")
print("-" * 85)
all_num_ok = True
for (p_val, u_val, v_val) in primes_data:
    assert u_val**2 + v_val**2 == p_val
    for (e, m) in [(2,1), (3,1), (3,2), (4,2), (5,3)]:
        l = e - m
        # Construct w = (1+i) π^(e-m) \bar{π}^(e+m)
        from sympy import I as Im
        pi_num = u_val + v_val * Im
        pi_bar_num = u_val - v_val * Im
        w_num = (1 + Im) * (pi_num ** (e - m)) * (pi_bar_num ** (e + m))
        w_num = sp.expand(w_num)
        U_num = int(sp.re(w_num))
        V_num = int(sp.im(w_num))
        Delta_num = U_num**2 - V_num**2
        
        # J_m = |Im(pi_bar^(4m))|
        pi_bar_4m_num = sp.expand(pi_bar_num ** (4*m))
        J_m = abs(int(sp.im(pi_bar_4m_num)))
        RHS_num = 2 * (p_val ** (2*l)) * J_m
        
        match = (abs(Delta_num) == RHS_num)
        all_num_ok = all_num_ok and match
        print(f"  {p_val:>3} ({u_val},{v_val})    ({e},{m})  {abs(Delta_num):>18} {RHS_num:>18}     {'OK' if match else 'FAIL'}")

print()
print(f"Numerical verification all pass: {all_num_ok}")
print()
print("=" * 80)
print(f"STAGE 1 RESULT: identity Delta_l = ±2 p^(2l) · Im(\\bar{{π}}^(4m)) verified symbolically and numerically")
print("=" * 80)
