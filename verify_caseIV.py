#!/usr/bin/env python3
"""
Verify Case IV of Proposition 5.5.

Setting: j₀ ∈ C. W_a, W_b involve t_{j₀}; W_c does not.
We show Im(W_a⁴+W_b⁴) is non-constant in t_{j₀} unless ε_a, ε_b
are complementary (which gives R_a+R_b=0, handled separately).

Key check: the t_{j₀}⁴ and t_{j₀}³ coefficients of Im(W_a⁴+W_b⁴)
can't both vanish for non-complementary ε_a, ε_b.
"""
from sympy import symbols, expand, I as symI, im as sym_im, Poly, factor
from itertools import product as iprod

t, t2, t3 = symbols('t t2 t3', real=True)

def Im_poly(e):
    return expand(sym_im(expand(e)))

print("=" * 60)
print("  Case IV: degree of Im(W_a⁴+W_b⁴) in t_{j₀}")
print("=" * 60)

# k=3 with j₀ = variable t (first position)
# Other active vars: t2, t3
# Test all (s_a, s_b) × (P, Q) combinations
scenarios = [
    ("same sign, diff P/Q",  1,  1, (t2+symI)*(t3+symI), (t2+symI)*(t3-symI)),
    ("diff sign, same P/Q",  1, -1, (t2+symI)*(t3+symI), (t2+symI)*(t3+symI)),
    ("diff sign, diff P/Q",  1, -1, (t2+symI)*(t3+symI), (t2+symI)*(t3-symI)),
    ("complementary",        1, -1, (t2+symI)*(t3+symI), (t2-symI)*(t3-symI)),
]

for name, sa, sb, P, Q in scenarios:
    Wa = (t + sa*symI) * P
    Wb = (t + sb*symI) * Q
    Im_sum = Im_poly(Wa**4 + Wb**4)
    p = Poly(Im_sum, t)
    deg = p.degree()

    # Extract leading coefficients
    c4 = p.nth(4) if deg >= 4 else 0
    c3 = p.nth(3) if deg >= 3 else 0

    print(f"\n  {name}:")
    print(f"    s_a={sa:+d}, s_b={sb:+d}")
    print(f"    deg(Im, t) = {deg}")
    if deg >= 0:
        print(f"    coeff t⁴ = {factor(c4)}")
        print(f"    coeff t³ = {factor(c3)}")
    if deg == 0 and Im_sum == 0:
        print(f"    Im ≡ 0 → complementary case (R_a+R_b=0)")
    elif deg > 0:
        print(f"    Non-constant → contradicts Im(W_c⁴)=const  ✓")

print(f"\n{'='*60}")
print(f"  Exhaustive check: all k=3 pairs (ε_a,ε_b), ε_a≠ε_b")
print(f"{'='*60}")

all_eps = list(iprod([0, 1], repeat=3))
n_nonconst = 0
n_complementary = 0

for ea in all_eps:
    for eb in all_eps:
        if ea == eb:
            continue
        is_comp = all(ea[j] + eb[j] == 1 for j in range(3))

        # j₀ = 0 (first coordinate), P/Q from coords 1,2
        sa = (-1)**ea[0]
        sb = (-1)**eb[0]
        P = (t2 + (-1)**ea[1]*symI) * (t3 + (-1)**ea[2]*symI)
        Q = (t2 + (-1)**eb[1]*symI) * (t3 + (-1)**eb[2]*symI)
        Wa = (t + sa*symI) * P
        Wb = (t + sb*symI) * Q
        Im_sum = Im_poly(Wa**4 + Wb**4)
        p = Poly(Im_sum, t)
        deg = p.degree()

        if is_comp:
            assert deg <= 0 and Im_sum == 0, f"Complementary but Im≠0: {ea},{eb}"
            n_complementary += 1
        else:
            assert deg > 0, f"Non-complementary but Im is constant: {ea},{eb}, deg={deg}"
            n_nonconst += 1

print(f"  Non-complementary pairs with Im non-constant: {n_nonconst}")
print(f"  Complementary pairs with Im ≡ 0: {n_complementary}")
print(f"  Total pairs: {n_nonconst + n_complementary} (= 56 = 8·7)")
print(f"\n  ALL non-complementary pairs have Im(W_a⁴+W_b⁴) non-constant")
print(f"  in t_{{j₀}}.  Case IV argument verified.  ✓")
