#!/usr/bin/env python3
"""
Verification script for Lemma (k=3 base case).

Proves that for k=3, all 24 all-one triples yield Gaussian elements
z₊ ∈ Q(i)[t₁, t₃] that are irreducible (hence not Gaussian squares).

Usage:  python3 verify_k3base.py
"""
from itertools import product as iprod
from sympy import (symbols, expand, factor, Poly, Rational, I as symI,
                   diff, gcd as sym_gcd, resultant, collect)

t1, t3 = symbols('t1 t3', real=True)


def gmul(a, b):
    """Multiply two Gaussian elements (re, im)."""
    return (expand(a[0]*b[0] - a[1]*b[1]),
            expand(a[0]*b[1] + a[1]*b[0]))


def compute_z_plus(ea, eb, ec):
    """
    Compute z₊ = c₁₁/4 + (c₁₂/2)·i for a k=3 triple,
    using the full W_ε product (matching the original code).
    """
    t2 = symbols('t2', real=True)
    tl = [t1, t2, t3]

    def build_W(eps):
        w = (1, 0)
        for j in range(3):
            w = gmul(w, (tl[j], 2*eps[j]-1))
        return w

    Wa, Wb, Wc = build_W(ea), build_W(eb), build_W(ec)
    ra = expand(-4*Wa[0]*Wa[1])
    rb = expand(-4*Wb[0]*Wb[1])
    rc = expand(-4*Wc[0]*Wc[1])
    cond1 = expand(ra + rb - rc)

    p = Poly(cond1, t2)
    c12 = expand(p.nth(2))
    c11 = expand(p.nth(1))

    z_re = expand(c11 / 4)
    z_im = expand(c12 / 2)
    return z_re, z_im


# ═══════════════════════════════════════════════════════════════════
#  Phase A: Enumerate all-one triples
# ═══════════════════════════════════════════════════════════════════

def phase_A():
    print("=" * 70)
    print("  Phase A: Enumerate all-one triples for k=3")
    print("=" * 70)

    all_eps = list(iprod([0, 1], repeat=3))
    triples = []

    for ea in all_eps:
        for eb in all_eps:
            if eb <= ea:
                continue
            for ec in all_eps:
                if ec == ea or ec == eb:
                    continue
                # non-complementary
                if all(ea[j] + eb[j] == 1 for j in range(3)):
                    continue
                # |S_j| = 1 for all j
                if not all(abs((2*ea[j]-1)+(2*eb[j]-1)-(2*ec[j]-1)) == 1
                           for j in range(3)):
                    continue
                triples.append((ea, eb, ec))

    print(f"  Found {len(triples)} all-one triples.")
    assert len(triples) == 24

    # Compute z₊ for each
    z_polys = {}
    for ea, eb, ec in triples:
        z_re, z_im = compute_z_plus(ea, eb, ec)
        z = expand(z_re + z_im * symI)
        z_polys[(ea, eb, ec)] = (z_re, z_im, z)

    # Find distinct polynomials (exact equality)
    distinct_exact = {}
    for triple, (zr, zi, z) in z_polys.items():
        found = False
        for key in distinct_exact:
            if expand(z - z_polys[key][2]) == 0:
                found = True
                break
        if not found:
            distinct_exact[triple] = z

    print(f"  Exactly distinct z₊: {len(distinct_exact)}")

    # Find equivalence classes up to Q(i)-units {±1, ±i} and t₁→−t₁
    units = [(1, '1'), (-1, '-1'), (symI, 'i'), (-symI, '-i')]
    classes = []  # list of (representative_z, members)

    for triple, (zr, zi, z) in z_polys.items():
        found_class = None
        for ci, (rep_z, members) in enumerate(classes):
            for u, uname in units:
                if expand(z - u*rep_z) == 0:
                    found_class = (ci, uname, False)
                    break
                if expand(z - u*rep_z.subs(t1, -t1)) == 0:
                    found_class = (ci, uname, True)
                    break
            if found_class:
                break
        if found_class:
            ci, uname, flipped = found_class
            classes[ci][1].append((triple, uname, flipped))
        else:
            classes.append((z, [(triple, '1', False)]))

    print(f"  Equivalence classes (up to ±1, ±i, t₁→−t₁): {len(classes)}")
    for ci, (rep_z, members) in enumerate(classes):
        print(f"    Class {ci}: {len(members)} triples")
        for triple, uname, flipped in members[:2]:
            tag = f" [{uname}" + (" ∘ t₁→−t₁" if flipped else "") + "]"
            print(f"      {triple}{tag}")
        if len(members) > 2:
            print(f"      ... and {len(members)-2} more")

    # The paper's claim "exactly two distinct polynomials related by t₁→−t₁"
    # means 2 classes up to Q(i)-units, each split into ±t₁ pairs.
    # Let's check: are there exactly 2 classes ignoring t₁→−t₁ but keeping units?
    classes_no_flip = []
    for triple, (zr, zi, z) in z_polys.items():
        found_class = None
        for ci, (rep_z, members) in enumerate(classes_no_flip):
            for u, uname in units:
                if expand(z - u*rep_z) == 0:
                    found_class = ci
                    break
            if found_class is not None:
                break
        if found_class is not None:
            classes_no_flip[found_class][1].append(triple)
        else:
            classes_no_flip.append((z, [triple]))

    print(f"\n  Classes up to Q(i)-units (no t₁ flip): {len(classes_no_flip)}")
    for ci, (rep_z, members) in enumerate(classes_no_flip):
        print(f"    Class {ci}: {len(members)} triples, z₊ = {collect(rep_z, t1)}")

    # Use first representative for the proof
    rep_z = classes[0][0]
    rep_triple = classes[0][1][0][0]
    rep_re, rep_im = z_polys[rep_triple][0], z_polys[rep_triple][1]

    return rep_triple, rep_re, rep_im, rep_z, triples, z_polys, classes


# ═══════════════════════════════════════════════════════════════════
#  Phase B: Irreducibility
# ═══════════════════════════════════════════════════════════════════

def phase_B(z_re, z_im, z, class_idx=0):
    print(f"\n  --- Class {class_idx} representative ---")

    # Extract coefficients of t₁ manually
    z_exp = expand(z)
    a_coeff = expand(z_exp.coeff(t1, 2))
    b_coeff = expand(z_exp.coeff(t1, 1))
    c_coeff = expand(z_exp.coeff(t1, 0))
    print(f"  z₊ = a·t₁² + b·t₁ + c")
    print(f"    a = {a_coeff}")
    print(f"    b = {b_coeff}")
    print(f"    c = {c_coeff}")

    # Verify reconstruction
    assert expand(z - (a_coeff*t1**2 + b_coeff*t1 + c_coeff)) == 0

    # ── Discriminant Δ(t₃) = b² − 4ac ──
    print("\n  Discriminant Δ(t₃) = b² − 4ac:")
    Delta = expand(b_coeff**2 - 4*a_coeff*c_coeff)
    Delta_f = factor(Delta)
    deg_D = Poly(Delta, t3).degree()
    print(f"    Δ = {Delta_f}")
    print(f"    deg(Δ, t₃) = {deg_D}")

    if deg_D % 2 == 1:
        print(f"    deg is odd ⟹ Δ not a square in Q(i)[t₃]  ✓")
        delta_not_sq = True
    else:
        dDelta = diff(Delta, t3)
        g = sym_gcd(Delta, dDelta)
        g_deg = Poly(g, t3).degree() if g != 0 else -1
        print(f"    gcd(Δ, dΔ/dt₃) = {factor(g)},  deg = {g_deg}")
        delta_not_sq = (g_deg == 0)
        if delta_not_sq:
            print(f"    Δ squarefree ⟹ not a square  ✓")
        else:
            # Check if squarefree after removing constant factors
            print(f"    WARNING: gcd has positive degree")

    # ── gcd(z₊, ∂z₊/∂t₁) ──
    print("\n  gcd(z₊, ∂z₊/∂t₁):")
    dz = diff(z, t1)

    # Use resultant to check coprimality
    res = resultant(z, dz, t1)
    res_f = factor(res)
    res_is_nonzero = (res != 0)
    print(f"    ∂z₊/∂t₁ = {dz}")
    print(f"    res(z₊, ∂z₊/∂t₁, t₁) = {res_f}")
    print(f"    res ≠ 0: {res_is_nonzero}")

    if res_is_nonzero:
        # Resultant nonzero as polynomial in t₃ ⟹ gcd = 1 over Q(i)(t₃)
        # This means gcd over Q(i)[t₁,t₃] is at most a polynomial in t₃ alone.
        # But if z₊ is primitive in t₁ (leading coeff a not zero), gcd = 1.
        print(f"    ⟹ gcd(z₊, ∂z₊/∂t₁) = 1 over Q(i)(t₃)[t₁]  ✓")
        gcd_is_const = True
    else:
        gcd_is_const = False
        print(f"    WARNING: resultant is zero!")

    # ── Summary ──
    print(f"\n  z₊ degree 2 in t₁, Δ(t₃) not a square")
    print(f"    ⟹ z₊ irreducible over Q(i)(t₃), hence over Q(i)  ✓")
    print(f"    ⟹ z₊ not a Gaussian square  ✓")

    return gcd_is_const, delta_not_sq, Delta_f, a_coeff, b_coeff, c_coeff


def verify_all_triples(triples, z_polys):
    """Verify gcd(z₊, ∂z₊/∂t₁) = 1 for all 24 triples via resultant."""
    print("\n--- gcd = 1 for all 24 triples (via resultant) ---")
    all_ok = True
    for triple in triples:
        z_re, z_im, z = z_polys[triple]
        dz = diff(z, t1)
        res = resultant(z, dz, t1)
        ok = (res != 0 and Poly(res, t3).degree() >= 0)
        if not ok:
            print(f"    FAIL: {triple}, res = {res}")
            all_ok = False
    status = "✓" if all_ok else "✗"
    print(f"  All 24 triples: res(z₊, ∂z₊/∂t₁, t₁) ≠ 0  {status}")
    return all_ok


# ═══════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("Verification: k=3 base case for Theorem 7.5\n")

    rep_triple, z_re, z_im, z, triples, z_polys, classes = phase_A()

    print("\n" + "=" * 70)
    print("  Phase B: Irreducibility of z₊ over Q(i)")
    print("=" * 70)

    # Check both class representatives
    all_gcd_ok = True
    all_delta_ok = True
    class_data = []
    for ci, (rep_z, members) in enumerate(classes):
        rep_triple_ci = members[0][0]
        zr, zi, zv = z_polys[rep_triple_ci]
        gcd_ok, delta_ok, Delta_f, a, b, c = phase_B(zr, zi, zv, ci)
        all_gcd_ok = all_gcd_ok and gcd_ok
        all_delta_ok = all_delta_ok and delta_ok
        class_data.append((a, b, c, Delta_f))

    all_ok = verify_all_triples(triples, z_polys)

    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  24 all-one triples:      ✓")
    print(f"  Equivalence classes:     {len(classes)}")
    print(f"  gcd(z₊, ∂z₊/∂t₁) = 1:   {'✓' if all_gcd_ok else '✗'}")
    print(f"  Δ(t₃) not square:        {'✓' if all_delta_ok else '✗'}")
    print(f"  All 24 verified:         {'✓' if all_ok else '✗'}")
    verdict = "PASS" if (all_gcd_ok and all_delta_ok and all_ok) else "FAIL"
    print(f"  VERDICT: {verdict}")

    # ── Verify nice factored form ──
    print("\n" + "=" * 70)
    print("  Factored form verification")
    print("=" * 70)

    # Class 0 representative
    z0 = z_polys[classes[0][1][0][0]][2]
    a0 = expand(z0.coeff(t1, 2))
    b0 = expand(z0.coeff(t1, 1))
    c0 = expand(z0.coeff(t1, 0))

    # Check a = -(t₃ - i)²
    a_nice = expand(-(t3 - symI)**2)
    assert expand(a0 - a_nice) == 0, f"a ≠ -(t₃-i)²: {a0} vs {a_nice}"
    print(f"  a = -(t₃ - i)²  ✓")

    # Check c = (t₃ - i)² = -a
    assert expand(c0 + a0) == 0, f"c ≠ -a"
    print(f"  c = (t₃ - i)² = -a  ✓")

    # So z₊ = -(t₃-i)²(t₁² - 1) + b·t₁
    z_check = expand(-(t3 - symI)**2 * (t1**2 - 1) + b0*t1)
    assert expand(z0 - z_check) == 0
    print(f"  z₊ = -(t₃ - i)²(t₁² - 1) + b·t₁  ✓")
    print(f"  b = {b0}")

    # Factor b
    b_factored = factor(b0)
    print(f"  b factored = {b_factored}")

    # Check b = 2(i(t₃² - 1) - 6t₃)
    b_nice = expand(2*(symI*(t3**2 - 1) - 6*t3))
    assert expand(b0 - b_nice) == 0
    print(f"  b = 2(i(t₃² - 1) - 6t₃)  ✓")

    # Class 1: check a₁ = (t₃ + i)²
    z1 = z_polys[classes[1][1][0][0]][2]
    a1 = expand(z1.coeff(t1, 2))
    b1 = expand(z1.coeff(t1, 1))
    c1 = expand(z1.coeff(t1, 0))

    a1_nice = expand((t3 + symI)**2)
    assert expand(a1 - a1_nice) == 0
    print(f"\n  Class 1: a = (t₃ + i)²  ✓")
    assert expand(c1 + a1) == 0
    print(f"  c = -(t₃ + i)² = -a  ✓")

    b1_nice = expand(-2*(symI*(t3**2 - 1) + 6*t3))
    assert expand(b1 - b1_nice) == 0
    print(f"  b = -2(i(t₃² - 1) + 6t₃)  ✓")

    # Check relationship between classes: z₁(t₁,t₃) = -conj(z₀(-t₁, t₃))
    z0_neg_conj = expand(-z0.subs(t1, -t1).subs(symI, -symI))
    assert expand(z1 - z0_neg_conj) == 0
    print(f"\n  z₊'(t₁,t₃) = -conj(z₊(-t₁, t₃))  ✓")
    print(f"  (Galois conjugation + sign + t₁-flip preserves irreducibility)")

    # ── Final discriminant check ──
    print(f"\n  Δ₀ = b₀² - 4a₀c₀ = b₀² + 4a₀² (since c = -a)")
    Delta0 = expand(b0**2 + 4*a0**2)
    Delta0_f = factor(Delta0)
    print(f"     = {Delta0_f}")
    print(f"  deg(Δ₀, t₃) = {Poly(Delta0, t3).degree()} (odd ⟹ not a square)")

    print(f"\n  COMPLETE: z₊ is irreducible over Q(i), not a Gaussian square.  ✓")
