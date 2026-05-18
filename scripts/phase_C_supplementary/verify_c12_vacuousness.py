"""
verify_c12_vacuousness.py
=========================

Verifies Proposition `prop:c12-vacuousness` (Session 21):

For every k >= 2 and every admissible rational specialization t_1 -> rho
(i.e., rho not in {0, +1, -1}), the family
    {Phi_S(rho, t_J\{1})}_{|S| odd, S ⊆ J}
remains Q-linearly independent in Q[t_J\{1}], where J = {1, 3, 4, ..., k}
and Phi_S = prod_{j∈S} t_j * prod_{l∈J\S} (t_l^2 - 1) (eq. chiPhi-def in §5.5).

This independence, combined with Step 1 of Lemma C3, implies that
c_{12}(rho, t_3, ..., t_k) ≢ 0 in Q[t_3, ..., t_k] for any admissible
rho, since the same FLT-type constraint that would force vanishing reduces,
at lambda_a = lambda_b = lambda_c = 1, to
    tau_a + tau_b - tau_c = 0  in {±1}^3,
which has no integer solution.

VERIFICATION STRATEGY (3 claims):
  Claim L (leading monomial structure):
    For each odd-card S ⊆ J, Phi_S(rho, t_rest) has a UNIQUE leading monomial
    of maximal total degree, with non-zero coefficient (rho if 1∈S, rho²-1
    if 1∉S).
  Claim B (bijection):
    The map S → S △ {1} restricted to odd-card S ⊆ J is a bijection onto
    the power set of J\{1}.
  Claim I (independence):
    The Q-rank of {Phi_S(rho, t_rest)}_{|S| odd} equals 2^{k-2} (full).

  Claim Z (Case (D) vacuousness, direct test):
    For each k=3 non-complementary triple, compute c_{12}(t_1, t_3) and
    verify that no rational rho makes c_{12}(rho, t_3) ≡ 0 in Q[t_3].

All claims tested for k ∈ {3, 4, 5, 6, 7} (and Claim Z exhaustively for k=3).

Author: Champagne project, Session 21 (2026-05-17)
"""

import sympy as sp
from itertools import combinations, product


def phi_S_specialized(S, J, t_dict, rho_val):
    """Phi_S(rho, t_rest) with t_1 -> rho_val.

    Phi_S = prod_{j∈S} t_j * prod_{l∈J\S} (t_l^2 - 1).
    """
    f = sp.Integer(1)
    for j in J:
        if j == 1:
            if j in S:
                f *= rho_val
            else:
                f *= (rho_val ** 2 - 1)
        else:
            if j in S:
                f *= t_dict[j]
            else:
                f *= (t_dict[j] ** 2 - 1)
    return sp.expand(f)


def leading_monomial_data(phi_spec, rest_syms):
    """Return (max total degree, list of (monom, coeff) at that degree, full Poly)."""
    if phi_spec == 0:
        return (0, [], None)
    p = sp.Poly(phi_spec, *rest_syms)
    monoms = p.monoms()
    if not monoms:
        return (0, [], p)
    total_degs = [sum(m) for m in monoms]
    max_deg = max(total_degs)
    max_monoms = [(m, c) for m, c in zip(monoms, p.coeffs()) if sum(m) == max_deg]
    return (max_deg, max_monoms, p)


def claim_L_test(k, rho_val):
    """Test Claim L for given (k, rho). Returns (all_unique, all_distinct, all_nonzero)."""
    J = [1] + list(range(3, k + 1))
    rest = [j for j in J if j != 1]
    t_syms = {j: sp.Symbol(f't{j}') for j in J}
    rest_syms = [t_syms[j] for j in rest]

    odd_subsets = []
    for r in range(1, len(J) + 1, 2):
        for S in combinations(J, r):
            odd_subsets.append(set(S))

    leading_data = []
    for S in odd_subsets:
        phi_spec = phi_S_specialized(S, J, t_syms, rho_val)
        deg, max_monoms, _ = leading_monomial_data(phi_spec, rest_syms)
        leading_data.append((S, deg, max_monoms, phi_spec))

    all_unique_lm = all(len(d[2]) == 1 for d in leading_data)
    if not all_unique_lm:
        return (False, False, False, leading_data)

    lms = [d[2][0][0] for d in leading_data]
    coeffs = [d[2][0][1] for d in leading_data]
    lms_distinct = len(set(lms)) == len(lms)
    all_nonzero = all(c != 0 for c in coeffs)

    return (all_unique_lm, lms_distinct, all_nonzero, leading_data)


def claim_B_test(k):
    """Verify the bijection S ↔ S△{1} from odd-card subsets of J to all subsets of J\{1}."""
    J = [1] + list(range(3, k + 1))
    rest = [j for j in J if j != 1]

    odd_subsets = []
    for r in range(1, len(J) + 1, 2):
        for S in combinations(J, r):
            odd_subsets.append(set(S))

    images = []
    for S in odd_subsets:
        if 1 in S:
            images.append(frozenset(S - {1}))
        else:
            images.append(frozenset(S))

    expected = set(
        frozenset(s) for r in range(0, len(rest) + 1) for s in combinations(rest, r)
    )
    actual = set(images)
    return (actual == expected) and (len(images) == len(actual))


def claim_I_test(k, rho_val, leading_data):
    """Compute Q-linear rank of {Phi_S(rho, t_rest)}."""
    J = [1] + list(range(3, k + 1))
    rest = [j for j in J if j != 1]
    t_syms = {j: sp.Symbol(f't{j}') for j in J}
    rest_syms = [t_syms[j] for j in rest]

    all_monoms = set()
    for S, deg, max_monoms, phi_spec in leading_data:
        if phi_spec == 0:
            continue
        p = sp.Poly(phi_spec, *rest_syms)
        for m in p.monoms():
            all_monoms.add(m)
    all_monoms = sorted(all_monoms)
    rows = []
    for S, deg, max_monoms, phi_spec in leading_data:
        if phi_spec == 0:
            rows.append([0] * len(all_monoms))
            continue
        p = sp.Poly(phi_spec, *rest_syms)
        cd = dict(zip(p.monoms(), p.coeffs()))
        rows.append([cd.get(m, 0) for m in all_monoms])
    return sp.Matrix(rows).rank(), len(leading_data)


# ----------------------------------------------------------------------
# Claim Z: direct exhaustive check at k=3
# ----------------------------------------------------------------------

def r_eps_k3(s1, s2, s3, t1, t2, t3):
    I = sp.I
    W = (t1 + s1 * I) * (t2 + s2 * I) * (t3 + s3 * I)
    return sp.expand(sp.im(W ** 2))


def get_c12_k3(eps_a, eps_b, eps_c, t1, t2, t3):
    sa = tuple((-1) ** b for b in eps_a)
    sb = tuple((-1) ** b for b in eps_b)
    sc = tuple((-1) ** b for b in eps_c)
    ra = r_eps_k3(*sa, t1, t2, t3)
    rb = r_eps_k3(*sb, t1, t2, t3)
    rc = r_eps_k3(*sc, t1, t2, t3)
    Cond1 = sp.expand(ra + rb - rc)
    Cond1_poly = sp.Poly(Cond1, t2)
    return sp.expand(Cond1_poly.coeff_monomial(t2 ** 2))


def claim_Z_test_k3():
    """For all non-complementary k=3 triples, verify Case (D) is empty."""
    t1, t2, t3 = sp.symbols('t1 t2 t3', real=True)
    bitstrings = list(product([0, 1], repeat=3))
    all_ones = (1, 1, 1)
    problematic = []
    total = 0
    for ea in bitstrings:
        for eb in bitstrings:
            for ec in bitstrings:
                xor_abc = tuple((a + b + c) % 2 for a, b, c in zip(ea, eb, ec))
                if xor_abc == all_ones:
                    continue
                total += 1
                c12 = get_c12_k3(ea, eb, ec, t1, t2, t3)
                if c12 == 0:
                    problematic.append((ea, eb, ec, "IDENTICALLY_ZERO"))
                    continue
                c12_in_t3 = sp.Poly(c12, t3)
                coeffs = c12_in_t3.all_coeffs()
                common = None
                for c in coeffs:
                    c_exp = sp.expand(c)
                    if c_exp == 0:
                        continue
                    if not c_exp.has(t1):
                        common = set()
                        break
                    roots = sp.solve(c_exp, t1)
                    rational_roots = set(r for r in roots if r.is_rational)
                    if common is None:
                        common = rational_roots
                    else:
                        common = common & rational_roots
                    if common == set():
                        break
                if common and len(common) > 0:
                    problematic.append((ea, eb, ec, common))
    return total, problematic


# ----------------------------------------------------------------------
# Main test runner
# ----------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 72)
    print(" verify_c12_vacuousness.py")
    print(" Verification of Proposition prop:c12-vacuousness (Session 21)")
    print("=" * 72)
    print()

    # Test rho values: representative non-boundary rationals
    test_rhos = [sp.Rational(2, 1), sp.Rational(3, 2), sp.Rational(7, 4),
                 sp.Rational(-2, 1), sp.Rational(-3, 2), sp.Rational(11, 5)]

    # ---- Claim L, B, I for k = 3, 4, 5, 6, 7 ----
    print("[Claim L, B, I]  Leading-monomial, bijection, and independence")
    print()
    all_pass = True
    for k in [3, 4, 5, 6, 7]:
        print(f"  k = {k} (|J| = {k - 1}, 2^(k-2) = {2 ** (k - 2)} odd subsets):")
        B_ok = claim_B_test(k)
        print(f"    Claim B (bijection S ↔ S△{{1}}):              {'PASS' if B_ok else 'FAIL'}")
        if not B_ok:
            all_pass = False
        for rho_val in test_rhos:
            unique, distinct, nonzero, ldata = claim_L_test(k, rho_val)
            rank, total = claim_I_test(k, rho_val, ldata)
            L_ok = unique and distinct and nonzero
            I_ok = (rank == total)
            print(f"    rho = {str(rho_val):>6s}: L = {'PASS' if L_ok else 'FAIL'}  "
                  f"(unique LM: {unique}, distinct: {distinct}, nonzero coeff: {nonzero});  "
                  f"I = {'PASS' if I_ok else 'FAIL'} (rank {rank}/{total})")
            if not (L_ok and I_ok):
                all_pass = False
        print()

    # ---- Claim Z: exhaustive k=3 ----
    print("[Claim Z]  Exhaustive Case (D) check for k=3 non-complementary triples")
    total, problematic = claim_Z_test_k3()
    print(f"    triples scanned: {total}")
    print(f"    problematic (Case D non-empty): {len(problematic)}")
    if problematic:
        all_pass = False
        for p in problematic[:5]:
            print(f"      {p}")
        print("    FAIL")
    else:
        print("    PASS: every non-complementary k=3 triple has empty Case (D).")
    print()

    print("=" * 72)
    if all_pass:
        print(" OVERALL: PASS")
        print(" Proposition prop:c12-vacuousness verified for k ∈ {3,4,5,6,7}")
        print(" and exhaustively for k=3 non-complementary triples.")
    else:
        print(" OVERALL: FAIL — see details above.")
    print("=" * 72)
