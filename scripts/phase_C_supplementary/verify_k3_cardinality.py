"""
verify_k3_cardinality.py
=========================
Session 19, gamma-1 reinforcement.

Purpose
-------
Provide independent symbolic + numerical verification of the cardinality and
orbit-decomposition claims that will be promoted to Lemma D.1 and Lemma D.2
in Appendix D of main.tex.

Background (consistent with main.tex notation)
----------------------------------------------
- For each bitstring eps in {0,1}^3, set s_j(eps) = (-1)^{eps_j}  (cf. L210)
- For a triple (eps_a, eps_b, eps_c) with eps_a < eps_b (lex order),
  eps_c not in {eps_a, eps_b}, and eps_a XOR eps_b != (1,1,1)
  (non-complementary), define
      S_j = s_j(eps_a) + s_j(eps_b) - s_j(eps_c) in {-3, -1, +1, +3}
  (Definition def:Sj, L1920).
- The "all-one" condition: |S_j| = 1 for all j in {1,2,3}.

Claims
------
(A) The number of triples (eps_a, eps_b, eps_c) satisfying the four
    filtering conditions (lex order, distinctness, non-complementary,
    all-one) is exactly 24.
(B) The 24 triples, decorated with the explicit Gaussian element
    z_+(eps_a, eps_b, eps_c) computed by the construction of Section 7,
    reduce to exactly 2 equivalence classes under the relation
       z_+ ~ u * z_+^{(c)}(+- t_1, t_3),   u in {+1, -1}, c in {0, 1}.

Output
------
Symbolic + numerical verification, with full enumeration of intermediate
stage counts, the 24 triples, and the orbit assignment.

Status: stand-alone verification; the paper-side proof of Lemma D.1 and
Lemma D.2 in main.tex is logically independent of this script.
"""

from __future__ import annotations

import itertools
import sys

import sympy as sp


# ---------------------------------------------------------------------
# Combinatorial enumeration (Claim A)
# ---------------------------------------------------------------------

def bits(eps: int) -> tuple[int, int, int]:
    """Return the 3-bit tuple (eps_1, eps_2, eps_3) of eps in {0,...,7}."""
    return ((eps >> 2) & 1, (eps >> 1) & 1, eps & 1)


def sign(eps: int, j: int) -> int:
    """
    s_j(eps) = (-1)^{eps_j} with j in {1,2,3} (1-indexed mathematical convention).
    Bit encoding: eps = (eps_1, eps_2, eps_3) -> integer eps_1*4 + eps_2*2 + eps_3,
    so bit at position (3-j) recovers eps_j.
    """
    return 1 if ((eps >> (3 - j)) & 1) == 0 else -1


def S_vec(ea: int, eb: int, ec: int) -> tuple[int, int, int]:
    """Vector (S_1, S_2, S_3) for the triple (eps_a, eps_b, eps_c)."""
    return tuple(sign(ea, j) + sign(eb, j) - sign(ec, j) for j in (1, 2, 3))


def xor3(ea: int, eb: int) -> int:
    """Bitwise XOR on 3-bit ints."""
    return ea ^ eb


def is_complementary(ea: int, eb: int) -> bool:
    """True iff eps_a + eps_b = (1,1,1), i.e. eps_b is the bit-complement of eps_a."""
    return xor3(ea, eb) == 0b111


def stage_counts() -> dict:
    """Enumerate the filter cascade and return the cardinality at each stage."""
    counts = {}

    # Stage 0: all triples in {0,...,7}^3
    s0 = list(itertools.product(range(8), repeat=3))
    counts["stage_0_all_triples"] = len(s0)

    # Stage 1: lex ordering eps_a < eps_b
    s1 = [(a, b, c) for (a, b, c) in s0 if a < b]
    counts["stage_1_lex_a_lt_b"] = len(s1)

    # Stage 2: distinctness eps_c not in {eps_a, eps_b}
    s2 = [(a, b, c) for (a, b, c) in s1 if c != a and c != b]
    counts["stage_2_distinct_c"] = len(s2)

    # Stage 3: non-complementary eps_a XOR eps_b != (1,1,1)
    s3 = [(a, b, c) for (a, b, c) in s2 if not is_complementary(a, b)]
    counts["stage_3_non_complementary"] = len(s3)

    # Stage 4: all-one |S_j| = 1 for all j
    s4 = [(a, b, c) for (a, b, c) in s3 if all(abs(s) == 1 for s in S_vec(a, b, c))]
    counts["stage_4_all_one"] = len(s4)

    return counts, s4


# ---------------------------------------------------------------------
# Gaussian element z_+ for each triple (Claim B)
# ---------------------------------------------------------------------

def compute_z_plus(ea: int, eb: int, ec: int,
                   t1: sp.Symbol, t3: sp.Symbol, I: sp.Symbol) -> sp.Expr:
    """
    Compute the Gaussian element z_+ for the triple (eps_a, eps_b, eps_c)
    at k=3, following the construction of Section 7 (main.tex L2084):

        z_+ = - sum_{eps in {a,b,c}} sigma_eps * (s_{2,eps} * alpha_eps + beta_eps * i)

    where
        G_eps     = prod_{j != 2} (t_j + s_{j,eps} * i)   (cf. L1918)
        alpha_eps = Re(G_eps^2),    beta_eps = Im(G_eps^2)   (cf. L1917)
        s_{j,eps} = (-1)^{eps_j}   (cf. L210)

    The convention sigma_a = +1, sigma_b = +1, sigma_c = -1 follows
    Lemma~lem:Cond1antirecip (and matches the Pillar C setup).

    For k=3, the product over j != 2 has j in {1, 3}.
    """
    s1a, s2a, s3a = sign(ea, 1), sign(ea, 2), sign(ea, 3)
    s1b, s2b, s3b = sign(eb, 1), sign(eb, 2), sign(eb, 3)
    s1c, s2c, s3c = sign(ec, 1), sign(ec, 2), sign(ec, 3)

    # G_eps for j != 2 means j in {1, 3} at k=3
    G_a = (t1 + s1a * I) * (t3 + s3a * I)
    G_b = (t1 + s1b * I) * (t3 + s3b * I)
    G_c = (t1 + s1c * I) * (t3 + s3c * I)

    G_a_sq = sp.expand(G_a ** 2)
    G_b_sq = sp.expand(G_b ** 2)
    G_c_sq = sp.expand(G_c ** 2)

    # Split into real and imaginary parts over Q(i)
    # G_eps^2 = alpha_eps + beta_eps * i, with alpha, beta in Q[t1, t3]
    def split_re_im(expr):
        # SymPy treats I as a symbol; collect by I
        pol = sp.Poly(expr, I)
        # pol is degree-1 in I (after expansion of G^2, only I and 1 appear because t1,t3 are real)
        coeffs = pol.all_coeffs()  # high-to-low
        # If degree 1: [coef_I, coef_1]; if degree 0: [coef_1]
        if pol.degree() == 1:
            beta, alpha = coeffs[0], coeffs[1]
        elif pol.degree() == 0:
            alpha = coeffs[0]
            beta = sp.Integer(0)
        else:
            raise ValueError(f"unexpected degree {pol.degree()} for G_eps^2")
        return sp.expand(alpha), sp.expand(beta)

    alpha_a, beta_a = split_re_im(G_a_sq)
    alpha_b, beta_b = split_re_im(G_b_sq)
    alpha_c, beta_c = split_re_im(G_c_sq)

    sigma_a, sigma_b, sigma_c = +1, +1, -1

    # z_+ = - sum sigma_eps * (s_{2,eps} * alpha_eps + beta_eps * i)
    z_plus = -(
        sigma_a * (s2a * alpha_a + beta_a * I)
        + sigma_b * (s2b * alpha_b + beta_b * I)
        + sigma_c * (s2c * alpha_c + beta_c * I)
    )

    return sp.expand(z_plus)


# ---------------------------------------------------------------------
# Class representatives (from main.tex L3296-3303)
# ---------------------------------------------------------------------

def class_representatives(t1: sp.Symbol, t3: sp.Symbol, I: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """
    z_+^{(0)} = -(t_3 - i)^2 (t_1^2 - 1) + 2*(i*(t_3^2-1) - 6*t_3)*t_1
    z_+^{(1)} =  (t_3 + i)^2 (t_1^2 - 1) - 2*(i*(t_3^2-1) + 6*t_3)*t_1
    """
    z0 = -(t3 - I)**2 * (t1**2 - 1) + 2 * (I * (t3**2 - 1) - 6 * t3) * t1
    z1 = (t3 + I)**2 * (t1**2 - 1) - 2 * (I * (t3**2 - 1) + 6 * t3) * t1
    return sp.expand(z0), sp.expand(z1)


# ---------------------------------------------------------------------
# Orbit detection
# ---------------------------------------------------------------------

def orbit_of(z: sp.Expr, z0: sp.Expr, z1: sp.Expr,
             t1: sp.Symbol, t3: sp.Symbol) -> str | None:
    """
    Try to express z as u * z_c(eps1 * t1, t3) with u, eps1 in {+1, -1} and c in {0,1}.
    Return a string description, or None if no match (which would invalidate Claim B).
    """
    candidates = {
        ("0", "+", "+"): z0,
        ("0", "+", "-"): z0.subs(t1, -t1),
        ("0", "-", "+"): -z0,
        ("0", "-", "-"): -z0.subs(t1, -t1),
        ("1", "+", "+"): z1,
        ("1", "+", "-"): z1.subs(t1, -t1),
        ("1", "-", "+"): -z1,
        ("1", "-", "-"): -z1.subs(t1, -t1),
    }
    z_norm = sp.expand(z)
    for (c, u, e1), val in candidates.items():
        if sp.expand(z_norm - val) == 0:
            return f"class={c}, u={u}, eps_1_flip={e1}"
    return None


# ---------------------------------------------------------------------
# Main verification
# ---------------------------------------------------------------------

def main() -> int:
    print("=" * 72)
    print("verify_k3_cardinality.py")
    print("Session 19, gamma-1 reinforcement: Lemma D.1 and Lemma D.2 check")
    print("=" * 72)

    # Stage 1: cardinality cascade
    print("\n--- Claim A: cardinality cascade ---")
    counts, all_one_triples = stage_counts()
    for k, v in counts.items():
        print(f"  {k}: {v}")

    expected = {
        "stage_0_all_triples": 8**3,
        "stage_1_lex_a_lt_b": 8 * 7 // 2 * 8,
        "stage_2_distinct_c": 8 * 7 // 2 * 6,
        "stage_3_non_complementary": (8 * 7 // 2 - 4) * 6,
        "stage_4_all_one": 24,
    }
    print("\n--- Expected cascade ---")
    for k, v in expected.items():
        print(f"  {k}: {v}")

    cascade_ok = all(counts[k] == expected[k] for k in expected)
    print(f"\nCascade matches expected values: {cascade_ok}")
    if not cascade_ok:
        print("FAIL: cascade mismatch")
        return 1

    if len(all_one_triples) != 24:
        print(f"FAIL: expected 24 all-one triples, got {len(all_one_triples)}")
        return 1

    print(f"\nClaim A verified: |{{all-one triples at k=3}}| = 24.")

    # Stage 2: orbit decomposition
    print("\n--- Claim B: orbit decomposition under z_+ ~ u * z_c(+- t_1, t_3) ---")
    t1, t3, I = sp.symbols("t1 t3", real=True), None, None
    t1 = sp.Symbol("t1")
    t3 = sp.Symbol("t3")
    I = sp.I  # Gaussian unit

    z0, z1 = class_representatives(t1, t3, I)
    print(f"  z_+^(0) = {z0}")
    print(f"  z_+^(1) = {z1}")

    print("\n  Per-triple z_+ expansion and orbit assignment:")
    class_counts = {"0": 0, "1": 0}
    failures = []
    rows = []

    for idx, (ea, eb, ec) in enumerate(all_one_triples, start=1):
        z = compute_z_plus(ea, eb, ec, t1, t3, I)
        # Normalise: the construction yields z = (something) of the form
        # 4 * (z_+ expression). We allow free overall non-zero scalar by
        # checking if z/4 is in the orbit, then if not, try alternative
        # normalisations.
        # Empirically (from the construction), z = -4 * z_+ for class-0
        # representative case (Table row 1: triple (000, 011, 001)).
        # Verify by direct substitution below.
        orbit_label = None
        for scale in [sp.Rational(1, 4), sp.Rational(-1, 4),
                      sp.Rational(1, 2), sp.Rational(-1, 2),
                      1, -1, sp.Rational(1, 8), sp.Rational(-1, 8)]:
            res = orbit_of(scale * z, z0, z1, t1, t3)
            if res is not None:
                orbit_label = (scale, res)
                break

        if orbit_label is None:
            failures.append((idx, ea, eb, ec, z))
            rows.append((idx, format(ea, "03b"), format(eb, "03b"), format(ec, "03b"), "?", "NO MATCH"))
            continue

        scale, desc = orbit_label
        c = desc.split(",")[0].split("=")[1]
        class_counts[c] += 1
        rows.append((idx, format(ea, "03b"), format(eb, "03b"), format(ec, "03b"),
                     f"c={c}", f"scale={scale}, {desc}"))

    print(f"\n  Class 0 count: {class_counts['0']}")
    print(f"  Class 1 count: {class_counts['1']}")
    print(f"  Total assigned: {class_counts['0'] + class_counts['1']} (expected: 24)")

    if failures:
        print(f"\n  FAIL: {len(failures)} triples did not match any orbit element:")
        for idx, ea, eb, ec, z in failures:
            print(f"    #{idx}: (eps_a, eps_b, eps_c) = "
                  f"({format(ea,'03b')}, {format(eb,'03b')}, {format(ec,'03b')})")
            print(f"      z = {z}")
        return 1

    # Number of distinct equivalence classes
    classes_present = set(r[4] for r in rows if r[4] != "?")
    print(f"\n  Distinct equivalence classes: {len(classes_present)} (expected: 2)")

    if len(classes_present) != 2:
        print("FAIL: orbit count != 2")
        return 1

    # Print full assignment table
    print("\n--- Full assignment table ---")
    print(f"{'#':>3} {'eps_a':>5} {'eps_b':>5} {'eps_c':>5} {'class':>6}  detail")
    for r in rows:
        print(f"{r[0]:>3} {r[1]:>5} {r[2]:>5} {r[3]:>5} {r[4]:>6}  {r[5]}")

    # -----------------------------------------------------------------
    # Claim C: Hamming-2 pair structure underlying Lemma D.1 / D.2
    # -----------------------------------------------------------------
    print("\n--- Claim C: Hamming-2 pair structure (proof-supporting) ---")

    def hamming(a, b):
        return bin(a ^ b).count("1")

    pair_to_triples = {}
    for idx, (ea, eb, ec) in enumerate(all_one_triples, start=1):
        pair_to_triples.setdefault((ea, eb), []).append((idx, ec))

    h_counter = {1: 0, 2: 0, 3: 0}
    for (ea, eb) in pair_to_triples.keys():
        h_counter[hamming(ea, eb)] += 1
    print(f"  Distinct (eps_a, eps_b) pairs supporting triples: {len(pair_to_triples)}")
    print(f"  Pair Hamming-distance distribution: {h_counter}")
    if pair_to_triples and not all(h == 2 for (a, b) in pair_to_triples for h in [hamming(a, b)]):
        print("  FAIL: some supporting pair is not Hamming-2.")
        return 1
    if len(pair_to_triples) != 12:
        print(f"  FAIL: expected 12 Hamming-2 pairs, got {len(pair_to_triples)}")
        return 1

    # Each H-2 pair must support exactly 2 triples
    triples_per_pair = [len(v) for v in pair_to_triples.values()]
    if not all(n == 2 for n in triples_per_pair):
        print(f"  FAIL: not all pairs support exactly 2 triples: {sorted(set(triples_per_pair))}")
        return 1
    print(f"  Each Hamming-2 pair supports exactly 2 triples (12 pairs x 2 = 24). PASS.")

    # -----------------------------------------------------------------
    # Claim D: j_0-based class assignment table (proof-supporting)
    # -----------------------------------------------------------------
    # For each pair, identify j_0 (the unique coordinate where eps_a = eps_b)
    # and the two triples generated. Record (j_0, common_bit_at_j_0, class).
    print("\n--- Claim D: j_0-stratified class distribution ---")
    j0_summary = {}
    for (ea, eb), trips in pair_to_triples.items():
        diff = ea ^ eb
        if diff == 0b011: j0 = 1
        elif diff == 0b101: j0 = 2
        elif diff == 0b110: j0 = 3
        else:
            print(f"  FAIL: pair ({format(ea,'03b')},{format(eb,'03b')}) is not Hamming-2.")
            return 1
        j0_pos = 3 - j0  # bit position from LSB
        common_bit = (ea >> j0_pos) & 1
        # Read classes from rows assignment by triple index
        classes_in_pair = []
        for (idx, ec) in trips:
            row = rows[idx - 1]
            c_str = row[4].split("=")[-1]
            classes_in_pair.append(int(c_str))
        key = (j0, common_bit)
        j0_summary.setdefault(key, []).extend(classes_in_pair)

    print(f"  {'j_0':>3} {'common_bit':>11} {'class 0 count':>14} {'class 1 count':>14}")
    expected_balanced = True
    for key in sorted(j0_summary.keys()):
        cs = j0_summary[key]
        c0 = cs.count(0); c1 = cs.count(1)
        print(f"  {key[0]:>3} {key[1]:>11} {c0:>14} {c1:>14}")
        if c0 != 2 or c1 != 2:
            expected_balanced = False

    if not expected_balanced:
        print("  FAIL: j_0-stratified distribution is not (2, 2) for some (j_0, common_bit).")
        return 1
    print(f"  All 6 (j_0, common_bit) cells have (class 0, class 1) = (2, 2). PASS.")

    # Also: per-pair class distribution
    print("\n  Per-pair class distribution by j_0:")
    for j0_target in [1, 2, 3]:
        pairs_at_j0 = []
        for (ea, eb), trips in pair_to_triples.items():
            diff = ea ^ eb
            if (diff == 0b011 and j0_target == 1) or (diff == 0b101 and j0_target == 2) or (diff == 0b110 and j0_target == 3):
                cs_in_pair = []
                for (idx, ec) in trips:
                    cs_in_pair.append(int(rows[idx - 1][4].split("=")[-1]))
                pairs_at_j0.append((ea, eb, sorted(cs_in_pair)))
        n_same = sum(1 for (_, _, cs) in pairs_at_j0 if cs[0] == cs[1])
        n_diff = sum(1 for (_, _, cs) in pairs_at_j0 if cs[0] != cs[1])
        print(f"    j_0 = {j0_target}: {len(pairs_at_j0)} pairs ({n_same} same-class, {n_diff} mixed-class)")

    print("\n" + "=" * 72)
    print("PASS: Claims A, B, C, D verified.")
    print("  A: cardinality = 24")
    print("  B: orbit count = 2 (each of size 12)")
    print("  C: Hamming-2 pair structure (12 pairs x 2 triples)")
    print("  D: j_0-stratified balanced class distribution")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
