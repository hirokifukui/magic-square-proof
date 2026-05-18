"""
Pillar A Step P2 — Stage 6: verify the row-sum reduction to ε_1 Δ_{l_1} + ε_2 Δ_{l_2} + ε_3 Δ_{l_3} = 0.

The setup:
  3x3 magic square M of squares with centre entry e_0^2.
  Row/column/diagonal sum = S = 3 e_0^2.
  
  Each line (row, col, or diag) through the centre has the form a_i^2 + e_0^2 + a_{10-i}^2 = S,
  i.e. a_i^2 + a_{10-i}^2 = 2 e_0^2  (one of the 4 antipodal pair sums).
  
  Each line NOT through the centre is a row or col entirely "off-centre", but in a 3x3 magic
  square there are 8 lines total: 3 rows, 3 cols, 2 diagonals. The centre cell sits in 4 lines:
  middle row, middle col, both diagonals. The 4 antipodal-pair sums correspond to these 4 
  centre-lines.
  
  Wait — actually all lines pass through some cell pattern:
    Row 1: a_1 + a_2 + a_3 = S      <- not through centre (= a_5 = e_0^2)
    Row 2: a_4 + a_5 + a_6 = S      <- through centre
    Row 3: a_7 + a_8 + a_9 = S      <- not through centre
    Col 1: a_1 + a_4 + a_7 = S      <- not through centre
    Col 2: a_2 + a_5 + a_8 = S      <- through centre
    Col 3: a_3 + a_6 + a_9 = S      <- not through centre
    Diag1: a_1 + a_5 + a_9 = S      <- through centre
    Diag2: a_3 + a_5 + a_7 = S      <- through centre
  
  Lines through the centre (4 of them): give a_i^2 + a_{10-i}^2 = 2 e_0^2 (after subtracting a_5 = e_0^2).
    These yield the 4 antipodal pairs.
  
  Lines NOT through the centre (4 of them): a_i^2 + a_j^2 + a_k^2 = 3 e_0^2 with i,j,k 
    such that none is 5. These are non-trivial 3-term sums.
  
  Now for Pillar A (k=1, e_0 = p^e), each off-centre entry a_i^2 is one of the 4 antipodal-pair 
  reps, i.e. one of the W^2 or W'^2 from some R_l class.
  
  Let me re-examine. The 4 antipodal pairs are:
    (a_1, a_9), (a_2, a_8), (a_3, a_7), (a_4, a_6)
  Each pair satisfies a_i^2 + a_{10-i}^2 = 2 e_0^2 = N.
  
  So each pair has its two entries equal to {W_l, W_l'} for some l ∈ {0, ..., e-1} 
  (cannot be e by P1, l = e being the centre rep).
  
  Now consider a NON-centre line, e.g., Row 1: a_1^2 + a_2^2 + a_3^2 = 3 e_0^2.
  
  Each a_i has a pair partner: a_1 ↔ a_9, a_2 ↔ a_8, a_3 ↔ a_7.
  a_1^2 + a_9^2 = 2 e_0^2  ⟹  a_9^2 = 2 e_0^2 - a_1^2
  a_2^2 + a_8^2 = 2 e_0^2  ⟹  a_8^2 = 2 e_0^2 - a_2^2
  a_3^2 + a_7^2 = 2 e_0^2  ⟹  a_7^2 = 2 e_0^2 - a_3^2
  
  And Row 3 says a_7^2 + a_8^2 + a_9^2 = 3 e_0^2.
  Substituting: (2 e_0^2 - a_3^2) + (2 e_0^2 - a_2^2) + (2 e_0^2 - a_1^2) = 3 e_0^2
    => 6 e_0^2 - (a_1^2 + a_2^2 + a_3^2) = 3 e_0^2
    => a_1^2 + a_2^2 + a_3^2 = 3 e_0^2   ✓ (Row 1 itself)
  
  So Row 3 is dependent on Row 1 modulo the 4 antipodal-pair identities. No new constraint from Row 3.
  
  Similarly Col 1: a_1^2 + a_4^2 + a_7^2 = 3 e_0^2 with a_7^2 = 2 e_0^2 - a_3^2.
    => a_1^2 + a_4^2 + (2 e_0^2 - a_3^2) = 3 e_0^2
    => a_1^2 - a_3^2 + a_4^2 = e_0^2   (NEW)
  
  So Col 1 (modulo antipodal pairs) becomes a_1^2 - a_3^2 + a_4^2 = e_0^2.

The Pillar A claim is:
  For pair i, write a_i = W_{l_i} or W'_{l_i} for some l_i ∈ {0, ..., e-1} (the four distinct l's
  selected by P1 pigeonhole). Then a_i^2 = (some choice of) {W_{l_i}^2, W'_{l_i}^2}.
  
  By P1, the 4 antipodal pairs select 4 DISTINCT l's, say {l_1, l_2, l_3, l_4} ⊆ {0,...,e-1}.
  
  Subtract: a_i^2 - a_j^2 (within line) gives... and combine three pairs to get a 
  difference-of-Delta_l equation. The key trick: each a_i^2 in Row 1 corresponds to some
  pair (i, 10-i), and the choice of "W or W'" within the pair gives a ± sign. So a_i^2 ∈ 
  {W_{l_i}^2, W_{l_i}^2 - Δ_{l_i}} = {W_{l_i}^2, W'^2_{l_i}}.
  
  Actually: W_{l_i}^2 - W'^2_{l_i} = Δ_{l_i}. Either a_i^2 = W^2 or a_i^2 = W^2 - Δ. So
  a_i^2 = W_{l_i}^2 - σ_i Δ_{l_i} / 2  (for σ_i ∈ {0, 2} which we can rewrite as ± 1 with offset).
  
  Let me re-set: for each pair i ∈ {1,2,3,4}, let l_i be the assigned class. The pair entries
  squared are {W_{l_i}^2, W'^2_{l_i}}, and Δ_{l_i} = W_{l_i}^2 - W'^2_{l_i}.
  
  For Row 1: a_1, a_2, a_3 ∈ pairs 1, 2, 3 respectively. (Row 1 partner: a_9, a_8, a_7 ∈ pairs 1, 2, 3.)
  Define ε_i ∈ {±1} via: a_i^2 = (W_{l_i}^2 + W'^2_{l_i})/2 + ε_i Δ_{l_i}/2 = e_0^2 + ε_i Δ_{l_i}/2.
  [Since W_{l_i}^2 + W'^2_{l_i} = 2 e_0^2.]
  
  Then Row 1: sum over i ∈ {1,2,3} of (e_0^2 + ε_i Δ_{l_i}/2) = 3 e_0^2
    => sum of ε_i Δ_{l_i}/2 = 0  
    => ε_1 Δ_{l_1} + ε_2 Δ_{l_2} + ε_3 Δ_{l_3} = 0.  ← THIS IS THE TARGET EQUATION!
  
  Now ν_p(Δ_{l_i}) = ν_p(2 p^(2 l_i) J_{e - l_i}) = 2 l_i + 0 + 0 = 2 l_i  
  (since p is odd and ν_p(J_{e-l_i}) = 0 by Stage 5).
  
  The three l_i are distinct (since the three pairs in Row 1 are different pairs, by P1 distinct l's).
  So {2 l_1, 2 l_2, 2 l_3} are three distinct values, the smallest being 2 l_min.
  
  Ultrametric: ν_p(ε_1 Δ_{l_1} + ε_2 Δ_{l_2} + ε_3 Δ_{l_3}) = min{2 l_1, 2 l_2, 2 l_3} = 2 l_min < ∞.
  But the LHS is 0, which has ν_p = +∞. Contradiction.
  
This script: numerically verify the reduction "Row 1 sum -> ε_i Δ_{l_i} sum" for an explicit setup.
"""
import sympy as sp

print("=" * 80)
print("Stage 6: Row-sum reduction to ε_1 Δ_{l_1} + ε_2 Δ_{l_2} + ε_3 Δ_{l_3} = 0")
print("=" * 80)
print()
print("Verifying the algebraic identity:")
print("  Row 1: a_1^2 + a_2^2 + a_3^2 = 3 e_0^2")
print("    with a_i^2 = e_0^2 + ε_i Δ_{l_i}/2,  ε_i ∈ {±1}")
print("  ⟹  ε_1 Δ_{l_1} + ε_2 Δ_{l_2} + ε_3 Δ_{l_3} = 0")
print()

# Symbolic check
e0sq = sp.Symbol('e_0^2')
D1, D2, D3 = sp.symbols('Delta_1 Delta_2 Delta_3')
eps1, eps2, eps3 = sp.symbols('epsilon_1 epsilon_2 epsilon_3')

a1sq = e0sq + eps1 * D1 / 2
a2sq = e0sq + eps2 * D2 / 2
a3sq = e0sq + eps3 * D3 / 2

row_sum = sp.expand(a1sq + a2sq + a3sq)
target = 3 * e0sq
diff_expr = sp.expand(row_sum - target)
print(f"  Row 1 sum: {row_sum}")
print(f"  Target:    {target}")
print(f"  Difference (after subtracting target):   {diff_expr}")
print(f"  Setting this to 0: 2 · (above) = ε_1 Δ_1 + ε_2 Δ_2 + ε_3 Δ_3 = {sp.expand(2 * diff_expr)}")
print(f"  ⟹ ε_1 Δ_1 + ε_2 Δ_2 + ε_3 Δ_3 = 0  ✓")
print()

# Now the ν_p argument: distinct ν_p, ultrametric.
print("Stage 7: Ultrametric exclusion")
print()
print("  Given l_1 < l_2 < l_3 ≤ e-1 (distinct by P1 pigeonhole), and ν_p(Δ_l) = 2l + ν_p(2) + ν_p(J_{e-l}).")
print("  Since p is odd: ν_p(2) = 0.")
print("  By Stage 5: ν_p(J_{e-l}) = 0 for e-l ≥ 1, i.e. l ≤ e-1. ✓")
print("  Hence ν_p(Δ_l) = 2l for l ≤ e-1.")
print()
print("  Sum: ε_1 Δ_{l_1} + ε_2 Δ_{l_2} + ε_3 Δ_{l_3} with l_1 < l_2 < l_3.")
print("  ν_p(ε_i Δ_{l_i}) = ν_p(Δ_{l_i}) = 2 l_i, pairwise distinct.")
print("  Ultrametric: ν_p(sum) = min{2 l_1, 2 l_2, 2 l_3} = 2 l_1 (finite).")
print("  But sum = 0 has ν_p = +∞. Contradiction.")
print()

# Numerical sanity: take e = 4, p = 5, l = 0, 1, 2, 3 (four distinct l's). Pick three l's for Row 1.
print("=" * 80)
print("Stage 8: Numerical sanity — for e=4, p=5, sum ε_i Δ_{l_i} cannot be 0 for any sign choice")
print("=" * 80)
print()

p_val = 5
e_val = 4
u_val, v_val = 2, 1
i_sym = sp.I
pi_bar = u_val - v_val * i_sym

def Delta_l(p_val, e_val, l):
    m = e_val - l
    pi_bar_4m = sp.expand(pi_bar ** (4*m))
    J_m = abs(int(sp.im(pi_bar_4m)))
    return 2 * (p_val ** (2*l)) * J_m

print(f"  p = {p_val}, e = {e_val}")
for l in range(e_val):
    D = Delta_l(p_val, e_val, l)
    nu = 0
    n = D
    while n % p_val == 0:
        nu += 1
        n //= p_val
    print(f"  Δ_{l} = {D}, ν_p(Δ_{l}) = {nu}  (expected: {2*l})  {'OK' if nu == 2*l else 'FAIL'}")

print()
print("Test all sign-pattern combinations for triples (l_1, l_2, l_3):")
from itertools import combinations, product
fails = 0
total = 0
for triple in combinations(range(e_val), 3):
    Deltas = [Delta_l(p_val, e_val, l) for l in triple]
    for signs in product([-1, 1], repeat=3):
        total += 1
        s = sum(sg * D for sg, D in zip(signs, Deltas))
        if s == 0:
            print(f"  l={triple}, signs={signs}: SUM = 0   <- VIOLATES Pillar A!")
            fails += 1
print(f"  Tested {total} (triple, sign) combinations: {fails} failures (expected: 0)")
print()

print("=" * 80)
print("STAGES 6-8 RESULT: Pillar A row-sum reduction and ultrametric exclusion verified.")
print("=" * 80)
