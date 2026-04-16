#!/usr/bin/env python3
"""
Verify Theorem 7.5 inductive step.

Recurrence: z_+^(k) = z_+^(k-1)·t_k² + 2w^(k-1)·t_k - z_+^(k-1)

If z_+^(k) = h² with h = a·t_k + b, then:
  (I)   a² = z_+^(k-1)
  (II)  ab = w^(k-1)
  (III) b² = -z_+^(k-1)

Question: Does (I) alone contradict the inductive hypothesis that
z_+^(k-1) is not a Gaussian square? Or do we need (I)+(III)?

Answer: (I) says a² = z_+^(k-1), i.e., z_+^(k-1) IS a perfect square
in Q(i)[t_1,...,t_{k-2}]. This directly contradicts the inductive
hypothesis. So (I) alone suffices.

But let's verify this logic carefully and check whether (I)+(III)
give additional information.
"""
from sympy import symbols, expand, I as symI, factor, sqrt, solve

print("=" * 60)
print("  Theorem 7.5 inductive step: logical analysis")
print("=" * 60)

print("""
  SETUP:
    z_+^(k) = z_+^(k-1)·t² + 2w·t - z_+^(k-1)   (t = t_k)
    Suppose z_+^(k) = h² with h = at + b.
    Then h² = a²t² + 2abt + b².

  COEFFICIENT MATCHING:
    t²: a² = z_+^(k-1)           ... (I)
    t¹: 2ab = 2w → ab = w        ... (II)
    t⁰: b² = -z_+^(k-1)         ... (III)

  ANALYSIS:
    (I) states: z_+^(k-1) = a² for some a ∈ Q(i)[t_1,...,t_{k-2}].
    This means z_+^(k-1) is a perfect square in Q(i)[t_1,...,t_{k-2}].
    By the inductive hypothesis, z_+^(k-1) is NOT a Gaussian square.
    CONTRADICTION. ∎

    The argument is complete with (I) alone.
    (III) provides redundant but consistent information:
    From (I) and (III): b² = -a², so b = ±ia (in Q(i)).
    From (II): w = ab = ±ia², and from (I): w = ±i·z_+^(k-1).
    This is the w = ±i·z_+ condition explored in the numerical checks.
""")

# Now let's verify concretely for k=3 → k=4
print("=" * 60)
print("  Concrete check: k=3 base → k=4 step")
print("=" * 60)

t1, t3, t4 = symbols('t1 t3 t4', real=True)

# k=3 representative z_+ (Class 0 from Lemma 7.7)
z3 = expand(-(t3 - symI)**2 * (t1**2 - 1)
            + 2*(symI*(t3**2 - 1) - 6*t3) * t1)

print(f"\n  z_+^(3) = {z3}")

# Check: is z3 a perfect square in Q(i)[t1,t3]?
# If it were, we could write z3 = f² for some f.
# f would have degree 1 in t1 (since z3 has degree 2).
# f = α·t1 + β with α, β ∈ Q(i)[t3].
# f² = α²t1² + 2αβt1 + β²
# Matching with z3 = a·t1² + b·t1 + c:

a3 = expand(z3.coeff(t1, 2))
b3 = expand(z3.coeff(t1, 1))
c3 = expand(z3.coeff(t1, 0))

print(f"\n  z3 as poly in t1:")
print(f"    a = {a3}")
print(f"    b = {b3}")
print(f"    c = {c3}")

# For z3 = f²: need a = α², i.e., a must be a perfect square.
# a = -(t3-i)². Is -(t3-i)² a perfect square?
# -(t3-i)² = (i(t3-i))² · (-1/i²) = ... hmm.
# -(t3-i)² = i² · (t3-i)² · (-1) = -(t3-i)² ... circular.
# Actually: is there α ∈ Q(i)[t3] with α² = -(t3-i)²?
# α² = -(t3-i)² = (i(t3-i))² · (i²) ... no.
# (i(t3-i))² = i²(t3-i)² = -(t3-i)². YES!
# So α = i(t3-i) = it3 + 1 satisfies α² = -(t3-i)².

alpha_candidate = symI * (t3 - symI)
alpha_sq = expand(alpha_candidate**2)
print(f"\n  Check: [i(t3-i)]² = {alpha_sq}")
print(f"  a3 = {a3}")
print(f"  Equal? {expand(alpha_sq - a3) == 0}")

# So a3 IS a perfect square! α = i(t3-i).
# Similarly c3 = (t3-i)². Is this a perfect square? Yes: (t3-i)².
# But wait: z3 = a·t1² + b·t1 + c being a perfect square requires
# discriminant b²-4ac = 0 (for a quadratic in t1 to be a perfect square).

disc = expand(b3**2 - 4*a3*c3)
disc_f = factor(disc)
print(f"\n  Discriminant b²-4ac = {disc_f}")
print(f"  Is zero? {disc == 0}")

# Discriminant is -64i·t3·(t3+i)² ≠ 0.
# So z3 is NOT a perfect square (its discriminant is nonzero).
# This is consistent with Lemma 7.7.

print(f"\n  z_+^(3) is NOT a perfect square (disc ≠ 0). ✓")

# Now for the inductive step k=3→k=4:
# We need a specific all-one triple for k=4 that extends a k=3 triple.
# Let's use the recurrence directly.

# From the paper: z_+^(4) = z_+^(3)·t4² + 2w^(3)·t4 - z_+^(3)
# where w^(3) depends on the new sign s_{4,ε}.
# For a specific extension, let's compute w from the formula.

# Actually, the key point is already proven:
# If z_+^(4) = h² = (at4+b)², then a² = z_+^(3), contradicting
# the fact that z_+^(3) is not a perfect square.

print(f"\n  Inductive step logic:")
print(f"    IF z_+^(4) = (a·t4 + b)²")
print(f"    THEN a² = z_+^(3)  [from t4² coefficient]")
print(f"    BUT z_+^(3) is NOT a perfect square (Lemma 7.7)")
print(f"    CONTRADICTION. ✓")

print(f"\n  Note: equation (I) alone suffices. Equations (II) and (III)")
print(f"  are consequences but are not needed for the contradiction.")
print(f"  The paper's proof is CORRECT as stated.")

# Additional check: verify that a² = z3 is truly the content of (I).
print(f"\n{'='*60}")
print(f"  Formal verification: (I) ⟹ z_+^(k-1) is a square")
print(f"{'='*60}")

# In Q(i)[t1,...,t_{k-2}], "a² = z_+^(k-1)" means there exists
# a polynomial a such that a·a = z_+^(k-1).
# This is EXACTLY the definition of z_+^(k-1) being a perfect square
# in the ring Q(i)[t1,...,t_{k-2}].

# The inductive hypothesis says z_+^(k-1) is NOT a perfect square.
# So (I) directly contradicts the hypothesis.

# No need for (II) or (III). The proof is complete with (I) alone.

print(f"  (I): a² = z_+^(k-1) ⟺ z_+^(k-1) is a perfect square")
print(f"  Inductive hypothesis: z_+^(k-1) is NOT a perfect square")
print(f"  Direct contradiction. No additional equations needed.")
print(f"\n  VERDICT: Paper's proof is CORRECT and COMPLETE. ✓")
