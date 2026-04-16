#!/usr/bin/env python3
"""
Verify Theorem 7.5 inductive step: domain of a, b in h = a·t_k + b.

The recurrence is:
  z_+^(k) = z_+^(k-1)·t_k² + 2w^(k-1)·t_k - z_+^(k-1)

with z_+^(k-1) and w^(k-1) in Q(i)[t_1,...,t_{k-2}].
Note: the SUBSCRIPT indices on t are 1, 3, 5, ... (odd, skipping t_2
which is the eliminated variable). For k=3, the active variables are
t_1, t_3. For k=4, the active variables are t_1, t_3, t_4 (t_4 is
the new variable added in the inductive step).

Question: if z_+^(k) = h² with h ∈ Q(i)[t_1,...,t_{k-1}], does h
necessarily have the form a·t_k + b with a, b independent of t_k?
And more importantly: can a, b depend on t_{k-1}?

Analysis:
  z_+^(k) is degree 2 in t_k, with coefficients in Q(i)[t_1,...,t_{k-2}].
  (Here t_{k-2} means the second-to-last active variable, not t_{k-2}
  literally.)

  The key: z_+^(k-1) and w^(k-1) do NOT involve t_k.
  They live in Q(i)[t_1, t_3, ..., t_{k-1}].

  Wait -- let me re-read the paper. The paper says:
  "w^(k-1) ∈ Q(i)[t_1,...,t_{k-2}]"
  and "is independent of t_k".

  The paper's variable convention: for k primes, the variables are
  t_1, t_2, ..., t_k, with t_2 eliminated. So the active variables
  for z_+^(k) are t_1, t_3, t_4, ..., t_k.

  The recurrence adds t_k as the new variable:
  z_+^(k) involves t_1, t_3, ..., t_{k-1}, t_k
  z_+^(k-1) involves t_1, t_3, ..., t_{k-1}  (all EXCEPT t_k)
  w^(k-1) involves t_1, t_3, ..., t_{k-1}     (same as z_+^(k-1))

  Actually, I need to be more careful. Let me re-read the paper:

  "w^(k-1) = ... ∈ Q(i)[t_1, ..., t_{k-2}]"

  The paper writes t_1, ..., t_{k-2} for the variables of w^(k-1).
  With k primes and t_2 eliminated, the variables for the (k-1)-prime
  problem are t_1, t_3, ..., t_{k-1}. That's k-2 variables.
  The paper writes these as t_1, ..., t_{k-2} (relabeled).

  So the paper's "t_1, ..., t_{k-2}" is a RELABELING of the actual
  variables t_1, t_3, ..., t_{k-1}. In particular, "t_{k-2}" in the
  paper's notation IS t_{k-1} in the original numbering.

  This means: a, b ∈ Q(i)[t_1,...,t_{k-2}] in the paper's notation
  is the same as a, b ∈ Q(i)[t_1, t_3, ..., t_{k-1}] in original
  numbering.

  The paper's claim is: h = a·t_k + b with a, b NOT involving t_k.
  This is correct because z_+^(k) is degree 2 in t_k, so h must
  be degree 1 in t_k, and the coefficients a, b can depend on
  ALL other active variables (t_1, t_3, ..., t_{k-1}).

  The NOTATION "a, b ∈ Q(i)[t_1,...,t_{k-2}]" is potentially
  misleading because it suggests a, b are in a ring with k-2
  variables, but the k-2 refers to the NUMBER of remaining
  variables after removing t_2 (eliminated) and t_k (the new one).
"""
from sympy import symbols, expand, I as symI, Poly

print("=" * 65)
print("  Theorem 7.5: variable domain analysis")
print("=" * 65)

# Let's work with the actual variables for k=4.
# Active variables: t_1, t_3, t_4 (t_2 is eliminated).
# z_+^(3) lives in Q(i)[t_1, t_3].
# z_+^(4) = z_+^(3)·t_4² + 2w·t_4 - z_+^(3) lives in Q(i)[t_1, t_3, t_4].
# w^(3) lives in Q(i)[t_1, t_3] (does NOT involve t_4).

t1, t3, t4 = symbols('t1 t3 t4', real=True)

# z_+^(3) from Lemma 7.7, Class 0 representative
z3 = expand(-(t3 - symI)**2 * (t1**2 - 1)
            + 2*(symI*(t3**2 - 1) - 6*t3) * t1)

print(f"\n  z_+^(3) ∈ Q(i)[t1, t3]:")
print(f"    = {z3}")
print(f"    Variables: t1, t3 (does NOT involve t4)")

# If z_+^(4) = h² with h ∈ Q(i)[t1, t3, t4]:
# Since z_+^(4) has degree 2 in t4, h has degree 1 in t4.
# h = a·t4 + b where a, b ∈ Q(i)[t1, t3].
# NOTE: a, b CAN involve t3! The paper writes "a, b ∈ Q(i)[t1,...,t_{k-2}]"
# which for k=4 means a, b ∈ Q(i)[t1, t3] (two variables). This is CORRECT.

print(f"\n  If z_+^(4) = h², then h = a·t4 + b with:")
print(f"    a, b ∈ Q(i)[t1, t3]  (the variables of z_+^(3))")
print(f"    Coefficient of t4²: a² = z_+^(3)")
print(f"    → z_+^(3) must be a perfect square in Q(i)[t1, t3]")

# Now the reviewer's question: the paper writes t_1,...,t_{k-2}.
# For k=4: t_1, t_2 (= t_3 after relabeling). Two variables.
# For k=5: t_1, t_2, t_3 (= t_1, t_3, t_4 after relabeling). Three variables.
#
# The paper's "t_1,...,t_{k-2}" is a shorthand for "all active variables
# except t_k (the new one)". The CLAIM is correct.

print(f"\n  The paper's notation 't_1,...,t_{{k-2}}' is a shorthand")
print(f"  for the k-2 active variables OTHER than t_k.")
print(f"  For k=4: these are t_1, t_3 (two variables).")
print(f"  For general k: t_1, t_3, t_4, ..., t_{{k-1}} (k-2 variables).")

# However, the paper literally writes "a, b ∈ Q(i)[t_1, ..., t_{k-2}]".
# A reviewer reading this as the LITERAL variables t_1, ..., t_{k-2}
# (in original numbering) would get t_1, t_2, ..., t_{k-2}, which
# INCLUDES t_2 (the eliminated variable) -- WRONG.
#
# And it EXCLUDES t_{k-1} -- which is one of the active variables!
#
# So the notation IS misleading. Let me check what the paper actually says.

print(f"\n{'='*65}")
print(f"  Checking: does the paper's notation match?")
print(f"{'='*65}")

# The paper says h ∈ Q(i)[t_1,...,t_{k-1}] and a, b ∈ Q(i)[t_1,...,t_{k-2}].
# The subscript convention in §7.3: the variables are t_1, t_3, ..., t_k
# (with t_2 eliminated). But the paper's equations use the literal
# numbering t_1, t_3, t_4, ..., t_k.
#
# When the paper writes "h ∈ Q(i)[t_1,...,t_{k-1}]", it means
# all variables up to t_{k-1}, which includes t_1, t_3, ..., t_{k-1}.
# (t_2 is not among the generators because it was eliminated.)
#
# When the paper writes "a, b ∈ Q(i)[t_1,...,t_{k-2}]", this is
# ambiguous: does it mean t_1, t_2, ..., t_{k-2} (literal)?
# Or t_1, t_3, ..., t_{k-2} (skipping t_2)?
#
# In the context of §7.3, t_2 is eliminated, so Q(i)[t_1,...,t_{k-2}]
# should be read as "the polynomial ring in t_1, t_3, ..., t_{k-2}",
# which has k-3 generators.
#
# But the correct ring for a, b is Q(i)[t_1, t_3, ..., t_{k-1}],
# which has k-2 generators.
#
# The issue: t_{k-1} IS an active variable (it appears in z_+^(k-1)
# and w^(k-1)), and a, b should be allowed to depend on it.
# The paper writes "t_1,...,t_{k-2}" which EXCLUDES t_{k-1}.
#
# This is a BUG in the paper's notation!

print(f"  Paper says:  h = a·t_k + b with a, b ∈ Q(i)[t_1,...,t_{{k-2}}]")
print(f"  Correct:     h = a·t_k + b with a, b ∈ Q(i)[t_1, t_3,...,t_{{k-1}}]")
print(f"  The paper's 't_{{k-2}}' should be 't_{{k-1}}'.")
print(f"")
print(f"  This is a NOTATIONAL error, not a logical error.")
print(f"  The proof's conclusion (a² = z_+^(k-1) → contradiction)")
print(f"  is valid regardless, since z_+^(k-1) lives in the SAME ring")
print(f"  Q(i)[t_1, t_3,...,t_{{k-1}}] as a and b.")
print(f"")
print(f"  Similarly, the UFD parenthetical should say")
print(f"  'Q(i)[t_1,...,t_{{k-1}}] is a UFD' (the ring containing h),")
print(f"  not 'Q(i)[t_1,...,t_{{k-1}}]' which the paper already writes.")

# Wait, let me re-read the paper EXACTLY.
# The paper says:
# "h = a t_k + b with a, b ∈ Q(i)[t_1,...,t_{k-2}]
#  (the factorisation is unique since Q(i)[t_1,...,t_{k-1}] is a UFD)"
#
# So the UFD remark says Q(i)[t_1,...,t_{k-1}], but a, b are claimed
# to be in Q(i)[t_1,...,t_{k-2}]. The mismatch IS the bug.

print(f"\n  CONCLUSION: The paper should say")
print(f"  'a, b ∈ Q(i)[t_1,...,t_{{k-1}}]' (not t_{{k-2}}).")
print(f"  The proof logic is unaffected; this is a notation fix.")
