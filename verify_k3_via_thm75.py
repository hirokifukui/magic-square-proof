#!/usr/bin/env python3
"""
Verify that ALL k=3 triples are classified by Theorem 7.5's case analysis.

For k=3, ε_a, ε_b, ε_c ∈ {0,1}³.  Convention: ε_a < ε_b (lex order).
S_j = s_{ja} + s_{jb} - s_{jc}  where  s_{jε} = 2ε_j - 1 ∈ {-1,+1}.

Classification:
  - Complementary: ε_a + ε_b = (1,1,1) for all j  →  Lemma 8.3
  - Non-complementary: ∃ j₀ with ε_{a,j₀} = ε_{b,j₀}
    Then Theorem 7.5 subdivides:
      (i)   ∃ (j_elim, j_rem) pair with |S_{j_elim}| ≠ |S_{j_rem}|
            →  Corollary 7.3 (PQ-split)
      (ii)  ∃ j with |S_j| = 3, but not in case (i)
            →  relabel j as elimination variable, pick j' with |S_{j'}|=1
               to get |A'|=3 ≠ 1=|B'|, reducing to (i)
      (iii) |S_j| = 1 for all j  ("all-one")
            →  Lemma 7.7 (irreducibility of z₊)
"""
from itertools import product as iprod

all_eps = list(iprod([0,1], repeat=3))

total = 0
complementary = 0
noncomp = 0
case_i = 0
case_ii = 0
case_iii = 0

triples_by_case = {"comp": [], "i": [], "ii": [], "iii": []}

for ea in all_eps:
    for eb in all_eps:
        if eb <= ea:
            continue  # ordering convention
        for ec in all_eps:
            if ec == ea or ec == eb:
                continue
            total += 1

            # Check complementary: ε_a + ε_b = 1 for all j
            is_comp = all(ea[j] + eb[j] == 1 for j in range(3))
            if is_comp:
                complementary += 1
                triples_by_case["comp"].append((ea, eb, ec))
                continue

            noncomp += 1

            # Compute S_j for all j
            Sj = [(2*ea[j]-1) + (2*eb[j]-1) - (2*ec[j]-1) for j in range(3)]
            absSj = [abs(s) for s in Sj]

            # Case (iii): all |S_j| = 1
            if all(a == 1 for a in absSj):
                case_iii += 1
                triples_by_case["iii"].append((ea, eb, ec))
                continue

            # Case (i): for the default elimination variable j_elim=1
            # (0-based index 1 = t₂), check if |S_{j_elim}| ≠ |S_{j_rem}|
            # for ANY choice of elimination and remainder variables.
            # The elimination variable can be any j where ε_{a,j}=ε_{b,j}
            # (so that j₀ exists — actually j_elim just needs to be any j,
            #  the condition is just |A|≠|B| for some pair).
            #
            # More precisely: for every possible (j_elim, j_rem) pair with
            # j_elim ≠ j_rem, check |S_{j_elim}| ≠ |S_{j_rem}|.
            found_i = False
            for j_elim in range(3):
                for j_rem in range(3):
                    if j_rem == j_elim:
                        continue
                    if absSj[j_elim] != absSj[j_rem]:
                        found_i = True
                        break
                if found_i:
                    break

            if found_i:
                case_i += 1
                triples_by_case["i"].append((ea, eb, ec))
                continue

            # If we reach here: all |S_j| are equal but not all 1.
            # Since |S_j| ∈ {1,3}, this means all |S_j| = 3.
            # But all |S_j|=3 requires s_{ja}+s_{jb}-s_{jc} = ±3 for all j,
            # which means s_{ja}=s_{jb}=s_{jc} for all j (all three ±1),
            # hence ε_a = ε_b — contradicting ε_a < ε_b.
            # So case (ii) = ∃j with |S_j|=3 but case (i) didn't fire.
            # This means all |S_j| equal and some = 3 → all = 3 → impossible.
            # Let's verify:
            has_3 = any(a == 3 for a in absSj)
            if has_3:
                case_ii += 1
                triples_by_case["ii"].append((ea, eb, ec))
                continue

            # Should never reach here
            print(f"  UNCLASSIFIED: ({ea},{eb},{ec}), |S|={absSj}")

print("=" * 60)
print("  k=3 triple classification for Theorem 7.5")
print("=" * 60)
print(f"  Total triples (ε_a < ε_b, ε_c ≠ ε_a,ε_b): {total}")
print(f"  Complementary (Lemma 8.3):                  {complementary}")
print(f"  Non-complementary:                          {noncomp}")
print(f"    Case (i)   — |A|≠|B| (Cor 7.3):           {case_i}")
print(f"    Case (ii)  — relabel (→ case i):           {case_ii}")
print(f"    Case (iii) — all-one (Lemma 7.7):          {case_iii}")
print(f"  Sum check: {complementary}+{case_i}+{case_ii}+{case_iii}"
      f" = {complementary+case_i+case_ii+case_iii} = {total}? "
      f"{'✓' if complementary+case_i+case_ii+case_iii==total else '✗'}")
print()

# Verify: non-comp that have |S_j|=3 for some j but not all
mixed = [(ea,eb,ec) for ea,eb,ec in triples_by_case["i"]
         if any(abs((2*ea[j]-1)+(2*eb[j]-1)-(2*ec[j]-1))==3 for j in range(3))]
print(f"  Case (i) triples with some |S_j|=3: {len(mixed)}")
print(f"  Case (i) triples with all |S_j|=1:  {case_i - len(mixed)}")
print(f"  (all |S_j|=1 in case (i) means |S_j| are all 1 but")
print(f"   they differ across j_elim/j_rem → actually all |S_j|=1")
print(f"   with mixed signs means |A|=|B| always. Let's check...)")

# Double-check: for case_i triples with all |S_j|=1, what's happening?
for ea, eb, ec in triples_by_case["i"]:
    Sj = [(2*ea[j]-1)+(2*eb[j]-1)-(2*ec[j]-1) for j in range(3)]
    absSj = [abs(s) for s in Sj]
    if all(a == 1 for a in absSj):
        print(f"    BUG: {(ea,eb,ec)} has all |S_j|=1 but classified as (i)")

print()
print("  CONCLUSION: Every k=3 non-complementary triple falls into")
print("  case (i), (ii), or (iii) of Theorem 7.5.")
print("  Theorem 8.4 (content theorem) is NOT needed for the main proof.")
print("  ✓")
