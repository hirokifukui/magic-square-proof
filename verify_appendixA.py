#!/usr/bin/env python3
"""
Verify Appendix A, Steps 1-2.

Parametrisation (eq. 35):
  [ e₀²+a,      e₀²+d,        e₀²-a-d      ]
  [ e₀²-2a-d,   e₀²,          e₀²+2a+d     ]
  [ e₀²+a+d,    e₀²-d,        e₀²-a        ]

Check: all 8 line-sums = 3e₀², and the 4 symmetric-pair identities hold.
"""
from sympy import symbols, expand

e, a, d = symbols('e a d')

M = [
    [e + a,       e + d,       e - a - d    ],
    [e - 2*a - d, e,           e + 2*a + d  ],
    [e + a + d,   e - d,       e - a        ],
]
S = 3 * e

print("=== 8 line-sums ===")
ok = True
for i in range(3):
    rs = expand(sum(M[i]))
    r = rs == S
    ok = ok and r
    print(f"  Row {i+1}: {rs} = 3e? {r}")
for j in range(3):
    cs = expand(sum(M[i][j] for i in range(3)))
    r = cs == S
    ok = ok and r
    print(f"  Col {j+1}: {cs} = 3e? {r}")
d1 = expand(M[0][0] + M[1][1] + M[2][2])
d2 = expand(M[0][2] + M[1][1] + M[2][0])
ok = ok and d1 == S and d2 == S
print(f"  Diag1: {d1} = 3e? {d1 == S}")
print(f"  Diag2: {d2} = 3e? {d2 == S}")
print(f"  All 8: {'✓' if ok else '✗'}")

print("\n=== 4 symmetric-pair identities ===")
entries = [M[i][j] for i in range(3) for j in range(3)]
pairs = [(0,8),(1,7),(2,6),(3,5)]
pok = True
for i, j in pairs:
    s = expand(entries[i] + entries[j])
    r = s == 2*e
    pok = pok and r
    print(f"  a_{i+1}+a_{j+1} = {s} = 2e? {r}")
print(f"  All 4: {'✓' if pok else '✗'}")

print(f"\n  VERDICT: {'PASS' if ok and pok else 'FAIL'}")
