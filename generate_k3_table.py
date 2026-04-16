#!/usr/bin/env python3
"""
Generate LaTeX table of the 24 all-one triples for k=3 (Lemma 7.7).
Outputs k3_triples_table.tex.
"""
from itertools import product as iprod
from sympy import symbols, expand, I as symI

t1, t3 = symbols('t1 t3', real=True)

def gmul(a, b):
    return (expand(a[0]*b[0]-a[1]*b[1]), expand(a[0]*b[1]+a[1]*b[0]))

def compute_z(ea, eb, ec):
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
    from sympy import Poly
    p = Poly(expand(ra+rb-rc), t2)
    c12 = expand(p.nth(2))
    c11 = expand(p.nth(1))
    return expand(c11/4 + c12/2 * symI)

# Class representatives
z_rep0 = expand(-(t3-symI)**2*(t1**2-1) + 2*(symI*(t3**2-1)-6*t3)*t1)
z_rep1 = expand((t3+symI)**2*(t1**2-1) - 2*(symI*(t3**2-1)+6*t3)*t1)

units = [(1,'1'), (-1,'-1'), (symI,'i'), (-symI,'-i')]

all_eps = list(iprod([0,1], repeat=3))
rows = []

for ea in all_eps:
    for eb in all_eps:
        if eb <= ea:
            continue
        for ec in all_eps:
            if ec == ea or ec == eb:
                continue
            if all(ea[j]+eb[j]==1 for j in range(3)):
                continue
            if not all(abs((2*ea[j]-1)+(2*eb[j]-1)-(2*ec[j]-1))==1 for j in range(3)):
                continue
            z = compute_z(ea, eb, ec)
            # Find class and relation
            found = False
            for cls, z_rep in [(0, z_rep0), (1, z_rep1)]:
                for u_val, u_name in units:
                    for flip, flip_name in [(False, ''), (True, '-')]:
                        zr = z_rep if not flip else z_rep.subs(t1, -t1)
                        if expand(z - u_val*zr) == 0:
                            sign_str = flip_name + 't_1'
                            rows.append((ea, eb, ec, cls, u_name, sign_str))
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if not found:
                print(f"UNMATCHED: {ea},{eb},{ec}")

assert len(rows) == 24, f"Expected 24, got {len(rows)}"

# Format LaTeX
def fmt_eps(e):
    return ''.join(str(x) for x in e)

lines = []
lines.append(r"\begin{table}[ht]")
lines.append(r"\centering")
lines.append(r"\footnotesize")
lines.append(r"\begin{tabular}{r@{\;}l@{\;}l@{\;}l@{\quad}c@{\quad}l}")
lines.append(r"\hline")
lines.append(r"\# & $\eps_a$ & $\eps_b$ & $\eps_c$ & Class & Relation \\")
lines.append(r"\hline")

for i, (ea, eb, ec, cls, u_name, sign_str) in enumerate(rows):
    n = i + 1
    ea_s = fmt_eps(ea)
    eb_s = fmt_eps(eb)
    ec_s = fmt_eps(ec)

    if u_name == '1' and sign_str == 't_1':
        rel = f"$z_+^{{({cls})}}$"
    else:
        u_part = '' if u_name == '1' else u_name + r'\cdot '
        rel = f"${u_part}z_+^{{({cls})}}({sign_str}, t_3)$"

    lines.append(f"{n:2d} & ${ea_s}$ & ${eb_s}$ & ${ec_s}$ & {cls} & {rel} \\\\")

lines.append(r"\hline")
lines.append(r"\end{tabular}")
lines.append(r"\caption{The $24$ all-one triples for $k=3$ and their equivalence class.}")
lines.append(r"\label{tab:k3triples}")
lines.append(r"\end{table}")

tex = '\n'.join(lines)
with open('k3_triples_table.tex', 'w') as f:
    f.write(tex)

print(tex)
print(f"\n  Written to k3_triples_table.tex ({len(rows)} rows)")
