"""
Correct rep counting using math.isqrt. Small N only.
"""
from math import isqrt

def count_reps_correct(N):
    count = 0
    reps = []
    V_max = isqrt(N // 2)
    for V in range(0, V_max + 1):
        U_sq = N - V*V
        U = isqrt(U_sq)
        if U*U == U_sq and U >= V:
            count += 1
            reps.append((U, V))
    return count, reps

print(f"{'p':>4} {'e':>3} {'N':>12} {'#reps':>8} {'e+1':>5}  {'match':>6}  {'reps'}")
print("-" * 80)
all_ok = True
for p in [5, 13, 17, 29]:
    for e in range(1, 4):
        N = 2 * p**(2*e)
        cnt, reps = count_reps_correct(N)
        ok = (cnt == e + 1)
        all_ok = all_ok and ok
        flag = 'OK' if ok else 'FAIL'
        print(f"  {p:>3} {e:>3} {N:>12} {cnt:>8} {e+1:>5}  {flag:>6}  {reps}")
print()
print(f"All cases match #reps = e+1: {all_ok}")
