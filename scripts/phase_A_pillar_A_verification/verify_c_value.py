"""
Verify the constant c = 2^(4m-1) * (-i) in Im(\bar{π}^(4m)) ≡ u^(4m) * c (mod π).

We claimed: Im(\bar{π}^(4m)) ≡ (2u)^(4m) / (2i) (mod π)
         = 2^(4m) u^(4m) / (2i)
         = 2^(4m - 1) u^(4m) / i
         = 2^(4m - 1) u^(4m) · (-i)        [since 1/i = -i]
         = u^(4m) · (2^(4m - 1) · (-i))
         = u^(4m) · c
where c = 2^(4m - 1) · (-i).

Numerical verification: compute LHS = Im(\bar{π}^(4m)) and RHS = u^(4m) · c, 
reduce both mod π, check they're equal in Z[i]/(π) ≅ F_p.

In Z[i]/(π) where π = u + vi, we have i ≡ -u/v (mod π) [from u + vi ≡ 0, so vi ≡ -u, so i ≡ -u/v].
But we work in F_p instead by reducing everything mod p and using the projection Z[i]/(π) → F_p
that sends a + b i to a + b · (i mod π) where i mod π is some specific element.

Actually simpler: any α ∈ Z[i] reduces mod π to a unique element of {0, 1, ..., p-1}.
Compute α mod π via norm trick: α mod π = (α · \bar{π}) mod p ... hmm, not quite.

Cleanest approach: 
  Z[i]/(π) ≅ F_p via the homomorphism that sends i to v_inv * (-u) mod p
  (since i ≡ -u/v ≡ -u · v^{-1} mod π).
  
Verify by computing both sides modulo p via this isomorphism, and checking equality.
"""
import sympy as sp
from sympy import I as iSym

def reduce_mod_pi(alpha_int_re, alpha_int_im, u, v, p):
    """
    α = alpha_int_re + alpha_int_im * i ∈ Z[i].
    Map to F_p via i ↦ -u * v^{-1} mod p.
    """
    v_inv = pow(v, p - 2, p)  # v^{-1} mod p
    i_image = (-u * v_inv) % p
    return (alpha_int_re + alpha_int_im * i_image) % p

# Test with several primes and m values
for p_val, u_val, v_val in [(5,2,1), (13,3,2), (17,4,1), (29,5,2), (37,6,1)]:
    print(f"--- p={p_val}, (u,v)=({u_val},{v_val}) ---")
    pi_bar = u_val - v_val * iSym
    for m in range(1, 6):
        # LHS: Im(\bar{π}^(4m))
        pi_bar_4m = sp.expand(pi_bar ** (4*m))
        Im_LHS = int(sp.im(pi_bar_4m))
        # Reduce to F_p
        LHS_mod_p = Im_LHS % p_val
        
        # RHS: u^(4m) · (2^(4m-1) · (-i)), reduced mod π
        u_4m = pow(u_val, 4*m, p_val)
        coef = pow(2, 4*m - 1, p_val)
        # RHS = u^(4m) · 2^(4m-1) · (-i), as element of Z[i]:
        # re part = 0, im part = -u^(4m) · 2^(4m-1)
        # But wait, we have it as: RHS_in_Zi = u^(4m) · 2^(4m-1) · (-i) = 0 + (-u^(4m) · 2^(4m-1)) i
        rhs_re = 0
        rhs_im = (-u_val**(4*m) * pow(2, 4*m - 1)) 
        RHS_mod_p = reduce_mod_pi(rhs_re, rhs_im, u_val, v_val, p_val)
        
        ok = (LHS_mod_p == RHS_mod_p)
        print(f"  m={m}: LHS mod p = {LHS_mod_p}, RHS mod p = {RHS_mod_p}, {'OK' if ok else 'FAIL'}")
