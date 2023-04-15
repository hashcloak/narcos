from sage.rings.integer_ring import Z, ZZ
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF

"""
This is a straightforward modification of the code in cado-nfs
See https://gforge.inria.fr/projects/cado-nfs/
and cado-nfs/polyselect/{alpha.sage,auxiliary.c}
See cado-nfs/sieve/makefb.sage for computing explicitly the roots mod p^k
"""

def number_of_roots(f,Fpx):
    """
    cado-nfs function, (c) cado-nfs/polyselect/alpha.sage
    Counts the roots of f mod p, without multiplicities. Projective roots
    are also counted (without multiplicities).
    """
    #fp= GF(p)['x'](f)
    fp = Fpx(f)
    if fp != 0:
        s=len([r[1] for r in fp.roots()])
        if (f.degree()>fp.degree()):
            s+=1
    else:
        # f reducing to zero mod p is a degenerate case.
        # we should never enter this branch because f = f0/cont(f0)
        print("Warning, counting roots of zero polynomial\n")
        s=f.degree()
    return s

def no_roots_mod_pk_rec(f, p, k, Fpx):
    v = f.content().valuation(p)
    fv = f // p**v
    if v >= k:
        return  p**k
    elif k == 1:
        return len(Fpx(fv).roots())
    else:
        n_pk = 0
        dfv = fv.derivative()
        dfvp = Fpx(dfv)
        for (r,e) in Fpx(fv).roots():
            if dfvp(r) != 0:
                n_pk += p**v
            else:
                s = ZZ(r)
                f2 = fv(s+p*f.parent().gen(0)) // p
                n_pk += p**v*no_roots_mod_pk_rec(f2, p, k-v-1, Fpx)
        return n_pk

def no_roots_mod_pk(f, p, k):
    """ Compute the number of zeros of f modulo p^k recursively from the roots modulo p
    ("roots" but f can have many more "roots" than deg(f) if p divides disc(f))
    This is just another way to rewrite alpha for NFS
    """
    disc_f = f.discriminant()
    Fpx = GF(p)['x']
    if ((disc_f % p) != 0):
        return number_of_roots(f,Fpx)
    else:
        n_pk_aff = no_roots_mod_pk_rec(f,p,k,Fpx)
        if ((f.leading_coefficient() % p) == 0):
            n_pk_pro = 1
            if k >= 2:
                n_pk_pro = no_roots_mod_pk_rec((f.reverse()(p*f.parent().gen(0)))//p,p,k-1,Fpx)
        else:
            n_pk_pro = 0
        return n_pk_aff + n_pk_pro

