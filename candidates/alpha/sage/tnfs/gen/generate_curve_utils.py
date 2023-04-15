from sage.functions.log import log
from sage.rings.real_mpfr import RR, RealField
from sage.rings.integer_ring import ZZ
from sage.functions.other import floor, ceil
from sage.arith.functions import lcm
from sage.arith.misc import GCD, gcd
from sage.arith.misc import valuation
from sage.misc.misc_c import prod
from sage.rings.fast_arith import prime_range
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF

from tnfs.simul.enumerate_sparse_T import bit_positions_2naf, bit_positions

def find_u_mod_m(poly, m=None):
    """
    Congruence conditions u mod m such that poly(m*x+u) is an integer

    INPUT:
    - `poly`: a univariate polynomial with rational coefficients
    RETURN: m and a list of pairs (ui, ci) where ci is the content of poly(m*x+ui)
    if the poly already has integer coefficients, it returns 1, [(0, content(poly))]
    """
    if m is None:
        m = lcm([pj.denom() for pj in poly.list()])
    res = []
    _x_ = poly.variables()[0]
    for i0 in range(m):
        pi = poly(m*_x_ +i0) # evaluate at x*m + i0
        lcm_den = lcm([pj.denom() for pj in pi.list()])
        if lcm_den == 1:
            gcd_num = gcd([pj.numer() for pj in pi.list()])
            res.append((i0, gcd_num))
    return res

def small_factors(N, prod_small):
    """
    Compute the B-smooth factor of N with successive GCD

    INPUT:
    - `N`: (large) integer
    - `prod_small`: large integer, product o small primes up to some bound B (e.g. 10^6)

    RETURN: n, R = N//n where n is B-smooth and N//n has no prime factor smaller than B
   """
    pri = gcd(N, prod_small)
    pr = pri
    n = N // pri
    while pri > 1:
        pri = gcd(n, pri)
        pr *= pri
        n = n // pri
    return pr, n

def find_min(poly_p, length_p, m=1, umodm=0):
    """
    find the minimum seed m0 such that poly_p(m0) has length_p bits
    """
    RRa = RealField(ceil(length_p*log(2) + 2))
    roots_P = (poly_p - 2**(length_p-1)).roots(RRa)
    u_min = []
    for r in roots_P:
        u = ZZ(floor(r[0]))
        if m > 1: # works also for umodm=0
            u = u - (u % m) + umodm
        if u > 0:
            for c in [100,10,1]:
                while ZZ(poly_p(u)).nbits() < length_p:
                    u += m*c
            # now poly(u) has length_p bits or more, decrease u
            for c in [100,10,1]:
                while (u-m*c > 0) and ZZ(poly_p(u-m*c)).nbits() >= length_p:
                    u -= m*c
        else:
            for c in [100,10,1]:
                while ZZ(poly_p(u)).nbits() < length_p:
                    u -= m*c
            # now poly(u) has length_p bits or more, increase u (decrease abs(u))
            for c in [100,10,1]:
                while (u+m*c < 0) and ZZ(poly_p(u+m*c)).nbits() >= length_p:
                    u += m*c
        u_min.append(u)
    return u_min

def find_max(poly_p, length_p, m=1, umodm=0):
    """ Find maximum in terms of absolute value"""
    RRa = RealField(ceil(length_p*log(2) + 2))
    roots_P = (poly_p - (2**length_p-1)).roots(RRa)
    u_max = []
    for r in roots_P:
        u = ceil(r[0])
        if m > 1: # works also for umod = 0
            u = u - (u % m) + umodm
        if u > 0:
            for c in [100,10,1]:
                while (u-m*c) > 0 and ZZ(poly_p(u)).nbits() > length_p:
                    u -= m*c
            # u is s.t. poly(u) is length_p bits long or less
            for c in [100,10,1]:
                while ZZ(poly_p(u+m*c)).nbits() < (length_p+1):
                    u += m*c
        else:
            for c in [100,10,1]:
                while (u+m*c) < 0 and ZZ(poly_p(u)).nbits() > length_p:
                    u += m*c
            # u is s.t. poly(u) is length_p bits long or less
            for c in [100,10,1]:
                while ZZ(poly_p(u-m*c)).nbits() < (length_p+1):
                    u -= m*c
        u_max.append(u)
    return u_max

def find_min_max_list_mod_m(poly_p, length_p_min, length_p_max, m=1, umodm=[0]):
    umin = find_min(poly_p, length_p_min, m, umodm[-1])
    umax = find_max(poly_p, length_p_max, m, umodm[0])
    # for negative u, min and max are swapped: this is in absolute value
    # so re-swap to have u_min_neg < u_max_neg < 0, and
    # poly_p(u_min_neg) -> length_p_max
    # poly_p(u_max_neg) -> length_p_min
    # because |u_min_neg| > |u_max_neg| > 0
    u_min_pos = ZZ(min([um for um in umin if um > 0]))
    u_max_neg = ZZ(max([um for um in umin if um < 0]))
    u_max_pos = ZZ(max([um for um in umax if um > 0]))
    u_min_neg = ZZ(min([um for um in umax if um < 0]))
    # adjust u_min_pos and u_max_pos
    i = len(umodm)-2
    ok = True
    while i >= 0 and ok:
        u_min = u_min_pos - (u_min_pos % m) + umodm[i]
        ok = u_min < u_min_pos and (ZZ(poly_p(ZZ(u_min)))).nbits() >= length_p_min
        if ok:
            u_min_pos = u_min
        i += 1
    i = 1
    ok = True
    while i < len(umodm) and ok:
        u_max = u_max_pos - (u_max_pos % m) + umodm[i]
        ok = u_max > u_max_pos and (ZZ(poly_p(ZZ(u_max)))).nbits() <= length_p_max
        if ok:
            u_max_pos = u_max
        i += 1
    # adjust u_min_neg and u_max_neg
    # careful: for negative values, u_min_neg <-> length_max and u_max_neg <-> length_min
    i = len(umodm)-2
    ok = True
    while i >= 0 and ok:
        u_min = u_min_neg - (u_min_neg % m) + umodm[i]
        ok = u_min < u_min_neg and (-ZZ(poly_p(ZZ(u_min)))).nbits() <= length_p_max
        if ok:
            u_min_neg = u_min
        i += 1
    i = 1
    ok = True
    while i < len(umodm) and ok:
        u_max = u_max_neg - (u_max_neg % m) + umodm[i]
        ok = u_max > u_max_neg and (-ZZ(poly_p(ZZ(u_max)))).nbits() >= length_p_min
        if ok:
            u_max_neg = u_max
        i += 1
    
    assert ZZ(poly_p(ZZ(u_min_pos))).nbits() == length_p_min
    assert ZZ(poly_p(ZZ(u_min_pos-m))).nbits() == length_p_min-1
    assert ZZ(poly_p(ZZ(u_max_pos))).nbits() == length_p_max
    assert ZZ(poly_p(ZZ(u_max_pos+m))).nbits() == length_p_max+1

    assert ZZ(poly_p(ZZ(u_min_neg))).nbits() == length_p_max
    assert ZZ(poly_p(ZZ(u_min_neg-m))).nbits() == length_p_max+1
    assert ZZ(poly_p(ZZ(u_max_neg))).nbits() == length_p_min
    assert ZZ(poly_p(ZZ(u_max_neg+m))).nbits() == length_p_min-1
    
    return u_min_pos, u_max_pos, u_min_neg, u_max_neg

def lift_simple_root_ZZpoly(six, ri, ui, mi, Li, l, Fl, L):
    """
    Lift a simple root with Hensel lifting

    INPUT:
    - `six`: s_i(x) a univariate polynomial in ZZ[x] (with integer coefficients)
    - `ri`: a root mod l of six lifted in NN, positive integer in {0,...,l-1} s.t. si(ri)=0 mod l
    - `ui`: a residue modulo mi
    - `mi`: a modulus, positive integer
    - `Li`: positive integer
    - `l`: a prime number
    - `Fl`: GF(l) the prime finite field of order l
    - `L`: upper bound on Li

    OUTPUT: (uj, mj, Lj) such that s_i(mj*x + uj) = 0 mod l^(Lj-Li)
    Initially, there is a polynomial s_0(x), and
    s_0(mj*x + uj) = 0 mod l^Lj
    """
    x = six.parent().gen(0)
    while Li < L:
        six = six(l*x + ri)//l
        ri = six.roots(Fl)[0][0]
        ri = ZZ(ri)
        ui = ui + ri * mi
        mi = mi * l
        Li = Li + 1
    return (ui, mi, Li)

def lift_multiple_root_ZZpoly(six, ri, ui, mi, Li, l, Fl, L, verbose=True):
    """
    Lift a multiple root with Hensel lifting

    INPUT:
    - `six`: s_i(x) a univariate polynomial in ZZ[x] (with integer coefficients)
    - `ri`: a root mod l of six lifted in NN, positive integer in {0,...,l-1} s.t. si(ri)=0 mod l
    - `ui`: a residue modulo mi
    - `mi`: a modulus, positive integer
    - `Li`: positive integer
    - `l`: a prime number
    - `Fl`: GF(l) the prime finite field of order l
    - `L`: upper bound on Li

    OUTPUT: set {(uj, mj, Lj)} such that s_i(mj*x + uj) = 0 mod l^(Lj-Li)
    The set can be empty, or it can have several elements
    Because for multiple roots, the lifting pattern varies
    Initially, there is a polynomial s_0(x), and
    s_0(mj*x + uj) = 0 mod l^Lj
    """
    x = six.parent().gen(0)
    S = [(six, ri, ui, mi, Li)]
    # because in Python there is no proper linked list,
    # here two lists are used and the new one erases the old one
    R = [] # for the results
    while len(S) > 0:
        new_S = []
        for (six, ri, ui, mi, Li) in S: # initially the step ui = u_{i-1} + ri*mi was already done before
            six = six(l*x + ri)//l
            cont = gcd([ZZ(sij) for sij in six.list()]) # content
            vv = valuation(cont, l)
            if vv > 0:
                six = six // l**vv
            Lj = Li+1+vv # very important to add vv here
            if Lj >= L:
                R.append((ui, mi, Lj))
                if verbose:
                    print("bad root that stayed bad, obtained uj = {}={:#x} mod mj = {}={:#x}, expected valuation Lj = {}".format(ui, ui, mi, mi, Lj))
            else:
                mj = mi*l # mi*l not mi*l^vv
                roots_six = six.roots(Fl)
                for (ri, ei) in roots_six: # several roots possible?
                    ri = ZZ(ri)
                    uj = ui + ri*mi
                    if ei > 1:
                        new_S.append((six, ri, uj, mj, Lj))
                    else:
                        (uk, mk, Lk) = lift_simple_root_ZZpoly(six, ri, uj, mj, Lj, l, Fl, L)
                        R.append((uk, mk, Lk))
                        if verbose:
                            print("bad root that stabilises, obtained uj={}={:#x} mod mj={}={:#x}, expected valuation Lj = {}".format(uk, uk, mk, mk, Lk))
        S = [ii for ii in new_S] # update the Python list
    return R

def find_seed_congruence_for_adicity_QQpoly(vx, u_mod_m, m, val_l, l, verbose=True):
    """
    compute the congruences of the seed ui mod l^n so that l^val_l | vx(h * l^n + ui)
    return a list of congruences, where l^val_l | vx(ui)

    INPUT:
    - `vx`: polynomial with coefficients in QQ
    - `u_mod_m`: a list of congruences ui mod m so that vx(m*x+ui) is an integer
    - `m`: integer, modulus (a divisor of the ppcm of the denominators of the coeffs of vx)
    - `val_l`: valuation at l for vx(ui)
    - `l`: prime

    RETURN: dict {mi: [ui]} so that vx(mi*x + ui) = 0 mod l^val_l

    It uses the functions lift_simple_root_ZZpoly and lift_multiple_root_ZZpoly
    Inspired by this work:
    Aurore Guillevic and Shashank Singh,
    On the alpha value of polynomials in the TNFS algorithm,
    Appendix A.
    Mathematical Cryptology, number 1.
    """
    l = ZZ(l)
    if not l.is_prime():
        return []
    Fl = GF(l)
    list_congruences = []

    x = vx.parent().gen(0)
    for ui in u_mod_m:
        vix = vx(m*x + ui) # so that vix has integer coefficients
        if len([1 for vi in vix.list() if vi.denom() != 1]) != 0:
            continue # polynomial has not integer coefficients
        c = gcd([ZZ(vi) for vi in vix.list()]) # content of the polynomial
        v = valuation(c, l)
        if v > 0:
            vix = vix // l**v
        # compute the roots mod l
        for (ri,ei) in vix.roots(Fl):
            ri = ZZ(ri)
            uj = ui + ri*m # mod m*l
            if ei == 1: # simple root, easy doing
                # ri is a root of vix mod l
                # now to have the initial vx to be integer and divisible by l, we need u = something mod m*l
                (uj, mj, Lj) = lift_simple_root_ZZpoly(vix, ri, uj, m*l, 1+v, l, Fl, val_l)
                list_congruences.append((uj, mj))
                if verbose:
                    v_uj = ZZ(vx(uj))
                    print("uj={}={:#x} mod mj={}={:#x}={} v(uj)={:#x} {} bits {}-valuation = {}, expected {}".format(uj, uj, mj, mj, mj.factor(), v_uj, v_uj.nbits(), l, v_uj.valuation(l), Lj))
                    vujx = vx(mj*x+uj)
                    cont = gcd([ZZ(vj) for vj in vujx.list()])
                    val_vx = cont.valuation(l)
                    print("{}-valuation of vx(mj*x+uj) = {}".format(l, val_vx))
            else: # multiple root, tricky
                list_bad_iroots_mod_l = lift_multiple_root_ZZpoly(vix, ri, uj, m*l, 1+v, l, Fl, val_l)
                list_congruences = list_congruences + [(uj, mj) for (uj, mj, Lj) in list_bad_iroots_mod_l]
                if verbose:
                    for (uj, mj, Lj) in list_bad_iroots_mod_l:
                        v_uj = ZZ(vx(uj))
                        print("bad root uj={}={:#x} mod mj={}={:#x}={} v(uj)={:#x} {} bits {}-valuation = {}, expected {}".format(uj, uj, mj, mj, mj.factor(), v_uj, v_uj.nbits(), l, v_uj.valuation(l), Lj))
    print("len(list_congruences) = {}".format(len(list_congruences)))
    dict_mj_uj = {}
    for (uj, mj) in list_congruences:
        if mj in dict_mj_uj:
            dict_mj_uj[mj].append(uj)
        else:
            dict_mj_uj[mj] = [uj]
    for mj in dict_mj_uj:
        (dict_mj_uj[mj]).sort()
    print("dict_mj_uj = {}".format(dict_mj_uj))
    return dict_mj_uj

def find_seed_congruence_for_adicity_ZZpoly(vx, val_l, l):
    """
    compute the congruences of the seed ui mod l^n so that l^val_l | vx(h * l^n + ui)
    return a list of congruences, where l^val_l | vx(ui), and N s.t. all u = ui mod N will satisfy the constrains

    INPUT:
    - `vx`: polynomial with coefficients in ZZ
    - `val_l`: valuation at l for vx(ui)
    - `l`: prime integer

    RETURN: dict {mi: [ui]} so that vx(mi*x + ui) = 0 mod l^val_l

    It uses the functions lift_simple_root_ZZpoly and lift_multiple_root_ZZpoly
    Inspired by this work:
    Aurore Guillevic and Shashank Singh,
    On the alpha value of polynomials in the TNFS algorithm,
    Appendix A.
    Mathematical Cryptology, number 1.
    """
    l = ZZ(l)
    if not l.is_prime():
        return {}
    if len([1 for vi in vx.list() if vi.denom() != 1]) != 0:
        return {}
    Fl = GF(l)
    list_congruences = []
    x = vx.parent().gen(0)

    c = gcd([ZZ(vi) for vi in vx.list()]) # content of the polynomial
    v = valuation(c, l)
    vix = vx // l**v
    # compute the roots mod l
    for (ri,ei) in vix.roots(Fl):
        ri = ZZ(ri)
        uj = ri # mod l
        if ei == 1: # simple root, easy doing
            (uj, mj, Lj) = lift_simple_root_ZZpoly(vix, ri, ri, l, 1+v, l, Fl, val_l)
            list_congruences.append((uj, mj))
        else:# multiple root, tricky
            list_bad_iroots_mod_l = lift_multiple_root_ZZpoly(vix, ri, ri, l, 1+v, l, Fl, val_l)
            list_congruences = list_congruences + [(uj, mj) for (uj, mj, Lj) in list_bad_iroots_mod_l]
    dict_mj_uj = {}
    for (uj, mj) in list_congruences:
        if mj in dict_mj_uj:
            dict_mj_uj[mj].append(uj)
        else:
            dict_mj_uj[mj] = [uj]
    for mj in dict_mj_uj:
        (dict_mj_uj[mj]).sort()
    print("dict_mj_uj = {}".format(dict_mj_uj))
    return dict_mj_uj

def str_py_binary_form(u):
    bp = bit_positions(u)
    bp.reverse()
    str_u = ""
    for ei in bp:
        if ei == 0:# this is +1
            str_u += "+1"
        elif ei == 1:
            str_u += "+2"
        else:
            str_u += "+2**{}".format(ei)
    return str_u, len(bp)

def str_binary_form(u):
    bp = bit_positions(u)
    bp.reverse()
    str_u = ""
    for ei in bp:
        if ei == 0:# this is +1
            str_u += "+1"
        elif ei == 1:# this is +2 (not +2^1)
            str_u += "+2"
        else:
            str_u += "+2^{}".format(ei)
    return str_u, len(bp)

def str_py_2NAF(u):
    bp = bit_positions_2naf(u)
    bp.reverse()
    str_u = ""
    for (ei,si) in bp:
        if ei == 0:
            str_u += "{:+d}".format(si)
        elif ei == 1:
            if si == 1:
                str_u += "+2"
            else:
                str_u += "-2"
        else:
            if si == 1:
                str_u += "+2**{}".format(ei)
            else:
                str_u += "-2**{}".format(ei)
    return str_u, len(bp)

def str_2NAF(u):
    bp = bit_positions_2naf(u)
    bp.reverse()
    str_u = ""
    for (ei,si) in bp:
        if ei == 0:
            str_u += "{:+d}".format(si)
        elif ei == 1:
            if si == 1:
                str_u += "+2"
            else:
                str_u += "-2"
        else:
            if si == 1:
                str_u += "+2^{}".format(ei)
            else:
                str_u += "-2^{}".format(ei)
    return str_u, len(bp)


def get_next_seed(u_start, u_stop, m, u_mod_m, poly_p, poly_r, rnbits, allow_cofact_r=False, polys_cofact_twist=[],verbose=False,allowed_cofactor=1,factor_r=False,k=None):
    """ computes the next seed u, starts from u_start and increase/decrease u so that (u % m) in u_m
    if u_start > u_stop, then u should decrease (u = u-m for example)
    find u s.t. p=poly_p(u), r=poly_r(u), tw=poly_r_twist(u) are prime, and u is within the interval of bounds u_start, u_stop
    if rx(u) has a cofactor, then the prime part should be larger than rnbits
    u_start is a valid integer, so that poly_p(u_start) and poly_r(u_start) are integers,
    i.e. (u_start) % m  in u_mod_m
    if valuation is not None and m == 1 then actually u should be multiplied by valuation
    """
    if allow_cofact_r:
        prod_primes = prod(prime_range(10**8))
    else:
        prod_primes = prod(prime_range(10**7))
    i = 0
    # save the gcd of parameters to detect possible systematic cofactor
    gcd_p = 0
    gcd_r = 0
    gcd_ci = [0 for j in range(len(polys_cofact_twist))]
    u0 = u_start
    u = u0
    cond = True
    while (cond) and ((u_start <= u_stop and u_start <= u0 and u0 <= u_stop) or (u_start >= u_stop and u_start >= u0 and u0 >= u_stop)):
        #print("u0 = {:#x} u = {:#x}".format(u0, u))
        p = ZZ(poly_p(u))
        r = ZZ(poly_r(u))
        gcd_pi, p1 = small_factors(p, prod_primes) #gcd(prod_primes, p)
        gcd_ri, r1 = small_factors(r, prod_primes) #gcd(prod_primes, r)
        gcd_p = gcd(gcd_p, gcd_pi)
        gcd_r = gcd(gcd_r, gcd_ri)
        cond = gcd_pi > 1 or (not allow_cofact_r and gcd_ri > 1) or (allow_cofact_r and r1.nbits() < rnbits)
        if not cond: # compute cofactors
            tw= [ZZ(c_i(u)) for c_i in polys_cofact_twist]
            j = 0
            for ci in tw:
                gcd_cii, cii1 = small_factors(ci, prod_primes) #gcd(prod_primes, ci)
                gcd_ci[j] = gcd(gcd_ci[j], gcd_cii)
                cond = cond or gcd_cii > allowed_cofactor
                j += 1
        else:
            tw = []
        cond = cond or not p.is_pseudoprime()
        if not cond:
            cond = cond or not ((not allow_cofact_r and r.is_pseudoprime()) or (allow_cofact_r and r1.is_pseudoprime()))
            if k is not None:
                p_mod_k = "(p={}%{}) ".format((p%k), k)
            else:
                p_mod_k = ""
            if verbose and factor_r:
                print("# factoring r1 of {} bits".format(r1.nbits()))
                fact_r1 = r1.factor()
                r1max = max([ri for (ri,ei) in fact_r1])
                print("u={}={:#x}, p is prime of {} bits {}but r = {} * {} ({} bits)".format(u,u,p.nbits(), p_mod_k, gcd_ri.factor(), fact_r1, r1max.nbits()))
                if r1.nbits() >= rnbits:
                    gcd_ri *= (r1//r1max)
                    r1 = r1max
                    cond = False
            elif verbose and allow_cofact_r:
                print("u={}={:#x}, p is prime of {} bits {}but r = {} * r1; r1 = {} #({} bits)".format(u,u,p.nbits(), p_mod_k, gcd_ri.factor(), r1, r1.nbits()))
        if not cond:
            for ci in tw:
                cond = cond or not ci.is_pseudoprime()
        i +=1
        if verbose and (i % 10**5) == 0:
            print(i)
            print("gcd all p: {} gcd all r: {} gcd all ci: {}".format(gcd_p, gcd_r, gcd_ci))
            # reset
            gcd_p = 0
            gcd_r = 0
            gcd_ci = [0 for j in range(len(polys_cofact_twist))]
        if cond:
            #u += u_step # u_step can be negative
            if m == 1 or len(u_mod_m) == 1:
                if u_start <= u_stop:
                    u0 = u0 + m
                else:
                    u0 = u0 - m
            else:
                u0_mod_m = u0 % m
                idx_u_mod_m = u_mod_m.index(u0_mod_m)
                if u_start <= u_stop:
                    if idx_u_mod_m == len(u_mod_m) - 1:
                        u0 = u0 - u0_mod_m + m + u_mod_m[0]
                    else:
                        u0 = u0 - u0_mod_m + u_mod_m[idx_u_mod_m+1]
                else:
                    if idx_u_mod_m == 0:
                        u0 = u0 - u0_mod_m - m + u_mod_m[-1]
                    else:
                        u0 = u0 - u0_mod_m + u_mod_m[idx_u_mod_m-1]
            u = u0

            
    # end of while, either a u was found or u reached the bound of the interval
    if verbose:
        print("gcd(pi) = {}\ngcd(ri) = {}".format(gcd_p,gcd_r))
    # if u is within the interval:
    if ((u_start <= u_stop and u_start <= u0 and u0 <= u_stop) or (u_start >= u_stop and u_start >= u0 and u0 >= u_stop)):
        #print("u0={}, u={}, i = {}".format(u0, u, i))
        if allow_cofact_r:
            return u, p, r1, [gcd_ri] + tw, i
        else:
            return u, p, r, tw, i
    else:
        return None, None, None, None, i
