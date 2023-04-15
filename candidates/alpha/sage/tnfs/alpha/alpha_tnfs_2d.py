"""
SageMath implementation of alpha for TNFS
function:
alpha_TNFS_2d(f,h,B)

Implementation: Aurore Guillevic, Inria Nancy

Directly inspired from cado-nfs/polyselect/alpha.sage
Written by Bai, Gaudry, Hanrot, Thom{\'e}, Zimmermann
See https://gforge.inria.fr/projects/cado-nfs/
"""
from sage.rings.real_mpfr import RR
from sage.rings.rational_field import Q, QQ
from sage.rings.integer_ring import Z, ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fast_arith import prime_range
from sage.functions.log import log
from sage.arith.functions import lcm
from sage.arith.misc import GCD, gcd
from sage.arith.misc import valuation
from sage.misc.flatten import flatten
import time

def number_of_roots_TNFS(f, FqZ, Fq):
    """Computes the number of roots of f, a polynomial in Kh[x], mod q
    1st, map f to Fq[z], then computes the number of distinct roots"""
    proj_I = False
    fq = FqZ([Fq(coeff_i) for coeff_i in f.list()])
    if fq != 0:
        s = len(fq.roots())
        if fq.degree() < f.degree():
            s += 1
            proj_I = True # the ideal divides the leading coeff
    else:
        # f reducing to zero mod I is a degenerate case. Not clear what
        # we should return...
        # note: this should never happen, otherwise it means that the content of f is not 1
        print("Warning, counting roots of zero polynomial")
        s = f.degree()
    return s, proj_I

def lift_Fl_ZZ(a):
    """lift from ZZ/lZZ to ZZ, minimize the absolute value."""
    l = a.parent().characteristic()
    lift_l = ZZ(a)
    if abs(lift_l) <= abs(lift_l - l):
        return lift_l
    else:
        return lift_l-l

def average_valuation_affine_TNFS_I_principal_rec(I, Norm_I, gam, Fq, FqZ, Kh, Oh, L):
    """ 
    I is a principal prime ideal of Oh
    gam is a generator of I
    L is a list of sub-polys of f (the f(r_i + x*gam)) to be processed
    v is a valuation (only used to return a value)
    this procedure gets the 1st element of L, process it, and
    add at the end of L the new polynomials f2 = f(r_i+gam*x) to be processed."""
    # remove and process 1st element of L
    f,j = L.pop(0)
    # the Content function in Magma does not compute a Gcd in ZZ, but
    # in SageMath it looks ok
    val_i = valuation(f.content_ideal(), I)
    if val_i > 0:
        fv = f // (gam**val_i) # maybe will not work. Will work only if f has coefficients in Oh, does it?
    else:
        fv = f
    val_i = QQ(val_i)
    dfv = fv.derivative()
    # now map coefficients of f and df to Oh/(I) = F_{p^e}
    # to get f and df in F_{p^e} [x]
    fv_FqZ = FqZ([Fq(coeff_i) for coeff_i in fv.list()])
    dfv_FqZ = FqZ([Fq(coeff_i) for coeff_i in dfv.list()])
    # When do we actually use f' mod I?
    for (r,_) in fv_FqZ.roots(): # easy to compute roots in a finite field
        if dfv_FqZ(r) != 0: # r is a simple root, with multiplicity one.
            val_i += QQ(1/(Norm_I-1))
        else:# r is a multiple root, meaning a root of f and f',
            # continue the lifting process until f'(r) !=0
            # if r is in a prime field, there is no r0.list(), r0 itself is it.
            if Fq.degree() == 1:
                lift_r = Kh(lift_Fl_ZZ(r))
            else:
                r0_list = r.polynomial().list()
                len_r0 = len(r0_list)
                coeffs_r0 = [lift_Fl_ZZ(coeff_i) for coeff_i in r0_list] + [ZZ(0) for i in range(len_r0, Kh.degree())]
                lift_r = Kh(coeffs_r0)
            f2 = fv(lift_r + Kh(gam)*fv.variables()[0])
            L.append((f2,j+1))
    v = QQ(val_i/(Norm_I**j))
    return v

def average_valuation_affine_TNFS_I_principal(f, I, Norm_I, gam, Fq, FqZ, Kh, Oh):
    """ f is a polynomial in Oh[x]
        I is a principal prime ideal """
    L = [(f,0)] # list of (f,j) where j is the recursive call depth
    v = QQ(0)
    while len(L) > 0:
        v += average_valuation_affine_TNFS_I_principal_rec(I, Norm_I, gam, Fq, FqZ, Kh, Oh, L)
    return QQ(v)

def map_Kh_Fq(Kh, Fq, c):
    """ c is in Kh, wants to map it to Fq. Kh = QQ[y]/(h(y)), Fq = Fp[z]/(h(z)), same polynomial h """
    # but sometimes, Fq = Fp[z]/m(z) where deg(m) < deg(h)
    if c in QQ:
        return Fq(c) # easy mapping
    #print("c = {}, c.parent() = {}".format(c, c.parent()))
    denom_coeffs = lcm([ZZ(ci.denom())  for ci in c.polynomial().list()])
    #print("gcd denom = {}".format(denom_coeffs))
    coeffs = [ZZ(denom_coeffs*ci) for ci in c.polynomial().list()]
    coeffs += [ZZ(0) for i in range(len(coeffs)+1, Fq.degree())]
    if len(coeffs) == Fq.degree():
        return Fq(coeffs)
    return Fq(sum([Fq.gen(0)**i*coeffs[i] for i in range(len(coeffs))]))

# needs a mapping function from Kh[x1,x2] to Fq[z1,z2]
def map_Kh_Fq_bivariate_pol_ring(Kh_x1x2, Kh, Fq_z1z2, Fq, fv):
    monomials_Fq = [map_Kh_Fq(Kh, Fq, fv.monomial_coefficient(m))* Fq_z1z2.gen(0)**m.degrees()[0]*Fq_z1z2.gen(1)**m.degrees()[1] for m in fv.monomials()]
    return sum(monomials_Fq)

def roots_bivariate_enum(fv, fv_z1z2, Fq, Fq_Z, Oh, Oh1, delta, gam, I):
    # pb: what to do with [0,0]?
    # --> if there is a projective root, then this function is called with Evaluate(Reverse(f), delta*x1+gamma*x2)
    # if val_I(ld(f)) > 1, then val_I
    roots_I = []
    roots_IJ = []
    roots_Fq = []
    diff_deg = Oh.degree() - Fq.degree()
    nb_roots = 0
    for i2 in Fq: # enumerate the values in [1,...,q-1] if q is a prime,
        # enumerate in Fq if q is a prime power
        # because Magma/SageMath does not compute roots of bivariate polynomials
        fi2 = fv_z1z2([Fq_Z.gen(0),i2]) # so that now fi2 is univariate in Z
        if fi2 in Fq: # this is a constant, no root
            list_i1 = []
        else: # fi2.degree() >= 1:
            list_i1 = fi2.roots() # Magma/SageMath knows how to compute roots of univariate polynomials
            
        nb_roots += len(list_i1)
        for (i1,_) in list_i1:
            if (i2 != 0) or (i1 != 0): # how do I do with [0,0] ???
                if Fq.degree() == 1:
                    Ri = lift_Fl_ZZ(i1)
                    Rj = lift_Fl_ZZ(i2)
                    rI = Oh.ideal(Ri*delta + Rj*gam) #assume that Fq is prime
                else:
                    # lift from Fq = Fp^i to ZZ[x]
                    # Magma/SageMath needs a list of exactly the good length, not smaller, so fill with zeros
                    Ri = Oh([lift_Fl_ZZ(i1e) for i1e in i1.polynomial().list()] + [ZZ(0) for k in range(diff_deg)])
                    Rj = Oh([lift_Fl_ZZ(i2e) for i2e in i2.polynomial().list()] + [ZZ(0) for k in range(diff_deg)])
                    rI = Oh.ideal(Ri*delta + Rj*gam)
                # test if the lift of the root is also a root in Kh
                fv_rI = fv(Ri,Rj)
                if (Oh.ideal(fv_rI) + I) == Oh1:
                    continue # actually this is not a root
                # test if the roots are "duplicates",
                # that is Oh.ideal(i1*delta + j1*gam) = Oh.ideal(i2*delta + j2*gam)
                j = 0
                while (j < len(roots_I)) and ((roots_Fq[j] == (0,0)) or (rI + roots_I[j]) == Oh1):# jump over [0,0]
                    j += 1
                if j >= len(roots_I):
                    # the root (i1,i2) = i2*delta + j2*gam is coprime to the other roots
                    roots_I.append(rI)
                    roots_IJ.append((Ri,Rj))
                    roots_Fq.append((i1,i2))
    # now [0,0]
    if (len(roots_I) < (Fq.cardinality()-1)) and (fv_z1z2((0,0)) == 0) and ((fv(0,0) + I) != Oh1):
        # otherwise, adding [0,0] as another root will make a full set of roots:
        # we are supposed to have #roots = #Fq <=> val_I(f) > 0 but this is not the case
        roots_I.append(Oh.ideal(0))
        roots_IJ.append((0,0))
        roots_Fq.append((Fq(0),Fq(0)))
    #print("computing roots by enumeration. Found {} roots, {} were duplicate roots.".format(nb_roots, nb_roots-len(roots_Fq)))
    return roots_Fq, roots_IJ

def average_valuation_affine_TNFS_I_not_principal_rec(I, Norm_I, delta, gam, Fq, FqZ, Fq_z1z2, Kh, Kh_x1x2, Oh, Oh1, L, Reverse_f=False):
    """I is a non-principal prime ideal of Oh -> there is no generator gamma of I
    we can try to deal with two generators <gamma,delta>
    L is a list of sub-polys of f (the f(r_i + x1*gamma+x2*delta)) to be processed
    the polynomials f are bivariate in x1,x2
    this procedure gets the 1st element of L, process it, and
    add at the end of L the new polynomials f2 = f(r_i+gamma*x1+delta*x2) to be processed.
    returns: a valuation v"""
    f,j = L.pop(0)
    val_i = min([valuation(Oh.ideal(fi), I) for fi in f.coefficients()])
    #val_i = valuation(f.content_ideal(), I) # -> not valid for multivariate polynomials
    if val_i > 0:
        fv = f // (gam**val_i)
        # now just re-multiply to get integer coefficients in ZZ instead of in QQ
        # re-multiply by Oh.ideal(gam).norm()/Norm_I?
        lcm_coeffs = lcm([ZZ(fv_ij.denominator()) for fv_ij in flatten([fv_i.list() for fv_i in fv.coefficients()])])
        fv = gcd(ZZ(lcm_coeffs),ZZ(abs(ZZ(Oh.ideal(gam).norm())/Norm_I))**val_i) *fv
    else:
        fv = f
    val_i = QQ(val_i)
    if len(fv.variables()) == 1: # f is actually univariate the 1st time
        dfv = fv.derivative()
        # now map coefficients of f and df to Oh/(I) = F_{p^e}
        # to get f and df in F_{p^e} [x]
        fv_FqZ = FqZ([Fq(coeff_i) for coeff_i in fv.list()])
        dfv_FqZ = FqZ([Fq(coeff_i) for coeff_i in dfv.list()])
        for (r,_) in fv_FqZ.roots(): # easy to compute roots in a finite field
            #if not (Evaluate(dfv, lift_r) in I) then // fv'(r) != 0
            if (dfv_FqZ(r) != 0): # fv'(r) != 0
                # r is a simple root, with multiplicity one.
                val_i += QQ(1/(Norm_I-1))
            else: # r is a multiple root, meaning a root of f and f',
                # continue the lifting process until f'(r) !=0
                # if r is in a prime field, there is no r0.list(), r0 itself is it.
                if Fq.degree() == 1:
                    lift_r = Kh(lift_Fl_ZZ(r))
                else:
                    r0_list = r.polynomial().list()
                    len_r0 = len(r0_list)
                    coeffs_r0 = [lift_Fl_ZZ(coeff_i) for coeff_i in r0_list] + [ZZ(0) for i in range(len_r0, Kh.degree())]
                    lift_r = Kh(coeffs_r0)
                f2 = fv(lift_r + Kh(delta)*Kh_x1x2.gen(0) + Kh(gam)*Kh_x1x2.gen(1))
                L.append((f2,j+1))
        v = QQ(val_i/(Norm_I**j))
        return v
    else:
        # if j > 0 then fv will be multivariate and derivative is with respect to a given variable, e.g.
        dfv1 = fv.derivative(f.variables()[0])
        dfv2 = fv.derivative(f.variables()[1])
        fv_Fq_z1z2 = map_Kh_Fq_bivariate_pol_ring(Kh_x1x2, Kh, Fq_z1z2, Fq, fv)
        dfv1_Fq_z1z2 = map_Kh_Fq_bivariate_pol_ring(Kh_x1x2, Kh, Fq_z1z2, Fq, dfv1)
        dfv2_Fq_z1z2 = map_Kh_Fq_bivariate_pol_ring(Kh_x1x2, Kh, Fq_z1z2, Fq, dfv2)
        # f is bivariate, the roots function is not defined.
        if fv_Fq_z1z2 in Fq: # actually it is a constant, no root
            v = QQ(val_i/(Norm_I**j))
        else: # roots of bivariate polynomial not implemented
            list_r, list_r_Oh = roots_bivariate_enum(fv, fv_Fq_z1z2, Fq, FqZ, Oh, Oh1, delta, gam, I)
            for i_r in range(len(list_r)): # r is [r1,r2]
                r = list_r[i_r] ; r0, r1 = r
                lift_r = list_r_Oh[i_r]
                if (dfv1_Fq_z1z2(r) != 0) or (dfv2_Fq_z1z2(r) != 0):
                    val_i += QQ(1/(Norm_I-1))
                else:
                    f2 = fv([lift_r[0] + Kh(delta)*Kh_x1x2.gen(0), lift_r[1] + Kh(gam)*Kh_x1x2.gen(1)])
                    L.append((f2,j+1))
            v = QQ(val_i/(Norm_I**j))
    return v
    
def average_valuation_affine_TNFS_I_not_principal(f, I, Norm_I, delta, gam, Fq, FqZ, Fq_z1z2, Kh, Kh_x1x2, Oh, Oh1):
    """ f is univariate or bivariate """
    L = [(f,0)]
    v = QQ(0)
    while len(L) > 0:
        v += average_valuation_affine_TNFS_I_not_principal_rec(I, Norm_I, delta, gam, Fq, FqZ, Fq_z1z2, Kh, Kh_x1x2, Oh, Oh1, L)
    return QQ(v)

def is_bad_ideal(I, disc, Oh1):
    return ((disc + I) != Oh1)

def average_valuation_homogeneous_coprime(f,I_disc_f,I,Norm_I,Kh,Oh,Oh1,Kh_x1x2,test_principal=True):
    """ return the average valuation and wether I is a bad prime or not """
    Fq = I.residue_field('z')
    #FqZ.<Z> = Fq[]
    FqZ = Fq['Z']; (Z,) = FqZ._first_ngens(1)
    t_good = RR(0.0)
    t_bad = RR(0.0)
    if ((I_disc_f + I) == Oh1): # good primes, exact formula
        t1 = time.time()
        n_roots, proj_ideal = number_of_roots_TNFS(f,FqZ,Fq)
        t2 = time.time()
        t_good += RR(t2-t1)
        return QQ(n_roots/(Norm_I-1)*Norm_I/(Norm_I+1)), False, proj_ideal, t_good # nr * p / (p^2-1)
    else:# bad prime ideal, recursive formula
        #print("    bad prime")
        t1 = time.time()        
        ld_I = Oh.ideal(f.leading_coefficient())
        proj_ideal = (ld_I + I) != Oh1
        #if proj_ideal:
        #    print("      projective ideal")
        #else:
        #    print("      integral ideal")
        #print("test_principal = {}".format(test_principal))
        # test if gam == 0 because in that case either I = (l) or actually I is principal
        delta, gam = I.gens_two() # delta=l in ZZ, gam with 'w'
        #print("l={}, I = {}, delta,gam = {}, {}".format(Norm_I,I,delta,gam))
        if gam == 0 or (gam.norm() == Norm_I) or (test_principal and I.is_principal()):
            #print("      is principal")
            if gam == 0:
                gam = delta
            #elif gam.norm() == Norm_I: then gam plays the role of generator
            elif (gam.norm() != Norm_I):
                gam = I.gens_reduced()[0]
            r = average_valuation_affine_TNFS_I_principal(f, I, Norm_I, gam, Fq, FqZ, Kh, Oh) * Norm_I
            # only if I is not coprime to the leading coefficient of f:
            if proj_ideal:
                f_1x = f.reverse()(gam*f.variables()[0])
                r_pr = average_valuation_affine_TNFS_I_principal(f_1x, I, Norm_I, gam, Fq, FqZ, Kh, Oh)
                r += r_pr
            r = r/(Norm_I+1)
        else:
            #print("      not principal")
            Fq_z1z2 = Fq['z1, z2']; (z1, z2,) = Fq_z1z2._first_ngens(2)
            r = average_valuation_affine_TNFS_I_not_principal(f, I, Norm_I, delta, gam, Fq, FqZ, Fq_z1z2, Kh, Kh_x1x2, Oh, Oh1) * Norm_I
            # only if I is not coprime to the leading coefficient of f:
            if proj_ideal:
                f_1x = (f.reverse())(delta*(Kh_x1x2.gen(0)) + gam*(Kh_x1x2.gen(1)))
                pr = average_valuation_affine_TNFS_I_not_principal(f_1x, I, Norm_I, delta, gam, Fq, FqZ, Fq_z1z2, Kh, Kh_x1x2, Oh, Oh1)
                r += pr
            r = r/(Norm_I+1)
        t2 = time.time()
        t_bad += RR(t2-t1)
        return QQ(r), True, proj_ideal, t_bad

def alpha_TNFS_l(f,I_disc_f,l,Kh,Oh,Oh1,Kh_x1x2,test_principal=True):
    sum_val_I = QQ(0)
    time_good = RR(0.0)
    time_bad = RR(0.0)
    good = 0
    bad = 0
    
    """
    Fl = GF(l)
    Flz = PolynomialRing(Fl, names=('z',))
    (z,) = Flz._first_ngens(1)
    list_ideals_l = [Oh.ideal([l, fac_l_i]) for (fac_l_i,_) in Flz(f).factor() if fac_l_i.degree() < f.degree()]
    if len(list_ideals_l) == 0:
        list_ideals_l = [Oh.ideal(l)]
    print("l={}, list_ideals_l = {}".format(l, list_ideals_l))
    for I in list_ideals_l:
        if not I.is_prime():
            continue
    """
    #print("l={}".format(l))
    for II in Oh.ideal(l).factor():
        I = II[0]
        #print("I = {}".format(I))
        #Norm_I = I.norm()
        #d_I = valuation(Norm_I, l)
        # or:
        d_I = I.residue_class_degree()
        Norm_I = l**d_I
        # Computes the contribution at prime ideal I of the alpha value of f
        r, is_bad, proj_ideal, time_I = average_valuation_homogeneous_coprime(f,I_disc_f,I,Norm_I,Kh,Oh,Oh1,Kh_x1x2,test_principal=test_principal)
        sum_val_I += d_I*r
        #if is_bad:
        #    print("    was bad")
        #if proj_ideal:
        #    print("    was proj")
        # in case there are ideals of different degrees above the same prime (it happens for h of degree 3 for instance),
        # multiply by d_I inside the loop, where d_I is the degree of the ideal, and log(|Norm(I)|) = d_I*Log(l)
        if is_bad:
            time_bad += time_I
            bad += 1
        else:
            time_good += time_I
            good += 1
    return RR( log(l) * ( 1/(l-1)-sum_val_I)), time_good, time_bad, good, bad

    

def alpha_TNFS_2d(f,h,B,test_principal=True,print_timings=False):
    """
    Computes the alpha value of f, up to prime bound B
    f: univariate polynomial (in ZZ[y][x], coefficients in ZZ[y])
    h: univariate polynomial in ZZ[y]
    B: bound on the primes to use, usually 1000 or 2000
    """
    # Kh.<w> = NumberField(h)
    Kh = NumberField(h, names=('w',)); (w,) = Kh._first_ngens(1)
    KhX = Kh['X']; (X,) = KhX._first_ngens(1)
    fh = KhX([Kh(fi) for fi in f.list()])
    # f has coefficients in variable 'y'
    # h is a polynomial in variable 'y'
    # Kh_x1x2.<x1,x2> = Kh[]
    Kh_x1x2 = Kh['x1, x2']; (x1, x2,) = Kh_x1x2._first_ngens(2)
    Oh = Kh.maximal_order()
    Oh1 = Oh.ideal(1)
    disc_f = Kh(fh.discriminant()) # if f in ZZ[y][x], disc in ZZ[y], if f in Kh[x], disc in Kh
    I_disc_f = Oh.ideal(disc_f)
    sum_alpha_l = RR(0.0)
    time_good = RR(0.0)
    time_bad = RR(0.0)
    good = 0 # counter
    bad = 0  # counter
    for l in prime_range(B+1):
        sum_alpha_l_i, time_good_i, time_bad_i, good_i, bad_i = alpha_TNFS_l(fh, I_disc_f, l, Kh, Oh, Oh1, Kh_x1x2,test_principal=test_principal)
        sum_alpha_l += sum_alpha_l_i
        time_good += time_good_i
        time_bad += time_bad_i
        good += good_i
        bad += bad_i
    if print_timings:
        print("good ideals: {:4d}, time: {:12.6f} ms, per ideal: {:12.6f} ms".format(good, float(time_good*1000), float(RR(time_good*1000/good))))
        print("bad  ideals: {:4d}, time: {:12.6f} ms, per ideal: {:12.6f} ms".format(bad, float(time_bad*1000), float(RR(time_bad*1000/bad))))
    return sum_alpha_l
    #return sum([alpha_TNFS_l(f, I_disc_f, l, Kh, Oh, Oh1, Kh_x1x2) for l in prime_range(B+1)])
