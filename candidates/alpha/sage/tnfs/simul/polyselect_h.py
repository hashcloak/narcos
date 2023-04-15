""" functions to precompute polynomials h"""

from sage.rings.integer_ring import Z, ZZ
from sage.rings.real_mpfr import RR
from sage.arith.misc import GCD, gcd
from sage.misc.prandom import randint
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field import NumberField

from sage.libs.pari import pari

from tnfs.simul.polyselect_utils import get_coeffs_from_counter
from tnfs.simul.polyselect_utils import pretty_print_coeffs_from_coeffs, pretty_print_poly_from_coeffs, write_tab_h
from tnfs.simul.polyselect_utils import discard_reverse, discard_palindrome, discard_f_x, discard_reverse_f_x, discard_duplicate, discard_root_1

print("Required: symbolic link to enumerate_sparse_T.py downloaded from")
print("https://gitlab.inria.fr/smasson/cocks-pinch-variant/blob/master/enumerate_sparse_T.py")
print("to generate polynomials of given degree")
from tnfs.simul.enumerate_sparse_T import get_sparse_T_HW_gen, get_sparse_T_uptoHW_gen, get_sparse_T_HW_NAF_gen, get_sparse_T_uptoHW_NAF_gen, bits_2naf, bit_positions

def enumerate_sparse_poly(non_zero_coeffs, deg):
    # a polynomial of degree deg has deg+1 coefficients, but since we force the constant coefficient
    # to be non-zero, it remains only (deg+1)-1 coefficients.
    # the leading coefficient is 1 and is counted in get_sparse_T_HW_gen()
    for support in get_sparse_T_HW_gen(deg,non_zero_coeffs-1):
        sb = [0] + [bi+1 for bi in bit_positions(support)]
        n = len(sb) # this should be equals to Hw
        # for each bit position except the one for deg, we have two choices: -1 or 1.
        # the leading coefficient is necessarily 1:
        # sb[n-1] == deg
        # the second leading coefficient is >= 0 (otherwise this is a duplicate)
        if sb[n-2] == deg-1:
            m = n-2
        else:
            m = n-1
        # the constant coefficient is non-zero
        for s in range(2**m):
            hc = [0 for i in range(deg+1)]
            hc[deg] = 1
            if sb[n-2] == deg-1:
                hc[deg-1] = 1
            si = s
            for i in range(m):
                if (si % 2) == 0:
                    hc[sb[i]] = 1
                else:
                    hc[sb[i]] = -1
                si //= 2
            if not discard_duplicate(hc, deg) and not discard_root_1(hc, deg):
                yield hc

def get_list_irr_polys_nonzero_coeffs(deg, nonzero_coeffs, compute_zeta=True, inv_zeta_Kh_min=0.0, simul_zeta=False, N_simul_zeta=10**6, output_file="", onthefly=False, verbose=False):
    """Returns (or compute on the fly) polynomials of degree deg and nonzero_coeffs 
    Require enumerate_sparse_T.py from gitlab.inria.fr
    :param             deg: degree
    :param  nonzero_coeffs: number of non-zero coefficients, including leading and constant coeffs
    :param    compute_zeta: for polynomials h, compute 1/zeta_Kh(2) where Kh=NumberField(h)
    :param inv_zeta_Kh_min: minimum value of 1/zeta_Kh(2) (between 0.5 and 1.0)
    :param      simul_zeta: if False, use pari-gp lfun, if True, computes an average of coprime ideals for N_simul_zeta(=10^6) samples. For deg > 16, pari does fail, but the other one is very very slow
    :param    N_simul_zeta: number of random samples for experimental value of coprime ideals in Kh
    :param     output_file: filename for output, no extension (will generate .py and .gp)
    :param        onthefly: do not store in a table, only write them on-the-fly to a file (useful for deg >= 12)
    :returns: list of polynomials h (list of coefficients) (formerly: in ZZ['y']), None if onthefly
    the coefficients are in {-1,0,1} and the polynomial is monic.
    """
    max_coeff = 1
    ZZy = PolynomialRing(ZZ, names=('y',)) # needed to check irreducibility
    (y,) = ZZy._first_ngens(1)
    tab_h = []
    non_detected_duplicates = 0
    number_irr = 0 # number of irreducible polynomials h
    number_inv_zeta = 0
    if onthefly and len(output_file) == 0:
        print("Error please provide a filename to write the data on-the-fly")
        return
    elif onthefly:
        out_file = open(output_file, 'w+')
        if not out_file:
            print("# error opening file"+str(output_file))
            return
        print("# results written in file "+str(output_file))
        out_file.write("tab_h_{}_{}_{} = [ \\\n".format(deg,max_coeff,nonzero_coeffs))

    length = deg+1
    for hc in enumerate_sparse_poly(nonzero_coeffs, deg):
        #if gcd(hc) != 1: # the content of the polynomial
        #    continue
        h = ZZy(hc)
        if not h.is_irreducible():
            continue
        number_irr += 1
        # write nicely hc
        hc_string = pretty_print_coeffs_from_coeffs(hc)
        h_string = pretty_print_poly_from_coeffs(hc)
        if compute_zeta:
            Kh = NumberField(h, 'ah')
            number_roots_unity = len(Kh.roots_of_unity())
            if verbose:
                print("h="+h_string)
            recomputed = True
            if not simul_zeta:
                inv_zeta_Kh = RR(1/RR(pari('lfun(' + str(h) + ', 2)')))
            elif inv_zeta_Kh_min == 0.0:
                inv_zeta_Kh = compute_experimental_zeta_Kh(hc, log2_cost=192, samples=N_simul_zeta)
            else:
                # first compute with lower precision,
                # then refine only if the value is likely to be >= inv_zeta_Kh_min
                inv_zeta_Kh = compute_experimental_zeta_Kh(hc, log2_cost=192, samples=N_simul_zeta//10)
                if inv_zeta_Kh >= inv_zeta_Kh_min - 0.1:
                    if verbose:
                        print("inv_zeta_Kh = {0:.8f}, recomputing with better precision".format(float(inv_zeta_Kh)))
                    inv_zeta_Kh = compute_experimental_zeta_Kh(hc, log2_cost=192, samples=N_simul_zeta)
                else:
                    recomputed = False
                        
            if (inv_zeta_Kh >= inv_zeta_Kh_min-0.02) and recomputed:
                number_inv_zeta += 1
                if not onthefly:
                    tab_h.append((inv_zeta_Kh,number_roots_unity,hc))
                elif len(output_file) > 0:
                    out_file.write("    ({0:.8f},{1:2d}, {2}), #{3}\n".format(float(inv_zeta_Kh), number_roots_unity, hc_string, h_string))
                if verbose:
                    print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), number_roots_unity, hc_string, h_string))
            elif verbose:
                print("inv_zeta_Kh = {0:.8f}".format(float(inv_zeta_Kh)))

        else:
            if not onthefly:
                tab_h.append(hc)
            elif len(output_file) > 0:
                out_file.write("    {},\n".format(hc_string))
            if verbose:
                print("    {}, # {}".format(hc_string, h))
    
    print("# There were {} irreducible polynomials".format(number_irr))
    if compute_zeta and inv_zeta_Kh_min > 0.0:
        print("# There were {} where inv_zeta_Kh >= {}".format(number_inv_zeta, inv_zeta_Kh_min))
    if not onthefly and len(output_file)>0:
        if compute_zeta:
            print("# sorting the polynomials by decreasing value of 1/zeta_Kh(2)")
            tab_h.sort(reverse=True)
            write_tab_h(tab_h, deg, "{}_{}".format(max_coeff,nonzero_coeffs), output_file, with_zeta=True)
        else:
            write_tab_h(tab_h, deg, "{}_{}".format(max_coeff,nonzero_coeffs), output_file, with_zeta=False)
    elif onthefly and len(output_file)>0:
        out_file.write("]\n\n")
        out_file.flush()
        out_file.close()
        return None
    return tab_h

def random_vec(A,size):
    return [randint(-A,A) for i in range(size)]

def random_vec_positive_lc(A,size):
    return [randint(-A,A) for i in range(size-1)] + [randint(1,A)]

def compute_experimental_zeta_Kh(hc, samples=10**6,log2_cost=192,verbose=False):
    """pick random pairs and check if they are coprime """
    ZZy = PolynomialRing(ZZ, names=('y',)) # needed to check irreducibility
    (y,) = ZZy._first_ngens(1)
    h = ZZy(hc)
    Kh = NumberField(h, 'ah')
    duplicates = 0
    deg_h = h.degree()
    A = max(10, round(RR(2**(log2_cost/(2*deg_h)))))
    si = 0
    while si < samples:
        si += 1
        #a=random_vec(A,deg_h)
        #b=random_vec(A,deg_h)
        a=random_vec_positive_lc(A,deg_h)
        b=random_vec_positive_lc(A,deg_h)
        coprime_ideal=(Kh.ideal(Kh(a))+Kh.ideal(Kh(b))==Kh.ideal(1))
        while not coprime_ideal:
            si += 1
            a=random_vec(A,deg_h)
            b=random_vec(A,deg_h)
            duplicates += 1
            coprime_ideal=(Kh.ideal(Kh(a))+Kh.ideal(Kh(b))==Kh.ideal(1))
    inv_zeta_Kh = float(RR((samples-duplicates)/samples))
    if verbose:
        print("    ({:.8f}, {}),".format(inv_zeta_Kh, h))
    return inv_zeta_Kh
