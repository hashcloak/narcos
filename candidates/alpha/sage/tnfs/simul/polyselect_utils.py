from sage.functions.generalized import sgn, sign
from sage.rings.integer_ring import Z, ZZ
from sage.rings.rational_field import Q, QQ
from sage.rings.real_mpfr import RR
from sage.arith.misc import GCD, gcd
from sage.arith.functions import lcm
from sage.misc.functional import minimal_polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field import NumberField

from sage.libs.pari import pari

def number_poly(deg, max_coeff, monic=True):
    if monic:
        return (2*max_coeff)*(2*max_coeff+1)**(deg-2)*(max_coeff+1)
    else:
        return (2*max_coeff)*(2*max_coeff+1)**(deg-2)*(max_coeff+1)*max_coeff

def pretty_print_coeffs_from_coeffs(hc):
    hc_string = "["
    for jj in range(len(hc)-2):
        hc_string += "{:2d},".format(hc[jj])
    hc_string += "{:1d},".format(hc[-2]) # second leading coeff is >= 0
    hc_string += "{:1d}]".format(hc[-1]) # leading coeff is > 0
    return hc_string

def pretty_print_poly_from_coeffs(hc):
    if hc[-1] == 1:
        h_string="y^{}".format(len(hc)-1)
    else:
        h_string="{}*y^{}".format(hc[-1], len(hc)-1)
    for jj in range(len(hc)-2,1,-1):
        if hc[jj] > 1 or hc[jj] < -1:
            h_string += "{:+d}*y^{}".format(hc[jj], jj) # print in format where there is always a sign + or -
        elif hc[jj] == -1:
            h_string += "-y^{}".format(jj) # print in format where there is always a sign + or -
        elif hc[jj] == 1:
            h_string += "+y^{}".format(jj) # print in format where there is always a sign + or -
    if hc[1] > 1 or hc[1] < -1:
        h_string += "{:+d}*y".format(hc[1])
    elif hc[1] == -1:
        h_string += "-y"
    elif hc[1] == 1:
        h_string += "+y"
    if hc[0] > 1 or hc[0] < -1:
        h_string += "{:+d}".format(hc[0])
    elif hc[0] == -1:
        h_string += "-1"
    elif hc[0] == 1:
        h_string += "+1"
    return h_string

def write_tab_h(tab_h, deg, max_coeff, filename, header="", with_zeta=False):
    out_file = open(filename, 'w+')
    if not out_file:
        print("# error opening file"+str(filename))
        return
    print("# results written in file "+str(filename))
    if len(header) > 0:
        out_file.write(header)
    out_file.write("tab_h_{}_{} = [ \\\n".format(deg,max_coeff))
    if with_zeta:
        for i in range(len(tab_h)):
            inv_z, w, hc = tab_h[i]
            hc_string = pretty_print_coeffs_from_coeffs(hc)
            h_string = pretty_print_poly_from_coeffs(hc)
            out_file.write("    ({:.8f},{:2d},{}),#{}\n".format(float(inv_z), int(w), hc_string, h_string))
        out_file.write("]\n\n")
    else:
        for i in range(len(tab_h)):
            hc = tab_h[i]
            hc_string = pretty_print_coeffs_from_coeffs(hc)
            h_string = pretty_print_poly_from_coeffs(hc)
            out_file.write("    {},#{}\n".format(hc_string, h_string))
        out_file.write("]\n\n")
    out_file.flush()
    out_file.close()

def is_palindrome(h_coeffs):
    sign_h0 = sgn(h_coeffs[0]) # because h[0] is an int
    deg = len(h_coeffs)-1
    j = 0
    while (2*j <= deg) and h_coeffs[deg-j] == h_coeffs[j]*sign_h0:
        j += 1
    return (2*j > deg)

def is_antipalindrome(h_coeffs):
    """True if h(-y).reverse() = (+/-) h or h.reverse()(-y) = (+/-) h"""
    # 1. h.reverse()(-y)
    sign_h0 = sgn(h_coeffs[0]) # because h[0] is an int
    reverse_h = [sign_h0*h_coeffs[i] for i in range(len(h_coeffs)-1, -1, -1)]
    # now reverse_h(-y):
    reverse_h_y = [hi for hi in reverse_h]
    for i in range(1, len(reverse_h_y), 2):
        reverse_h_y[i] = -reverse_h_y[i]
    if reverse_h_y[-1] < 0:
        reverse_h_y = [-hi for hi in reverse_h_y]
    i=0
    while i < len(reverse_h_y) and reverse_h_y[i] == h_coeffs[i]:
        i += 1
    if (i >= len(reverse_h_y)):
        return True
    # 2. h(-y).reverse()
    h_y = [hi for hi in h_coeffs]
    for i in range(1, len(h_y), 2):
        h_y[i] = -h_y[i]
    h_y.reverse()
    if h_y[-1] < 0:
        h_y = [-hi for hi in h_y]
    i=0
    while i < len(h_y) and h_y[i] == h_coeffs[i]:
        i += 1
    return (i >= len(h_y))

def is_even(h_coeffs):
    deg = len(h_coeffs)-1
    if (deg % 2) == 1:
        return False
    j = 1
    for j in range(1, deg, 2):
        if h_coeffs[j] != 0:
            return False
    return True

def automorphism_factor(h_coeffs):
    if len(h_coeffs) == 3: # h deg 2, x -> -x-h1 is an automorphism, where h1 is tiny
        return 2
    if len(h_coeffs) == 2:
        return 1
    if len(h_coeffs) <= 1:
        return 0
    aut_h = 1
    if is_antipalindrome(h_coeffs): # x -> -1/x is an automorphism
        aut_h = 2
    if is_palindrome(h_coeffs): # x -> 1/x is an automorphism
        aut_h = 2
    if is_even(h_coeffs): # x -> -x is an automorphism
        aut_h *= 2
    return aut_h

def divide_degree_by_d(coeffs_poly, d, u):
    """ input: an array of coefficients of a univariate polynomial p
    returns a list of coefficients of a polynomial P of degree deg/d such that P(u^d) = p(u)
    """
    m = len(coeffs_poly)
    cP = [0]*(((m-1) // d) +1)
    j=0
    while d*j < m:
        i=0
        while (i < d) and (d*j+i < m):
            cP[j] += coeffs_poly[d*j+i]*u**i
            i+=1
        j += 1
    return cP, [-u**d, 1]

def divide_degree_by_two_even(coeffs_poly, u):
    if not is_even(coeffs_poly):
        raise ValueError("Error the polynomial is not even, there is no automorphism x -> -x")
    return [coeffs_poly[i] for i in range(0, len(coeffs_poly),2)], [-u**2, 1]

def divide_degree_by_two_palindrome(coeffs_poly, u):
    if not is_palindrome(coeffs_poly):
        raise ValueError("Error the polynomial is not a palindrome, there is no automorphism x -> 1/x")
    # px palindrome, computes minimal_polynomial(a+1/a) in QQ[a]/(px(a))
    Rx = QQ['x']
    Kp = NumberField(Rx(coeffs_poly), 'v')
    v = Kp.gens()[0]
    mx = minimal_polynomial(v+1/v)
    mx_denom = lcm([mi.denom() for mi in mx.list()])
    poly_P = mx_denom*mx
    cP = [ZZ(pi) for pi in poly_P.list()]
    #u_sq = u+1/u # x - (u+1/u) <=>  u*x - (u**2+1)
    # (u+1/u) is root of mx
    poly_u = [-(u**2+1), u]
    return cP, poly_u

def discard_reverse(f,deg):
    """Return True if f should be discarded because reverse(f) > f w.r.t. ordering"""
    # since f and the reversed polynomial are equivalent, we can assume
    # |f[deg]| < |f[0]| or (|f[deg]| = |f[0]| and |f[deg-1]| < |f[1]|) or ...
    # (c) Paul Zimmermann, cado-nfs
    j = 0
    while (2*j < deg) and abs(f[deg-j]) == abs(f[j]):
        j += 1
    if ((2*j < deg) and (abs(f[deg-j]) > abs(f[j]))):
        return True
    return False
    # if 2j == deg then f is a palindrome w.r.t. absolute value of coefficients

def discard_palindrome(f, deg):    
    """return True if f should be discarded because palindrome(f) > f w.r.t. ordering """
    j = 0
    while (2*j < deg) and abs(f[deg-j]) == abs(f[j]):
        j += 1
    is_palindrome = (j*2 >= deg)
    if is_palindrome:# and f[deg-1]>0: # to avoid overlap with next test
        # this one is a palindrome for the absolute value of coefficients.
        # there might be duplicates: -reverse(f), +/-reverse(f)(-x), +/-reverse(f(-x))
        # TODO
        # duplicate: f.reverse()*sgn(f[0])
        # order the signs: arbitrarily choose to keep the polys with positive high degree
        # issue: conflict with the next function that discards (-1)^deg*f(-x)
        sign_f0 = sgn(f[0]) # because f[0] is an int
        j = 0
        while (2*j <= deg) and f[deg-j] == f[j]*sign_f0:
            j += 1
        if ((2*j <= deg) and f[deg-j] < 0): # what about f[deg-j] == 0 ? # or (2*j == deg and f[deg-j] < 0 and sign_f0 == -1):
            #print("h = {:48s}, discarding. stopped at j={}, f[deg-{}] = {}, f[{}]*({}) = {} ".format(f,j,j,f[deg-j],j,sign_f0,sign_f0*f[j]))
            #rev_f = f ; rev_f.reverse()
            #if sign_f0 < 0:
            #    rev_f = [-ri for ri in rev_f]
            #print("r = {:48s}".format(rev_f))
            return True
        #if (j*2 > deg): # h == h.reverse()*sign_h0 -> it is not possible with two distinct counters, there is no duplicate
    return False

def discard_f_x(f,deg):
    """ Return True if f should be discarded because f(-x) > f(x) w.r.t. ordering """
    # if (deg%2) == 0 we compare to f(-x) and the even coefficients are the same, the odd coefficients changed
    # f[deg-1], f[deg-3] ... changed (deg is even)
    #else we compare to -f(-x) and the even coefficients are changed, the odd coefficients changed
    # f[deg-1], f[deg-3] ... changed (deg is odd)
    #
    # detect duplicate (Paul Zimmermann)
    # Since f(x) is equivalent to (-1)^deg*f(-x), if f[deg-1] = 0, then the
    # largest i = deg-3, deg-5, ..., such that f[i] != 0 should have f[i] > 0.
    # Out of the 44070 remaining polynomials for degree 4 and bound = 6, this test
    # discards 3276, i.e., about 7.4%
    if f[deg-1] == 0:
        j = deg-3
        while j >= 0 and f[j] == 0:
            j -= 2
        if (j >= 0) and f[j] < 0:
            return True
    return False
    
def discard_reverse_f_x(f,deg):
    """ Return False if f should be discarded because reverse(f(-x)) > f(x) w.r.t. ordering """
    # compare to reverse(f(-x))
    f_x = [f[i] if (i%2) == 0 else -f[i] for i in range(deg+1)]
    if f_x[0] < 0:
        rf_x = [-f_x[deg-i] for i in range(deg+1)]
    else:
        rf_x = [f_x[deg-i] for i in range(deg+1)]

    # abort if reverse(f(-x)) appears not to be in the list already
    if rf_x[deg-1] < 0:
        return False # can continue with this f, reverse(f(-x)) is not valid and would be discarded
    if discard_reverse(rf_x,deg):
        return False
    # condition f(-x) but for reverse(f(-x))
    if discard_f_x(rf_x,deg):
        return False # can continue with this f, reverse(f(-x)) would be discarded
    if discard_palindrome(rf_x, deg):
        return False # can continue with this f, reverse(f(-x)) is not valid and would be discarded
    # now decide of an order between f and rf_x and see if we keep f
    j = deg
    while j >= 0 and f[j] == rf_x[j]:
        j -= 1
    # choose the one whose diff coeff of highest degree is zero, or of min abs value, or positive.
    return (j >= 0) and (rf_x[j] == 0 or abs(rf_x[j]) < abs(f[j]) or rf_x[j] > f[j])

def discard_duplicate(fc, deg):
    return discard_reverse(fc,deg) \
        or discard_f_x(fc,deg) \
        or discard_palindrome(fc,deg) \
        or discard_reverse_f_x(fc,deg)

def discard_root_1(fc, deg):
    f_even = sum([fc[i] for i in range(0,deg+1,2)])
    f_odd = sum([fc[i] for i in range(1,deg+1,2)])
    return (f_even + f_odd == 0) or (f_even - f_odd == 0)

def get_coeffs_from_counter(counter, deg, max_coeff=1, monic=True, jump_duplicates=False):
    """ return a list of coefficients as a decoding of ``counter`` in basis ``max_coeff``
    :param   counter: integer
    :param       deg: degree of the corresponding polynomial, there are deg+1 coeffs
    :param max_coeff: bound (included) on the coefficients
    :param     monic: the corresponding polynomial f is monic (f has leading coefficient 1)

    last coeff is strictly positive (poly will have strictly positive leading coefficient)
    source: cado-nfs/polyselect/dlpolyselect.c:polygen_JL_f
    and cado-nfs/utils/mpz_poly.cpp:3930:mpz_poly_setcoeffs_counter() from Paul Zimmermann
    """
    next_counter = counter
    
    f = [int(0) for i in range(deg+1)]
    if monic:
        f[deg] = 1
    else:
        f[deg] = 1 + (counter % max_coeff) # f[deg] > 0
        counter = counter // max_coeff
    f[deg-1] = counter % (max_coeff+1)     # f[deg-1] >= 0
    counter = counter // (max_coeff+1)
    for i in range(deg-2,0,-1):
        f[i] = (counter % (2 * max_coeff + 1)) - max_coeff
        counter = counter // (2*max_coeff + 1)
    # non-zero constant coefficient
    f0 = (counter % (2 * max_coeff)) - max_coeff
    if f0 < 0:
        f[0] = f0
    else:
        f[0] = f0+1

    # computes the next valid counter (jump duplicates)
    counter = next_counter
    next_counter += 1

    # does not work properly
    #if jump_duplicates and not monic and abs(f[deg]) >= abs(f[0]) and max_coeff > 1:
    #    next_counter = counter - (counter % max_coeff) + max_coeff
    #    if(abs(f[0]) == 1):
    #        if deg > 2 and abs(f[deg-1]+1) > abs(f[1]):
    #            idx = next_counter
    #            next_counter = idx - (idx % (max_coeff*(max_coeff+1))) + (idx % max_coeff)
    #            next_counter += max_coeff*(max_coeff+1)
    #            j = 2
    #            mod_j = max_coeff*(max_coeff+1)*(2*max_coeff+1)
    #            mod_k = mod_j*(2*max_coeff+1)
    #            while ((deg-j > j) and (f[j-1] == 0)):
    #                if (abs(f[deg-j]+1) > abs(f[j])):
    #                    if ((f[deg-j]+1) < f[j]):
    #                        idx = next_counter
    #                        next_counter = idx - (idx % mod_j) + (-abs(f[j])+max_coeff)*mod_j
    #                    else:
    #                        idx = next_counter
    #                        next_counter = idx - (idx % mod_j) + (-abs(f[j])+max_coeff)*mod_j + mod_k
    #                j += 1
    #                mod_j = mod_k
    #                mod_k *= (2*max_coeff+1)
    ## tentative
    #elif jump_duplicates and monic and 1 == abs(f[0]) and max_coeff > 1:
    #    # jump the polynomials x^d + m*x^(d-1) + ... + a*x +/- 1 where m > a
    #    if deg > 3 and abs(f[1]) == 1:
    #        next_counter = counter - (counter % (max_coeff+1)) + max_coeff+1
    #        if(abs(f[deg-2]+1) > abs(f[2])): # there will be a +1 in f[deg-2] because of counter + max_coeff+1
    #            idx = next_counter
    #            next_counter = idx - (idx % (max_coeff*(max_coeff+1))) + (idx % max_coeff)
    #            next_counter += max_coeff*(max_coeff+1)
    #            j = 3
    #            mod_j = max_coeff*(max_coeff+1)*(2*max_coeff+1)
    #            mod_k = mod_j*(2*max_coeff+1)
    #            while ((deg-j > j) and (f[j-1] == 0)):
    #                if (abs(f[deg-j]+1) > abs(f[j])):
    #                    if ((f[deg-j]+1) < f[j]):
    #                        idx = next_counter
    #                        next_counter = idx - (idx % mod_j) + (-abs(f[j])+max_coeff)*mod_j
    #                    else:
    #                        idx = next_counter
    #                        next_counter = idx - (idx % mod_j) + (-abs(f[j])+max_coeff)*mod_j + mod_k
    #                j += 1
    #                mod_j = mod_k
    #                mod_k *= (2*max_coeff+1)

    if f[0] == 0 or f[-1] == 0:
        return None, next_counter
    # detect duplicates
    if discard_duplicate(f,deg):
        return None, next_counter
    # check that 1, -1 are not roots
    if discard_root_1(f, deg):
        return None, next_counter
    # content of the polynomial
    if gcd(f) != 1:
        return None, next_counter

    return f, next_counter

def get_list_irr_poly(deg, max_coeff=1, monic=True, compute_zeta=False, inv_zeta_Kh_min=0.0, simul_zeta=False, N_simul_zeta=10**6, output_file="", onthefly=False, start_counter=None, stop_counter=None, verbose=False):
    """ get the list of (monic) irreducible univariate polynomials of degree ``deg`` and coefficients bounded by ``max_coeff``
    :param             deg: degree
    :param       max_coeff: bound on the coefficients, inclusive
    :param           monic: monic poly (True/False)
    :param    compute_zeta: for polynomials h, compute 1/zeta_Kh(2) where Kh=NumberField(h)
    :param inv_zeta_Kh_min: minimum value of 1/zeta_Kh(2) (between 0.5 and 1.0)
    :param      simul_zeta: if False, use pari-gp lfun, if True, computes an average of coprime ideals for N_simul_zeta(=10^6) samples. For deg > 16, pari does fail, but the other one is very very slow
    :param    N_simul_zeta: number of random samples for experimental value of coprime ideals in Kh
    :param   start_counter: for parallel running of this function
    :param    stop_counter: for parallel running of this function
    :param     output_file: filename for output, no extension (will generate .py and .gp)
    :param        onthefly: do not store in a table, only write them on-the-fly to a file (useful for deg >= 12)
    :returns: list of polynomials h (list of coefficients) (formerly: in ZZ['y']), None if onthefly

    len(tab_h_6) = 84
    len(tab_h_7) = 229
    len(tab_h_8) = 731
    len(tab_h_9) = 2165
    len(tab_h_10) = 6736
    len(tab_h_11) = 20396
    len(tab_h_12) = 62457
    len(tab_h_13) = 189719
    len(tab_h_14) = 578504
    len(tab_h_15) = 1758123
    len(tab_h_16) = 5344377

    len(tab_h_18) = 49075276
    """
    ZZy = PolynomialRing(ZZ, names=('y',)) # needed to check irreducibility
    (y,) = ZZy._first_ngens(1)
    tab_h = []
    non_detected_duplicates = 0
    number_irr = 0 # number of irreducible polynomials h
    number_inv_zeta = 0
    counter = 0
    max_counter = number_poly(deg, max_coeff, monic=monic)
    if start_counter != None and start_counter >= 0 and start_counter < max_counter:
        counter = start_counter
    if stop_counter != None and stop_counter >= 0 and stop_counter < max_counter:
        max_counter = stop_counter
    min_counter = counter
    
    if onthefly and len(output_file) == 0:
        print("Error please provide a filename to write the data on-the-fly")
        return
    elif onthefly:
        out_file = open(output_file, 'w+')
        if not out_file:
            print("# error opening file"+str(output_file))
            return
        print("# results written in file "+str(output_file))
        out_file.write("tab_h_{}_{} = [ \\\n".format(deg,max_coeff))
    
    while counter < max_counter:
        hc, next_counter = get_coeffs_from_counter(counter, deg, max_coeff, monic=monic)
        counter = next_counter
        if hc == None:
            continue
        mc = max([abs(ai) for ai in hc])
        # do not consider the polynomials with all abs(coeffs) < max_coeff
        if mc != max_coeff:
            continue

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
                elif verbose:
                    print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), number_roots_unity, hc_string, h_string))

        else:
            if not onthefly:
                tab_h.append(hc)
            elif len(output_file) > 0:
                out_file.write("    {},\n".format(hc_string))
            if verbose:
                print("    {}, # {}".format(hc_string, h))

    print("# There were {} entries for counter from {} to {}".format(number_irr, min_counter, max_counter))
    if compute_zeta and inv_zeta_Kh_min > 0.0:
        print("# There were {} where inv_zeta_Kh >= {}".format(number_inv_zeta, inv_zeta_Kh_min))
    if not onthefly and len(output_file) > 0:
        if compute_zeta:
            print("# sorting the polynomials by decreasing value of 1/zeta_Kh(2)")
            tab_h.sort(reverse=True)
            write_tab_h(tab_h, deg, max_coeff, output_file, with_zeta=True)
        else:
            write_tab_h(tab_h, deg, max_coeff, output_file, with_zeta=False)
    elif onthefly and len(output_file) > 0:
        out_file.write("]\n\n")
        out_file.flush()
        out_file.close()
        return None
    return tab_h

