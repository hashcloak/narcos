import sage

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
from operator import itemgetter
from sage.functions.log import log
from sage.functions.other import ceil
from sage.rings.integer import Integer
from sage.rings.integer_ring import Z, ZZ
from sage.rings.real_mpfr import RR
from sage.arith.misc import GCD, gcd
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field import NumberField
from sage.matrix.constructor import Matrix
from sage.misc.prandom import randint

from sage.libs.pari import pari

from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d
from tnfs.alpha.alpha2d import alpha2d
from tnfs.alpha.alpha3d import alpha3d

from tnfs.simul.enumerate_sparse_T import get_sparse_T_HW_gen, get_sparse_T_uptoHW_gen, get_sparse_T_HW_NAF_gen, get_sparse_T_uptoHW_NAF_gen, bits_2naf
# file downloaded from
# https://gitlab.inria.fr/smasson/cocks-pinch-variant/blob/master/enumerate_sparse_T.py

from tnfs.simul.polyselect_utils import get_coeffs_from_counter, number_poly, pretty_print_poly_from_coeffs


from tnfs.poly_h.tab_h_2_1 import tab_h_2_1
from tnfs.poly_h.tab_h_2_2_all import tab_h_2_2
from tnfs.poly_h.tab_h_2_3_all import tab_h_2_3
from tnfs.poly_h.tab_h_2_4_all import tab_h_2_4
from tnfs.poly_h.tab_h_2_5_all import tab_h_2_5

from tnfs.poly_h.tab_h_3_1 import tab_h_3_1
from tnfs.poly_h.tab_h_3_2_all import tab_h_3_2
from tnfs.poly_h.tab_h_3_3_all import tab_h_3_3
from tnfs.poly_h.tab_h_3_4_all import tab_h_3_4

from tnfs.poly_h.tab_h_4_1 import tab_h_4_1
from tnfs.poly_h.tab_h_4_2_all import tab_h_4_2
from tnfs.poly_h.tab_h_4_3_all import tab_h_4_3

from tnfs.poly_h.tab_h_5_1 import tab_h_5_1
from tnfs.poly_h.tab_h_5_2_all import tab_h_5_2

from tnfs.poly_h.tab_h_6_1 import tab_h_6_1
from tnfs.poly_h.tab_h_6_2 import tab_h_6_2

from tnfs.poly_h.tab_h_7_1 import tab_h_7_1
from tnfs.poly_h.tab_h_8_1 import tab_h_8_1
from tnfs.poly_h.tab_h_9_1 import tab_h_9_1

from tnfs.poly_h.tab_h_10_1_sparse_3_7 import tab_h_10_1_sparse_3_7
from tnfs.poly_h.tab_h_11_1_sparse_3_7 import tab_h_11_1_sparse_3_7
from tnfs.poly_h.tab_h_12_1_sparse_3_6 import tab_h_12_1_sparse_3_6
from tnfs.poly_h.tab_h_13_1_sparse_3_6 import tab_h_13_1_sparse_3_6
from tnfs.poly_h.tab_h_14_1_sparse_3_6 import tab_h_14_1_sparse_3_6

from tnfs.poly_h.tab_h_16_1_sparse_2_5_min_08 import tab_h_16_1_sparse_2_5
from tnfs.poly_h.tab_h_18_1_3_5_min_inv_zeta_08_prec_10e3 import tab_h_18_1
from tnfs.poly_h.tab_h_20_1_3_zeta_10e4 import tab_h_20_1_3
from tnfs.poly_h.tab_h_20_1_5_zeta_10e4 import tab_h_20_1_5
from tnfs.poly_h.tab_h_24_1_sparse_3_5_prec_10e4_min_09 import tab_h_24_1_sparse_3_5

## global parameter
tab_precomputed_h = {2:tab_h_2_1+tab_h_2_2+tab_h_2_3+tab_h_2_4+tab_h_2_5,
                     3:tab_h_3_1+tab_h_3_2+tab_h_3_3+tab_h_3_4,
                     4:tab_h_4_1+tab_h_4_2[:19]+tab_h_4_3[:38],
                     5:tab_h_5_1+tab_h_5_2[:120],
                     6:tab_h_6_1+tab_h_6_2,
                     7:tab_h_7_1,
                     8:tab_h_8_1[:80],
                     9:tab_h_9_1,
                     10:tab_h_10_1_sparse_3_7,  #  730 polynomials
                     11:tab_h_11_1_sparse_3_7,  # 1362 polynomials
                     12:tab_h_12_1_sparse_3_6,  #  602 polynomials
                     13:tab_h_13_1_sparse_3_6,  #  761 polynomials
                     14:tab_h_14_1_sparse_3_6,  #  979 polynomials
                     16:tab_h_16_1_sparse_2_5,
                     18:tab_h_18_1,
                     20:tab_h_20_1_3+tab_h_20_1_5[:1000],
                     24:tab_h_24_1_sparse_3_5[:100],
}
const_log_2_f = log(float(2.0))
const_2_f = float(2.0)

class Polyselect():

    def __init__(self, E=None, p=None, k=None, deg_h=None):
        if E is None and (p is None or k is None):
            raise ValueError("Polyselect init: please provide class pairing-friendly curve E or prime integer p, extension degree k >= 1")
            
        self._E = E # a pairing-friendly curve E (BN, BLS...)
        if E is None:
            self._k = k
            self._p = p
            self._Fp = GF(p, proof=False)
            self._Fpz = self._Fp['z']; (self._z,) = self._Fpz._first_ngens(1)
        else:
            self._k = E.k()
            self._p = E.p()
            self._Fp = E.Fp()
            self._Fpz, self._z = self._E.Fpz()
        self._pk = self._p**self._k
        self._list_h = {}
        if (deg_h is not None) and deg_h > 1 and self._k % deg_h == 0:
            self._deg_h = deg_h
        else:
            self._deg_h = None
        #self.compute_h(deg_h)
        self._alpha_not_implemented = float(0.5)

    def get_Fp(self):
        return self._Fp
    def get_Fpz(self):
        return self._Fpz, self._z
    
    def compute_h(self, deg_h=None, tab_h=None):
        """ list of possible h for each possible degree """
        if deg_h is None: # any deg_h dividing the embedding degree k
            range_deg_h = [i for i in range(2,self._k+1) if ((self._k % i) == 0)]
        elif deg_h >= 2 and (self._k % deg_h) == 0:
            range_deg_h = [deg_h]
        else:
            raise ValueError("Error: deg_h = {} invalid, it should be 2 <= deg_h <= k and deg_h | k, where k={} is the embedding degree".format(deg_h, self._k))
        if tab_h is not None:
            tab_precomp_h = tab_h
        else:
            tab_precomp_h = tab_precomputed_h
        for deg_h in range_deg_h:
            if not deg_h in tab_precomp_h:
                print("Warning there isn't a precomputed table of irreducible poly h of degree {}".format(deg_h))
                print("        Please compute it before running Polyselect")
                raise ValueError("Error deg_h = {}, unable to find precomputed table of h, please compute it separately with get_list_irr_poly_nonzero_coeffs()".format(deg_h))
                
            print("computing list h of degree {} from {} irreducible polys over Q".format(deg_h, len(tab_precomp_h[deg_h])))
            self._list_h[deg_h] = [tup for tup in tab_precomp_h[deg_h] if self._Fpz(tup[2]).is_irreducible() ]
            print("{} found".format(len(self._list_h[deg_h])))

    def get_h(self, deg_h=None):
        if deg_h is None:
            return self._list_h # all
        if deg_h > 1 and self._k % deg_h == 0 and not deg_h in self._list_h:
            self.compute_h(deg_h)
        if deg_h in self._list_h:
            return self._list_h[deg_h]
        raise ValueError("Error deg_h = {}, unable to find precomputed h or to compute them".format(deg_h))
            
    def print_h(self, deg_h=None):
        if deg_h is None:
            range_deg_h = [i for i in range(2,self._k+1) if ((self._k % i) == 0)]
        else:
            range_deg_h = [deg_h]
        for deg_h in range_deg_h:
            if not deg_h in self._list_h:
                print("Error there isn't a precomputed table of irreducible poly h of degree {}".format(deg_h))
            print("#deg_h = {}".format(deg_h))
            for tup in self._list_h[deg_h]:
                print(tup)

    def Base_m(self, deg_f):
        """
        :param deg_f: degree of polynomial f (highest degree)
        :returns: polynomials f and g s.t. f.degree() is deg_f, g.degree() is 1, (f.resultant(g) % p) = 0, and f ang g are list of coefficients
        
        Computes polynomials f and g with the base-m method.
        f and g are monic and share a common root modulo p.
        """
        p = self._p
        s0 = ceil(RR(p**(1/(deg_f+1))))
        f_coeffs = (p.digits(s0))
        if len(f_coeffs) == deg_f+2:
            f_coeffs[deg_f] += f_coeffs[deg_f+1]*s0
            #f_coeffs[deg_f+1] = 0
            f_coeffs = f_coeffs[:-1] # remove last element which is 0
        g_coeffs = [-s0, 1]
        return f_coeffs, g_coeffs

    
    def JouxLercier(self, deg_f, max_coeff, monic=False):
        """ Joux--Lercier polynomial selection [MathComp2003]
        See http://www.ams.org/journals/mcom/2003-72-242/S0025-5718-02-01482-5
        :param     deg_f: highest degree of polynomials
        :param max_coeff: bound on the coefficients of f: in [-max_coeff, max_coeff]
        :returns        : a pair f,g of polynomials and max_gi = float(log_2(abs(largest coefficient of g)))
        
        f ang g share a root mod p, g has degree deg_f-1, resultant(f,g) % p == 0
        """
        # re-use methods to enumerate f
        counter = 0
        ZZx = PolynomialRing(ZZ, names=('x',)) # needed to check irreducibility
        (x,) = ZZx._first_ngens(1)
        p = self._p
        Fp = self._Fp
        Fpz,z = self._Fpz,self._z
        list_fg_polys = []
        max_counter = number_poly(deg_f, max_coeff, monic=monic)
        while counter < max_counter :
            f_coeffs, next_counter = get_coeffs_from_counter(counter, deg_f, max_coeff, monic=monic)
            counter = next_counter
            if f_coeffs is None:
                continue
            
            # is it a duplicate? see cado-nfs code
            #print("f_coeffs = {}".format(f_coeffs))
            fp = Fpz(f_coeffs)
            rr = fp.roots()
            if len(rr) > 0:
                f = ZZx(f_coeffs)
                if f.is_irreducible():
                    # f is irreducible over Q but has at least one root mod p
                    # computes g
                    if deg_f == 2:# try to detect duplicates
                        all_rows_r = []
                    for ri in rr:
                        r0 = ZZ(ri[0])
                        # row i=0
                        row = [p] + [0 for i in range(1,deg_f)]
                        # row i=1, 2, 3... deg_f-1
                        for i in range(1, deg_f):
                            rowi = [0 for j in range(0,i-1)] + [-r0, 1] + [0 for j in range(i+1, deg_f)]
                            row = row + rowi
                        M = Matrix(deg_f,deg_f, row)
                        #print("M = {}".format(M))
                        M = M.LLL(delta=0.9999,eta=0.50001)
                        #print("LLL(M) = {}".format(M))
                        for i in range(deg_f):
                            g_coeffs = M[i].coefficients() # or .list() ?
                            #print("g_coeffs = {}".format(g_coeffs))
                            c = gcd(g_coeffs)
                            if c > 1:
                                g_coeffs = [ci // c for ci in g_coeffs]
                            if g_coeffs[-1] < 0:
                                g_coeffs = [-ci for ci in g_coeffs]
                            # when deg_g = 1, there are only two coeffs, and sometimes the second row is a rotation of the 1st row
                            if deg_f == 2:
                                row0 = [g_i for g_i in g_coeffs]
                                if len(all_rows_r) == 0:
                                    all_rows_r.append(row0)
                                else:# compare
                                    cond_continue = False
                                    i_r=0
                                    while i_r < len(all_rows_r) and not cond_continue:
                                        row_ = all_rows_r[i_r]
                                        cond_continue = (abs(row_[0]) == abs(row0[1]) and abs(row_[1]) == abs(row0[0])) or \
                                                        (abs(row_[0]) == abs(row0[0]) and abs(row_[1]) == abs(row0[1]))
                                        i_r += 1
                                    if cond_continue:
                                        continue
                                    else:
                                        all_rows_r.append(row0)
                            gi = ZZx(g_coeffs)
                            if gi.is_irreducible(): # irreducibility test over Q
                                g_max_coeff = max([log(float(abs(ci)), const_2_f) for ci in g_coeffs])
                                list_fg_polys.append((f, gi, g_max_coeff))
                                # or with list of coefficients instead of polynomials:
                                # list_fg_polys.append((f_coeffs, g_coeffs, g_max_coeff))
        list_fg_polys.sort(key=itemgetter(2)) # sort according to 3rd value (max coeff of g)
        return list_fg_polys

    def Conjugation(self, max_coeff=4, monic=False):
        """ returns two polynomials, one quadratic and with small coefficients, the second linear,
        so that they share a common root modulo p.
        [BarbulescuGaudryGuillevicMorain15] https://eprint.iacr.org/2016/605
        :param max_coeff: bound on the coefficients of aux1
        :returns: univariate polynomials aux1 and aux0
        aux1 has degree 2
        aux0 = v*x - u where u/v is a rational reconstruction of a root of aux1 mod p
        """
        # 1. enumerate the (non-monic) polynomials of degree 2
        # 2. For those which have a root mod p, computes a rational reconstruction of the roots
        # 3. returns the one with smallest coefficients, or returns the list
        return self.JouxLercier(deg_f=2, max_coeff=max_coeff, monic=monic)

    def GeneralizedJouxLercier(self, deg_f, max_coeff, deg_phi, monic=False, sieving_dim=2, compute_alpha=True, B1_alpha=2000, max_test_polys=500):
        """ Generalized Joux--Lercier method for NFS in extension fields GF(p^n)
        See [Barbulescu PhD] and [BGGM EC'15]
        :param     deg_f: the high degree of polynomials
        :param max_coeff: bound on the coefficients of f: in [-max_coeff, max_coeff]
        :param   deg_phi: the degree of the irreducible factor of f mod p
        :returns        : f,g, max_gi where f, g are (coefficients of) polynomials
        of degree deg_f, deg_f-1,
        resultant(f,g) % p^n == 0 and max_gi = float(log_2(abs(largest coefficient of g)))
        
        For now, works only for univariate polynomials.
        TODO: write a TNFS-compatible version, where f and g have coefficients in Kh
        """
        # re-use methods to enumerate f
        counter = 0
        ZZx = PolynomialRing(ZZ, names=('x',)) # needed to check irreducibility
        (x,) = ZZx._first_ngens(1)
        p = self._p
        Fp = self._Fp
        Fpz,z = self._Fpz, self._z
        list_fg_polys = []
        tests_polys = 0
        
        score_min = None
        max_counter = number_poly(deg_f, max_coeff, monic=monic)
        while counter < max_counter and tests_polys < max_test_polys:
            #print("counter = {}".format(counter))
            f_coeffs, next_counter = get_coeffs_from_counter(counter, deg_f, max_coeff, monic=monic)
            #print("f_coeffs = {}, next_counter = {}".format(f_coeffs, next_counter))
            counter = next_counter
            if f_coeffs is None:
                continue
            
            #print("f = {}".format(f_coeffs))
            fp = Fpz(f_coeffs)
            fa = fp.factor()
            if deg_phi in [fai[0].degree() for fai in fa]:
                f = ZZx(f_coeffs)
                if f.is_irreducible():
                    for fai in [fai[0] for fai in fa if fai[0].degree() == deg_phi]:
                        #print("f = {}, \nirr factor mod p = {}".format(f, fai))
                        # build matrix
                        #fai.coefficients(sparse=false)
                        phi_coeffs = [ZZ(ai) for ai in fai.list()]

                        # buid a deg_f * deg_f matrix -> will correspond to (deg_f-1) polynomials
                        row = [p] + [0 for i in range(1,deg_f)]
                        # row i=1, 2, 3... deg_f-2
                        # if deg_phi = deg_f-1, there is only one row with phi_coeffs
                        for i in range(1, deg_phi):
                            row_i = [0 for j in range(0,i)] + [p] + [0 for j in range(i+1, deg_f)]
                            row = row + row_i
                        row_d = phi_coeffs + [0 for i in range(deg_phi+1, deg_f)]
                        row = row + row_d
                        for i in range(deg_phi+1, deg_f):
                            row_i = [0 for j in range(deg_phi,i)] + phi_coeffs + [0 for i in range(i+1, deg_f)]
                            row = row + row_i
                
                        M = Matrix(deg_f, deg_f, row)
                        M = M.LLL()
                        for i in range(deg_f):
                            g_coeffs = M[i].coefficients()
                            c = gcd(g_coeffs)
                            if c > 1:
                                g_coeffs = [ci // c for ci in g_coeffs]
                            if g_coeffs[-1] < 0:
                                g_coeffs = [-ci for ci in g_coeffs]
                            g = ZZx(g_coeffs)
                            if g.is_irreducible(): # irreducibility test over Q
                                tests_polys += 1
                                g_max_coeff = max([log(float(abs(ci)),const_2_f) for ci in g_coeffs])
                                f_max_coeff = max(f_coeffs)
                                score = (log(float(f_max_coeff),const_2_f) + g_max_coeff)*sieving_dim
                                if score_min is None or score < score_min:
                                    score_min = score
                                    ff = pretty_print_poly_from_coeffs(f_coeffs)
                                    print("f={:32s} |g| = {:.2f}b score={:.4f}".format(ff, g_max_coeff, score))
                                list_fg_polys.append((f, g, f_max_coeff, g_max_coeff, score))
                                #print("f = {}".format(f))
                                #print("g = {}".format(g))
                                #print("max g_i = {}".format(g_max_coeff))
                                # or with list of coefficients instead of polynomials:
                                # list_fg_polys.append((f_coeffs, g_coeffs, g_max_coeff))
        print("found a list of {} pairs of polys".format(len(list_fg_polys)))
        list_fg_polys.sort(key=itemgetter(4)) # sort according to 4rd value (max coeff of g)
        if compute_alpha:
            score_min = None
            list_fg_polys = list_fg_polys[:min(len(list_fg_polys)//2, 250)]
            i = 0
            for item in list_fg_polys:
                f, g, f_max_coeff, g_max_coeff, score = item
                if sieving_dim == 2:
                    alpha_f = alpha2d(f,B1_alpha)
                    alpha_g = alpha2d(g,B1_alpha)
                elif sieving_dim == 3:
                    alpha_f = alpha3d(f,B1_alpha)
                    alpha_g = alpha3d(g,B1_alpha)
                else:
                    alpha_f = self._alpha_not_implemented
                    alpha_g = self._alpha_not_implemented
                    
                score = (log(float(f_max_coeff), const_2_f) + float(g_max_coeff))*sieving_dim + float(alpha_f+alpha_g)/const_log_2_f
                if score_min is None or score < score_min:
                    score_min = score
                    ff = pretty_print_poly_from_coeffs(f.coefficients(sparse=False))
                    print("f={:32s} |g| = {:.2f}b alpha_f={:7.4f} alpha_g={:7.4f} score={:.4f}".format(ff, g_max_coeff, alpha_f, alpha_g, score))
                list_fg_polys[i] = (f, g, f_max_coeff, g_max_coeff, alpha_f, alpha_g, score)
                i += 1
            list_fg_polys.sort(key=itemgetter(6)) # sort according to score with alpha
            # print the best one
            f, g, f_max_coeff, g_max_coeff, alpha_f, alpha_g, score = list_fg_polys[0]
            ff = pretty_print_poly_from_coeffs(f.coefficients(sparse=False))
            print("f={:32s} |g| = {:.2f}b alpha_f={:7.4f} alpha_g={:7.4f} score={:.4f}".format(ff, g_max_coeff, alpha_f, alpha_g, score))
        else:
            # print the best one
            f, g, f_max_coeff, g_max_coeff, score = list_fg_polys[0]
            ff = pretty_print_poly_from_coeffs(f.coefficients(sparse=False))
            print("f={:32s} |g| = {:.2f}b score={:.4f}".format(ff, g_max_coeff, score))
            
        return list_fg_polys

    def Resultant_univariate(self, deg_g, aux0, aux1, test_g_mod_p_irreducible=True):
        """
        For usual NFS.
        :param deg_g:
        :param aux0: a list of coefficients (e.g. f given by JouxLercier or Conjugation)
        :param aux1: a list of coefficients (e.g. g given by JouxLercier or Conjugation)
        """
        Rx = PolynomialRing(ZZ, names=('x',)) # needed to check irreducibility
        (x,) = Rx._first_ngens(1)
        RxU = PolynomialRing(Rx, names=('U',))
        (U,) = RxU._first_ngens(1)
        list_cyclic_g_t = {1: x-U,
                           #2: x**2-U, will never work with MNT6 curves because p=X^2+1
                           #2: x**2+x-U,
                           2: x**2-U*x+1,
                           3: x**3 - U*x**2 - (U+3)*x - 1,
                           4: x**4 -U*x**3 -6*x**2 + U*x + 1,
                           6: x**6 -2*U*x**5 - (5*U+15)*x**4 -20*x**3 +5*U*x**2 + (2*U+6)*x + 1}
        if deg_g in list_cyclic_g_t:
            gt = list_cyclic_g_t[deg_g]
        else:
            # no cyclic solution
            gt = x**deg_g + x + U + 1
            while not gt.is_irreducible() or (not test_g_mod_p_irreducible or not (self._Fpz(Rx(gt.resultant(Rx(aux1)(U))))).is_irreducible()):
                gt += 1
        f = Rx(gt.resultant(Rx(aux0)(U))) # self._E.polynomial_p
        g = Rx(gt.resultant(Rx(aux1)(U))) # self._E.u()-U
        # positive leading coeff
        if f.leading_coefficient() < 0:
            f = -f
        if g.leading_coefficient() < 0:
            g = -g
        # content-free
        f = f // f.content()
        g = g // g.content()
        max_fi = max([abs(fi) for fi in f.list()])
        max_gi = max([abs(gi) for gi in g.list()])

        if (f.resultant(g)).mod(self._pk) != 0:
            print("error in Resultant_univariate() res(f, g) mod p^k != 0")
            return
        # test if g is irreducible mod p --> but for Sarkar--Singh, it will not be the case.
        # for Sarkar Singh, it should give an irreducible factor of degree k
        if test_g_mod_p_irreducible:
            gp = self._Fpz(g)
            if not gp.is_irreducible():
                print("g = res(aux0={}, gt={}) is not irreducible mod p".format(aux0, gt))
                return
        if not f.is_irreducible():
            return
        if not g.is_irreducible():
            return
        # the automorphism order is deg_g
        return f, g, max_fi, max_gi, deg_g
        
    def Resultant(self, deg_g, aux0, aux1, h=None, with_y=False, mult_y=True, MNT6=False):
        """
        :param    h: a list of coefficients
        :param aux0: a list of coefficients (e.g. f given by JouxLercier or Conjugation)
        :param aux1: a list of coefficients (e.g. g given by JouxLercier or Conjugation)
        """
        ZZy = PolynomialRing(ZZ, names=('y',)) # needed to check irreducibility
        (y,) = ZZy._first_ngens(1)
        Rxy = PolynomialRing(ZZy, names=('x',))
        (x,) = Rxy._first_ngens(1)
        ZZyxU = PolynomialRing(Rxy, names=('U',))
        (U,) = ZZyxU._first_ngens(1)
        #print("deg_g = {}, aux0 = {}, aux1 = {}, h = {}".format(deg_g, aux0, aux1, h))
        list_cyclic_g_t = {1: x-U,
                           #2: x**2-U, will never work with MNT6 curves because p=X^2+1
                           #2: x**2+x-U,
                           2: x**2-U*x+1,
                           3: x**3 - U*x**2 - (U+3)*x - 1,
                           4: x**4 -U*x**3 -6*x**2 + U*x + 1,
                           6: x**6 -2*U*x**5 - (5*U+15)*x**4 -20*x**3 +5*U*x**2 + (2*U+6)*x + 1}
        if MNT6:
            list_cyclic_g_t[2] = x**2+x-U
        list_cyclic_g_ty_ = {1: x-U-y,
                            2: x**2-y-U,
                            3: x**3 - (U+y)*x**2 - (U+y+3)*x - 1,
                            4: x**4 -(U+y)*x**3 -6*x**2 + (U+y)*x + 1,
                            6: x**6 -2*(U+y)*x**5 - (5*(U+y)+15)*x**4 -20*x**3 +5*(U+y)*x**2 + (2*(U+y)+6)*x + 1}
        list_cyclic_g_ty = {1: x-U*y,
                            2: x**2-U*y,
                            3: x**3 - (U*y)*x**2 - (U*y+3)*x - 1,
                            4: x**4 -(U*y)*x**3 -6*x**2 + (U*y)*x + 1,
                            6: x**6 -2*(U*y)*x**5 - (5*(U*y)+15)*x**4 -20*x**3 +5*(U*y)*x**2 + (2*(U*y)+6)*x + 1}
        if not with_y:
            if deg_g in list_cyclic_g_t:
                gtU = list_cyclic_g_t[deg_g]
                g = Rxy(gtU.resultant(ZZy(aux1)(U))) # self._E.u()-U
                gp = self._Fpz(g)
                uu = U
                gt = gtU
                while not gt.is_irreducible() or not gp.is_irreducible():
                    uu = uu+1
                    gt = gtU(uu)
                    g = Rxy(gt.resultant(ZZy(aux1)(U))) # self._E.u()-U
                    gp = self._Fpz(g)
            else:
                # no cyclic solution
                gt = x**deg_g + x + U + 1
                g = Rxy(gt.resultant(ZZy(aux1)(U))) # self._E.u()-U
                gp = self._Fpz(g)
                while not gt.is_irreducible() or not gp.is_irreducible():
                    gt += 1
                    g = Rxy(gt.resultant(ZZy(aux1)(U))) # self._E.u()-U
                    gp = self._Fpz(g)
            f = Rxy(gt.resultant(ZZy(aux0)(U))) # self._E.polynomial_p
            g = Rxy(gt.resultant(ZZy(aux1)(U))) # self._E.u()-U
        else:
            if mult_y and deg_g in list_cyclic_g_ty:
                gt = list_cyclic_g_ty[deg_g]
            elif (not mult_y) and deg_g in list_cyclic_g_ty_:
                gt = list_cyclic_g_ty_[deg_g]
            else:
                # no cyclic solution
                gt = x**deg_g + x + U + y + 1
                while not gt.is_irreducible():
                    gt += 1
            
            f = Rxy(gt.resultant(ZZy(aux0)(U))) # self._E.polynomial_p
            g = Rxy(gt.resultant(ZZy(aux1)(U))) # self._E.u()-U
            # now reduce mod h(y) the coefficients in y of f, g
            if with_y and h is not None:
                hy = ZZy(h)
                f = Rxy([fi % hy for fi in f.list()])
                g = Rxy([gi % hy for gi in g.list()])
        # positive leading coeff
        if (f.list()[-1]).list()[-1] < 0:
            f = -f
        if (g.list()[-1]).list()[-1] < 0:
            g = -g
        # content-free
        gcd_fi = gcd([gcd([fij for fij in fi.list()]) for fi in f.list() if len(fi.list()) > 0])
        if gcd_fi > 1:
            #f = Rxy((1/Integer(gcd_fi))*f)
            f = f // Integer(gcd_fi)
        gcd_gi = gcd([gcd([gij for gij in gi.list()]) for gi in g.list() if len(gi.list()) > 0])
        if gcd_gi > 1:
            g = g // Integer(gcd_gi)
        max_fij = max([max([abs(fij) for fij in fi.list()]) for fi in f.list() if len(fi.list()) > 0])
        max_gij = max([max([abs(gij) for gij in gi.list()]) for gi in g.list() if len(gi.list()) > 0])

        Res_fgh = ZZ((f.resultant(g)).resultant(h))
        if h is not None and Res_fgh.mod(self._pk) != 0:
            kk = self._k-1
            while kk >= 1 and Res_fgh.mod(self._p**kk) != 0:
                kk -= 1
            if kk > 0:
                print("Polyselect error Res(Res(f,g),h) mod p^k != 0, this is 0 mod p^{}".format(kk))
            else:
                print("Polyselect error Res(Res(f,g),h) mod p^k != 0")
                print("f={}; g={}; h={}".format(f,g,h))
            return
        if h is None and ZZ((f.resultant(g)).mod(self._p) != 0):
            print("Polyselect error Res(f,g) mod p != 0")
            return
        # test if g is irreducible mod p
        #print("g = {}".format(g))
        #print("type(g) = {}".format(type(g)))
        # if g is univariate
        if g.degree() > 1 and not with_y:
            gp = self._Fpz(g)
            if not gp.is_irreducible():
                print("Error Resultant: g ={} is not irreducible mod p: {}".format(g, gp.factor()))
                return
        if not f.is_irreducible():
            print("Polyselect error f is not irreducible")
            return
        if not g.is_irreducible():
            print("Polyselect error g is not irreducible")
            return
        # the automorphism order is deg_g
        return f, g, max_fij, max_gij, deg_g

    def NFS_JLSV1(self, deg_g, sieving_dim=2, compute_alpha=True, verbose=False, B0_alpha=800):
        """
        compute two auxiliary polynomials x-m and v*x-u s.t. m ~ sqrt(p), u/v = m mod p and u,v ~ sqrt(p)
        INPUT:
        - `deg_g`: degree of the polynomials, irreducible mod p
        - `sieving_dim`: 2 or 3, others not implemented
        - `compute_alpha`: boolean, False: set to default value
        - `verbose`: print steps
        - `B0_alpha`: upper bound on the small primes in computing alpha

        Polynomial selection method from [Joux Lercier Smart Vecauteren CRYPTO'2006]
        https://www.iacr.org/archive/crypto2006/41170323/41170323.pdf
        """
        p = self._p
        m0 = ceil(RR(p).sqrt())
        M = Matrix(2,2, [p,0,m0,1])
        M = M.LLL(delta=0.9999,eta=0.50001)
        aux1 = [m0,1]
        result = []
        if (M[0][0] == m0 and M[0][1] == 1) or (M[0][0] == -m0 and M[0][1] == -1) \
           or (M[0][1] == m0 and M[0][0] == 1) or (M[0][1] == -m0 and M[0][0] == -1):
            aux0 = [M[1][0], M[1][1]]
        else:
            aux0 = [M[0][0], M[0][1]]
        ans = self.Resultant_univariate(deg_g, aux0, aux1)
        if ans is None:
            return None
        f, g, max_fi, max_gi, aut_g = ans
        if compute_alpha:
            if verbose == 2:
                print("computing alpha_f")
            if sieving_dim == 2:
                alpha_f = float(alpha2d(f, B0_alpha))
            elif sieving_dim == 3:
                alpha_f = float(alpha3d(f, B0_alpha))
            else:
                alpha_f = self._alpha_not_implemented
                print("alpha not implemented for sieving_dim = {}, setting to default value {}".format(sieving_dim, self._alpha_not_implemented))
            if verbose == 2:
                print("computing alpha_g")
            if sieving_dim == 2:
                alpha_g = float(alpha2d(g, B0_alpha))
            elif sieving_dim == 3:
                alpha_g = float(alpha3d(g, B0_alpha))
            else:
                alpha_g = self._alpha_not_implemented
            sum_alpha = alpha_f+alpha_g
        else:
            alpha_f = float(0.0); alpha_g = float(0.0); sum_alpha = float(0.0)
        score = (log(float(max_fi), const_2_f)+log(float(max_gi), const_2_f))*sieving_dim + (alpha_f+alpha_g)/const_log_2_f
        result.append((f, g, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, score))
        return result[0]

    def NFS_Conjugation(self, deg_g, sieving_dim=2, max_coeff=4, monic=False, compute_alpha=True, verbose=False, number_results=5, B0_alpha=800, B1_alpha=2000):
        """Conjugation method for GF(p^k), deg(f) = 2*deg(g), and deg_g = k
        INPUT:
        - `deg_g`: degree of the extension. g will be irreducible mod p
        - `sieving_dim`: sieving dimension for computing alpha, for now 2 and 3 are possible
        - `max_coeff`: max absolute value of the coeffs of the auxiliary polynomial of degree 2
        - `monic`: if f has to be monic
        - `compute_alpha`: True/False
        - `verbose`: for printing intermediate data
        - `number_results`: output the list of best possibilities
        - `B0_alpha`: bound for computing alpha the first time
        - `B1_alpha`: bound for computing alpha with more accuracy the second time

        Polynomial selection from [Barbulescu-Gaudry-Guillevic-Morain EUROCRYPT'15](https://eprint.iacr.org/2016/605)
        """
        list_aux = self.Conjugation(max_coeff=max_coeff, monic=monic)
        result = []
        print("{} auxiliary polynomials to test".format(len(list_aux)))
        for i in range(len(list_aux)):
            aux0, aux1, max_lg_gi = list_aux[i] # aux0 has degree 2, aux1 has degree 1
            ans = self.Resultant_univariate(deg_g, aux0, aux1)
            if ans is None:
                continue
            f, g, max_fi, max_gi, aut_g = ans
            if compute_alpha:
                if verbose == 2:
                    print("computing alpha_f")
                if sieving_dim == 2:
                    alpha_f = float(alpha2d(f,B0_alpha))
                elif sieving_dim == 3:
                    alpha_f = float(alpha3d(f,B0_alpha))
                else:
                    alpha_f = self._alpha_not_implemented
                    print("alpha not implemented for sieving_dim = {}, setting to default value {}".format(sieving_dim, self._alpha_not_implemented))
                if verbose == 2:
                    print("computing alpha_g")
                if sieving_dim == 2:
                    alpha_g = float(alpha2d(g,B0_alpha))
                elif sieving_dim == 3:
                    alpha_g = float(alpha3d(g,B0_alpha))
                else:
                    alpha_g = self._alpha_not_implemented
                sum_alpha = alpha_f+alpha_g
            else:
                alpha_f = 0.0; alpha_g = 0.0; sum_alpha = 0.0
            # score: (log(|f|^sieving_dim) + alpha_f + log(|g|^sieving_dim) + alpha_g)/log(2)
            score = (log(float(max_fi), const_2_f)+log(float(max_gi), const_2_f))*sieving_dim + (alpha_f+alpha_g)/const_log_2_f
            result.append((f, g, aux0, aux1, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, float(score)))
            if verbose==2:
                print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} aux0={} aux1={} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, float(log(max_fi, const_2_f)), float(log(max_gi, const_2_f)), float(score), aux0.list(), aux1.list(), f, g))
            elif verbose==1:
                print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, float(log(max_fi, const_2_f)), float(log(max_gi, const_2_f)), float(score)))
            result.sort(key=itemgetter(10))
        if compute_alpha:
            result = result[:10] # select the best 5 selected polynomials
            print("recompute alpha for the best pairs (f,g)")
            for i in range(min(10,len(result))):
                f, g, aux0, aux1, max_fi, max_gi, aut_g = result[i][:7]
                if sieving_dim == 2:
                    alpha_f = float(alpha2d(f,B1_alpha))
                    alpha_g = float(alpha2d(g,B1_alpha))
                elif sieving_dim == 3:
                    alpha_f = float(alpha3d(f,B1_alpha))
                    alpha_g = float(alpha3d(g,B1_alpha))
                else:
                    alpha_f = self._alpha_not_implemented
                    alpha_g = self._alpha_not_implemented
                sum_alpha = alpha_f+alpha_g
                score = (log(float(max_fi), const_2_f)+log(float(max_gi), const_2_f))*sieving_dim + (alpha_f+alpha_g)/const_log_2_f
                print("{}: {: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} aux0={} aux1={} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, float(log(max_fi, const_2_f)), float(log(max_gi, const_2_f)), float(score), aux0, aux1, f, g))
                result[i] = (f, g, aux0, aux1, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, score)
            result.sort(key=itemgetter(10))
        print("final list:")
        i=0
        for vi in result:
            f, g, aux0, aux1, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, score = vi
            print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} aux0={} aux1={} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, float(log(max_fi, const_2_f)), float(log(max_gi, const_2_f)), float(score), aux0, aux1, f, g))
            i+=1

        return result # will return a list of at most 10 pairs of polynomials

    def NFS_Special(self, deg_g, poly_p=None, u=None):
        """ uses the special form of parameter p. This is like Joux-Pierrot but with resultants."""
        if poly_p is None:
            aux0 = self._E.polynomial_p # will fail if E is None
        else:
            aux0 = poly_p
        if u is None:
            aux1 = [-self._E.u(), 1] # U-u
        else:
            aux1 = [-u, 1] # U-u
        return self.Resultant_univariate(deg_g, aux0, aux1)
        
    def TNFS_Special(self, deg_g, h=None, poly_p=None, u=None, with_y=False, MNT6=False):
        """ uses the special form of parameter p. This is like Joux-Pierrot but with resultants."""
        if poly_p is None:
            aux0 = self._E.polynomial_p
        else:
            aux0 = poly_p
        if u is None:
            aux1 = [-self._E.u(), 1] # U-u
        else:
            aux1 = [-u, 1] # U-u
        #print("deg_g = {}, aux0 = {}, aux1 = {}, h = {}".format(deg_g, aux0, aux1, h))
        return self.Resultant(deg_g, aux0, aux1, h=h, with_y=with_y, mult_y=(not MNT6))
    
    def TNFS_Conjugation(self, deg_g, h, with_y=False, max_coeff=4, monic=False,compute_alpha=False,alpha_test_principal=True,verbose=False,number_results=1,B0_alpha=800,B1_alpha=1200):
        """Conjugation method. """
        list_aux = self.Conjugation(max_coeff=max_coeff, monic=monic)
        result = []
        deg_h = h.degree()
        print("{} auxiliary polynomials to test".format(len(list_aux)))
        for i in range(len(list_aux)):
            aux0, aux1, max_gi = list_aux[i]
            ans = self.Resultant(deg_g, aux0, aux1, h=h, with_y=with_y)
            if ans is None:
                continue
            f, g, max_fij, max_gij, aut_g = ans
            if compute_alpha:
                if verbose == 2:
                    print("computing alpha_f")
                alpha_f = float(alpha_TNFS_2d(f,h,B0_alpha,test_principal=alpha_test_principal))
                if verbose == 2:
                    print("computing alpha_g")
                alpha_g = float(alpha_TNFS_2d(g,h,B0_alpha,test_principal=alpha_test_principal))
                sum_alpha = alpha_f+alpha_g
            else:
                alpha_f = 0.0; alpha_g = 0.0; sum_alpha = 0.0
            #score = RR(log(max_fij,2.0))*(deg_h-1) + RR(alpha_f/log(2.0)) + RR(log(max_gij,2.0))*(deg_h-1) + RR(alpha_g/log(2.0))
            # score: (log(|f|^deg_h) + alpha_f + log(|g|^deg_h) + alpha_g)/log(2.0)
            score = log(float(max_fij), const_2_f)*deg_h + float(alpha_f)/const_log_2_f + log(float(max_gij), const_2_f)*deg_h + float(alpha_g)/const_log_2_f
            result.append((f, g, max_fij, max_gij, aut_g, alpha_f, alpha_g, sum_alpha, score))
            if verbose==2:
                print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} aux0{}aux1{} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, log(float(max_fij), const_2_f), log(float(max_gij), const_2_f), score, aux0.list(), aux1.list(), f, g))
            elif verbose==1:
                print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, log(float(max_fij), const_2_f), log(float(max_gij), const_2_f), score))
            result.sort(key=itemgetter(8))
        if compute_alpha:
            result = result[:5]
            print("recompute alpha for the best pairs (f,g)")
            for i in range(min(5,len(result))):
                f, g, max_fij, max_gij, aut_g = result[i][:5]
                alpha_f = float(alpha_TNFS_2d(f,h,1000,test_principal=alpha_test_principal))
                alpha_g = float(alpha_TNFS_2d(g,h,1000,test_principal=alpha_test_principal))
                sum_alpha = alpha_f+alpha_g
                #score = RR(log(max_fij,2.0))*(deg_h-1) + RR(alpha_f*log(2.0)) + RR(log(max_gij,2.0))*(deg_h-1) + RR(alpha_g*log(2.0))
                score = log(float(max_fij), const_2_f)*deg_h + float(alpha_f)/const_log_2_f + log(float(max_gij), const_2_f)*deg_h + float(alpha_g)/const_log_2_f
                #print("{}: {: .4f} {: .4f} {: .4f}, score {:.4f}".format(i, float(alpha_f*log(2)), float(alpha_g*log(2)), float(sum_alpha*log(2)), float(score)))
                print("{}: {: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, log(float(max_fij), const_2_f), log(float(max_gij), const_2_f), score, f, g))
                result[i] = (f, g, max_fij, max_gij, aut_g, alpha_f, alpha_g, sum_alpha, score)
            result.sort(key=itemgetter(8))
        print("final list:")
        i=0
        for vi in result:
            f, g, max_fij, max_gij, aut_g, alpha_f, alpha_g, sum_alpha, score = vi
            print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, log(float(max_fij), const_2_f), log(float(max_gij), const_2_f), score, f, g))
            i+=1
        if number_results<=1:
            return result[0]
        return result[:number_results]

    def TNFS_SarkarSingh_JL(self, deg_aux1, deg_g, h=None, with_y=False, max_coeff=4, monic=False):
        """ First variant of Sarkar Singh. The second variant is with GJL. """
        list_aux = self.JouxLercier(deg_f=deg_aux1, max_coeff=max_coeff, monic=monic)
        aux0, aux1, max_gi = list_aux[0]
        return self.Resultant(deg_g, aux0, aux1, h=h, with_y=with_y)

    def TNFS_GJL(self, deg_phi, deg_f, h, with_y=False, max_coeff=4, monic=False, compute_alpha=True, alpha_test_principal=True, B0_alpha=800, B1_alpha=1200, number_results=1, verbose=2):
        list_fg_polys = self.GeneralizedJouxLercier(deg_f, max_coeff, deg_phi, monic, sieving_dim=2, compute_alpha=False, max_test_polys=500)
        # do not compute alpha(f), alpha(g) as it does not mean anything -> we are going to compute alpha_tnfs_2d(f, h) and alpha_tnfs_2d(g, h)
        # if compute_alpha=True: the output looks like
        # f, g, f_max_coeff, g_max_coeff, alpha_f, alpha_g, score
        # otherwise it is:
        # f, g, f_max_coeff, g_max_coeff, score
        result = []
        deg_h = h.degree()
        aut_g = 1 # no automorphim with GJL
        i=0
        threshold_max = 30
        if len(list_fg_polys) > threshold_max:# to a first selection of those with smallest coefficients before computing alpha
            list_fg_polys.sort(key=itemgetter(3))
            list_fg_polys = list_fg_polys[:threshold_max]
        for f, g, max_fi, max_gi, score in list_fg_polys:
            # max_gi is already the log in basis 2 of the absolute value of coefficients
            if compute_alpha:
                # compute alpha with h as base field, under f and g. Assume that deg(h) coprime to deg(phi)
                alpha_f = float(alpha_TNFS_2d(f, h, B0_alpha, test_principal=alpha_test_principal))
                alpha_g = float(alpha_TNFS_2d(g, h, B0_alpha, test_principal=alpha_test_principal))
                sum_alpha = alpha_f+alpha_g
            else:
                alpha_f = 0.0; alpha_g = 0.0; sum_alpha = 0.0
            score = log(float(max_fi), const_2_f)*deg_h + float(alpha_f)/const_log_2_f + max_gi*deg_h + float(alpha_g)/const_log_2_f
            result.append((f, g, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, float(score)))
            if verbose==2:
                print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, log(float(max_fi), const_2_f), max_gi, score, f, g))
            elif verbose==1:
                print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} f={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, log(float(max_fi), const_2_f), max_gi, score, f))
            result.sort(key=itemgetter(8))
            i += 1
        if compute_alpha:
            result = result[:5]
            print("recompute alpha for the best pairs (f,g)")
            for i in range(min(5,len(result))):
                f, g, max_fi, max_gi, aut_g = result[i][:5]
                alpha_f = float(alpha_TNFS_2d(f, h, B1_alpha, test_principal=alpha_test_principal))
                alpha_g = float(alpha_TNFS_2d(g, h, B1_alpha, test_principal=alpha_test_principal))
                sum_alpha = alpha_f+alpha_g
                score = log(float(max_fi), const_2_f)*deg_h + float(alpha_f)/const_log_2_f + max_gi*deg_h + float(alpha_g)/const_log_2_f
                print("{}: {: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, log(float(max_fi), const_2_f), max_gi, score, f, g))
                result[i] = (f, g, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, score)
            result.sort(key=itemgetter(8))# sort according to the score
        print("final list:")
        i=0
        for vi in result:
            f, g, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, score = vi
            print("{:2d}:{: .4f} {: .4f} {: .4f}, fi:{:.2f}, gi:{:.2f} score {:.4f} f={} g={}".format(i, float(alpha_f)/const_log_2_f, float(alpha_g)/const_log_2_f, float(sum_alpha)/const_log_2_f, log(float(max_fi), const_2_f), max_gi, score, f, g))
            i+=1
        if number_results<=1:
            return result[0]
        return result[:number_results]

    def NFS_SarkarSingh_JL(self, deg_phi_top, deg_aux1, deg_phi_base, max_coeff=2, monic=False, compute_alpha=True, sieving_dim=2, max_test_poly_aux1=500, B1_alpha=2000, verbose=False):
        """
        Sarkar-Singh polynomial selection https://eprint.iacr.org/2015/944
        INPUT:
        - `deg_phi_top`: degree of the GCD of aux1 and aux0
        - `deg_aux1`: degree of aux1 >= deg_phi_top+1 where f=Resultant(aux1, phi_base), g=Resultant(aux0, phi_base)
        - `deg_phi_base`: degree of the irreducible polynomial phi_base "at the bottom"
        - `max_coeff`: bound on the absolute value of the coefficients of aux1
        - `monic`: if aux1 has to be monic
        - `compute_alpha`
        phi2 is irreducible of degree deg_phi_top,
        phi1 is irreducible of degree deg_phi_base and deg_phi_top*deg_phi_base = k
        deg_aux1 >= deg_phi_top+1
        aux1,aux0 = GJL(deg_f=deg_aux1, deg_phi=deg_phi_top) and deg(aux0) = deg(aux1)-1 by design
        phi1 is obtained as the irreducible GCD polynomial of aux1 and aux0
        return f=Resultant(aux1,phi2), g=Resultant(aux0,phi2)
        """
        if deg_phi_top == 1:
            if verbose:
                print("deg_phi_top = 1, Joux-Lercier technique with deg_f = {}, max_coeff={}, monic = {}".format(deg_aux1, max_coeff, monic))
            list_aux = self.JouxLercier(deg_f=deg_aux1, max_coeff=max_coeff, monic=monic)
            #aux0, aux1, max_gi = list_aux[0]
        else:
            if verbose:
                print("deg_phi_top = {}, generalized Joux-Lercier technique with deg_f = {}, max_coeff={}, monic = {}".format(deg_phi_top, deg_aux1, max_coeff, monic))
            list_aux = self.GeneralizedJouxLercier(deg_f=deg_aux1, max_coeff=max_coeff, deg_phi=deg_phi_top, monic=monic, sieving_dim=sieving_dim, compute_alpha=False, max_test_polys=max_test_poly_aux1)
            # it would be meaningless to compute alpha there
            #aux0, aux1, max_aux0i, max_aux1i, score = list_aux[0]
        if verbose:
            print("obtained a list of {} auxiliary polynomials".format(len(list_aux)))
        if len(list_aux) == 0:
            raise ValueError("function NFS_SarkarSingh_JL(): Empty list of polynomials. Try increasing max_coeff={} to a higher value.".format(max_coeff))
        # 1. compute SSingh polynomials and sort according to a score made of the coefficient size of the polynomials
        score_min = None
        list_fg_polys = [None for i in range(len(list_aux))]
        i = 0
        for item in list_aux:
            aux0, aux1 = item[:2]
            res = self.Resultant_univariate(deg_phi_base, aux0, aux1, test_g_mod_p_irreducible=False)
            if res is None:
                continue
            f, g, max_fi, max_gi, aut_fg = res
            max_fi = max([abs(fi) for fi in f.coefficients()])
            max_gi = max([log(float(abs(gi)), const_2_f) for gi in g.coefficients()])
            
            score = (log(float(max_fi), const_2_f) + max_gi)*sieving_dim
            if score_min is None or score < score_min:
                score_min = score
                ff = pretty_print_poly_from_coeffs(f.coefficients(sparse=False))
                if verbose:
                    print("f={:40s} |g| = {:.2f}b score={:.4f}".format(ff, max_gi, score))
            list_fg_polys[i] = (f, g, max_fi, max_gi, aut_fg, score)
            i += 1
        list_fg_polys = list_fg_polys[:i] # to suppress the None at the end
        list_fg_polys.sort(key=itemgetter(5)) # sort according to score based on coefficient size

        # recompute the score with alpha
        if compute_alpha:
            score_min = None
            i = 0
            for item in list_fg_polys:
                f, g, max_fi, max_gi, aut_fg, score = item
                if sieving_dim == 2:
                    alpha_f = alpha2d(f, B1_alpha)
                    alpha_g = alpha2d(g, B1_alpha)
                elif sieving_dim == 3:
                    alpha_f = alpha3d(f, B1_alpha)
                    alpha_g = alpha3d(g, B1_alpha)
                else:
                    alpha_f = self._alpha_not_implemented
                    alpha_g = self._alpha_not_implemented
                score = (log(float(max_fi), const_2_f) + max_gi)*sieving_dim + float(alpha_f+alpha_g)/const_log_2_f
                if score_min is None or score < score_min:
                    score_min = score
                    ff = pretty_print_poly_from_coeffs(f.coefficients(sparse=False))
                    if verbose:
                        print("f={:40s} |g| = {:.2f}b alpha_f={:7.4f} alpha_g={:7.4f} score={:.4f}".format(ff, max_gi, alpha_f, alpha_g, score))
                list_fg_polys[i] = (f, g, max_fi, max_gi, aut_fg, alpha_f, alpha_g, score)
                i += 1
            list_fg_polys.sort(key=itemgetter(7)) # sort according to score with alpha
        return list_fg_polys[0]
