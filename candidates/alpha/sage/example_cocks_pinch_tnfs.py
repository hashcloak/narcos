from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.rings.number_field.number_field import NumberField
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.structure.proof.all import arithmetic
from sage.arith.misc import GCD, gcd
from sage.arith.functions import lcm
from sage.functions.log import log

# imports from TNFS package
import tnfs
from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d

from tnfs.simul.polyselect import Polyselect
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.simul.simulation_nfs import Simulation_NFS

from tnfs.simul.polyselect_utils import automorphism_factor
from tnfs.simul.polyselect_utils import pretty_print_coeffs_from_coeffs, pretty_print_poly_from_coeffs
from tnfs.gen.generate_curve_utils import str_binary_form, str_py_binary_form

# polynomial rings
ZZy = ZZ['y']; (y,) = ZZy._first_ngens(1)
Rxy = PolynomialRing(ZZy, names=('x',)); (x,) = Rxy._first_ngens(1)

def test_finite_field_tnfs(q, r, k, cost, samples=100000, special=False, conj=False, sarkarsingh=False, jouxlercier=False, all_deg_h=None, qx=None, u=None, max_coeff=2, deg_f=None, deg_phi_base=None, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False):
    """
    run TNFS for the field GF(q^k) with r a prime divisor of the cyclotomic subgroup
    INPUT:
    - `q`: prime integer, characteristic of the target field
    - `r`: prime integer, order of subgroup of GF(q^k)
    - `k`: extension degree
    - `cost`: expected cost of DL
    - `samples`: number of samples in the Monte Carlo simulation
    - `special`: Special-TNFS, requires q of special form given with qx, u, such that q = qx(u)
    - `conj`: Conjugation-TNFS
    - `sarkarsingh`: Sarkar-Singh-TNFS, requires k composite
    - `jouxlercier`: Joux-Lercier-TNFS
    - `all_deg_h`: a list of targeted degrees of subfield in TNFS. If k is prime, it is deg_h=k by default.
    - `qx`: polynomial for the special form of q if special is chosen
    - `u`: seed such that q = qx(u) for the special STNFS method
    - `max_coeff`: for Conj or Sarkar-Singh or Joux-Lercier methods
    - `deg_f`: for JL and GJL (jouxlercier) deg_f > k/deg_h, or Sarkar-Singh: deg_f >= (deg_phi_top+1)*deg_phi_base
    - `deg_phi_base`: for Sarkar-Singh: deg_phi_top*deg_phi_base = k/deg_h, gcd(aux0, aux1) = phi_top, f=Resultant(aux1, phi_base), g = Resultant(aux0, phi_base)
    """
    poly_init = Polyselect(p=q, k=k)
    # 1. TNFS
    if all_deg_h is None:# test all possible degrees of subfields, so that it divides k
        # if k is prime, the only possible degree is k. deg = 1 corresponds to NFS (below)
        all_deg_h = [i for i in range(k,1,-1) if (k % i) == 0]
    for deg_h in all_deg_h:
        deg_g = k // deg_h
        poly_init.compute_h(deg_h=deg_h)
        list_h = poly_init.get_h()[deg_h]
        print("polynomials h")
        for item in list_h:
            inv_zeta_Kh, w, hc = item
            hc_str = pretty_print_coeffs_from_coeffs(hc)
            h_str = pretty_print_poly_from_coeffs(hc)
            print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
        print("")
        for item in list_h[:5]:
            inv_zeta_Kh, w, hc = item
            hc_str = pretty_print_coeffs_from_coeffs(hc)
            h_str = pretty_print_poly_from_coeffs(hc)
            print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
            h = ZZy(hc)

            # actually with D=3, one can have aux = x^2+x+1, and with D=1, aux=x^2+1, but that's all.
            if special:
                print("special-TNFS")
                res_poly = poly_init.TNFS_Special(deg_g=deg_g, h=h, poly_p=qx, u=u, with_y=(gcd(deg_g, deg_h)==1))
            elif conj:
                print("conj-TNFS")
                res_poly = poly_init.TNFS_Conjugation(deg_g=deg_g, h=h, with_y=(gcd(deg_g, deg_h) > 1), max_coeff=max_coeff, monic=False, compute_alpha=compute_alpha, alpha_test_principal=alpha_test_principal,verbose=2, number_results=1, B0_alpha=B0_alpha, B1_alpha=B1_alpha)
            elif sarkarsingh:
                print("Sarkar-Singh-TNFS")
                res_poly = poly_init.TNFS_SarkarSingh_JL(deg_phi_top=k//(deg_h*deg_phi_base), deg_aux1=deg_f//deg_phi_base, deg_g=deg_g, h=h, with_y=(gcd(deg_g,deg_h) > 1), max_coeff=max_coeff, B0_alpha=B0_alpha, B1_alpha=B1_alpha)
            elif jouxlercier:
                print("Joux-Lercier-TNFS")
                res_poly = poly_init.TNFS_GJL(deg_phi=k//deg_h, deg_f=deg_f, h=h, with_y=False, max_coeff=max_coeff, monic=False, compute_alpha=compute_alpha, alpha_test_principal=alpha_test_principal, B0_alpha=B0_alpha, B1_alpha=B1_alpha, number_results=1)

            if res_poly == None:
                raise ValueError("Error in Polynomial selection")
            f, g, max_fij, max_gij, aut_fg = res_poly[:5]
            Kh = NumberField(h, names=('ah',)); (ah,) = Kh._first_ngens(1)
            print("h = {} # {}".format(h, h.list()))
            print("inv_zeta_Kh, w = {:.6f},{:2d}".format(inv_zeta_Kh, w))
            assert (ZZ(h.resultant(f.resultant(g))) % q**k) == 0
            if f.parent() is not Rxy:
                f = f(x)
                g = g(x)

            if (gcd(deg_g,deg_h) == 1):
                aut_h = automorphism_factor(hc)
            else:
                aut_h = 1
            aut = aut_h*aut_fg
            print("f = {}".format(f))
            print("g = {}".format(g))
            print("max_fi = {:6d}, log_2 max_gi = {:6.2f}".format(max_fij, float(log(max_gij,2))))
            print("aut = {}".format(aut))

            # computing alpha, this takes at least a few seconds

            if f.degree() <= 12:# and gcd(deg_g,deg_h) == 1 and gcd(f.degree(),deg_h) == 1:
                alpha_f = float(alpha_TNFS_2d(f,h,1000))
                alpha_g = float(alpha_TNFS_2d(g,h,1000))
            else:
                alpha_f = 0.5
                alpha_g = 0.5
            sum_alpha = alpha_f+alpha_g
            print("alpha_f = {:.4f} alpha_g = {:.4f} sum_alpha = {:.4f}".format(alpha_f,alpha_g,sum_alpha))
            print("    ({:.8f},{:2d}, {:40s}, {:.4f}, {:.4f}, {:.4f}),".format(float(inv_zeta_Kh),int(w),str(h),float(alpha_f),float(alpha_g),float(sum_alpha)))

            #initialisation of data
            Fp = poly_init.get_Fp()
            Fpz = poly_init.get_Fpz()
            simul = Simulation_TNFS(q,r,Fp,Fpz,h,f,g,Rxy,cost,aut,inv_zeta_Kh,count_sieving=True,alpha_f=alpha_f,alpha_g=alpha_g)

            simul.print_params()
            simul.simulation(samples=samples) #takes few seconds for 10^4, mins for 10^5, up to 20 min for 10^6
            simul.print_results()
            print("#::::::::::::::")
            # if there is not enough relations of there are too many relations, re-run with the same polynomials but with a higher/smaller cost
            simul.adjust_cost(samples=samples)
            print("############")

def test_cocks_pinch_k5_gmt():
    """
    Run the TNFS simulation tool on the parameters of a Cocks-Pinch k=5 curve
    Pairing-friendly elliptic curve from https://hal.inria.fr/hal-02305051v2
    Cocks-Pinch curves of embedding degrees five to eight and optimal ate pairing computation
    Aurore Guillevic, Simon Masson, Emmanuel Thome, Design, Codes and Cryptography, 2020
    https://dx.doi.org/10.1007/s10623-020-00727-w
    """
    u = ZZ(2**64 - 2**61 + 2**15)
    D = 10**10 + 147
    i = 1
    h_t = 3
    h_y = 0x11e36418c7c8b454
    p = 0x40000138cd26ab94b86e1b2f7482785fa18f877591d2a4476b4760217f860bfe8674e2a4610d669328bda13044c030e8cc836a5b363f2d4c8abcab71b12091356bb4695c5626bc319d38bf65768c5695f9ad97
    r = 0x9610000000015700ab80000126012600c4007000a800e000f000200040008001
    k = 5
    assert r == cyclotomic_polynomial(k)(u)
    tr0 = u+1
    Fr = GF(r)
    sqrt_D = sqrt(Fr(-D))
    yr = (tr0-2)/sqrt_D
    y0 = ZZ(yr)
    if (r-y0) < y0:
        y0 = r-y0
    assert p == ((tr0 + h_t*r)**2 + D*(y0 + h_y*r)**2)//4
    tr = tr0 + h_t*r
    y = y0 + h_y*r
    assert ((p+1-tr0) % r) == 0
    assert ((p+1-tr) % r) == 0
    assert (cyclotomic_polynomial(5)(p) % r) == 0
    Fp = GF(p)
    b5 = 0x3dd2d2b0b2e68770bf01b41946ab867390cf9ecc4a858004fc769c278f079574677c7db3e7201c938b099f85eb6e85f200b95a80b24fdbdf584098d690c6b91b21d00f52cc79473a11123b08ab2a616b4a4fbf
    E5 = EllipticCurve([Fp(-3), Fp(b5)])
    P = E5.random_element()
    assert (p+1-tr)*P == E5(0)

    test_finite_field_tnfs(p, r, k, cost=180, samples=100000, conj=True, all_deg_h=[5], max_coeff=2, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=150, samples=100000, jouxlercier=True, all_deg_h=[5], max_coeff=1, deg_f=3, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=136, samples=100000, jouxlercier=True, all_deg_h=[5], max_coeff=1, deg_f=4, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=129, samples=100000, jouxlercier=True, all_deg_h=[5], max_coeff=1, deg_f=5, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=129, samples=100000, jouxlercier=True, all_deg_h=[5], max_coeff=1, deg_f=6, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=131, samples=100000, jouxlercier=True, all_deg_h=[5], max_coeff=1, deg_f=7, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)

def test_cocks_pinch_k8_gmt():
    """
    Run the NFS simulation tool on the parameters of a Cocks-Pinch k=8 curve
    Pairing-friendly elliptic curve from https://hal.inria.fr/hal-02305051v2
    Cocks-Pinch curves of embedding degrees five to eight and optimal ate pairing computation
    Aurore Guillevic, Simon Masson, Emmanuel Thome, Design, Codes and Cryptography, 2020
    https://dx.doi.org/10.1007/s10623-020-00727-w
    k = 8, D = 4, r = cyclotomic_polynomial(8)(u) = u^4+1, t = (u^i+1) mod r + ht*r, y = (tr0-1)*u^2/2 + hy*r = (u^i-1)*u^2/2 mod r + hy*r
    """
    u = ZZ(0xffc00020fffffffc)
    D = 4
    i = 1
    ht = 1
    hy = 0xdc04
    k = 8
    tr0 = u**i+1
    r = ZZ(cyclotomic_polynomial(8)(u))
    tr = tr0 + ht*r
    y0 = (((u**i-1)*u**2) % r)//2
    y_ = y0 + hy*r
    p4 = (tr**2 + D*y_**2)
    assert (p4 % 4) == 0
    p = p4 // 4

    assert u == 2**64 - 2**54 + 2**37 + 2**32 - 2**2
    q = ZZ(0xbb9dfd549299f1c803ddd5d7c05e7cc0373d9b1ac15b47aa5aa84626f33e58fe66943943049031ae4ca1d2719b3a84fa363bcd2539a5cd02c6f4b6b645a58c1085e14411)
    s = ZZ(0xff0060739e18d7594a978b0ab6ae4ce3dbfd52a9d00197603fffdf0000000101)
    assert q == p
    assert r == s
    assert (q % s) == u**i
    assert s == u**4 + 1
    assert p == ((tr0 + ht*r)**2 + D*(y0 + hy*r)**2)//4
    # obtain the special form for px
    X = y # because of ZZy and Rxy
    rx = X**4 + 1
    t0x = X + 1
    tx = t0x + ht*rx
    y0x = (t0x-2) * X**2/2
    yx = y0x + hy*rx
    px = (tx**2 + D*yx**2)
    assert px(u)//4 == p
    print("px={}".format(px))
    Fp = GF(p)
    E = EllipticCurve([Fp(2), Fp(0)])
    P = E.random_element()
    assert (p+1-tr)*P == E(0)

    # in addition to the TNFS polynomial selections as for the k=5 curve, there is the Special technique and the SarkarSingh technique.
    # Special
    test_finite_field_tnfs(p, r, k, cost=363, samples=100000, special=True, all_deg_h=[2], qx=px, u=u, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=165, samples=100000, special=True, all_deg_h=[4], qx=px, u=u, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=136, samples=100000, special=True, all_deg_h=[8], qx=px, u=u, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    # Conj
    test_finite_field_tnfs(p, r, k, cost=131, samples=100000, conj=True, all_deg_h=[2], max_coeff=2, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=153, samples=100000, conj=True, all_deg_h=[4], max_coeff=2, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    test_finite_field_tnfs(p, r, k, cost=209, samples=100000, conj=True, all_deg_h=[8], max_coeff=2, B0_alpha=800, B1_alpha=1200, compute_alpha=True, alpha_test_principal=False)
    # Sarkar-Singh
    # todo
    # Joux-Lercier, check what's happening when gcd(deg(h), k/deg(h)) > 1
    # todo

if __name__ == "main":
    test_cocks_pinch_k5_gmt()
    test_cocks_pinch_k8_gmt()
