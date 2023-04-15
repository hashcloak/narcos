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
from tnfs.alpha.alpha2d import alpha2d
from tnfs.alpha.alpha3d import alpha3d

from tnfs.simul.polyselect import Polyselect
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.simul.simulation_nfs import Simulation_NFS

from tnfs.simul.polyselect_utils import automorphism_factor
from tnfs.simul.polyselect_utils import pretty_print_coeffs_from_coeffs, pretty_print_poly_from_coeffs
from tnfs.gen.generate_curve_utils import str_binary_form, str_py_binary_form

# polynomial rings
Rx = ZZ['x']; (x,) = Rx._first_ngens(1)

def test_finite_field_nfs(q, r, k, cost, sieving_dim=2, samples=100000, special=False, conj=False, sarkarsingh=False, jouxlercier=False, qx=None, u=None, max_coeff=2, deg_f=None, deg_phi_base=None, B0_alpha=800, B1_alpha=2000, compute_alpha=True, verbose=False):
    """
    run NFS for the field GF(q^k) with r a prime divisor of the cyclotomic subgroup
    INPUT:
    - `q`: prime integer, characteristic of the target field
    - `r`: prime integer, order of subgroup of GF(q^k)
    - `k`: extension degree
    - `cost`: expected cost of DL
    - `sieving_dim`: dimension of sieving
    - `samples`: number of samples in the Monte Carlo simulation
    - `special`: Special-NFS
    - `conj`: Conjugation-NFS
    - `sarkarsingh`: Sarkar-Singh-NFS
    - `jouxlercier`: Joux-Lercier-NFS
    - `qx`: polynomial for the special form of q if special is chosen
    - `u`: seed such that q = qx(u) for the special SNFS method
    - `max_coeff`: for Conj or Sarkar-Singh or Joux-Lercier methods
    - `deg_f`: for JL and GJL (jouxlercier) deg_f > k, or Sarkar-Singh: deg_f >= (deg_phi_top+1)*deg_phi_base
    - `deg_phi_base`: for Sarkar-Singh: deg_phi_top*deg_phi_base = k, gcd(aux0, aux1) = phi_top, f=Resultant(aux1, phi_base), g = Resultant(aux0, phi_base)
    - `B0_alpha`: bound on small primes for a first approximation on alpha
    - `B1_alpha`: bound on small primes for a more precise approximation of alpha
    - `compute_alpha`: True/False. If False, set alpha to 0.5
    """
    inv_zeta = float(1.0/zeta(float(sieving_dim)))
    if sieving_dim > 3:
        compute_alpha = False

    poly_init = Polyselect(p=q, k=k)
    # 2. NFS
    if special:
        print("special-NFS")
        res_poly = poly_init.NFS_Special(deg_g=k, poly_p=qx, u=u)
    elif conj:
        print("conj-NFS")
        res_poly = poly_init.NFS_Conjugation(deg_g=k, sieving_dim=sieving_dim, max_coeff=max_coeff, monic=False, compute_alpha=compute_alpha, verbose=2, number_results=1, B0_alpha=B0_alpha, B1_alpha=B1_alpha)
        # f, g, aux0, aux1, max_fi, max_gi, aut_g, alpha_f, alpha_g, sum_alpha, score
    elif sarkarsingh:
        print("Sarkar-Singh-NFS")
        res_poly = poly_init.NFS_SarkarSingh_JL(deg_phi_top=k//deg_phi_base, deg_aux1=deg_f//deg_phi_base, deg_phi_base=deg_phi_base, max_coeff=max_coeff, monic=False, compute_alpha=compute_alpha, sieving_dim=sieving_dim, max_test_poly_aux1=100, B1_alpha=B1_alpha, verbose=verbose)
    elif jouxlercier:
        print("(Generalized)Joux-Lercier-NFS")
        assert deg_f > k
        res_poly = poly_init.GeneralizedJouxLercier(deg_f=deg_f, max_coeff=max_coeff, deg_phi=k, monic=False, sieving_dim=sieving_dim, compute_alpha=compute_alpha, B1_alpha=B1_alpha, max_test_polys=100)

    if res_poly is None or len(res_poly) == 0:
        raise ValueError("Error in Polynomial selection")
    if conj:
        print("res_poly[0] = {}".format(res_poly[0]))
        f, g, aux0, aux1, max_fi, max_gi, aut_g = res_poly[0][:7]
        aut = aut_g
    elif jouxlercier:
        f, g, max_fi, max_gi = res_poly[0][:4]
        aut = 1
    else:
        f, g, max_fi, max_gi, aut_fg = res_poly[:5]
        aut = aut_fg

    print("f = {}".format(f))
    print("g = {}".format(g))
    print("max_fi = {:6d}, log_2 max_gi = {:6.2f}".format(max_fi, log(float(max_gi) ,2.0)))
    print("aut = {}".format(aut))

    # computing alpha, this takes at least a few seconds

    if sieving_dim == 2:
        alpha_f = float(alpha2d(f, B1_alpha))
        alpha_g = float(alpha2d(g, B1_alpha))
    elif sieving_dim == 3:
        alpha_f = float(alpha3d(f, B1_alpha))
        alpha_g = float(alpha3d(g, B1_alpha))
    else:
        alpha_f = 0.5
        alpha_g = 0.5
    sum_alpha = alpha_f+alpha_g
    print("alpha_f = {:.4f} alpha_g = {:.4f} sum_alpha = {:.4f}".format(alpha_f, alpha_g, sum_alpha))
    print("    ({:.8f}, {:.4f}, {:.4f}, {:.4f}),".format(float(inv_zeta), float(alpha_f),float(alpha_g),float(sum_alpha)))

    #initialisation of data
    Fp = poly_init.get_Fp()
    Fpz = poly_init.get_Fpz()
    simul = Simulation_NFS(q,r,Fp,Fpz,sieving_dim,f,g,Rx,cost,aut,count_sieving=True,alpha_f=alpha_f,alpha_g=alpha_g)

    simul.print_params()
    simul.simulation(samples=samples) #takes few seconds for 10^4, mins for 10^5, up to 20 min for 10^6
    simul.print_results()
    print("#::::::::::::::")
    # if there is not enough relations of there are too many relations, re-run with the same polynomials but with a higher/smaller cost
    simul.adjust_cost(samples=samples)
    print("############")

def test_cocks_pinch_k5_gmt():
    """
    Run the NFS simulation tool on the parameters of a Cocks-Pinch k=5 curve
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

    # with Conjugation, deg(f) = 10 and ||f|| = 1, deg(g) = 5 and ||g||=p^(1/2)
    test_finite_field_nfs(p, r, k, cost=127, sieving_dim=2, samples=100000, conj=True, max_coeff=3, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=131, sieving_dim=3, samples=100000, conj=True, max_coeff=3, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Sarkar-Singh, deg(f) = 15 and ||f|| = 1, deg(g) = 10 and ||g||=p^(1/3)
    test_finite_field_nfs(p, r, k, cost=202, sieving_dim=2, samples=100000, sarkarsingh=True, max_coeff=2, deg_f=15, deg_phi_base=5, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=151, sieving_dim=3, samples=100000, sarkarsingh=True, max_coeff=2, deg_f=15, deg_phi_base=5, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Joux-Lercier, deg(f) = 6 and ||f|| = 1, deg(g) = 5 and ||g||=p^(5/6)
    test_finite_field_nfs(p, r, k, cost=129, sieving_dim=2, samples=100000, jouxlercier=True, max_coeff=1, deg_f=6, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=160, sieving_dim=3, samples=100000, jouxlercier=True, max_coeff=1, deg_f=6, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Joux-Lercier, deg(f) = 7 and ||f|| = 1, deg(g) = 6 and ||g||=p^(5/7)
    test_finite_field_nfs(p, r, k, cost=132, sieving_dim=2, samples=100000, jouxlercier=True, max_coeff=1, deg_f=7, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=153, sieving_dim=3, samples=100000, jouxlercier=True, max_coeff=1, deg_f=7, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

def test_cocks_pinch_k7_gmt():
    """
    Run the NFS simulation tool on the parameters of a Cocks-Pinch k=7 curve
    Pairing-friendly elliptic curve from https://hal.inria.fr/hal-02305051v2
    Cocks-Pinch curves of embedding degrees five to eight and optimal ate pairing computation
    Aurore Guillevic, Simon Masson, Emmanuel Thome, Design, Codes and Cryptography, 2020
    https://dx.doi.org/10.1007/s10623-020-00727-w
    cyclotomic_polynomial(7) = x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
    """
    k = 7
    D = 20
    i = 6
    ht = -2
    hy = 0
    u = ZZ(2**43 - 2**41 - 0x47dfdb8)
    r = cyclotomic_polynomial(k)(u)
    tr0 = -(u**5 + u**4 + u**3 + u**2 + u) # (u**i + 1) % r
    tr = tr0 + ht*r
    # -20 is not a square in CyclotomicField(7), get sqrt(-20) modulo r
    Fr = GF(r)
    sqrt_D = ZZ(sqrt(Fr(-D)))
    y0 = ZZ(Fr((tr0-2)*sqrt_D/D))
    if r-y0 < y0:
        y0 = r-y0
    y = y0 + hy*r
    p4 = (tr**2 + D*y**2)
    assert (p4 % 4) == 0
    p = p4 // 4
    s = ZZ(0xb63ccd541c3aa13c7b7098feb312eecf5648fd215c0d2916714b429d14e8f889)
    q = ZZ(0x8f591a9876a6d2344ae66dd7540ea2fd28174755d16c4ae5c5cd5c1d208e639271b48c8ba7453c95a2a9be6434f2455504d419f13e35062aa5ebbc49ecfd30f9)
    z = ZZ(0x1a6823d30bb191e5e47729d6721be17804ee9c32ad1be1f9c3f77999cce79234)
    assert r == s
    assert p == q
    assert y == z
    Fp = GF(p)
    a0 = 11
    b7 = ZZ(0x15d384c76889d377dd63600fbe42628e0c386a3e87915790188d944845aab2b649964f386dc90b3a9b6120af5da9a2aaead5e415dd958c5cfa80ea61aac268b0)
    E = EllipticCurve([Fp(-3*a0**2), Fp(b7*a0**3)])
    P = E.random_element()
    assert (q+1-tr)*P == E(0)

    # with Conjugation, deg(f) = 14 and ||f|| = 1, deg(g) = 7 and ||g||=p^(1/2)
    test_finite_field_nfs(p, r, k, cost=168, sieving_dim=2, samples=100000, conj=True, max_coeff=3, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=139, sieving_dim=3, samples=100000, conj=True, max_coeff=3, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=144, sieving_dim=4, samples=100000, conj=True, max_coeff=3, B0_alpha=800, B1_alpha=2000, compute_alpha=False) # alpha_4d not implemented

    # with Sarkar-Singh, deg(f) = (7)*3 and ||f|| = 1, deg(g) = (7)*2 and ||g||=p^(1/3)
    test_finite_field_nfs(p, r, k, cost=202, sieving_dim=2, samples=10000, sarkarsingh=True, max_coeff=2, deg_f=21, deg_phi_base=7, B0_alpha=800, B1_alpha=2000, compute_alpha=False)
    test_finite_field_nfs(p, r, k, cost=151, sieving_dim=3, samples=10000, sarkarsingh=True, max_coeff=2, deg_f=21, deg_phi_base=7, B0_alpha=800, B1_alpha=2000, compute_alpha=False)

    # with Joux-Lercier, deg(f) = 8 and ||f|| = 1, deg(g) = 7 and ||g||=p^(7/8)
    test_finite_field_nfs(p, r, k, cost=143, sieving_dim=2, samples=100000, jouxlercier=True, max_coeff=1, deg_f=8, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=155, sieving_dim=3, samples=100000, jouxlercier=True, max_coeff=1, deg_f=8, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Joux-Lercier, deg(f) = 9 and ||f|| = 1, deg(g) = 8 and ||g||=p^(7/9)
    test_finite_field_nfs(p, r, k, cost=151, sieving_dim=2, samples=100000, jouxlercier=True, max_coeff=1, deg_f=9, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=154, sieving_dim=3, samples=100000, jouxlercier=True, max_coeff=1, deg_f=9, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

def test_cocks_pinch_k6_gmt():
    """
    Run the NFS simulation tool on the parameters of a Cocks-Pinch k=6 curve
    Pairing-friendly elliptic curve from https://hal.inria.fr/hal-02305051v2
    Cocks-Pinch curves of embedding degrees five to eight and optimal ate pairing computation
    Aurore Guillevic, Simon Masson, Emmanuel Thome, Design, Codes and Cryptography, 2020
    https://dx.doi.org/10.1007/s10623-020-00727-w
    """
    k = 6
    D = 3
    u = ZZ(2**128 - 2**124 - 2**69)
    r = ZZ(cyclotomic_polynomial(k)(u))
    i = 1
    ht = -1
    hy = ZZ(2**80 - 2**70 - 2**66 - 0x3fe0)
    tr0 = u**i + 1
    tr = tr0 + ht*r
    #y0 = (tr0-2)*(2*u-1)//3
    y0 = (u**2 + 2)//3
    y = y0 + hy*r
    p4 = tr**2 + D*y**2
    assert (p4 % 4) == 0
    p = p4//4
    q = ZZ(0x9401ff90f28bffb0c610fb10bf9e0fefd59211629a7991563c5e468d43ec9cfe1549fd59c20ab5b9a7cda7f27a0067b8303eeb4b31555cf4f24050ed155555cd7fa7a5f8aaaaaaad47ede1a6aaaaaaaab69e6dcb)
    s = ZZ(0xe0ffffffffffffc400000000000003ff10000000000000200000000000000001)
    z = ZZ(0xe0c43bffffffffc3d7cc6b000000040cf8abbfffffffff20b4b755555555554e59115555555555551576)
    assert r == s
    assert y == z
    assert p == q
    assert tr**2 - 4*p == -D*y**2
    Fp = GF(p)
    E = EllipticCurve([Fp(0), Fp(-1)])
    P = E.random_element()
    assert (p+1-tr)*P == E(0)
    # obtain the special form for px
    rx = Rx(cyclotomic_polynomial(k))
    t0x = x**i + 1
    tx = t0x + ht*rx
    y0x = (x**2 + 2)/3
    yx = y0x + hy*rx
    px = (tx**2 + D*yx**2)
    assert px(u)//4 == p
    print("px={}".format(px))

    # in addition to the NFS polynomial selections as for the k=5 curve, there is the Special technique and the SarkarSingh technique.
    # but what works the best is the TNFS technique.
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=2, samples=10000, special=True, qx=px, u=u, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=3, samples=10000, special=True, qx=px, u=u, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Conjugation, deg(f) = 12 and ||f|| = 1, deg(g) = 6 and ||g||=p^(1/2)
    test_finite_field_nfs(p, r, k, cost=128, sieving_dim=2, samples=10000, conj=True, max_coeff=2, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=3, samples=10000, conj=True, max_coeff=2, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Sarkar-Singh, deg(f) = (3)*3 and ||f|| = 1, deg(g) = (3)*2 and ||g||=p^(2/3)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=2, samples=10000, sarkarsingh=True, max_coeff=2, deg_f=9, deg_phi_base=3, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=3, samples=10000, sarkarsingh=True, max_coeff=2, deg_f=9, deg_phi_base=3, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    # with Sarkar-Singh, deg(f) = (2)*4 and ||f|| = 1, deg(g) = (2)*3 and ||g||=p^(3/4)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=2, samples=10000, sarkarsingh=True, max_coeff=2, deg_f=8, deg_phi_base=2, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=3, samples=10000, sarkarsingh=True, max_coeff=2, deg_f=8, deg_phi_base=2, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Joux-Lercier, deg(f) = 7 and ||f|| = 1, deg(g) = 6 and ||g||=p^(6/7)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=2, samples=10000, jouxlercier=True, max_coeff=1, deg_f=7, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=3, samples=10000, jouxlercier=True, max_coeff=1, deg_f=7, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Joux-Lercier, deg(f) = 8 and ||f|| = 1, deg(g) = 7 and ||g||=p^(6/8) = p^(3/4)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=2, samples=10000, jouxlercier=True, max_coeff=1, deg_f=8, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=200, sieving_dim=3, samples=10000, jouxlercier=True, max_coeff=1, deg_f=8, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

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
    r = cyclotomic_polynomial(8)(u)
    tr = tr0 + ht*r
    y0 = (((u**i-1)*u**2) % r)//2
    y = y0 + hy*r
    p4 = (tr**2 + D*y**2)
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
    rx = x**4 + 1
    t0x = x + 1
    tx = t0x + ht*rx
    y0x = (t0x-2) * x**2/2
    yx = y0x + hy*rx
    px = (tx**2 + D*yx**2)
    assert px(u)//4 == p
    Fp = GF(p)
    E = EllipticCurve([Fp(2), Fp(0)])
    P = E.random_element()
    assert (p+1-tr)*P == E(0)

    # in addition to the NFS polynomial selections as for the k=5 curve, there is the Special technique and the SarkarSingh technique.
    # but what works the best is the TNFS technique.
    test_finite_field_nfs(p, r, k, cost=904, sieving_dim=2, samples=100000, special=True, qx=px, u=u, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=534, sieving_dim=3, samples=100000, special=True, qx=px, u=u, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Conjugation, deg(f) = 16 and ||f|| = 1, deg(g) = 8 and ||g||=p^(1/2)
    test_finite_field_nfs(p, r, k, cost=199, sieving_dim=2, samples=100000, conj=True, max_coeff=2, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=156, sieving_dim=3, samples=100000, conj=True, max_coeff=2, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Sarkar-Singh, deg(f) = 4*3 and ||f|| = 1, deg(g) = 4*2 and ||g||=p^(2/3)
    test_finite_field_nfs(p, r, k, cost=169, sieving_dim=2, samples=100000, sarkarsingh=True, max_coeff=2, deg_f=12, deg_phi_base=4, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=156, sieving_dim=3, samples=100000, sarkarsingh=True, max_coeff=2, deg_f=12, deg_phi_base=4, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    # deg(f) = 2*5 and ||f|| = 1, deg(g) = 2*4 and ||g||=p^(4/5)
    test_finite_field_nfs(p, r, k, cost=162, sieving_dim=2, samples=100000, sarkarsingh=True, max_coeff=2, deg_f=10, deg_phi_base=2, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=163, sieving_dim=3, samples=100000, sarkarsingh=True, max_coeff=2, deg_f=10, deg_phi_base=2, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Joux-Lercier, deg(f) = 9 and ||f|| = 1, deg(g) = 8 and ||g||=p^(8/9)
    test_finite_field_nfs(p, r, k, cost=161, sieving_dim=2, samples=100000, jouxlercier=True, max_coeff=1, deg_f=9, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=169, sieving_dim=3, samples=100000, jouxlercier=True, max_coeff=1, deg_f=9, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

    # with Joux-Lercier, deg(f) = 10 and ||f|| = 1, deg(g) = 9 and ||g||=p^(8/10) = p^(4/5)
    test_finite_field_nfs(p, r, k, cost=171, sieving_dim=2, samples=100000, jouxlercier=True, max_coeff=1, deg_f=10, B0_alpha=800, B1_alpha=2000, compute_alpha=True)
    test_finite_field_nfs(p, r, k, cost=168, sieving_dim=3, samples=100000, jouxlercier=True, max_coeff=1, deg_f=10, B0_alpha=800, B1_alpha=2000, compute_alpha=True)

if __name__ == "main":
    test_cocks_pinch_k5_gmt()
    test_cocks_pinch_k7_gmt()
    test_cocks_pinch_k6_gmt()
    test_cocks_pinch_k8_gmt()
