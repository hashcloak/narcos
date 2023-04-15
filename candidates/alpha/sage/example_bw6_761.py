"""
Example: template file to run the code
Assume you have a pairing-friendly curve, and you have a polynomial form of the parameters, 
_p_(_x_), _r_(_x_) stand for the characteristic and the (prime subgroup) order of the curve.
You would like to use this code to estimate the security of the finite field GF(_p<sup>k</sup>_)
for some embedding degree _k_ and some integer seed _u_ such that _p_=_p_(_u_).

The polynomial selection is implemented for various parameters, but the 
automatic run of the code for all subfields is not (yet) implemented.
You would need to run it manually (that is, adjust the estimated cost for different deg(h)).

# this example file is given with parameters from https://eprint.iacr.org/2020/351
# Optimized and secure pairing-friendly elliptic curves suitable for one layer proof composition
# Youssef El Housni and Aurore Guillevic


"""
# imports from SageMath
import sage
import sage.all

from sage.rings.integer_ring import Z
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ, Q
from sage.structure.proof.all import arithmetic
from sage.arith.misc import GCD, gcd
from sage.arith.functions import lcm
from sage.functions.log import log
from sage.rings.number_field.number_field import NumberField

# imports from TNFS package
import tnfs
from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d

from tnfs.simul.polyselect import Polyselect
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.simul.simulation_nfs import Simulation_NFS

from tnfs.simul.polyselect_utils import automorphism_factor
from tnfs.simul.polyselect_utils import pretty_print_coeffs_from_coeffs, pretty_print_poly_from_coeffs

# Embedding degree
k = 6
# seed
u = 0x8508c00000000001
# polynomial form of the characteristic
QQs = QQ['s']; (s,) = QQs._first_ngens(1)
ZZy = ZZ['y']; (y,) = ZZy._first_ngens(1)

Px = (QQs(103*s**12 - 379*s**11 + 250*s**10 + 691*s**9 - 911*s**8 - 79*s**7 + 623*s**6 - 640*s**5 + 274*s**4 + 763*s**3 + 73*s**2 + 254*s + 229))/9
p = Px(ZZ(u))
assert p in ZZ
p = ZZ(p)
assert p.is_prime()
Rx = (QQs(s**6 - 2*s**5 + 2*s**3 + s + 1))/3
ell = Rx(ZZ(u))
assert ell in ZZ
ell = ZZ(ell) # cast in the ring of integers
assert ell.is_prime()

arithmetic(False) # or: proof.arithmetic(False) in Sage
# set a polynomial px of integer coefficients, there is no need to change the seed u,
# we will have px(u) = 0 mod p
Px_den = lcm([pi.denom() for pi in Px.list()])
px = ZZy(Px_den * Px)

poly_init = Polyselect(p=p, k=k)
Fp = poly_init.get_Fp()
Fpz,z = poly_init.get_Fpz()

####################################
cost = 127
####################################
# You need to provide a first guess on the cost of DL computation.
# if you have no idea, put 128 (this is in bits) and the code will adjust and try to converge itself.
# However it might take much more time.
# Update this value when you are more sure.
# for example, this is roughly 127 for deg_h = 6, 235 for deg_h = 3, for deg_h = 2.
# it varies few bits for two different polynomials h of same degree
# (e.g. from 127 to 131 for deg_h to be 6)

# 1. Special-TNFS family of DL computation
for deg_h in [i for i in range(k,1,-1) if (k % i) == 0]:
    # for k==6, we have deg_h in [2,3,6]
    # you might consider all possible divisors of k, including k, greater than 1
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

    for item in list_h:
        inv_zeta_Kh, w, hc = item
        hc_str = pretty_print_coeffs_from_coeffs(hc)
        h_str = pretty_print_poly_from_coeffs(hc)
        print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
        h = ZZy(hc)

        res_poly = poly_init.TNFS_Special(deg_g=deg_g, h=h, poly_p=px, u=u, with_y=(gcd(deg_g,deg_h)==1))
        if res_poly == None:
            raise ValueError("Error in Polynomial selection")
        f, g, max_fij, max_gij, aut_fg = res_poly[:5]
        Kh = NumberField(h, names=('ah',)); (ah,) = Kh._first_ngens(1)
        print("h = {} # {}".format(h, h.list()))
        print("inv_zeta_Kh, w = {:.6f},{:2d}".format(inv_zeta_Kh, w))
        assert (ZZ(h.resultant(f.resultant(g))) % p**k) == 0
        Rxy = f.parent()

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

        alpha_f = float(alpha_TNFS_2d(f,h,1000))
        alpha_g = float(alpha_TNFS_2d(g,h,1000))
        sum_alpha = alpha_f+alpha_g
        print("alpha_f = {:.4f} alpha_g = {:.4f} sum_alpha = {:.4f}".format(alpha_f,alpha_g,sum_alpha))
        print("    ({:.8f},{:2d}, {:40s}, {:.4f}, {:.4f}, {:.4f}),".format(float(inv_zeta_Kh),int(w),str(h),float(alpha_f),float(alpha_g),float(sum_alpha)))
        
        #initialisation of data
        simul = Simulation_TNFS(p,ell,Fp,Fpz,h,f,g,Rxy,cost,aut,inv_zeta_Kh,count_sieving=True,alpha_f=alpha_f,alpha_g=alpha_g)

        simul.print_params()
        simul.simulation(samples=100000) #takes few seconds for 10^4, mins for 10^5, up to 20 min for 10^6
        simul.print_results()
        print("#::::::::::::::")
        # if there is not enough relations of there are too many relations, re-run with the same polynomials but with a higher/smaller cost
        simul.adjust_cost(samples=100000)
        print("############")

    
# 2. Special-NFS family of DL computation (would correspond to deg h = 1)
# TODO template

