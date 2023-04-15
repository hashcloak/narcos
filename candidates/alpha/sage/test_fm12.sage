import sys
import re
import os
import pprint

import tnfs
from tnfs.curve.cyclo_kDe import Cyclo_kDe
from tnfs.curve.cyclo_kDe import polynomial_params
import tnfs.curve.fotiadis_martindale
from tnfs.simul.polyselect import Polyselect
from tnfs.simul.polyselect_utils import automorphism_factor
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.simul.simulation_nfs import Simulation_NFS
from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d
from tnfs.alpha.alpha2d import alpha2d
from tnfs.alpha.alpha3d import alpha3d

from tnfs.param.testvector_sparseseed import test_vector_fm12 as test_vector

args=sys.argv
print("usage: sage {} <idx> [<samples>]".format(args[0]))

if len(args)>=2:
    idx=Integer(args[1])
else:
    idx = 0
if len(args)>=3:
    samples=Integer(args[2])
else:
    samples=10^5

# the two curves in this set are
u399 = -2**64-2**63-2**11-2**10
u446 = -2**72-2**71-2**36
E1 = Cyclo_kDe(k=12,D=3,e0=17,u=u399,b=-2)
E2 = Cyclo_kDe(k=12,D=3,e0=17,u=u446,b=-2)

k=12
D=3
px, rx, tx, cx, yx, betax, lambx, D = tnfs.curve.fotiadis_martindale.polynomial_params_from_code(17)
print(" i    p   p^k cost dh snfs dim (all p, pk, cost in bits)")
i = idx-2
for v in test_vector[max(0,idx-2):min(len(test_vector),idx+3)]:
    u = v['u']
    p = ZZ(px(u))
    lnpk = (p^k).nbits()
    print("{:2d} {:1d} {:4d} {:5d}  {:3d} {:1d}".format(i, k, v['pnbits'], lnpk, v['cost_S'], v['deg_h_S']))
    i += 1

vec = test_vector[idx]
u, b, pnbits, rnbits, cost, deg_h, label = vec['u'],vec['b'],vec['pnbits'],vec['rnbits'],vec['cost_S'],vec['deg_h_S'],vec['label']

if pnbits >= 256:
    proof.arithmetic(False)

print("k={}, p {} bits, deg_h = {}, cost = {}".format(k, pnbits, deg_h, cost))

Special = True ; Conj = False; SSingh = False; GJL=False; polyselect_label = "Special"
#Special = False; Conj = True ; SSingh = False; GJL=False; polyselect_label = "Conj"
#Special = False; Conj = False; SSingh = True; GJL=False; polyselect_label = "SarkarSingh"
#Special = False; Conj = False; SSingh = False; GJL=True; polyselect_label = "GJL"

print("FM12-{} {}".format(k, pnbits, label))

E = Cyclo_kDe(k=k,D=D,e0=17,u=u,b=b)
E.print_parameters()
print("")

ZZy.<y> = PolynomialRing(ZZ)
Rxy.<x> = PolynomialRing(ZZy)
k = E.k()
p = E.p()
ell = E.r()
cofactor = E.c()
tr = E.tr()
Fp = E.Fp()
Fpz,z = E.Fpz()

ps = Polyselect(E)
poly_p = ZZy(E.poly_p())
poly_u = y
# change of variables
poly_p = ZZy(108*poly_p((y-2)/6))
poly_u = 6*y+2
print("poly_p = {}\npoly_u = {}\npoly_u(u) = {}".format(poly_p, poly_u, poly_u(u)))
print("poly_p(({})(u)) = {} mod p".format(poly_u, poly_p(poly_u(u)) % p ))

ps.compute_h(deg_h)
E_list_h = ps.get_h()
print("polynomials h")
ps.print_h(deg_h=deg_h)
print(polyselect_label+" TNFS")

deg_g=(k//deg_h)

for tup_h in E_list_h[deg_h]:
    inv_zeta_Kh, w, hc = tup_h
    inv_zeta_Kh = float(inv_zeta_Kh)
    
    if (gcd(deg_g,deg_h) ==1):
        aut_h = automorphism_factor(hc)
    else:
        aut_h = 1
    h = ZZy(hc)

    if Special:
        res = ps.TNFS_Special(deg_g=deg_g,h=h,poly_p=poly_p,u=ZZ(poly_u(u)),with_y=(gcd(deg_g,deg_h) > 1))
    elif Conj:
        res = ps.TNFS_Conjugation(deg_g=deg_g,h=h,with_y=(gcd(deg_g,deg_h) > 1),max_coeff=2)
    elif SSingh:
        res = ps.TNFS_SarkarSingh_JL(deg_aux1=deg_aux1, deg_g=deg_g, h=h, with_y=(gcd(deg_g,deg_h) > 1), max_coeff=max_coeff_aux1)

    if res == None:
        print("h = {}, TNFS_special returned nothing".format(h))
        continue
    f, g, max_fi, max_gi, aut_fg = res[:5]

    alpha_f = float(alpha_TNFS_2d(f,h,1000,test_principal=True))
    alpha_g = float(alpha_TNFS_2d(g,h,1000,test_principal=True))
    sum_alpha = alpha_f+alpha_g
    aut = aut_h*aut_fg
    print("inv_zeta_Kh, w = {:.6f},{:2d}".format(inv_zeta_Kh, w))
    Kh.<ah> = NumberField(h)
    clKh = Kh.class_number()
    print("h = {} # {} class number Kh = {}".format(h, h.list(), clKh))
    print("f = {}".format(f))
    print("g = {}".format(g))
    print("max_fi = {:6d}, log_2 max_gi = {:6.2f}".format(max_fi, float(log(max_gi,2))))
    print("aut = {}".format(aut))
    print("alpha_f = {:.4f} alpha_g = {:.4f}".format(alpha_f,alpha_g))
    print("    ({:.8f},{:2d}, {:40s}, {:.4f}, {:.4f}, {:.4f}),".format(float(inv_zeta_Kh),int(w),str(h),float(alpha_f),float(alpha_g),float(sum_alpha)))
    simul = Simulation_TNFS(p,ell,Fp,Fpz,h,f,g,Rxy,cost,aut,inv_zeta_Kh,count_sieving=True,alpha_f=alpha_f, alpha_g=alpha_g)
    simul.print_params()
    simul.simulation(samples=samples)
    simul.print_results()
    print("#::::::::::::::")
    # if there is not enough relations of there are too many relations, re-run with the same polynomials but with a higher/smaller cost
    simul.adjust_cost(samples)
    print("###############")
