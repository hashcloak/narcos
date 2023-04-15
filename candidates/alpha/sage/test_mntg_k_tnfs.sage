import sys
import re
import os
import pprint
##
## MNTG params (generalized MNT)
##
from tnfs.curve.mntg import MNTG
from tnfs.param.testvector_mnt_k import test_vector_MNTG, test_vector_MNT6_SB06
testvector = test_vector_MNTG + test_vector_MNT6_SB06

from tnfs.simul.polyselect import Polyselect
from tnfs.simul.polyselect_utils import automorphism_factor
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d

args=sys.argv
print("usage: {} [<idx>] [<samples>]".format(args[0]))
# k is given in the testvector
if len(args)>=2:
    idx=Integer(args[1])
else:
    idx = 0
if len(args)>=3:
    samples=Integer(args[2])
else:
    samples=10^5

k = testvector[idx]['k']
d = testvector[idx]['d']
c = testvector[idx]['c']

print("                       cost d cost d deg      s      s       d s    d d s")
print("  i k  c  d    p   p^k   ST h  SST h aux SNFS d Conj d   GJL f d SS f g d    (all p, pk, cost in bits)")
i = max(0,idx-2)
#i=0
#for vec in testvector:
for vec in testvector[max(0,idx-2):min(len(testvector),idx+3)]:
    k_ = vec['k']
    if k_ != k:
        i += 1
        continue
    u = vec['u']
    D = vec['D']
    c = vec['c']
    d = vec['d']
    pnbits = vec['pnbits']
    lnpk = k*pnbits
    try:
        cost_STNFS, deg_h_STNFS = vec['cost_S'], vec['deg_h_S']
        cost_SS_TNFS, deg_h_SS_TNFS, deg_aux_SS = vec['cost_SS'], vec['deg_h_SS'], vec['deg_aux_SS']
        cost_SNFS, d_SNFS = vec['cost_SNFS'], vec['sieving_dim_SNFS']
        cost_Conj, d_Conj = vec['cost_NFS_Conj'], vec['sieving_dim_Conj']
        cost_GJL, deg_f_GJL, d_GJL = vec['cost_NFS_GJL'], vec['GJL_deg_f'], vec['sieving_dim_GJL']
        cost_SS, deg_f_SS, deg_g_SS, d_SS = vec['cost_NFS_SS'], vec['SS_deg_f'], vec['SS_deg_g'], vec['sieving_dim_NFS_SS']

        print("{:3d} {:1d} {:2d} {:2d} {:4d} {:5d}  {:3d} {:1d}  {:3d} {:1d} {:1d}    {:3d} {:1d}  {:3d} {:1d}   {:3d} {:1d} {:1d}   {:3d} {:1d} {:1d} {:1d}".format(i, k_, c, d, pnbits, lnpk, cost_STNFS, deg_h_STNFS, cost_SS_TNFS, deg_h_SS_TNFS, deg_aux_SS, cost_SNFS, d_SNFS, cost_Conj, d_Conj,cost_GJL, deg_f_GJL, d_GJL,cost_SS, deg_f_SS, deg_g_SS, d_SS))
    except:
        print("{:3d} {:1d} {:2d} {:2d} {:4d} {:5d}".format(i, k_, c, d, pnbits, lnpk))
    i += 1

vec = testvector[idx]
print("vec = {}".format(vec))
k, u, D, c, d, a, b, pnbits, label = vec['k'], vec['u'], vec['D'], vec['c'], vec['d'], vec['a'], vec['b'], vec['pnbits'], vec['label']

try:
    cost_S, deg_h_S = vec['cost_S'], vec['deg_h_S']
    cost_SS, deg_h_SS, deg_aux_SS = vec['cost_SS'], vec['deg_h_SS'], vec['deg_aux_SS']
except:
    raise ValueError("missing parameters cost_S, deg_h_S, cost_SS, deg_h_SS, deg_aux_SS")

if float(log(abs(u),2)) >= 256.0:
    proof.arithmetic(False)

#Special = True ; Conj = False; SSingh = False; GJL=False; polyselect_label = "Special"
Special = False; Conj = True ; SSingh = False; GJL=False; polyselect_label = "Conj"
#Special = False; Conj = False; SSingh = True; GJL=False; polyselect_label = "SarkarSingh"
#Special = False; Conj = False; SSingh = False; GJL=True; polyselect_label = "GJL"

deg_aux1 = 3
if Special:
    cost, deg_h = cost_S, deg_h_S
elif Conj:
    cost, deg_h = cost_S-1, deg_h_S
elif GJL:
    raise ValueError("TNFS-GJL: estimated cost not given")
elif SSingh:
    cost, deg_h, deg_aux1 = cost_SS, deg_h_SS, deg_aux_SS

print("k={}, p {} bits, deg_h = {}, cost = {}".format(k, pnbits, deg_h, cost))
print("MNTG{}-{} {}".format(k, pnbits, label))

E = MNTG(k, u, D, c, d, a, b)
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
poly_p = ZZy(E.px())
print("poly_p = {}".format(poly_p))
poly_u = y
# for small values, this is better to keep the original poly
# for larger values, the following changes of variables helps (for k=6)
#poly_p = ZZy(poly_p(y/2))
#poly_u = 2*y
ps.compute_h(deg_h)
E_list_h = ps.get_h()
print("polynomials h")
ps.print_h(deg_h=deg_h)
print(polyselect_label+" TNFS")

deg_g=(k//deg_h)
if deg_aux1 <= 3:
    max_coeff_aux1 = 3
elif deg_aux1 <= 4:
    max_coeff_aux1 = 2
else:
    max_coeff_aux1 = 1
for tup_h in E_list_h[deg_h]:
    inv_zeta_Kh, w, hc = tup_h
    inv_zeta_Kh = float(inv_zeta_Kh)

    if (gcd(deg_g,deg_h) ==1):
        aut_h = automorphism_factor(hc)
    else:
        aut_h = 1
    h = ZZy(hc)

    if Special:
        res = ps.TNFS_Special(deg_g=deg_g,h=h,poly_p=poly_p,u=ZZ(poly_u(u)),with_y=(gcd(deg_g,deg_h) > 1), MNT6=True)
    elif Conj:
        res = ps.TNFS_Conjugation(deg_g=deg_g,h=h,with_y=(gcd(deg_g,deg_h) > 1),max_coeff=2)
    elif SSingh:
        res = ps.TNFS_SarkarSingh_JL(deg_aux1=deg_aux1, deg_g=deg_g, h=h, with_y=(gcd(deg_g,deg_h) > 1), max_coeff=max_coeff_aux1)

    if res == None:
        print("h = {}, TNFS_{} returned nothing".format(h, polyselect_label))
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

# it could be time consuming to adjust the parameters for each h.
# The strategy is first to find the best polynomial h that provides the larger number of relations for a fixed cost (even if there are not enough or too many relations),
# then for this particular h, adjust the cost.
"""


"""
