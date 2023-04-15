import sys
import re
import os
import pprint
##
## MNT params
##
from tnfs.curve.mnt_k import MNT_k
#from tnfs.param.testvector_mnt_k import testvector_MNT_k
#from tnfs.param.testvector_mnt_k import testvector_MNT6_with_c as testvector_MNT_k
from tnfs.param.testvector_mnt_k_cycle_D_1e9 import testvector_MNT_k_cycle_D_947_999632827 as testvector_MNT_k
from tnfs.simul.polyselect import Polyselect
from tnfs.simul.polyselect_utils import automorphism_factor
from tnfs.simul.simulation_nfs import Simulation_NFS
from tnfs.alpha.alpha2d import alpha2d
from tnfs.alpha.alpha3d import alpha3d

args=sys.argv
print("usage: {} [<idx>] [<samples>]".format(args[0]))
# k is given in the testvector
if len(args)>=2:
    idx=Integer(args[1])
else:
    idx = 30 # for cost = 191, log2 p = 1514
if len(args)>=3:
    samples=Integer(args[2])
else:
    samples=10^5

k = testvector_MNT_k[idx]['k']

print("                 cost d cost d deg      s      s       d s    d d s")
print("  i k    p   p^k   ST h  SST h aux SNFS d Conj d   GJL f d SS f g d    (all p, pk, cost in bits)")
i = max(0,idx-2)
#for vec in testvector_MNT_k:
for vec in testvector_MNT_k[max(0,idx-2):min(len(testvector_MNT_k),idx+3)]:
    k_ = vec['k']
    if k_ != k:
        i += 1
        continue
    u = vec['u']
    D = vec['D']
    c = vec['c']
    pnbits = vec['pnbits']
    if k_ == 6:
        p = 4*Integer(u)**2 + 1
    elif k_ == 4:
        p = Integer(u)**2 + Integer(u) + 1
    else:
        p = 12*Integer(u)**2 - 1
    lnpk = ceil((k_)*RR(log(p,2)))
    try:
        cost_STNFS, deg_h_STNFS = vec['cost_S'], vec['deg_h_S']
        cost_SS_TNFS, deg_h_SS_TNFS, deg_aux_SS = vec['cost_SS'], vec['deg_h_SS'], vec['deg_aux_SS']
        cost_SNFS, d_SNFS = vec['cost_SNFS'], vec['sieving_dim_SNFS']
        cost_Conj, d_Conj = vec['cost_NFS_Conj'], vec['sieving_dim_Conj']
        cost_GJL, deg_f_GJL, d_GJL = vec['cost_NFS_GJL'], vec['GJL_deg_f'], vec['sieving_dim_GJL']
        cost_SS, deg_f_SS, deg_g_SS, d_SS = vec['cost_NFS_SS'], vec['SS_deg_f'], vec['SS_deg_g'], vec['sieving_dim_NFS_SS']
        
        print("{:3d} {:1d} {:4d} {:5d}  {:3d} {:1d}  {:3d} {:1d} {:1d}    {:3d} {:1d}  {:3d} {:1d}   {:3d} {:1d} {:1d}   {:3d} {:1d} {:1d} {:1d}".format(i, k_, pnbits, lnpk, cost_STNFS, deg_h_STNFS, cost_SS_TNFS, deg_h_SS_TNFS, deg_aux_SS, cost_SNFS, d_SNFS, cost_Conj, d_Conj,cost_GJL, deg_f_GJL, d_GJL,cost_SS, deg_f_SS, deg_g_SS, d_SS))
    except:
        print("{:3d} {:1d} {:4d} {:5d}".format(i, k_, pnbits, lnpk))
    i += 1

vec = testvector_MNT_k[idx]
k, u, D, c, a, b, pnbits, label = vec['k'], vec['u'], vec['D'], vec['c'], vec['a'], vec['b'], vec['pnbits'], vec['label']

try:
    cost_STNFS, deg_h_STNFS = vec['cost_S'], vec['deg_h_S']
    cost_SS_TNFS, deg_h_SS_TNFS, deg_aux_SS = vec['cost_SS'], vec['deg_h_SS'], vec['deg_aux_SS']
    cost_SNFS, d_SNFS = vec['cost_SNFS'], vec['sieving_dim_SNFS']
    cost_Conj, d_Conj = vec['cost_NFS_Conj'], vec['sieving_dim_Conj']
    cost_GJL, deg_f_GJL, d_GJL = vec['cost_NFS_GJL'], vec['GJL_deg_f'], vec['sieving_dim_GJL']
    cost_SS, deg_f_SS, deg_g_SS, d_SS = vec['cost_NFS_SS'], vec['SS_deg_f'], vec['SS_deg_g'], vec['sieving_dim_NFS_SS']
except:
    raise ValueError("missing parameters of cost for SNFS and others")
    
if float(log(abs(u),2)) >= 256.0:
    proof.arithmetic(False)

#Special = True ; Conj = False; SSingh = False; GJL=False; polyselect_label = "Special"
Special = False; Conj = True ; SSingh = False; GJL=False; polyselect_label = "Conj"
#Special = False; Conj = False; SSingh = True; GJL=False; polyselect_label = "SarkarSingh"
#Special = False; Conj = False; SSingh = False; GJL=True; polyselect_label = "GJL"

sieving_dim = 2 # default value
if Special:
    cost, sieving_dim = cost_SNFS, d_SNFS
elif Conj:
    cost, sieving_dim = cost_Conj, d_Conj
elif GJL:
    cost, deg_f, sieving_dim = cost_GJL, deg_f_GJL, d_GJL
elif SSingh:
    cost, deg_f, deg_g, sieving_dim = cost_SS, deg_f_SS, deg_g_SS, d_SS
    deg_phi2 = gcd(deg_f,deg_g)
    deg_phi1 = k//deg_phi2
    deg_aux1 = deg_f // deg_phi2
if len(args)>=4:
    sieving_dim=Integer(args[3])

inv_zeta = float(1/zeta(sieving_dim))

print("MNT{}-{} {}".format(k, pnbits, label))
print("k={}, p {} bits, cost = {}".format(k, pnbits, cost))


E = MNT_k(k, u, a=a,b=b,c=c,D=D)
E.print_parameters()
print("")

ZZx.<x> = PolynomialRing(ZZ)
k = E.k()
p = E.p()
ell = E.r()
cofactor = E.c()
tr = E.tr()
Fp = E.Fp()
Fpz,z = E.Fpz()

ps = Polyselect(E)
poly_p = ZZx(E.px())
if k == 6:
    poly_p = ZZx(poly_p(x/2))
    poly_u = 2*x
else:
    poly_u = x

if Special:
    deg_g = k
    print("p {} bits, Special, deg_f={}, deg_g={}, sieving_dim={}, cost={}".format(pnbits, poly_p.degree()*k, k, sieving_dim, cost))
    res = ps.NFS_Special(deg_g=deg_g,poly_p=poly_p,u=poly_u(E.u()))
elif Conj:
    deg_g = k
    print("p {} bits, Conj, deg_f={}, deg_g={}, sieving_dim={}, cost={}".format(pnbits, 2*k, k, sieving_dim, cost))
    res = ps.NFS_Conjugation(deg_g=deg_g,max_coeff=6)
elif SSingh:
    max_coeff_aux1=1+max(0,5-deg_aux1)
    print("p {} bits, SSingh, deg_f={}, deg_g={}, sieving_dim={}, cost={}".format(pnbits, deg_phi2*deg_aux1, deg_phi2*(deg_aux1-1), sieving_dim, cost))
    res = ps.NFS_SarkarSingh_JL(deg_phi_top=deg_phi2, deg_aux1=deg_aux1, deg_phi_base=deg_phi1, max_coeff=max_coeff_aux1,max_test_poly_aux1=500)
elif GJL:
    print("p {} bits, GJL, deg_f = {}, sieving_dim = {}, cost = {}".format(pnbits, deg_f, sieving_dim, cost))
    max_coeff_f=1+max(0,6-deg_f)
    res = ps.GeneralizedJouxLercier(deg_f=deg_f, max_coeff=max_coeff_f, deg_phi=k, monic=False,compute_alpha=True,sieving_dim=sieving_dim)
    
if res == None:
    print("NFS_"+polyselect_label+" returned nothing")
    
if not (GJL or Conj):
    res_list = [res]
else:
    res_list = res
for res_i in res_list:
    if res_i == None:
        continue
    if GJL:
        f, g, max_fi, max_gi = res_i[:4]
        aut = 1
        #max_fi = max([abs(fi) for fi in f.coefficients()])
    elif Conj:
        f, g, aux0, aux1, max_fi, max_gi, aut = res_i[:7]
    else:
        f, g, max_fi, max_gi, aut = res_i[:5]
    if sieving_dim == 2:
        alpha_f = float(alpha2d(f,2000))
        alpha_g = float(alpha2d(g,2000))
    else:
        alpha_f = float(alpha3d(f,2000))
        alpha_g = float(alpha3d(g,2000))
    sum_alpha = alpha_f+alpha_g
    print("sieving_dim = {:1d}, inv_zeta({}) = {:.5f}".format(sieving_dim,sieving_dim,float(inv_zeta)))
    print("aux0 = {}".format(aux0))
    print("aux1 = {}".format(aux1))
    print("f = {}".format(f))
    print("g = {}".format(g))
    if GJL or SSingh: # max_gi is already a logarithm
        print("max_fi = {:6d}, log_2 max_gi = {:6.2f}".format(max_fi, float(max_gi)))
    else:
        print("max_fi = {:6d}, log_2 max_gi = {:6.2f}".format(max_fi, float(log(max_gi,2))))
    print("aut = {}".format(aut))
    print("alpha_f = {:.4f} alpha_g = {:.4f}, sum = {:.4f}".format(alpha_f,alpha_g, alpha_f+alpha_g))
    simul = Simulation_NFS(p,ell,Fp,Fpz,sieving_dim,f,g,ZZx,cost,aut,count_sieving=True,alpha_f=alpha_f, alpha_g=alpha_g)
    simul.print_params()
    simul.simulation(samples=samples)
    simul.print_results()
    print("#::::::::::::::")
    # if there is not enough relations of there are too many relations, re-run with the same polynomials but with a higher/smaller cost
    simul.adjust_cost(samples)
    print("###############")

"""
sort results
i=0
file=test_MNT_6_with_c_NFS_Conj_${i}.res ; grep $file -e '^cost' | sort -u ; cost=`grep $file -e '^cost' | sort -u | head -1 | cut -d ' ' -f3` ; cmd="grep $file -e '^cost' -e 'LogInt' -e '\* proba' | grep -e 'cost = $cost' -A 2 | sort -u" ; echo $cmd


"""
