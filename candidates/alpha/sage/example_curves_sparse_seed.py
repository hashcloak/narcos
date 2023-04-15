from sage.all_cmdline import *   # import sage library

import sage
import tnfs
from tnfs.curve.bn import BN
from tnfs.curve.bls import BLS
from tnfs.curve.kss import KSS
from tnfs.curve.aurifeuillean import Aurifeuillean
from tnfs.curve.fotiadis_martindale import FotiadisMartindale
from tnfs.curve.pairing_friendly_curve import BrezingWeng


from tnfs.simul.polyselect import Polyselect
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d
from tnfs.simul.polyselect_utils import automorphism_factor
from tnfs.simul.polyselect_utils import pretty_print_coeffs_from_coeffs, pretty_print_poly_from_coeffs
from tnfs.param.testvector_sparseseed import test_vector_bn
from tnfs.param.testvector_sparseseed import test_vector_bls12, test_vector_bls24, test_vector_bls24_cln, test_vector_bls48
from tnfs.param.testvector_sparseseed import test_vector_fm12, test_vector_fk12d3, test_vector_fm12d3
from tnfs.param.testvector_sparseseed import test_vector_kss16, test_vector_kss18, test_vector_kss32, test_vector_kss36, test_vector_kss40, test_vector_kss54

# sage example_curves_sparse_seed.py --bn --idx 12 --samples 100000
# sage example_curves_sparse_seed.py --bls -k 12 --idx 0 --samples 100000 ### the curve BLS12-377
# sage example_curves_sparse_seed.py --bls -k 12 --idx 1 --samples 100000 ### the curve BLS12-381
# sage example_curves_sparse_seed.py --bls -k 24 --idx 5 --samples 100000

allowed_k_bn = [12]
allowed_k_bls = [12,24,48]
allowed_k_kss = [16,18,32,36,40,54]
allowed_k_fk = [12]
allowed_k_fm = [12]

curve = None # KSS or BLS or BN or Aurifeuillean
k = None
D = None
e0 = None
a_auri = None
idx = 0
code = None

args=sys.argv
filename_stem = '.'.join((args[0].split('/')[-1]).split('.')[:-1])
print("usage: sage {} [--kss | --bls | --bn | --auri | --fm | --fk] [-k <k>] --idx <i> [-D <D>] [-e0 <e0>] [-a <auri_a>] [--code <fk or fm code>] [--samples <samples>] [--special|--conj|--sarkarsingh|--glv] [--nfs|--tnfs]".format(filename_stem))
i=1
cmd_options = ["--kss", "--bls", "--bn", "--auri", "--fm", "--fk", "-k", "--idx", "-D", "-e0", "-a"]

while i < len(args):
    if args[i] == "--bls":
        curve = BLS; str_curve = "BLS"; str_curve_lowcase = "bls"
        curve_path = tnfs.curve.bls
        allowed_k = allowed_k_bls
    elif args[i] == "--kss":
        curve = KSS; str_curve = "KSS"; str_curve_lowcase = "kss"
        curve_path = tnfs.curve.kss
        allowed_k = allowed_k_kss
    elif args[i] == "--bn":
        curve = BN; str_curve = "BN"; str_curve_lowcase = "bn"
        curve_path = tnfs.curve.bn
        allowed_k = allowed_k_bn
        k = 12
    elif args[i] == "--auri":
        curve = Aurifeuillean; str_curve = "Auri"; str_curve_lowcase = "auri"
        curve_path = tnfs.curve.aurifeuillean
        allowed_k = allowed_k_auri
    elif args[i] == "--fk":
        curve = FotiadisMartindale; str_curve = "FK"; str_curve_lowcase = "fk"
        curve_path = tnfs.curve.fotiadis_martindale
        allowed_k = allowed_k_fk
    #elif args[i] == "--fm":
    #    curve = FM; str_curve = "FM"; str_curve_lowcase = "fm"
    #    curve_path = tnfs.curve.fotiadis_martindale
    #    allowed_k = allowed_k_fm
    #    k = 12
    elif args[i] == "-k" and (i+1) < len(args) and args[i+1] not in cmd_options:
        k = Integer(args[i+1])
        i += 1
    elif args[i] == "--idx" and (i+1) < len(args) and args[i+1] not in cmd_options:
        idx = Integer(args[i+1])
        i += 1
    elif args[i] == "-D" and (i+1) < len(args) and args[i+1] not in cmd_options:
        D = Integer(args[i+1])
        i += 1
    elif args[i] == "-e0" and (i+1) < len(args) and args[i+1] not in cmd_options:
        e0 = Integer(args[i+1])
        i += 1
    elif args[i] == "-a" and (i+1) < len(args) and args[i+1] not in cmd_options:
        a_auri = Integer(args[i+1])
        i += 1
    elif args[i] == "--code" and (i+1) < len(args) and args[i+1] not in cmd_options:
        code = Integer(args[i+1])
        i += 1
    i += 1

if curve is None:
    raise ValueError("missing curve type: one of --kss or --bls or --bn or --auri")

if k is None:
    raise ValueError("missing embedding degree k, possible values with {}: {}".format(str_curve, allowed_k))
elif k not in allowed_k:
    raise ValueError("embedding degree k={} but possible values with {} are {}".format(k, str_curve, allowed_k))

proof.arithmetic(False)

QQx = QQ['x']; (x,) = QQx._first_ngens(1)
if curve == Aurifeuillean:
    px, rx, tx, cx, yx, betax, lambx, D, a_auri, exp_tr, automorphisms, m, u_m, cofactor_r = curve_path.polynomial_params(k, D, a_auri, e0)
else:
    px, rx, tx, cx, yx, betax, lambx, D = curve_path.polynomial_params(k)

px_denom = lcm([ci.denom() for ci in px.list()])
if px_denom > 1:
    print("px = ({})/{}".format(px*px_denom, px_denom))
else:
    print("px = {}".format(px))

if curve is BN:
    test_vector = test_vector_bn
elif curve is BLS:
    if k == 12:
        test_vector = test_vector_bls12
    elif k == 24:
        test_vector = test_vector_bls24+test_vector_bls24_cln
    elif k == 48:
        test_vector = test_vector_bls48
elif curve is KSS:
    if k == 16:
        test_vector = test_vector_kss16
    elif k == 18:
        test_vector = test_vector_kss18
    elif k == 32:
        test_vector = test_vector_kss32
    elif k == 36:
        test_vector = test_vector_kss36
    elif k == 40:
        test_vector = test_vector_kss40
    elif k == 54:
        test_vector = test_vector_kss54
elif curve is FotiadisMartindale:
    test_vector = test_vector_fm12d3
elif curve is FotiadisKonstantinou:
        test_vector = test_vector_fm12+test_vector_fk12d3

tv = test_vector[idx]

seed = ZZ(tv['u'])
label = tv['label']
deg_h = tv['deg_h_S']
deg_g = k // deg_h
cost =  tv['cost_S']
if D == 3:
    if curve is BN:
        E = curve(u=seed, b=tv['b'])
    else:
        E = curve(u=seed, b=tv['b'], k=k)
elif D == 1 or D==4:
    if curve is KSS:
        E = curve(k=k, u=seed, a=tv['a'])
    E = curve(u=seed, a=tv['a'], k=k)
else:
    E = curve(u=seed, a=tv['a'], b=tv['b'], k=k)

print(label)
print(E)
E.print_parameters()

# create an instance of Polyslect class
poly_init = Polyselect(E, deg_h=deg_h)
poly_init.compute_h(deg_h=deg_h)
list_h = poly_init.get_h()[deg_h]
for item in list_h:
    inv_zeta_Kh, w, hc = item
    hc_str = pretty_print_coeffs_from_coeffs(hc)
    h_str = pretty_print_poly_from_coeffs(hc)
    print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
print("")

ZZy = ZZ['y']; (y,) = ZZy._first_ngens(1)
poly_p = ZZy(px_denom*px)
poly_u = [-seed, 1]

for item in list_h:
    inv_zeta_Kh, w, hc = item
    hc_str = pretty_print_coeffs_from_coeffs(hc)
    h_str = pretty_print_poly_from_coeffs(hc)
    print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
    h = ZZy(hc)
    # Special Joux-Pierrot construction
    # since deg_h = 6, then deg_g = 12/6 = 2
    # since gcd(deg_h, deg_g) = 2 > 1, g should have algebraic coefficients: with_y=True
    res_poly = poly_init.TNFS_Special(deg_g=deg_g, h=h, poly_p=poly_p, u=seed, with_y=True)
    if res_poly == None:
        raise ValueError("Error in Polynomial selection")
    f, g, max_fij, max_gij, aut_fg = res_poly
    print(f)
    print(g)
    print(h)

    p = E.p()
    ell = E.r()
    Fp = E.Fp()
    Fpz,z = E.Fpz()
    Rxy = f.parent()

    assert (ZZ(h.resultant(f.resultant(g))) % p**k) == 0

    if (gcd(deg_g,deg_h) == 1):
        aut_h = automorphism_factor(hc)
    else:
        aut_h = 1

    aut = aut_h*aut_fg

    # computing alpha, this takes at least a few seconds

    alpha_f = float(alpha_TNFS_2d(f,h,1000))
    alpha_g = float(alpha_TNFS_2d(g,h,1000))
    sum_alpha = alpha_f+alpha_g
    print("alpha_f = {:.4f} alpha_g = {:.4f} sum_alpha = {:.4f}".format(alpha_f,alpha_g,sum_alpha))
    
    #initialisation of data
    simul = Simulation_TNFS(p,ell,Fp,Fpz,h,f,g,Rxy,cost,aut,inv_zeta_Kh,count_sieving=True,alpha_f=alpha_f,alpha_g=alpha_g)

    simul.print_params()
    simul.simulation(samples=100000) #takes few seconds for 10^4, mins for 10^5, up to 20 min for 10^6
    simul.print_results()
    print("#::::::::::::::")
    # if there is not enough relations of there are too many relations, re-run with the same polynomials but with a higher/smaller cost
    simul.adjust_cost(samples=100000)
    print("############")
