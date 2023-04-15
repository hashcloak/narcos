from sage.all_cmdline import *   # import sage library

import sage
import tnfs
from tnfs.curve.kss16 import KSS16
from tnfs.simul.polyselect import Polyselect
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d
from tnfs.simul.polyselect_utils import automorphism_factor
from tnfs.simul.polyselect_utils import pretty_print_coeffs_from_coeffs, pretty_print_poly_from_coeffs
from tnfs.param.testvector_sparseseed import test_vector_kss16

# create curve
k=16
KSS16_data = test_vector_kss16[0]
seed = KSS16_data['u'] #seed = -2**34+2**30-2**20-2**17-2**14+2**12+1
label = KSS16_data['label']
deg_h = KSS16_data['deg_h_S']
deg_g = k // deg_h
cost =  KSS16_data['cost_S']
# For KSS16 curves, actually Px(x-1) has smaller coefficients
# this is Px_Special
px, rx, tx, cx, yx, betax, lambx, D = tnfs.curve.kss16.polynomial_params()
x = px.parent().gens()[0]
px_special = px(x-1)
den = lcm([pi.denom() for pi in px_special.list()])
Px_Special = den*px_special

E = KSS16(u=seed, a=KSS16_data['a'])
print(label)
print(E)
E.print_parameters()
# create an instance of Polyslect class
poly_init = Polyselect(E, deg_h=deg_h)
# this takes some time because it tests a list of 880 polynomials for irreducibility mod p
#poly_init.compute_h(deg_h=deg_h)
#list_h = poly_init.get_h()[deg_h]
list_h = [
    (0.94663020, 2, [-1, 0, 0,-1, 0,-1, 0, 1, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^7-y^5-y^3-1
    (0.93916420, 2, [-1, 0, 0, 0, 0,-1, 0, 0, 1, 1, 0, 0, 0, 0, 0,0,1]), #y^16+y^9+y^8-y^5-1
    (0.93794293, 2, [-1,-1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,0,1]), #y^16+y^10+y^9-y-1
    (0.93684499, 2, [-1, 0, 0, 1, 0, 0, 0,-1, 0, 0, 0, 0, 0, 1, 0,0,1]), #y^16+y^13-y^7+y^3-1
    (0.93014791, 2, [-1,-1, 0,-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^7-y^3-y-1
    (0.92823560, 2, [-1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,0,1]), #y^16+y^11-y^3+y^2-1
    (0.92703257, 2, [-1, 1,-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^6-y^2+y-1
    (0.91709528, 2, [-1, 1, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,1]), #y^16+y^15-y^5+y-1
    (0.91556930, 2, [-1,-1, 0, 0, 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0,0,1]), #y^16-y^10+y^7-y-1
    (0.91406871, 2, [ 1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,0,1]), #y^16+y^9-y^2-y+1
    (0.91331619, 2, [-1, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 1, 0, 0, 0,0,1]), #y^16+y^11-y^7+y^5-1
    (0.91255236, 2, [ 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,0,1]), #y^16+y^14+y^3+y+1
    (0.90281171, 2, [-1, 0, 1,-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^5-y^3+y^2-1
    (0.89887313, 2, [-1, 0, 1, 0, 0, 0, 0, 0, 0,-1, 0, 1, 0, 0, 0,0,1]), #y^16+y^11-y^9+y^2-1
    (0.89877117, 2, [ 1,-1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1,0,1]), #y^16-y^14+y^9-y+1
    (0.89861046, 2, [-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,1]), #y^16+y^13+y^2-y-1
    (0.89749057, 2, [-1, 0,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,1]), #y^16+y^13-y^5-y^2-1
    (0.89603531, 2, [ 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0,-1,0,1]), #y^16-y^14-y^10+y+1
    (0.89486638, 2, [-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,0,1]), #y^16-y^12+y^11-y^3-1
    (0.89031776, 2, [-1, 0, 0, 0, 1, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0,0,1]), #y^16-y^10+y^7+y^4-1
    (0.88924763, 2, [ 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0,0,1]), #y^16+y^9+y^7+y^6+1
    (0.87854237, 2, [ 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,0,1]), #y^16+y^11+y^9+y^4+1
    (0.87774978, 2, [ 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,0,1]), #y^16+y^13+y^6+y^3+1
    (0.87336442, 2, [ 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,1]), #y^16+y^15-y^2-y+1
    (0.86527811, 2, [ 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^7+y^5+y^3+1
    (0.86100439, 2, [ 1, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 1, 0, 0,0,1]), #y^16+y^12-y^6+y+1
    (0.85978406, 2, [-1, 1, 0, 0, 0, 0,-1, 0, 1, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^8-y^6+y-1
    (0.85711700, 2, [ 1, 0,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16-y^4+y^3-y^2+1
    (0.85692724, 2, [-1, 0, 0, 1, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16-y^6+y^4+y^3-1
    (0.85596371, 2, [-1, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,1]), #y^16+y^13-y^4+y^2-1
    (0.85386879, 2, [ 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^5+y^4+y+1
    (0.85301754, 2, [ 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0,0,1]), #y^16-y^12+y^7+y^2+1
    (0.84860823, 2, [ 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^5+y^3+y+1
    (0.84435402, 2, [ 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0,-1,0,1]), #y^16-y^14-y^12+y+1
    (0.84313754, 2, [ 1, 1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0,0,1]), #y^16-y^8+y^2+y+1
    (0.84260145, 2, [ 1, 0, 0, 0, 0, 1, 0, 0, 1, 0,-1, 0, 0, 0, 0,0,1]), #y^16-y^10+y^8+y^5+1
    (0.84209705, 2, [ 1, 0, 0, 1, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0,0,1]), #y^16-y^10+y^5+y^3+1
    (0.84205437, 2, [-1, 0, 0,-1, 0,-1, 0, 0, 0, 0, 0, 1, 0, 0, 0,0,1]), #y^16+y^11-y^5-y^3-1
    (0.84155297, 2, [-1, 0, 0, 1, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0,0,1]), #y^16-y^8+y^6+y^3-1
    (0.83947403, 2, [-1,-1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,0,1]), #y^16+y^8+y^7-y-1
    (0.83834228, 2, [-1, 0, 0, 0, 0, 1, 0, 0,-1, 0, 1, 0, 0, 0, 0,0,1]), #y^16+y^10-y^8+y^5-1
    (0.83425256, 2, [-1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,1]), #y^16+y^12+y^3+y-1
    (0.82879405, 2, [ 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,-1, 0, 1, 0, 0,0,1]), #y^16+y^12-y^10+y^3+1
    (0.82455662, 2, [-1, 0, 0, 0, 0, 1, 1, 0,-1, 0, 0, 0, 0, 0, 0,0,1]), #y^16-y^8+y^6+y^5-1
    (0.82139783, 2, [ 1,-1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0,0,1]), #y^16-y^6+y^5-y+1
    (0.81758989, 2, [ 1, 0, 0, 0, 0, 1, 1, 0,-1, 0, 0, 0, 0, 0, 0,0,1]), #y^16-y^8+y^6+y^5+1
    (0.80730150, 2, [ 1, 0,-1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,0,1]), #y^16+y^9+y^4-y^2+1
    (0.80555197, 2, [-1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,0,1]), #y^16+y^14+y^6+y-1
]
for item in list_h:
    inv_zeta_Kh, w, hc = item
    hc_str = pretty_print_coeffs_from_coeffs(hc)
    h_str = pretty_print_poly_from_coeffs(hc)
    print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
print("")

for item in [list_h[i] for i in [14,39,25]]: # the best polynomials are at index 14, 39, 25
    inv_zeta_Kh, w, hc = item
    hc_str = pretty_print_coeffs_from_coeffs(hc)
    h_str = pretty_print_poly_from_coeffs(hc)
    print("    ({0:.8f},{1:2d}, {2}), #{3}".format(float(inv_zeta_Kh), w, hc_str, h_str))
    ZZy = ZZ['y']; (y,) = ZZy._first_ngens(1)
    h = ZZy(hc)
    # Special Joux-Pierrot construction
    # deg_g = k/deg_h = 1
    # since gcd(deg_h, deg_g) = 1, g should not have algebraic coefficients: with_y=False
    res_poly = poly_init.TNFS_Special(deg_g=deg_g, h=h, poly_p=Px_Special, u=seed+1, with_y=False)
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
    print("############")
