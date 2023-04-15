"""
Simulation of the (Special) Number Field Sieve
No tower

Author: Aurore Guillevic, Inria Nancy
File created: May 22, 2019

Method: a simulation of the elements in the relation collection
Monte-Carlo simulation with at least 10^5 samples taken at random
Computation of the error

"""

import sage
from sage.libs.pari import pari

from sys import version_info
if version_info[0] < 3:
    from exceptions import ValueError
#from operator import itemgetter
from sage.functions.log import log
from sage.functions.other import ceil, floor, sqrt
from sage.arith.misc import GCD, gcd
from sage.rings.integer import Integer
from sage.rings.integer_ring import Z, ZZ
from sage.rings.real_mpfr import RR
from sage.misc.prandom import randint

from sage.functions.transcendental import dickman_rho, zeta
from sage.symbolic.constants import euler_gamma
from sage.functions.exp_integral import li, log_integral

from sage.rings.number_field.number_field import NumberField
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from tnfs.simul.polyselect import Polyselect
from tnfs.simul.simulation_tnfs import time_linalg, find_log2_B, random_vec, random_vec_positive_lc, smoothness_proba, lg2_smoothness_proba

default_weight_per_row=200
# min and max according to appendix A in https://eprint.iacr.org/2019/431
cst_filtering_min = 9
cst_filtering_max = 382
default_cst_filtering = 20

def volume_sieving_space_NFS(A,sieving_dim):
    """
    Number of elements a_0 + a_1*x + ... + a_d*x^d where d = sieving_dim-1 and -A <= a_i <= A, and strictly positive leading term 0 < a_d <= A.

    :param A: bound (included) on the coefficients: -A <= a_i <= A
    :param sieving_dim: Elements sieved have sieving_dim coefficients bounded by A.
    :returns: (2*A+1)^(sieving_dim-1)*A
    """
    return (2*A+1)**(sieving_dim)*A

def lg2_volume_sieving_space_NFS(A,sieving_dim):
    """
    Number of elements a_0 + a_1*x + ... + a_d*x^d where d = sieving_dim-1 and -A <= a_i <= A, and strictly positive leading term 0 < a_d <= A.
    :param A: bound (included) on the coefficients: -A <= a_i <= A
    :param sieving_dim: Elements sieved have sieving_dim coefficients bounded by A.
    :returns: log_2((2*A+1)^(sieving_dim-1)*A)
    """
    return float((sieving_dim-1)*log(2*A+1,2) + log(A,2))

def core_volume_sieving_space_NFS(A, sieving_dim, inv_zeta=None):
    """
    :param A: bound (included) on the coefficients -A <= a_i <= A 
    :param sieving_dim: sieving dimension
    :returns: volume_sieving_space_NFS(A,sieving_dim)/zeta(sieving_dim)
    """
    if inv_zeta == None:
        return RR(volume_sieving_space_NFS(A,sieving_dim)/zeta(RR(sieving_dim)))
    else:
        return RR(volume_sieving_space_NFS(A,sieving_dim)/inv_zeta)

def lg2_core_volume_sieving_space_NFS(A, sieving_dim, inv_zeta=None):
    """
    :param A: bound (included) on the coefficients -A <= a_i <= A 
    :param sieving_dim: sieving dimension
    :returns: log_2(volume_sieving_space_NFS(A,sieving_dim)/zeta(sieving_dim))
    """
    if inv_zeta == None:
        return float(lg2_volume_sieving_space_NFS(A,sieving_dim) - log(zeta(RR(sieving_dim)),2))
    else:
        return float(lg2_volume_sieving_space_NFS(A,sieving_dim) - log(inv_zeta,2))

def find_A_NFS(sieving_dim, sieving_cost, start_A,log2_B=0):
    """ 
    Binary search on the value of A such that given h and an expected sieving 
    cost, the sieving volume will be minimal.
    :param deg_h: degree of polynomial h. Elements are considered modulo h. 
    """
    print("find_A, log2_B={}".format(log2_B))
    A = start_A
    v = lg2_volume_sieving_space_NFS(A, sieving_dim)
    if log2_B > 1:
        corr_log2_B = float(log(log(log2_B*log(2)),2))
    else:
        corr_log2_B = float(0.0)
    s = float(v + corr_log2_B)
    A_min = A
    A_max = A
    while s < sieving_cost: # increase A
        A = max(A+1, floor(1.1*A_min)) # if A < 9 it will loop if there is no max
        v = lg2_volume_sieving_space_NFS(A, sieving_dim)
        s = float(v + corr_log2_B)
        if s < sieving_cost:
            A_min = A
        else:
            A_max = A
    while s > sieving_cost: # decrease A
        A = min(A-1, floor(0.9*A_max))
        v = lg2_volume_sieving_space_NFS(A, sieving_dim)
        s = float(v + corr_log2_B)
        if s > sieving_cost:
            A_max = A
        else:
            A_min = A
    while ((A_max - A_min) > 1):
        A_mid = (A_min + A_max)//2
        v = lg2_volume_sieving_space_NFS(A_mid, sieving_dim)
        s = float(v + corr_log2_B)
        if s < sieving_cost:
            A_min = A_mid
        else:
            A_max = A_mid
    v1 = lg2_volume_sieving_space_NFS(A_min, sieving_dim)
    s1 = float(v1 + corr_log2_B)
    v2 = lg2_volume_sieving_space_NFS(A_max, sieving_dim)
    s2 = float(v2 + corr_log2_B)
    return A_min, s1, A_max, s2

def get_parameters_A_log2B_NFS(cost, ell_wordsize, aut, sieving_dim, cst_filtering=default_cst_filtering, weight_per_row=default_weight_per_row, count_sieving_cost=True):
    """ total expected cost in bits """
    cost = float(cost)
    prec = 0.0001
    start_log2B = float(cost/2.0 - log(cost/2.0))
    # time_linalg is a function
    # whose optional parameters are
    # cst_filtering, weight_per_row
    # maybe pu these functions inside the class...
    B_min, t1, B_max, t2 = find_log2_B(cost-1, start_log2B, prec, time_linalg, ell_wordsize, cst_filtering=cst_filtering, weight_per_row=weight_per_row)
    log2_B = B_max
    start_A = floor(2**(RR((cost-1)/sieving_dim)))
    if count_sieving_cost:
        A_min, s1, A_max, s2 = find_A_NFS(sieving_dim, float(cost-1+log(aut,2)), start_A, log2_B=log2_B)
    else:
        A_min, s1, A_max, s2 = find_A_NFS(sieving_dim, float(cost-1+log(aut,2)), start_A)
    # we should have s1 <= cost <= s2 and A_min + 1 = A_max
    if abs(cost-1+log(aut,2)-s1) < 0.5 and abs(cost-1+log(aut,2)-s2) >= 0.5:
        A = A_min
    else:
        A = A_max
    return A, log2_B




def statistics_NFS(f,g,sieving_dim,Rx,A,logB,samples,alpha_f=0.0,alpha_g=0.0):
    print("f={}\ng={}\nsieving_dim={}\nRx={}\nA={}\nlogB={}\n".format(f,g,sieving_dim,Rx,A,logB))
    
    log_Nf =[0.0 for i in range(samples)]
    log_Ng =[0.0 for i in range(samples)]
    proba_f = [0.0 for i in range(samples)]
    proba_g = [0.0 for i in range(samples)]

    duplicates = 0
    for i in range(samples):
        a=random_vec_positive_lc(A,sieving_dim)
        #a=random_vec(A,sieving_dim)
        #print("a = {}".format(a))
        while not gcd(a) == 1:
            a=random_vec(A,sieving_dim)
            duplicates += 1
        ax = Rx(a)
        log_Nf[i] = log(RR(abs(f.resultant(ax))))
        log_Ng[i] = log(RR(abs(g.resultant(ax))))
        proba_f[i] = smoothness_proba(log_Nf[i], logB, alpha_f)
        proba_g[i] = smoothness_proba(log_Ng[i], logB, alpha_g)

    # end for loop -> now there are norms for all samples
    # average
    avrg_log_Nf = RR(sum(log_Nf) / samples)
    avrg_log_Ng = RR(sum(log_Ng) / samples)
    avrg_proba_f = RR(sum(proba_f) / samples)
    avrg_proba_g = RR(sum(proba_g) / samples)
    prod_proba = [RR(proba_f[i]*proba_g[i]) for i in range(samples)]
    avrg_proba = RR(sum(prod_proba)/samples)
    # variance: sum square of diff between values and mean value    
    variance_norm_f = RR(sum([(log_Nf[i]-avrg_log_Nf)**2 for i in range(samples)]))
    variance_norm_g = RR(sum([(log_Ng[i]-avrg_log_Ng)**2 for i in range(samples)]))
    variance_proba_f = RR(sum([(proba_f[i]-avrg_proba_f)**2 for i in range(samples)]))
    variance_proba_g = RR(sum([(proba_g[i]-avrg_proba_g)**2 for i in range(samples)]))
    variance_proba = RR(sum([(prod_proba[i]-avrg_proba)**2 for i in range(samples)]))
    # standard deviation (ecart type)
    std_dev_norm_f = RR(sqrt(variance_norm_f/(samples-1)))
    std_dev_norm_g = RR(sqrt(variance_norm_g/(samples-1)))
    std_dev_proba_f = RR(sqrt(variance_proba_f/(samples-1)))
    std_dev_proba_g = RR(sqrt(variance_proba_g/(samples-1)))
    std_dev_proba = RR(sqrt(variance_proba/(samples-1)))
    
    # error
    err_norm_f = float(100*std_dev_norm_f/sqrt(samples)/avrg_log_Nf)
    err_norm_g = float(100*std_dev_norm_g/sqrt(samples)/avrg_log_Ng)
    err_proba_f = float(100*std_dev_proba_f/sqrt(samples)/avrg_proba_f)
    err_proba_g = float(100*std_dev_proba_g/sqrt(samples)/avrg_proba_g)
    err_proba = float(100*std_dev_proba/sqrt(samples)/avrg_proba)
    
    # cost
    avrg_log_Nf = float(avrg_log_Nf)
    avrg_log_Ng = float(avrg_log_Ng)
    avrg_proba_f = float(avrg_proba_f)
    avrg_proba_g = float(avrg_proba_g)
    avrg_proba = float(avrg_proba)
    
    std_dev_norm_f = float(std_dev_norm_f)
    std_dev_norm_g = float(std_dev_norm_g)
    std_dev_proba_f = float(std_dev_proba_f)
    std_dev_proba_g = float(std_dev_proba_g)
    std_dev_proba = float(std_dev_proba)
    
    return avrg_log_Nf, avrg_log_Ng, avrg_proba_f, avrg_proba_g, avrg_proba, err_norm_f, err_norm_g, err_proba_f, err_proba_g, err_proba, std_dev_norm_f, std_dev_norm_g, std_dev_proba_f, std_dev_proba_g, std_dev_proba, float(RR(RR(samples)/RR(samples+duplicates))), min(log_Nf), max(log_Nf), min(log_Ng), max(log_Ng)

class Simulation_NFS():

    def __init__(self, p, ell, Fp, Fpz, sieving_dim, f, g, Rx, log2_cost, automorphism,
                 alpha_f=None, alpha_g=None,
                 weight_per_row=None, cst_filtering=None,
                 count_sieving=None, label=None):
        """
        :param              p: prime integer defining a prime field
        :param            ell: prime integer, order of the mult. subgroup where to compute DL
                               (linear algebra is done modulo ell)
        :param             Fp: finite field of characteristic p
        :param            Fpz: polynomial ring of Fp
        :param    sieving_dim: in case of NFS-HD, sieve over a=a0+a1*x + ... + a_{dim-1}*x^(dim-1)
        :param              f: univariate irreducible polynomial
        :param              g: univariate irreducible polynomial
        :param             Rx: univariate polynomial ring over Q (or Z)
        :param      log2_cost: the cost in basis 2 of (S)NFS
        :param   automorphism: number of relations obtained for free thanks to automorphisms of Kf, Kg
                                (these are not extra rels, this is only a time speed-up)
        :param        alpha_f: alpha value of f
        :param        alpha_g: alpha value of g
        :param weight_per_row: default 200 (weight per row of the matrix)
        :param  cst_filtering: default 20, reduction factor of the matrix
        :param  count_sieving: True/False, whether to count the cost of sieving or not
        :param          label: for printing (such as name of the pairing-friendly curve)

        put somewhere else:
        :param  sieving_dim: sieving dimension, between 2 and min(deg(f), deg(g))
        :param            A: bound on the coefficients a_i in the relation collection -A <= a_i <= A
        :param       log2_B: B is the smoothness bound (relations involve primes <= B)
        :param      samples: number of random coprime tuples in the relation collection to compute an average smoothness probability, uses at least 10^6 (takes ~20 min to run)
        """
        self.p = p     #
        self.ell = ell
        self.ell_wordsize = ceil(RR(log(ell,2))/64)
        self.Fp = Fp   # there 3 params are in classes PFCurve (BN, BLS12, KSS18, BLS24)
        self.Fpz = Fpz #
        self.sieving_dim = sieving_dim
        self.f = f
        self.g = g
        self.Rx=Rx
        self.aut = automorphism
        self.inv_zeta = float(RR(1/zeta(sieving_dim)))
        self.list_cost_tested = list()
        
        if alpha_f != None and alpha_g != None:
            self.alpha_f = alpha_f
            self.alpha_g = alpha_g
            self.with_alpha = True
        else:
            self.alpha_f = 0.0
            self.alpha_g = 0.0
            self.with_alpha = False

        if weight_per_row == None or weight_per_row <= 0:
            self.weight_per_row = default_weight_per_row
        else:
            self.weight_per_row = weight_per_row
        
        if cst_filtering == None or cst_filtering <= 0:
            self.cst_filtering = default_cst_filtering
        else:
            self.cst_filtering = cst_filtering
        if count_sieving == None:
            self.count_sieving = True
        else:
            self.count_sieving = count_sieving
        self.label = label

        self.set_cost(log2_cost)

    def set_cost(self, cost):
        if cost <= 0:
            raise ValueError("Error cost: should be a strictly positive number but received {}".format(cost))        
        self.cost = cost
        self.A, self.log2_B = get_parameters_A_log2B_NFS(self.cost, self.ell_wordsize, self.aut, self.sieving_dim, self.cst_filtering, self.weight_per_row, count_sieving_cost=self.count_sieving)
        self.log2_core_vol_sieving_space = lg2_core_volume_sieving_space_NFS(self.A, self.sieving_dim, inv_zeta=self.inv_zeta)
        self.log2_vol_sieving_space = lg2_volume_sieving_space_NFS(self.A, self.sieving_dim)
        self.logB = float(self.log2_B*log(2))
        self.loglogB = float(log(self.logB)) # log(log(B)) where log=ln is natural logarithm (log(e) = 1)
        self.log2_time_relcol = float(self.log2_vol_sieving_space - log(self.aut,2))
        self.log2_time_relcol_sieving_cost = self.log2_time_relcol + float(log(self.loglogB,2))

        self.log2_time_linalg, self.log2_size_FB, self.log2_size_matrix, self.logB = time_linalg(self.log2_B, self.ell_wordsize)
        self.log2_total_time = float(log(2**self.log2_time_relcol + 2**self.log2_time_linalg ,2))
        self.log2_total_time_sieving_cost = float(log(2**self.log2_time_relcol_sieving_cost + 2**self.log2_time_linalg ,2))

        self.log2_size_FB_BD = float(1+self.log2_B - log(self.logB,2))
        self.log2_size_matrix_BD = float(self.log2_size_FB_BD - log(self.log2_B,2))
        self.log2_time_linalg_BD = float(5 + 2*self.log2_size_matrix_BD)
        self.log2_total_time_BD = float(log(2**self.log2_time_relcol + 2**self.log2_time_linalg_BD ,2))

    def print_params(self):
        print("cost = {} bits, sieving_dim = {}, A = {} = 2^{:.2f}, ell_wordsize = {}, B = 2^{:.3f} \n".format(self.cost, self.sieving_dim, self.A, float(log(self.A,2)), self.ell_wordsize, float(self.log2_B)))
        print("vol sieving space = 2^{0:.2f}, core vol sieving = 2^{1:.2f}, 1/zeta(2) = {2:.6f}".format(float(self.log2_vol_sieving_space), float(self.log2_core_vol_sieving_space), float(self.inv_zeta)))
        if self.with_alpha:
            print("alpha_f ={0:7.4f} (basis e), alpha_g ={1:7.4f} (basis e), sum ={2:7.4f}".format(float(self.alpha_f), float(self.alpha_g), float(self.alpha_f+self.alpha_g)))
            print("alpha_f ={0:7.4f} (basis 2), alpha_g ={1:7.4f} (basis 2), sum ={2:7.4f}".format(float(self.alpha_f/log(2)), float(self.alpha_g/log(2)), float(self.alpha_f/log(2)+self.alpha_g/log(2))))

    def simulation(self, samples):
        self.samples = samples
        self.avrg_log_Nf, self.avrg_log_Ng, self.avrg_smooth_proba_f, self.avrg_smooth_proba_g, \
            self.avrg_proba, self.err_norm_f, self.err_norm_g, self.err_proba_f, self.err_proba_g, self.err_proba, \
            self.std_dev_norm_f, self.std_dev_norm_g, self.std_dev_proba_f, self.std_dev_proba_g, self.std_dev_proba, \
            self.inv_zeta_exp, self.min_log_Nf, self.max_log_Nf, self.min_log_Ng, self.max_log_Ng = \
                statistics_NFS(self.f,self.g,self.sieving_dim,self.Rx,self.A,self.logB,self.samples,self.alpha_f,self.alpha_g)
        self.log2_proba_rough = float(lg2_smoothness_proba(self.avrg_log_Nf,self.logB,self.alpha_f) + lg2_smoothness_proba(self.avrg_log_Ng,self.logB,self.alpha_g))
        self.log2_avrg_proba = float(log(self.avrg_proba,2))
        self.log2_rels = float(self.log2_core_vol_sieving_space + self.log2_avrg_proba)
        self.log2_rels_rough = float(self.log2_core_vol_sieving_space + self.log2_proba_rough)


    def print_results(self):
        print("experimental 1/zeta(2) = {}".format(self.inv_zeta_exp))
        print("theoretical  1/zeta(2) = {}".format(self.inv_zeta))
        print("err_norm_f  = {}\nerr_norm_g  = {}\nerr_proba_f = {}\nerr_proba_g = {}\nerr_proba   = {}\n".format(self.err_norm_f, self.err_norm_g, self.err_proba_f, self.err_proba_g, self.err_proba))
        print("std_dev_norm_f = {}\nstd_dev_norm_g = {}\nstd_dev_proba_f = {}\nstd_dev_proba_g = {}\nstd_dev_proba   = {}\n".format(self.std_dev_norm_f, self.std_dev_norm_g, self.std_dev_proba_f, self.std_dev_proba_g, self.std_dev_proba))
        
        print("A = {0} = 2^{1:.3f}, B = 2^{2:.3f}, sieving_dim = {3}, samples = {4}".format(self.A,float(log(self.A,2)),float(self.log2_B),self.sieving_dim,self.samples))
        print("vol_sieving_space = (2*A+1)^{}*A = 2^{:.3f}".format(self.sieving_dim-1, float(self.log2_vol_sieving_space)))
        print("avrg_log2_Nf = {0:.4f} (bits) avrg_log2_Ng = {1:.4f} (bits) sum = {1:.4f} (bits)".format(float(self.avrg_log_Nf/log(2)), float(self.avrg_log_Ng/log(2)), float((self.avrg_log_Nf+self.avrg_log_Ng)/log(2))))
        print(" min_log2_Nf = {0:.4f} (bits)  min_log2_Ng = {1:.4f} (bits)".format(float(self.min_log_Nf/log(2)), float(self.min_log_Ng/log(2))))
        print(" max_log2_Nf = {0:.4f} (bits)  max_log2_Ng = {1:.4f} (bits)".format(float(self.max_log_Nf/log(2)), float(self.max_log_Ng/log(2))))
        
        print("smooth_proba(avrg_Nf,logB,alpha_f) = smooth_proba({0:.2f},2^{1:.2f},{2:.4f}) = 2^{3:.4f}".format(float(self.avrg_log_Nf*log(2)),float(self.log2_B),float(self.alpha_f/log(2)),float(lg2_smoothness_proba(self.avrg_log_Nf,self.logB,self.alpha_f))))
        print("smooth_proba(avrg_Ng,logB,alpha_g) = smooth_proba({0:.2f},2^{1:.2f},{2:.4f}) = 2^{3:.4f}".format(float(self.avrg_log_Ng*log(2)),float(self.log2_B),float(self.alpha_g/log(2)),float(lg2_smoothness_proba(self.avrg_log_Ng,self.logB,self.alpha_g))))
        print("prod of smooth proba of avrg norms: = 2^{0:.4f}".format(float(self.log2_proba_rough)))
        print("avrg_proba = 2^{0:.4f}".format(float(self.log2_avrg_proba)))

        print("relations obtained: core_vol_sieving_space * proba = 2^{0:.4f}".format(float(self.log2_rels)))
        print("relations obtained: core_vol_sieving_space * (rough proba) = 2^{0:.4f}".format(float(self.log2_rels_rough)))
        print("#{{factor base}} =~ 2LogIntegral(B) = 2^{0:.4f}".format(self.log2_size_FB))
        print("#{{factor base}} ~ 2*B/log(B)       = 2^{0:.4f}".format(self.log2_size_FB_BD))
        print("size matrix = size factor base / {0} = 2^{1:.4f} (worst case)".format(self.cst_filtering,self.log2_size_matrix))
        if not self.count_sieving:
            print("time relcol = 2^{0:.4f}".format(self.log2_time_relcol))
        else:
            print("time relcol = 2^{0:.4f} (counting sieving cost as ln(ln(B)))".format(self.log2_time_relcol_sieving_cost))
        print("time linalg = 2^{0:.4f}".format(self.log2_time_linalg))
        if not self.count_sieving:
            print("total time  = 2^{0:.4f}".format(self.log2_total_time))
        else:
            print("total time  = 2^{0:.4f} (counting sieving cost as ln(ln(B)))".format(self.log2_total_time_sieving_cost))
        print("BD18:")
        print("time linalg = 2^{0:.4f}".format(self.log2_time_linalg_BD))
        print("total time  = 2^{0:.4f}".format(self.log2_total_time_BD))
        
    def adjust_cost(self,samples):
        self.list_cost_tested.append(self.cost)
        while (self.log2_rels < self.log2_size_FB) or (self.log2_rels > 0.6+self.log2_size_FB):
            cost_old = self.cost
            if self.log2_rels < self.log2_size_FB:
                # compute adjusted cost at least for relation collection
                if self.log2_rels >= self.log2_size_FB - 0.8:
                    cost = self.cost+1
                else:
                    log2_cost_sieving = RR(log(self.loglogB,2))
                    log2_minimal_vol_sieving_space = RR(self.log2_size_FB - log(self.inv_zeta,2) - log(self.avrg_proba,2))
                    log2_minimal_time_relcol = RR(log2_minimal_vol_sieving_space - RR(log(self.aut,2)) + log2_cost_sieving)
                    cost = ceil(log2_minimal_time_relcol+1)
                if cost == cost_old:
                    cost = cost+1
                if cost in self.list_cost_tested:
                    break # there is a loop here
                print("not enough rels with cost={}, re-run with cost={}".format(cost_old, cost))
            else: #elif self.log2_rels > 0.6 + self.log2_size_FB:
                if self.log2_rels <= self.log2_size_FB + 1.5:
                    cost = self.cost-1
                else:
                    log2_cost_sieving = RR(log(self.loglogB,2))
                    log2_minimal_vol_sieving_space = RR(self.log2_size_FB - log(self.inv_zeta,2) - log(self.avrg_proba,2))
                    log2_minimal_time_relcol = RR(log2_minimal_vol_sieving_space - RR(log(self.aut,2)) + log2_cost_sieving)
                    cost = ceil(log2_minimal_time_relcol+1)
                if cost == cost_old:
                    cost = cost-1
                if cost in self.list_cost_tested:
                    break # there is a loop here
                print("too many rels with cost={}, re-run with cost={}".format(cost_old, cost))
            self.list_cost_tested.append(cost)
            self.set_cost(cost)
            self.print_params()
            self.simulation(samples=samples)
            self.print_results()
            print("#::::::::::::::")
