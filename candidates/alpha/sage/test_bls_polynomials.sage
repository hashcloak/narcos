import tnfs.curve.bls
from tnfs.curve.bls import BLS

import tnfs.curve.cyclo_kDe
import tnfs.curve.fst66
import tnfs.curve.bls12
import tnfs.curve.bls24

from tnfs.param.testvector_sparseseed import test_vector_bls12, test_vector_bls24

ZZx.<x> = ZZ[]

def test_k_0_mod_6():
    list_k_u0 = []
    list_k_cont_rx = []
    for l in range(1, 17):
        k = 6*l
        # check that all BLS, FST66 and Cyclo_kDe(k,3,1) match
        if (k % 18) == 0:
            continue
        print("k={}".format(k))
        params_bls = tnfs.curve.bls.polynomial_params(k)
        params_fst66 = tnfs.curve.fst66.polynomial_params(k)
        params_cyclo_kde = tnfs.curve.cyclo_kDe.polynomial_params(k,3,1)

        px = params_bls[0] # has usually denominator 3
        rx = params_bls[1] # has integer coefficients
        den_p = lcm([fi.denom() for fi in px.list()])
        if (den_p % 2) == 1:
            den_p = den_p * 2
        list_u0 = []
        Rx = px.parent()
        for u0 in range(den_p):
            p0x = px(Rx([u0, den_p])) # evaluate at den_p*x + u0
            den_p0x = lcm([fi.denom() for fi in p0x.list()])
            if den_p0x == 1 and (ZZx(p0x)).content() == 1:
                list_u0.append(u0)
        for u0 in list_u0:
            p0x = ZZx(px(Rx([u0, den_p]))) # evaluate at den_p*x + u0
            r0x = ZZx(rx(Rx([u0, den_p]))) # evaluate at den_p*x + u0
            cont_p0x = p0x.content()
            cont_r0x = r0x.content()
            if cont_r0x != 1:
                list_k_cont_rx.append((k, u0, cont_r0x))
            print("u={} px({}*x+{}) = {} * ({})".format(u0, den_p, u0, cont_p0x, p0x // cont_p0x))
            print("     rx({}*x+{}) = {} * ({})".format(den_p, u0, cont_r0x, r0x // cont_r0x))
        list_k_u0.append((k,list_u0))
        
        #px, rx, tx, cx, yx, betax, lambx, D
        label = ["px", "rx", "tx", "cx", "yx", "betax", "lambx", "D"]
        for j in [0,1,2,3,7]:
            if params_bls[j] != params_fst66[j] or params_bls[j] != params_cyclo_kde[j]:
                print("Error k = {} {} bls {} fst66 {} cyclo_kde {}".format(k, label[j], params_bls[j], params_fst66[j], params_cyclo_kde[j]))
        
        # yx: can be +/-
        j=4
        if (params_bls[j] != params_fst66[j] and params_bls[j] != -params_fst66[j]) or (params_bls[j] != params_cyclo_kde[j] and params_bls[j] != -params_cyclo_kde[j]):
            print("Error k = {} {} bls {} fst66 {} cyclo_kde {}".format(k, label[j], params_bls[j], params_fst66[j], params_cyclo_kde[j]))
        
        # betax, lambx: -bx-1, -lx-1
        for j in [5,6]:
            if (params_bls[j] != params_fst66[j] and params_bls[j] != (-params_fst66[j]-1)) or (params_bls[j] != params_cyclo_kde[j] and params_bls[j] != (-params_cyclo_kde[j]-1)):
                print("Error k = {} {} bls {} fst66 {} cyclo_kde {}".format(k, label[j], params_bls[j], params_fst66[j], params_cyclo_kde[j]))

    print("px is irreducible of content 1 for these k and u0 mod 6:\n{}".format(list_k_u0))
    print("rx has non-trivial content for these k and u0 mod 6: {}".format(list_k_cont_rx))

def test_k_3_mod_6():
    list_k_u0 = []
    list_k_cont_rx = []
    for k in [_ for _ in range(6,100) if (_ % 6) == 3]: # k is multiple of 3 but not 6
        print("\nk={}".format(k))
        params_bls = tnfs.curve.bls.polynomial_params(k)
        params_fst66 = tnfs.curve.fst66.polynomial_params(k)
        params_cyclo_kde = tnfs.curve.cyclo_kDe.polynomial_params(k,3,1)

        # check which value of x mod den is needed to have integer and prime px and rx
        px = params_bls[0] # has usually denominator 3
        rx = params_bls[1] # has integer coefficients
        den_p = lcm([fi.denom() for fi in px.list()])
        if (den_p % 2) == 1:
            den_p = den_p * 2
        list_u0 = []
        Rx = px.parent()
        for u0 in range(den_p):
            p0x = px(Rx([u0, den_p])) # evaluate at den_p*x + u0
            den_p0x = lcm([fi.denom() for fi in p0x.list()])
            if den_p0x == 1 and (ZZx(p0x)).content() == 1:
                list_u0.append(u0)
        for u0 in list_u0:
            p0x = ZZx(px(Rx([u0, den_p]))) # evaluate at den_p*x + u0
            r0x = ZZx(rx(Rx([u0, den_p]))) # evaluate at den_p*x + u0
            cont_p0x = p0x.content()
            cont_r0x = r0x.content()
            if cont_r0x != 1:
                list_k_cont_rx.append((k, u0, cont_r0x))
            print("u={} px({}*x+{}) = {} * ({})".format(u0, den_p, u0, cont_p0x, p0x // cont_p0x))
            print("     rx({}*x+{}) = {} * ({})".format(den_p, u0, cont_r0x, r0x // cont_r0x))
        list_k_u0.append((k,list_u0))
        
        #px, rx, tx, cx, yx, betax, lambx, D
        label = ["px", "rx", "tx", "cx", "yx", "betax", "lambx", "D"]
        
        for j in [0,1,2,3]:
            if params_bls[j] != params_fst66[j] and params_bls[j] != params_cyclo_kde[j]:
                den_bls = lcm([fi.denom() for fi in params_bls[j].list()])
                den_fst66 = lcm([fi.denom() for fi in params_fst66[j].list()])
                den_cyclo = lcm([fi.denom() for fi in params_cyclo_kde[j].list()])
                if den_bls == 1 and den_fst66 == 1 and den_cyclo == 1:
                    print("{} bls {}\n cyclo {}\n fst66 {}".format(label[j], params_bls[j], params_cyclo_kde[j], params_fst66[j]))
                else:
                    print("{} bls ({})/{}\n cyclo ({})/{}\n fst66 ({})/{}".format(label[j], params_bls[j]*den_bls, den_bls, params_cyclo_kde[j]*den_cyclo, den_cyclo, params_fst66[j]*den_fst66, den_fst66))
            elif params_bls[j] != params_fst66[j] and params_bls[j] == params_cyclo_kde[j]:
                den_bls = lcm([fi.denom() for fi in params_bls[j].list()])
                den_fst66 = lcm([fi.denom() for fi in params_fst66[j].list()])
                if den_bls == 1 and den_fst66 == 1:
                    print("{} bls, cyclo {}\n         fst66 {}".format(label[j], params_bls[j], params_fst66[j]))
                else:
                    print("{} bls, cyclo ({})/{}\n        fst66 ({})/{}".format(label[j], params_bls[j]*den_bls, den_bls, params_fst66[j]*den_fst66, den_fst66))

    print("px is irreducible of content 1 for these k and u0 mod 6:\n{}".format(list_k_u0))
    print("rx has non-trivial content for these k and u0 mod 6: {}".format(list_k_cont_rx))

def test_init_bls12():
    e = 0
    for vec in test_vector_bls12:
        u = vec['u']
        b = vec['b']
        try:
            E = BLS(12,u,b,1)
            print(E)
        except:
            raise
        #    e += 1
        #    print(err)
        #    continue

    print("{}/{} errors".format(e, len(test_vector_bls12)))

def test_init_bls24():
    e = 0
    for vec in test_vector_bls24:
        u = vec['u']
        b = vec['b']
        try:
            E = BLS(24,u,b,1)
            print(E)
        except:
            e += 1
            continue

    print("{}/{} errors".format(e, len(test_vector_bls24)))

if __name__ == "__main__":
    test_k_0_mod_6()
    test_k_3_mod_6()
    #test_init_bls12()
    #test_init_bls24()
