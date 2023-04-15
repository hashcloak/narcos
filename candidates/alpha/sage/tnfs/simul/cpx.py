from sage.functions.log import log
from sage.rings.integer import Integer

# classical L_N(alpha,c) function
def L_N_alpha_c(log2N, alpha, c):
    lnN = float(log(2))*float(log2N)
    return float(c) * float(log2N)**float(alpha) * float(log(lnN,2))**(float(1-alpha))

def guess_cost_STNFS_BLS(k, pknbits):
    alpha = float(Integer(1)/Integer(3))
    #c=(32/9)**(1/3)
    c = float((Integer(32)/Integer(9))**(Integer(1)/Integer(3)))
    cost_S = float(L_N_alpha_c(pknbits, alpha, c))
    if k==6:
        # this is L_{p^k}(1/3, (32/9)^(1/3)) and deg_h = 3
        deg_h_S = 3
    elif k==9:
        deg_h_S = 9
        cost_S -= 20
    elif k==12:
        deg_h_S = 6
    elif k==15:
        deg_h_S = 5
    elif k==27:
        deg_h_S = 27
    elif k==24:
        deg_h_S = 24
    return deg_h_S, cost_S

def guess_cost_ConjTNFS_BLS(k, pknbits):
    alpha = float(Integer(1)/Integer(3))
    c = float((Integer(32)/Integer(9))**(Integer(1)/Integer(3)))
    if k==6:
        # this is L_{p^k}(1/3, (48/9)^(1/3)) and deg_h = 2
        deg_h_C = 2
        #c=(64/9)**(1/3)
        c = float((Integer(64)/Integer(9))**(Integer(1)/Integer(3)))
        cost_C = float(L_N_alpha_c(pknbits, alpha, c)-21.0)
    elif k==9: # no
        deg_h_C = 3
        #c=(48/9)**(1/3)
        c = float((Integer(48)/Integer(9))**(Integer(1)/Integer(3)))
        cost_C = float(L_N_alpha_c(pknbits, alpha, c))
    elif k==12:
        deg_h_C = 6
        #c=(48/9)**(1/3)
        c =float((Integer(48)/Integer(9))**(Integer(1)/Integer(3)))
        cost_C = float(L_N_alpha_c(pknbits, alpha, c))
    elif k==15:
        deg_h_C = 3
        #c=(48/9)**(1/3)
        c =float((Integer(48)/Integer(9))**(Integer(1)/Integer(3)))
        cost_C = float(L_N_alpha_c(pknbits, alpha, c))
    elif k==24:
        deg_h_C = 24
        #c=(48/9)**(1/3)
        c = float((Integer(48)/Integer(9))**(Integer(1)/Integer(3)))
        cost_C = float(L_N_alpha_c(pknbits, alpha, c)-15.0)
    elif k==27:
        deg_h_C = 9
        #c=(48/9)**(1/3)
        c =float((Integer(48)/Integer(9))**(Integer(1)/Integer(3)))
        cost_C = float(L_N_alpha_c(pknbits, alpha, c))
    return deg_h_C, cost_C


