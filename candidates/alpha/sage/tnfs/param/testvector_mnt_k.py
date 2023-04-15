# https://crypto.stanford.edu/pbc/mnt.html
# to generate MNT6 parameters:
# ./gen/listmnt <D_start>
# ./gen/gendparam D
#
# ./gen/listmnt >list_mnt_params.txt
# for D in `grep list_mnt_params.txt -e ' [4-9][0-9][0-9],' | cut -d ',' -f1` ; do ./gen/gendparam $D ;done
# tail -n+2 list_mnt_params.txt | while read a b c ; do if [ $b = $c, ] ; then echo $b $c, $a ; fi; done | sort -n
#
# for D in `cat list_mnt_params_small_sorted.txt | cut -d ',' -f3` ; do echo "D=$D"; echo "D=$D">>list_mnt_params_ab_small_p_sorted.res; ./gen/gendparam $D >>list_mnt_params_ab_small_p_sorted.res ;done

# other parameters:
# 2006 DCC Scott Barreto Generating more MNT curves 

# (u, D, c, d, a, b, pnbits, label)
# c such that r*c = p+1-t and r is prime
# d such that r = Phi_k(u)/d and r is prime

# MNT4_298 curve BCTV (Ben-Sasson Chiesa Tromer Virza) CRYPTO'14 http://eprint.iacr.org/2014/595
#D = 614144978799019
#A4 = 2
#B4 = 0x3545a27639415585ea4d523234fc3edd2a2070a085c7b980f4e9cd21a515d4b0ef528ec0fd5
#q4 = 0x3bcf7bcd473a266249da7b0548ecaeec9635d1330ea41a9e35e51200e12c90cd65a71660001
#A6 = 11
#B6 = 0xd68c7b1dc5dd042e957b71c44d3d6c24e683fc09b420b1a2d263fde47ddba59463d0c65282
#q6 = 0x3bcf7bcd473a266249da7b0548ecaeec9635cf44194fb494c07925d6ad3bb4334a400000001
#r4 = q6
#r6 = q4
#t4 = q4+1-r4
#u4 = t4-1 # 0x1eef5546609756bec2a33f0dc9a1b671660000
#t6 = q6+1-r6
#u6 = (t6-1)//2 # -0xf77aaa3304bab5f61519f86e4d0db38b30000

test_vector_MNT_k = [
    {'k':4,
     'u':-0x72a74a7ba34d61b0,
     'D':684811,
     'c':1,
     'a':-3,
     'b':0xc4eb66bd690a0b02bdbf8f38bf33148,
     'pnbits':126,
     'rnbits':126,
     'cost_S':48, #h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':58,#h = y^2 + y + 2 # [2, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':54,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':54,# 54 for d=3, 58 for d=4, 63 for d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':57, # 57 with sieving_dim==2 and 3, 64 with d==4
     'SS_deg_f':6,
     'SS_deg_g':4,
     'sieving_dim_NFS_SS':2,
     'label':"MNT4_126"},# 0
    {'k':6,
     'u':0x3953a53dd1a6b0d8,
     'D':684811,
     'c':1,
     'a':-3,
     'b':0xe461ae1ea18b585a84b4cc13a54d28,
     'pnbits':126,
     'rnbits':126,
     'cost_S':59,
     'deg_h_S':2,
     'cost_SS':88, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':71,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':71,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':75,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':72,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_126"},# 1
    #STNFS (dh=2,cost=59) SNFS (sd=3,cost=77), (sd=4,cost=71) (sd=5,cost=74) GJL (df=7,sd=2,cost>=79) (df=7,sd=3,cost=75) (df=7,sd=4,cost>=80) (df=8,sd=3,cost>=77) (df=8,sd=2,cost>=84) (df=8,sd=4,cost>=80) SSingh (df=9,dg=6,dim=2,cost=90) (df=9,dg=6,dim=3,cost=72) (df=9,dg=6,dim=4,cost=73) (df=8,dg=6,dim=2,cost=85) (df=8,dg=6,dim=3,cost=73) (df=8,dg=6,dim=4,cost=77)
    #inv_zeta_Kh, w = 0.860829, 2
    #h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 - 4130826968157827288*x^2 - 4130826968157827291*x - 1
    #max_fi =     37, log_2 max_gi =  61.84
    #aut = 6
    #alpha_f = 0.3924 alpha_g = 2.7560, sum = 3.1485 (basis e)
    {'k':4,
     'u':-0x9576684660bded1c,
     'D':41950339,
     'c':1,
     'a':-3,
     'b':0x37e81bdbc93beef0fd967f2a9f98d56,
     'pnbits':127,
     'rnbits':127,
     'cost_S':48, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':57,# h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':54,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':54, # 54 with d=3, 64 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_127"},# 2
    {'k':6,
     'u':0x4abb3423305ef68e,
     'D':41950339,
     'c':1,
     'a':-3,
     'b':0x28eae4ef8038001e94454a5b33b80d66,
     'pnbits':127,
     'rnbits':127,
     'cost_S':59,
     'deg_h_S':2,
     'cost_SS':88, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':71,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':71,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':75,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':72,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_127"},# 3
    #STNFS (dh=2,cost=59,Px=4*x^2+1) SNFS (sd=4,cost=71) GJL (df=7,sd=3,cost=75) SSingh (df=9,dg=6,dim=3,cost=72)
    #inv_zeta_Kh, w = 0.860829, 2
    #h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 - 5384955105191589518*x^2 - 5384955105191589521*x - 1
    #max_fi =     37, log_2 max_gi =  62.22
    #aut = 6
    #alpha_f = 0.3924 alpha_g = 2.3875 sum = 2.7799
    {'k':4,
     'u':-0x932844866434d2d1a,
     'D':17851507,
     'c':1,
     'a':-3,
     'b':0x3aca3369a72f387b48cd12d161de383aec,
     'pnbits':135,
     'rnbits':135,
     'cost_S':50, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':59,# h = y^2 + 2*y - 2 # [-2, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':56,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':56, # 56 with d=3, 64 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_135"},# 4
    {'k':6,
     'u':0x49942243321a6968d,
     'D':17851507,
     'c':1,
     'a':-3,
     'b':0x1af25d60b5769a694efaaac7af240c2b2e,
     'pnbits':135,
     'rnbits':135,
     'cost_S':61,
     'deg_h_S':2,
     'cost_SS':91,# h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':74,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':74,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':78,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':74,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_135"},# 5
    #STNFS (dh=2,cost=61,Px=4*x^2+1) SNFS (sd=4,cost=74) GJL (df=7,sd=3,cost=78) SSingh (df=9,dg=6,dim=3,cost=74)
    #inv_zeta_Kh, w = 0.696878, 2
    #h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 - 84830405333130581645*x^2 - 84830405333130581648*x - 1
    #max_fi =     37, log_2 max_gi =  66.20
    #aut = 6
    #alpha_f = -0.5586 alpha_g = 2.2504 sum = 1.6918
    {'k':4,
     'u':0xa5fd1306547376302,
     'D':1430907,
     'c':1,
     'a':-3,
     'b':0xc90b53fe4e52907d1101b2d25178bd42,
     'pnbits':135,
     'rnbits':135,
     'cost_S':50, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':60,# h = y^2 + 1 # [1, 0, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':55,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':55, # 55 with d=3, 64 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_135"},# 6
    {'k':6,
     'u':-0x52fe89832a39bb181,
     'D':1430907,
     'c':1,
     'a':-3,
     'b':0x55481637e71613f9bdc20b4430fd51d9ce,
     'pnbits':135,
     'rnbits':135,
     'cost_S':62,
     'deg_h_S':2,
     'cost_SS':91,# h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':73,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':73,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':78,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':74,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_135"},# 7
    #STNFS (dh=2,cost=62,Px=4*x^2+1) SNFS (sd=4,cost=73) GJL (df=7,sd=3,cost=78) SSingh (df=9,dg=6,dim=3,cost=74)
    #inv_zeta_Kh, w = 0.696878, 2
    #h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 + 95685896826187919745*x^2 + 95685896826187919742*x - 1
    #max_fi =     37, log_2 max_gi =  66.37
    #aut = 6
    #alpha_f = -0.5586 alpha_g = 2.7527 sum = 2.1942
    {'k':4,
     'u':-0x5633623df28c77cf66d0,
     'D':8911723,
     'c':1,
     'a':-3,
     'b':0x1d04c3a9636757cdb16104a6c364a99c301750e0,
     'pnbits':157,
     'rnbits':157,
     'cost_S':54, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':61,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':60,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':60, # 60 with d=3, 67 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_157"},# 8
    {'k':6,
     'u':0x2b19b11ef9463be7b368,
     'D':8911723,
     'c':1,
     'a':-3,
     'b':0x270c8a1599bdcf6f5a03da2d43e9d2875e0b8a2,
     'pnbits':157,
     'rnbits':157,
     'cost_S':64,
     'deg_h_S':2,
     'cost_SS':94,# h = y^2 + 5*y + 5 # [5, 5, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':78,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':78,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':84,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':79,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_157"},# 9
    #STNFS (dh=2,cost=64,Px=4*x^2+1) SNFS (d=4,cost=78) GJL (df=7,d=2,cost=86) (df=7,d=3,cost=84) SSingh (df=9,dg=6,dim=3,cost=79)
    #inv_zeta_Kh, w = 0.860829, 2
    #h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 - 203535690277711545611112*x^2 - 203535690277711545611115*x - 1
    #max_fi =     37, log_2 max_gi =  77.43
    #aut = 6
    #alpha_f = 0.3924 alpha_g = -1.1539 sum =-0.7614
    {'k':4,
     'u':0x127e2e742407ae54d59ce,
     'D':12574563,
     'c':1,
     'a':-3,
     'b':0x28ef3cb4e288b2f78317fd2e711d8fa770a77753,
     'pnbits':161,
     'rnbits':161,
     'cost_S':55, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':63,# h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':61,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':61, # 61 with d=3, 67 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_161"},#10
    {'k':6,
     'u':-0x93f173a1203d72a6ace7,
     'D':12574563,
     'c':1,
     'a':-3,
     'b':0xf9d3b5b07822e0587a73606c314295bffe31a8fe,
     'pnbits':161,
     'rnbits':161,
     'cost_S':66,
     'deg_h_S':2,
     'cost_SS':95, # h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':79,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':79,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':85,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':79,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_161"},#11
    #STNFS (dh=2,cost=66,Px=4*x^2+1) SNFS (d=4,cost=79) GJL (df=7,d=2,cost=87) (df=7,d=3,cost=85) SSingh (df=9,dg=6,dim=3,cost=79)
    #inv_zeta_Kh, w = 0.696878, 2
    #h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 + 698641870279835749100775*x^2 + 698641870279835749100772*x - 1
    #max_fi =     37, log_2 max_gi =  79.21
    #aut = 6
    #alpha_f = -0.5586 alpha_g = 2.1540 sum = 1.5955
    {'k':4,
     'u':-0x2181f3c2869cc7173779c,
     'D':1807467,
     'c':1,
     'a':-3,
     'b':0x4398902c959ef541841bf8d2d705ef36b2e81fc4e,
     'pnbits':163,
     'rnbits':163,
     'cost_S':55, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':62,# h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':61,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':61, # 61 with d=3, 68 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_163"},#12
    {'k':6,
     'u':0x10c0f9e1434e638b9bbce,
     'D':1807467,
     'c':1,
     'a':-3,
     'b':0x26d24b29fa7828bfb1f08893aebf93084cbfa84e0,
     'pnbits':163,
     'rnbits':163,
     'cost_S':66,
     'deg_h_S':2,
     'cost_SS':95, # h = y^2 + 3*y - 1 # [-1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':79,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':79,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':85,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':80,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_163"},#13
    #STNFS (dh=2,cost=66,Px=4*x^2+1) SNFS (d=4,cost=79) (d=3,cost=83) GJL (df=7,d=2,cost=87) (df=7,d=3,cost=85) SSingh (df=9,dg=6,dim=3,cost=80)
    #inv_zeta_Kh, w = 0.860829, 2
    #h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 - 1265882309357691486190542*x^2 - 1265882309357691486190545*x - 1
    #max_fi =     37, log_2 max_gi =  80.07
    #aut = 6
    #alpha_f = 0.3924 alpha_g = 1.7802 sum = 2.1727
    {'k':4,
     'u':-0x16673a97a34cd60ce3d83dc,
     'D':28894627,
     'c':1,
     'a':-3,
     'b':0x1c094b7c78c115aef84a9e6d0cebdb94ba0b7a8ef13e8,
     'pnbits':177,
     'rnbits':177,
     'cost_S':58, # h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':64,# h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':63,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':63, # 63 with d=3, 69 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_177"},#14
    {'k':6,
     'u':0xb339d4bd1a66b0671ec1ee,
     'D':28894627,
     'c':1,
     'a':-3,
     'b':0xa9f951970f1dfcaa90a7ff24685017adc436bf442194,
     'pnbits':177,
     'rnbits':177,
     'cost_S':70,
     'deg_h_S':2,
     'cost_SS':96,# h = y^2 + 2*y - 2 # [-2, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':83,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':83,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':89,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':83,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_177"},#15
    #STNFS (dh=2,cost=70,Px=4*x^2+1) SNFS (d=4,cost=83) GJL (df=7,d=2,cost=89) (df=7,d=3,cost=89) SSingh (df=9,dg=6,dim=3,cost=83)
    #inv_zeta_Kh, w = 0.640123, 2
    #h = y^2 + 2*y - 2 # [-2, 2, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 - 216670820936602348867731950*x^2 - 216670820936602348867731953*x - 1
    #max_fi =     37, log_2 max_gi =  87.49
    #aut = 6
    #alpha_f = -0.3100 alpha_g = 3.4278 sum = 3.1178
    {'k':4,
     'u':0x2ad077548e6c85835006663ee4,
     'D':1060147,
     'c':1,
     'a':-3,
     'b':0x58530a2c040c0ae92373b08a4a1edd5cb0662f357acde77703c,
     'pnbits':203,
     'rnbits':203,
     'cost_S':62, # h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':68,# h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':67,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':67, # 67 with d=3, 72 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_203"},#16
    {'k':6,
     'u':-0x15683baa473642c1a803331f72,
     'D':1060147,
     'c':1,
     'a':-3,
     'b':0x322aba54eb77138eff3d82f38b9a126deecdde180fc998d2e0,
     'pnbits':203,
     'rnbits':203,
     'cost_S':74,
     'deg_h_S':2,
     'cost_SS':100,# h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':99,# h = y^2 + 3*y - 3 # [-3, 3, 1] class number Kh = 1
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':88,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':93,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':89,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_203"},#17
    #STNFS (dh=2,cost=74,Px=4*x^2+1) SNFS (d=4,cost=88) (d=3,cost=90) GJL (df=7,d=2,cost=93) (df=7,d=3,cost=94) SSingh (df=9,dg=6,dim=3,cost=89) (df=9,dg=6,dim=2,cost=100)
    #inv_zeta_Kh, w = 0.640123, 2
    #h = y^2 + 2*y - 2 # [-2, 2, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 + 1696049984560259057160035704690*x^2 + 1696049984560259057160035704687*x - 1
    #max_fi =     37, log_2 max_gi = 100.42
    #aut = 6
    #alpha_f = -0.3100 alpha_g = 2.4653 sum = 2.1553
    {'k':4,
     'u':0x5fcc6cf8c6616380f2a00a2d88,
     'D':9877443,
     'c':1,
     'a':-3,
     'b':0x278918bf852af9d0a6d281adbd995f82c664cbec6d85e6ebcce,
     'pnbits':206,
     'rnbits':206,
     'cost_S':63, # h = y^2 + 5*y - 1 # [-1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':70,# h = y^2 + 5*y - 1 # [-1, 5, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':68,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':68, # 68 with d=3, 73 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_206"},#18
    {'k':6,
     'u':-0x2fe6367c6330b1c079500516c4,
     'D':9877443,
     'c':1,
     'a':-3,
     'b':0x461b98fcdddd36f82fd435c029d21bed4ea6dee5f8551e02faa,
     'pnbits':206,
     'rnbits':206,
     'cost_S':75,
     'deg_h_S':2,
     'cost_SS':101,# h = y^2 + 5*y - 1 # [-1, 5, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':89,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':89,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':94,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':89,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_206"},#19
    #STNFS (dh=2,cost=75,Px=4*x^2+1) SNFS (d=4,cost=89) (d=3,cost=92) GJL (df=7,d=2,cost=94) (df=7,d=3,cost=95) SSingh (df=9,dg=6,dim=3,cost=89)
    #inv_zeta_Kh, w = 0.539679, 2
    #h = y^2 + 3*y - 2 # [-2, 3, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 + 3794971059826772759023915964100*x^2 + 3794971059826772759023915964097*x - 1
    #max_fi =     37, log_2 max_gi = 101.58
    #aut = 6
    #alpha_f = -1.4678 alpha_g = 2.6869 sum = 1.2191
    {'k':4,
     'u':-0xbf22d3eb2dc505cb471d89b1ed04,
     'D':496659,
     'c':1,
     'a':-3,
     'b':0x1f75699153d9292f46c13c0c0d7af14c2340d1f9af4e04877a20c636,
     'pnbits':224,
     'rnbits':224,
     'cost_S':66, # h = y^2 + 3*y - 1 # [-1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':71,# h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':71,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':71, # 71 with d=3, 74 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_224"},#20
    {'k':6,
     'u':0x5f9169f596e282e5a38ec4d8f682,
     'D':496659,
     'c':1,
     'a':-3,
     'b':0x12fd915c918ded9efe7a02e50cdedc0f6d312b211d1fb6e42cad676e,
     'pnbits':224,
     'rnbits':224,
     'cost_S':76,
     'deg_h_S':2,
     'cost_SS':101,# h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':92,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':92,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':96,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':92,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_224"},#21
    #STNFS (dh=2,cost=76,Px=4*x^2+1) SNFS (d=4,cost=92) (d=3,cost=93) (more rels with d=4) GJL (df=7,d=2,cost=96) SSingh (df=9,dg=6,dim=3,cost=92)
    #inv_zeta_Kh, w = 0.860829, 2
    #h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 - 1938349788736867434758349971322498*x^2 - 1938349788736867434758349971322501*x - 1
    #max_fi =     37, log_2 max_gi = 110.58
    #aut = 6
    #alpha_f = 0.3924 alpha_g = 2.5788 sum = 2.9712
    {'k':4,
     'u':0x1eb446a11547dfe43d2c45489ee86,
     'D':16460547,
     'c':1,
     'a':-3,
     'b':0x2878f4737dfaa1045f20de35c16e5ac2160d9a0df6b1d5cf7c721636f,
     'pnbits':226,
     'rnbits':226,
     'cost_S':65, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':70,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':72,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':72, # 72 with d=3, 74 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_226"},#22
    {'k':6,
     'u':-0xf5a23508aa3eff21e9622a44f743,
     'D':16460547,
     'c':1,
     'a':-3,
     'b':0x27170fa900c7098eec7ee38fcb99774aac7e67d1358dc290d34b2ea5c,
     'pnbits':226,
     'rnbits':226,
     'cost_S':76,
     'deg_h_S':2,
     'cost_SS':101, # h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':93,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':93,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':96,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':92,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_226"},#23
    #STNFS (dh=2,cost=76,Px=4*x^2+1) SNFS (d=4,cost=93) (d=3,cost=93) (more rels with d=4) GJL (df=7,d=2,cost=96) SSingh (df=9,dg=6,dim=3,cost=92)
    {'k':4,
     'u':-0x1103ba2e08d44a2299becf85b234e8,
     'D':34794363,
     'c':1,
     'a':-3,
     'b':0x7e195e922beede2d1485f751de25f84372f1803b73ee83fb3d294bada8,
     'pnbits':233,
     'rnbits':233,
     'cost_S':67, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':72,# h = y^2 + 3*y - 1 # [-1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':73,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':73, # 73 with d=3, 75 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_233"},#24
    {'k':6,
     'u':0x881dd17046a25114cdf67c2d91a74,
     'D':34794363,
     'c':1,
     'a':-3,
     'b':0x1006567cf3a70c153dd034084e51a4f819de8eae69b70a3c1be49a8e3e2,
     'pnbits':233,
     'rnbits':233,
     'cost_S':79,
     'deg_h_S':2,
     'cost_SS':103,# h = y^2 + y - 3 # [-3, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':94,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':94,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':97,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':94,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_233"},#25
    #STNFS (dh=2,cost=79,Px=4*x^2+1) SNFS (d=4,cost=94) (d=3,cost=94) (more rels with d=4) GJL (df=7,d=2,cost=97) SSingh (df=9,dg=6,dim=3,cost=94)
    {'k':4,
     'u':-0x485d3744c79714c5a38a140a2f3d6587a,
     'D':15496387,
     'c':1,
     'a':-3,
     'b':0x545d796af9e7b22d43de83ff3975768e57d5642ec064f2146c0d48b2a0d3fb175,
     'pnbits':261,
     'rnbits':261,
     'cost_S':71, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':75,# h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':76,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':76, # 76 with d=3, 77 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_261"},#26
    {'k':6,
     'u':0x242e9ba263cb8a62d1c50a05179eb2c3d,
     'D':15496387,
     'c':1,
     'a':-3,
     'b':0x392fe9201b8ec4985a99995286bb1b8f77500545abe035f9a218901e4fe440160,
     'pnbits':261,
     'rnbits':261,
     'cost_S':82,
     'deg_h_S':2,
     'cost_SS':106, # h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':98,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':98,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':101,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':99,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_261"},#27
    #STNFS (dh=2,cost=82,Px=4*x^2+1) SNFS (d=3,cost=98) GJL (df=7,d=2,cost=101) SSingh (df=9,dg=6,dim=3,cost=99)
    {'k':4,
     'u':-0x8a49c5efe7cef083b4feb8a01feeb0f70,
     'D':17960923,
     'c':1,
     'a':-3,
     'b':0xb387ac1274ae3d4523282ff63215b6711f489636c0822cfed9c9007d082bd4e02,
     'pnbits':263,
     'rnbits':263,
     'cost_S':71, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':77,# h = y^2 + 3*y - 3 # [-3, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':77,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':77, # 77 with d=3, 78 with d=2
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_263"},#28
    {'k':6,
     'u':0x4524e2f7f3e77841da7f5c500ff7587b8,
     'D':17960923,
     'c':1,
     'a':-3,
     'b':0x38de09175bbf755e5e4660c59d14690a268eee811b04a2b5d45bb06c641e3dfb34,
     'pnbits':263,
     'rnbits':263,
     'cost_S':83,
     'deg_h_S':2,
     'cost_SS':107, # h = y^2 + 3*y - 1 # [-1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':98,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':98,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':107,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':99,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_263"},#29
    #STNFS (dh=2,cost=83,Px=4*x^2+1) SNFS (d=3,cost=98) (d=4,cost=100) GJL (df=7,d=2,cost=101) SSingh (df=9,dg=6,dim=3,cost=99)
    {'k':4,
     'u':-0x96799dc5ea6686f885ab9e31dabf671375d0,
     'D':41614147,
     'c':1,
     'a':-3,
     'b':0x394ec7738f113baafcad03179ecb04cd53e00b9bc81c86e603cb5cbbb6256234f786c0b0,
     'pnbits':287,
     'rnbits':287,
     'cost_S':75, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':80,# h = y^2 + y + 2 # [2, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':79,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':79, # 80 with d=3, 79 with d=2
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_287"},#30
    {'k':6,
     'u':0x4b3ccee2f533437c42d5cf18ed5fb389bae8,
     'D':41614147,
     'c':1,
     'a':-3,
     'b':0x5748b9416f9a4e1a44705fc91e15922758eb82b6a6155ee0db15a6d13eae0e19491e88e2,
     'pnbits':287,
     'rnbits':287,
     'cost_S':86,
     'deg_h_S':2,
     'cost_SS':110,# h = y^2 + 3*y - 1 # [-1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':102,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':102,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':106,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':105,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_287"},#31
    #STNFS (dh=2,cost=86,Px=4*x^2+1) SNFS (d=3,cost=102) GJL (df=7,d=2,cost=104) SSingh (df=9,dg=6,dim=3,cost=103) (df=9,dg=6,dim=2,cost=108)
    {'k':4,
     'u':-0xe3a50a36e2c20e6aa1e5b4045ef8c6ec5d258,
     'D':1695003,
     'c':1,
     'a':-3,
     'b':0xa51144a05fa22552e42eacc7bd473d5d6ccb59f05ccd26a01edce2c9a352dcc9a34e169000,
     'pnbits':296,
     'rnbits':296,
     'cost_S':75, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':78,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':81,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':81, # 81 with d=3, 81 with d=2
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_296"},#32
    {'k':6,
     'u':0x71d2851b7161073550f2da022f7c63762e92c,
     'D':1695003,
     'c':1,
     'a':-3,
     'b':0x737cf382714499a73ca62b1655657bdea94304df533c2239e4421fea83b21e6b5e6e414f0a,
     'pnbits':296,
     'rnbits':296,
     'cost_S':86,
     'deg_h_S':2,
     'cost_SS':109,# h = y^2 + 5*y + 5 # [5, 5, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':102,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':102,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':106,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':105,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_296"},#33
    #STNFS (dh=2,cost=86,Px=4*x^2+1) SNFS (d=3,cost=102) GJL (df=7,d=2,cost=106) SSingh (df=9,dg=6,dim=3,cost=105)
    # BCTV14 Zexe
    # u=0x1eef5546609756bec2a33f0dc9a1b671660000
    # D=614144978799019
    # a = 2
    # b = 0x3545a27639415585ea4d523234fc3edd2a2070a085c7b980f4e9cd21a515d4b0ef528ec0fd5
    {'k':4,
     'u':0x1eef5546609756bec2a33f0dc9a1b671660000,
     'D':614144978799019,
     'c':1,
     'a':2,
     'b':0x3545a27639415585ea4d523234fc3edd2a2070a085c7b980f4e9cd21a515d4b0ef528ec0fd5,
     'pnbits':298,
     'cost_S':77, # h = y^2 + 5*y - 1 # [-1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':81,# h = y^2 + 5*y - 1 # [-1, 5, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':80,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':80, # 81 with d=3, 80 with d=2
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':79, # 79 with d=2, (f,g)=(6,4)
     'SS_deg_f':6,
     'SS_deg_g':4,
     'sieving_dim_NFS_SS':2,
     'label':"MNT4_298 sk-SNARK"},#34
    {'k':6,
     'u':-0xf77aaa3304bab5f61519f86e4d0db38b30000,
     'D':614144978799019,
     'c':1,
     'a':11,
     'b':0xd68c7b1dc5dd042e957b71c44d3d6c24e683fc09b420b1a2d263fde47ddba59463d0c65282,
     'pnbits':298,
     'cost_S':87, #with deg_h 2
     'deg_h_S':2,
     'cost_SS':109,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':103,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':102,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':106,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':105,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_298 sk-SNARK"},#35
    {'k':4,
     'u':-0x2b0d50fb869ca85baf29ed975b5a45834edb08,
     'D':6555651,
     'c':1,
     'a':-3,
     'b':0x5653cf1a4f5c76a9dd33792a2988f45726b70c634747fcb51e3ded5053ccd99655fb67fb97e,
     'pnbits':299,
     'rnbits':299,
     'cost_S':77, # h = y^2 + 3*y - 1 # [-1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':80,# h = y^2 + y - 3 # [-3, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':81,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':81,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_299"},#36
    {'k':6,
     'u':0x1586a87dc34e542dd794f6cbadad22c1a76d84,
     'D':6555651,
     'c':1,
     'a':-3,
     'b':0x5302210c210150e3d7d08ef2c1fe51b32c30b454eb046de6810b8c49a18e186efc43fbf612a,
     'pnbits':299,
     'rnbits':299,
     'cost_S':87,
     'deg_h_S':2,
     'cost_SS':109,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':103,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':103,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':106,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':105,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_299"},#37
    #STNFS (dh=2,cost=87,Px=4*x^2+1) SNFS (d=3,cost=103) GJL (df=7,d=2,cost=106)  (df=7,d=3,cost=114) SSingh (df=9,dg=6,dim=2,cost=109) (df=9,dg=6,dim=3,cost=105)
    {'k':4,
     'u':-0x271bb1ced62d1f4c983ae7f693e47c5f4fb3e56,
     'D':12121323,
     'c':1,
     'a':-3,
     'b':0x1d39c2b4dc4c2ad6ce97f225e9620add6b3b29f70b361d725d05ec38aa2ea44d8abf420df7b99,
     'pnbits':307,
     'rnbits':307,
     'cost_S':77, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':79,# h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':82,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':82,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_307"},#38
    {'k':6,
     'u':0x138dd8e76b168fa64c1d73fb49f23e2fa7d9f2b,
     'D':12121323,
     'c':1,
     'a':-3,
     'b':0x2432eed160d1c11dced0cf34c1b59c6db489e2c0333ef2d9161b7522c28834fac09b5eb7aa040,
     'pnbits':307,
     'rnbits':307,
     'cost_S':88,
     'deg_h_S':2,
     'cost_SS':109, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':104,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':104,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':107,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':106,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_307"},#39
    #STNFS (dh=2,cost=88,Px=4*x^2+1) SNFS (d=3,cost=104) GJL (df=7,d=2,cost=107) SSingh (df=9,dg=6,dim=3,cost=106)
    {'k':4,
     'u':0x11fc9a3244dddd91a88916b6b4b2f0fbe5f3f84e,
     'D':1475251,
     'c':1,
     'a':-3,
     'b':0x58b0f11b68e5b9c39c7375e2a258280115c27926f4a385eed54cc87fdd4a90f3deeefde34c2a64,
     'pnbits':313,
     'rnbits':313,
     'cost_S':78, # h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':81,# h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':82,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':82,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_313"},#40
    {'k':6,
     'u':-0x8fe4d19226eeec8d4448b5b5a59787df2f9fc27,
     'D':1475251,
     'c':1,
     'a':-3,
     'b':0x10279c4d5ef88311aee3d1344d21b662fa3945c0fc432e52956eb902275e8567ddc1b8c5b7769e6,
     'pnbits':313,
     'rnbits':313,
     'cost_S':89,
     'deg_h_S':2,
     'cost_SS':109,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':105,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':105,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':107,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':107,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_313"},#41
    #STNFS (dh=2,cost=89,Px=4*x^2+1) SNFS (d=3,cost=105) GJL (df=7,d=2,cost=107) SSingh (df=9,dg=6,dim=3,cost=107)
    {'k':4,
     'u':-0xa5acdc4a29d420f9f051248d169b6b21afc2a45362dd658,
     'D':1535827,
     'c':1,
     'a':-3,
     'b':0x3fccf15a9f970dd190dc17e8db1c8f36be3c410e8931d34247ab13177c0b0f23e77b986ccc6139fccfec7d5d6c957c,
     'pnbits':375,
     'rnbits':375,
     'cost_S':86, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':86,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':87,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':87,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_375"},#42
    {'k':6,
     'u':0x52d66e2514ea107cf82892468b4db590d7e15229b16eb2c,
     'D':1535827,
     'c':1,
     'a':-3,
     'b':0x2ce1a3b760100055fd4cedc791326388cfe8ed83e425d5f17de7206c368a16dc4877af87951c6f46941dde5cd589a2,
     'pnbits':375,
     'rnbits':375,
     'cost_S':96,
     'deg_h_S':2,
     'cost_SS':116,# h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':113,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':113,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':115,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':117,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_375"},#43
    #STNFS (dh=2,cost=96,Px=4*x^2+1) SNFS (d=3,cost=113) GJL (df=8,d=2,cost=122) (df=7,d=2,cost=115) SSingh (df=9,dg=6,dim=2,cost=117) (df=9,dg=6,dim=3,cost=117)
    {'k':4,
     'u':0x46d36fedb25f5674e4a95030d8cce1b56037abec3adb80dcb0087bca4,
     'D':27397843,
     'c':1,
     'a':-3,
     'b':0x7499aa6f3ec9512cd86708319e055641c445e6cafa0639ee90b5718def5dea6f6916cb1c864945738af9d6fcb078397682875016917dbaf58,
     'pnbits':453,
     'rnbits':453,
     'cost_S':95, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':94,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':93,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':93,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_453"},#44
    {'k':6,
     'u':-0x2369b7f6d92fab3a7254a8186c6670dab01bd5f61d6dc06e58043de52,
     'D':27397843,
     'c':1,
     'a':-3,
     'b':0xf5d83eebef21a3a281806a9fb3a833ae2737b165e4bfca40386cf33c134c8acda2ef271ca92d30f2be52afb0fccb156aa09fbb93608235b7a,
     'pnbits':453,
     'rnbits':453,
     'cost_S':106,
     'deg_h_S':2,
     'cost_SS':122,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':122,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':122,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':123,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':124,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_453"},#45
    #STNFS (dh=2,cost=106,Px=4*x^2+1) SNFS (d=3,cost=122) GJL (df=8,d=2,cost=130) (df=7,d=2,cost=123) SSingh (df=9,dg=6,dim=2,cost=124)
    {'k':4,
     'u':-0xab27ede9bc74898142916870b95c64fdc4d23f9356b19f67aa83350fea,
     'D':15411979,
     'c':1,
     'a':-3,
     'b':0x6bd7697e68e8ba844972fbc2a3f95d76c1b9a52bd2b62dc190c209b9d489ce67104c45be5e93421783f0b69285d6caf960f0131a51a8a5faa0a6,
     'pnbits':463,
     'rnbits':463,
     'cost_S':96, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':96,# h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':94,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':94,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_463"},#46
    {'k':6,
     'u':0x5593f6f4de3a44c0a148b4385cae327ee2691fc9ab58cfb3d5419a87f5,
     'D':15411979,
     'c':1,
     'a':-3,
     'b':0x4d2238f2bb4e20ce97011f8d68088907fd2362bf4ffb68b3d4b8d0976736d6fa546a94a3e980bfb2c2e66c6aae803d03914ee0356a2e895356b4,
     'pnbits':463,
     'rnbits':463,
     'cost_S':107,
     'deg_h_S':2,
     'cost_SS':124,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':123,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':123,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':124,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':124,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_463"},#47
    #STNFS (dh=2,cost=107,Px=4*x^2+1) SNFS (d=3,cost=123) GJL (df=8,d=2,cost=131) (df=7,d=2,cost=124) SSingh (df=9,dg=6,dim=2,cost=124)
    {'k':4,
     'u':-0x578ec1c9cfe43a7593b33c3093af30ced34a0b97d6b7a84aa4dee23ea3dd41de,
     'D':1974603,
     'c':1,
     'a':-3,
     'b':0xa4eac4e49b3f76275defa3ee105f2755aa70cfcd3601dfbbec4278ce3b951257dfaa9f8c3aa8aea6cb1ac64cbfeef73c1b1fc494a3878d86369613471bdaf2f,
     'pnbits':509,
     'rnbits':509,
     'cost_S':101, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':99, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':98,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':98,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_509"},#48
    {'k':6,
     'u':0x2bc760e4e7f21d3ac9d99e1849d7986769a505cbeb5bd425526f711f51eea0ef,
     'D':1974603,
     'c':1,
     'a':-3,
     'b':0x65d107b2cf02af186db37fe26b8fcc30dad7840273b93d606b346aa4a8248e1707abdfef50df4e2002cc3a0c0b6b5bb831f3cf6a8b5645776cdebe59ffe5628,
     'pnbits':509,
     'rnbits':509,
     'cost_S':112,
     'deg_h_S':2,
     'cost_SS':127,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':128,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':128,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':129,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':128,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_509"},#49
    #STNFS (dh=2,cost=112,Px=4*x^2+1) SNFS (d=3,cost=128) GJL (df=8,d=2,cost=136) (df=7,d=2,cost=129) SSingh (df=8,dg=6,sieving_dim=2,cost=128) (df=9,dg=6,sieving_dim=2,cost=128) (df=9,dg=6,sieving_dim=3,cost=135)
    {'k':4,
     'u':0x96d4f2a424f6e9f9c3cc0974885b79ff5bd3368bea7e3436577af8c4ac9a4f7cee902,
     'D':6347203,
     'c':1,
     'a':-3,
     'b':0x4cb80a2dc6460e8535240c0e36c52dc1824527eff1e4346a8dfff3f73875b1071df3c5225fb269d089b0875db64672a6c948a1774519f5968340a79dd292c84fdab3352146,
     'pnbits':551,
     'rnbits':551,
     'cost_S':105, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':103,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':100,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':100,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_551"},# 50
    {'k':6,
     'u':-0x4b6a7952127b74fce1e604ba442dbcffade99b45f53f1a1b2bbd7c62564d27be77481,
     'D':6347203,
     'c':1,
     'a':-3,
     'b':0x1e421b78380a143fa19af19a903825a48068325445077333d1310f34aa8caa2c4381f7a172de93c6d3c89f54c50d74d1e74ed3b4a35a388920a0e2889fd1d2c02d5d52812,
     'pnbits':551,
     'rnbits':551,
     'cost_S':116,
     'deg_h_S':2,
     'cost_SS':131,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':132,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':132,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':133,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':131,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_551"},# 51
    #STNFS (dh=2,cost=116,Px=4*x^2+1) SNFS (d=3,cost=132) GJL (df=7,d=2,cost=133) SSingh (df=9,dg=6,sieving_dim=2,cost=131)
    {'k':4,
     'u':-0x802606f2a5f718c9499e02250f0ccd025b8587915fa327bab50b26af25752f95f1b5b0,
     'D':2257627,
     'c':1,
     'a':-3,
     'b':0x207e59d01023b0b8d90b45c498455ccc8f971541d2f00956ac1ab65889aafbfebf938a67b4b8512510bc6e1925e2ff8ee0bc4f07331f86ab59906ef34809e17233efc37ea148,
     'pnbits':559,
     'rnbits':559,
     'cost_S':107, # h = y^2 + 3*y - 2 # [-2, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':105,# h = y^2 + 3*y - 2 # [-2, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':100,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':100,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_559"},# 52
    {'k':6,
     'u':0x4013037952fb8c64a4cf0112878666812dc2c3c8afd193dd5a85935792ba97caf8dad8,
     'D':2257627,
     'c':1,
     'a':-3,
     'b':0x2d7b5538e23c0d93f101f5f1d6bc54852daa31683586a12bb0a409a808101eb54f87a934cfbc24fbc768433be0f8f02af2140b024a602ef9a7390dac9c47cd74a4ca9f7b1c32,
     'pnbits':559,
     'rnbits':559,
     'cost_S':118,
     'deg_h_S':2,
     'cost_SS':133,# h = y^2 + 2*y - 2 # [-2, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':133,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':133,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':133,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':132,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_559"},# 53
    #STNFS (dh=2,cost=118,Px=4*x^2+1) SNFS (d=3,cost=133) GJL (df=7,d=2,cost=133) SSingh (df=9,dg=6,sieving_dim=2,cost=132)
    {'k':4,
     'u':-0x2f8dc35838ab0205ed01ae28a53c0cbe4361db8c0bfa73d71a60b804bae22d69e2ea446,
     'D':19353667,
     'c':1,
     'a':-3,
     'b':0x5913bf8c4e7d5cb3d980e307d094654a4ae1fd9bb785274715ca5cb7c30b8c3560e7c8be9456f109c813a9666cf3aa37343357815cf6c09d229e882c072d2f29e6e5b05681ae4,
     'pnbits':564,
     'rnbits':564,
     'cost_S':106, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':104,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':102,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':102,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_564"},# 54
    {'k':6,
     'u':0x17c6e1ac1c558102f680d714529e065f21b0edc605fd39eb8d305c025d7116b4f175223,
     'D':19353667,
     'c':1,
     'a':-3,
     'b':0x76278d7cca9afbc870a8ce59076fb5cbd029fc19035b8077f7ed62d90baa7e2e7cbc3388bc3787a6e336db5599a2f3b80099f5281d27f0a7fad0cb95351b3dc5024c331eb08dc,
     'pnbits':564,
     'rnbits':564,
     'cost_S':117,
     'deg_h_S':2,
     'cost_SS':131,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':133,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':133,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':134,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':132,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_564"},# 55
    #STNFS (dh=2,cost=117,Px=4*x^2+1) SNFS (d=3,cost=133) TNFS: cost=117 or 118? GJL (df=7,d=2,cost=134) SSingh (df=9,dg=6,sieving_dim=2,cost=132)
    {'k':4,
     'u':0x37bfbf1abe3f252d60d7e985c3cfdff1e07e64d512b77f29fe63e271022d8fc29f2b203bcc,
     'D':10128667,
     'c':1,
     'a':-3,
     'b':0x5ea1b568bb8ef31858f0183ce40088a7e533ee41812a1f75ee40d6b4a292c085fe1248012e7b8e0cf657a4c58a2e7d945a753c8620f2dae9344b4711aef6b0584aaf1b3e674e8678252,
     'pnbits':588,
     'rnbits':588,
     'cost_S':109, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':105,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':103,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':103,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_588"},# 56
    {'k':6,
     'u':-0x1bdfdf8d5f1f9296b06bf4c2e1e7eff8f03f326a895bbf94ff31f1388116c7e14f95901de6,
     'D':10128667,
     'c':1,
     'a':-3,
     'b':0x3319c54417c4f1f501005e12b595ef89b5e64f9a32b0f0c78acc77c75c4c87ff3b0e457a75812cda9bd672ddbaefb0f5f9146b24bf1977316ba1f4edbb7baf0ffc46854b0330aea3e1c,
     'pnbits':588,
     'rnbits':588,
     'cost_S':120,
     'deg_h_S':2,
     'cost_SS':134,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':136,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':136,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':136,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':134,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_588"},# 57
    #STNFS (dh=2,cost=120,Px=4*x^2+1) SNFS (d=3,cost=136) GJL (df=7,d=2,cost=136) SSingh (df=9,dg=6,sieving_dim=2,cost=134)
    {'k':4,
     'u':0x47f32e086e3d32db40d62685e13a5cd87342640d75e43faa22a617aeb7c916261c8423fc2a6,
     'D':7701427,
     'c':1,
     'a':-3,
     'b':0x136e3b01ea5baac5f7d5f9a8f768ff6b7aa9c0827a240382c5b29c29959d2dd52a8604c1c21ae80d48b424fb4725802e28fdcbf3d484d33a73de038de4bc82768de71a86030dd36d3bdc52,
     'pnbits':597,
     'rnbits':597,
     'cost_S':110, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':107,# h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':103,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':103,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_597"},# 58
    {'k':6,
     'u':-0x23f99704371e996da06b1342f09d2e6c39a13206baf21fd511530bd75be48b130e4211fe153,
     'D':7701427,
     'c':1,
     'a':-3,
     'b':0x75e8746129bfe38c441e7d1cf5e18ccbac249ee657d4cb28f4aed73f80dce02eb397ad2e6727c0e964c28762f50eceb873692d77068e7afcf62e3c29bdf44dfd0e24d9e6be63ac88eec4a,
     'pnbits':597,
     'rnbits':597,
     'cost_S':121,
     'deg_h_S':2,
     'cost_SS':135,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':137,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':137,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':137,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':135,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_597"},# 59
    #STNFS (dh=2,cost=121,Px=4*x^2+1) SNFS (d=2,cost=147), (d=3,cost=137) (d=4,cost=148) where d is sieving dimension GJL (df=7,d=2,cost=137) SSingh (df=9,dg=6,sieving_dim=2,cost=135)
    {'k':4,
     'u':-0x33219295b67a20cdffee07c4e745f918e8c7278010205d7e949f976529a81b2c524289fd6e0ab71ff78,
     'D':8317003,
     'c':1,
     'a':-3,
     'b':0xa046f925ef7659d3287ba75e8aa5cd007d6e936bff291c19826eb100b6001800674a6cbaabd8fa88f17c7d693998155af09cba7b45050428aecc7dfd74ae6272150aae3a52d4587094da9088c400de4a03308,
     'pnbits':660,
     'rnbits':660,
     'cost_S':116, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':112,# h = y^2 + y + 2 # [2, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':107,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':107,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_660"},#60
    {'k':6,
     'u':0x1990c94adb3d1066fff703e273a2fc8c746393c008102ebf4a4fcbb294d40d96292144feb7055b8ffbc,
     'D':8317003,
     'c':1,
     'a':-3,
     'b':0x8345302bce44caaf24cc80c44a41a076a1318561d52412d56a95ff3d8e0a9e8ebec0c39eb106d8d2c012680b8c58bc00a5940edec117960d37dee84d260a8976202ac158717a7a53b10662b28fe10d2ef710c,
     'pnbits':660,
     'rnbits':660,
     'cost_S':128,
     'deg_h_S':2,
     'cost_SS':141,# h = y^2 + 2*y - 2 # [-2, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':143,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':143,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':143,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':140,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_660"},#61
    #STNFS (dh=2,cost=128,Px=4*x^2+1) SNFS (d=2,cost=150), (d=3,cost=143) GJL (df=7,d=2,cost=143) SSingh (df=9,dg=6,sieving_dim=2,cost=140)
    #
    # CODA: MNT4 MNT6 753 bits
    # https://coinlist.co/build/coda/pages/MNT4753
    #q=0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB117E776F218059DB80F0DA5CB537E38685ACCE9767254A4638810719AC425F0E39D54522CDD119F5E9063DE245E8001
    #r=0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB26C5C28C859A99B3EEBCA9429212636B9DFF97634993AA4D6C381BC3F0057974EA099170FA13A4FD90776E240000001
    #a=2
    #b=0x01373684A8C9DCAE7A016AC5D7748D3313CD8E39051C596560835DF0C9E50A5B59B882A92C78DC537E51A16703EC9855C77FC3D8BB21C8D68BB8CFB9DB4B8C8FBA773111C36C8B1B4E8F1ECE940EF9EAAD265458E06372009C9A0491678EF4
    #tr = q-r+1
    #t=-0x15474b1d641a3fd86dcbcee5dcda7fe51852c8cbe26e600733b714aa43c31a66b0344c4e2c428b07a7713041ba17fff
    #u = t-1
    #assert q == u**2 + u + 1 and r == u**2 + 1
    #y = 0x141dc36427749fc20ad54657b76e627e5fab26c09f8c02c6958f6f61023f4132e92ef9990e0c67f0c9f0e8fb
    {'k':4,
     'u':-0x15474b1d641a3fd86dcbcee5dcda7fe51852c8cbe26e600733b714aa43c31a66b0344c4e2c428b07a7713041ba18000,
     'D':241873351932854907,
     'c':1,
     'a':2,
     'b':0x01373684a8c9dcae7a016ac5d7748d3313cd8e39051c596560835df0c9e50a5b59b882a92c78dc537e51a16703ec9855c77fc3d8bb21c8d68bb8cfb9db4b8c8fba773111c36c8b1b4e8f1ece940ef9eaad265458e06372009c9a0491678ef4,
     'pnbits':753,
     'rnbits':753,
     'cost_S':125, # h = y^2 + 3*y - 1 # [-1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':119,# h = y^2 + 3*y - 1 # [-1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':113,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':113,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':125, # with df=5
     'GJL_deg_f':6,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':118,
     'SS_deg_f':6,
     'SS_deg_g':4,
     'sieving_dim_NFS_SS':2,
     'label':"MNT4_753 CODA"},#62
    # CODA: MNT4 MNT6 753 bits
    # https://coinlist.co/build/coda/pages/MNT4753
    # r=0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB117E776F218059DB80F0DA5CB537E38685ACCE9767254A4638810719AC425F0E39D54522CDD119F5E9063DE245E8001
    # q=0x01C4C62D92C41110229022EEE2CDADB7F997505B8FAFED5EB7E8F96C97D87307FDB925E8A0ED8D99D124D9A15AF79DB26C5C28C859A99B3EEBCA9429212636B9DFF97634993AA4D6C381BC3F0057974EA099170FA13A4FD90776E240000001
    # a=11
    # b=0x7da285e70863c79d56446237ce2e1468d14ae9bb64b2bb01b10e60a5d5dfe0a25714b7985993f62f03b22a9a3c737a1a1e0fcf2c43d7bf847957c34cca1e3585f9a80a95f401867c4e80f4747fde5aba7505ba6fcf2485540b13dfc8468a
    # tr = q-r+1
    # u= (tr-1)//2
    # t = 0x15474b1d641a3fd86dcbcee5dcda7fe51852c8cbe26e600733b714aa43c31a66b0344c4e2c428b07a7713041ba18001
    # u = 0xaa3a58eb20d1fec36e5e772ee6d3ff28c296465f137300399db8a5521e18d33581a262716214583d3b89820dd0c000
    # assert q == 4*u**2 + 1 and  r == 4*u**2 -2*u + 1
    # same D,y
    # y=0x141dc36427749fc20ad54657b76e627e5fab26c09f8c02c6958f6f61023f4132e92ef9990e0c67f0c9f0e8fb
    # D=241873351932854907
    {'k':6,
     'u':0xaa3a58eb20d1fec36e5e772ee6d3ff28c296465f137300399db8a5521e18d33581a262716214583d3b89820dd0c000,
     'D':241873351932854907,
     'c':1,
     'a':11,
     'b':0x7da285e70863c79d56446237ce2e1468d14ae9bb64b2bb01b10e60a5d5dfe0a25714b7985993f62f03b22a9a3c737a1a1e0fcf2c43d7bf847957c34cca1e3585f9a80a95f401867c4e80f4747fde5aba7505ba6fcf2485540b13dfc8468a,
     'pnbits':753,
     'rnbits':753,
     'cost_S':137,
     'deg_h_S':2, # 153,3, 152,7,2, 148,9,6,2
     'cost_SS':148,# h = y^2 + 3*y - 2 # [-2, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':151,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':150,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':151,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':147,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_753 CODA"},#63
    {'k':4,
     'u':0x50b0f8ce495bd90daf2e28d4053fb1ff4cea5febbd8d91c9fe7d74aa0d99f90677aca2e2d6293053fee347f3cb048a67e,
     'D':4329347,
     'c':1,
     'a':-3,
     'b':0x126956b7faf419d8ab4614c3da993ab40145187b76f14265c6df220efbe68f228ecfec6cf5f9755cdcc5774c727ecfd4671d0495c35999d9d143bcb82f18da1f65a650b1d2890aee9ab6a270f9557b65e212d5476a4c3daa39699dc8be5fe5f8c9,
     'pnbits':773,
     'rnbits':773,
     'cost_S':126, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':120,# h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':114,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':114,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_773"},#64
    {'k':6,
     'u':-0x28587c6724adec86d797146a029fd8ffa6752ff5dec6c8e4ff3eba5506ccfc833bd651716b149829ff71a3f9e5824533f,
     'D':4329347,
     'c':1,
     'a':-3,
     'b':0xb37df0a30666263404ae9e7e950762b70c1c88847a96ae43f4a482023db263cf82743b2009561a83d1e6faf03b87c7a631e47eb6273aaa9a6b32631ff4602b5fbd06866309fb38f17a4fdec49070b88258a16d5bdc95a712c5391943806a93386,
     'pnbits':773,
     'rnbits':773,
     'cost_S':137,
     'deg_h_S':2,
     'cost_SS':147,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':153,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':153,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':152,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':148,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_773"},#65
    #STNFS (dh=2,cost=137,Px=4*x^2+1) SNFS (d=2,cost=157) (d=3,cost=153) GJL (df=7,d=2,cost=152) SSingh (df=9,dg=6,sieving_dim=2,cost=148)
    {'k':4,
     'u':0x27c80b709d6cacfa66e092fb568fc6ec3e5600831f390a12a2e0988bce228237ed2f8fd649eb538992d2a693b54cc9916609dc94a8717d3d4daa,
     'D':7983123,
     'c':1,
     'a':-3,
     'b':0x17052616aeea07833b314e42ac0b59f8740f677af160941b335f3dee19e66c1afa2d2d079c2f97a9dfb9f187816968963c0a640a1855c6b25d6aafc78e2f9db2b1745c976cef75cd7bd23301dea67cc5ae010915cb1b6aab3c7613022af03b580e6d407066ab14d4a8fe4c9b0937deef96a54d0,
     'pnbits':923,
     'rnbits':923,
     'cost_S':139, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':130,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':123,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':123,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_923"},#66
    {'k':6,
     'u':-0x13e405b84eb6567d3370497dab47e3761f2b00418f9c850951704c45e711411bf697c7eb24f5a9c4c9695349daa664c8b304ee4a5438be9ea6d5,
     'D':7983123,
     'c':1,
     'a':-3,
     'b':0x340851bc0067daaad8a5af86986373ff508e116e037fe90a6576c3a16e8a29b5bbcc84ed798f64e40f90f1abb9867de00d5709614b553b45d78bcaff531b44ddfdc97ad7e026eaac0582d73f7956e55916a7ad0568a11fe05ba2a6569cd15941f860f9ee0565302a766b8d5fcc32d876b8a3148,
     'pnbits':923,
     'rnbits':923,
     'cost_S':150,
     'deg_h_S':2,
     'cost_SS':157,# h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':164,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':164,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':164,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':157,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_923"},#67
    #STNFS (dh=2,cost=150,Px=4*x^2+1) SNFS (d=2,cost=164) (d=3,cost=166) GJL (df=7,d=2,cost=164) SSingh (df=9,dg=6,sieving_dim=2,cost=157)
    {'k':4,
     'u':0x319bb60a801abbf048b63730b25eb5acf20cbcc9de8dea73fcffb1f15af28a85ef0044aa3a8600d149db1b97251a2f02663910b6cdfe68fce6b4f2f42b764,
     'D':9493747,
     'c':1,
     'a':-3,
     'b':0x2c2bb84f339c31c85792922fe5eec610f5e0f70b80cf856e3cf2941d8a74e87f20674854443b5c319714dd5f6e72993f8975c077c16d6d31c451ddb550d08adf28ee230b18eedb47f6519f855072e5fae354b3b4b6bde385bce24aef9f5466bffef1e604a686d410f1b72c4976b96a4dacb53434aec584ad30591303a,
     'pnbits':996,
     'rnbits':996,
     'cost_S':145, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':135,#h = y^2 + 2*y - 1 # [-1, 2, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':127,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':127,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_996"},#68
    {'k':6,
     'u':-0x18cddb05400d5df8245b1b98592f5ad679065e64ef46f539fe7fd8f8ad794542f78022551d430068a4ed8dcb928d1781331c885b66ff347e735a797a15bb2,
     'D':9493747,
     'c':1,
     'a':-3,
     'b':0x7d2f7ee527518d6ba6cadf80b3f6ebd21516bc542dc3935c6eb1d148afbd310ad9cba5988bc7dd45c3e71607a336efbbc1a8504cad69eca91931ac46058be179e2a178d3759ca6bc8db86a370c96c6b4f2fe497d3838190e93b8c91ebf5a02a693a5faa1ddbb590084c8b45860c90444b269d51a219dd6b423036602c,
     'pnbits':996,
     'rnbits':996,
     'cost_S':156,
     'deg_h_S':2,
     'cost_SS':163,# h = y^2 + 3*y - 3 # [-3, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':168,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':168,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':169,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':162,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_996"},#69
    #STNFS (dh=2,cost=157(try156),Px=4*x^2+1) SNFS (d=2,cost=168) (d=3,cost=171) GJL (df=7,d=2,cost=169) SSingh (df=9,dg=6,sieving_dim=2,cost=162)
    {'k':4,
     'u':-0xbad96595098c3cafa0c09001b21d96dcbaa1c976f23813dfba36d10878302e0db7c316aae60c9e360ca6ffe309c31b5b82a9c40c5c5b4551a524dcfd8da3517ea44e1afbe967ecd1da7ccb44096,
     'D':29171923,
     'c':1,
     'a':-3,
     'b':0x5eac9dbc28e999c6764a95aefc4dba5119b540ec6e93c0b8994b18eeeaee6d2355ba1f663fae59d64a579189a77d02169109c2427f1ceec42d7f7cbd64328a509ada411314ef45bc30e5e84f0e23f5cba2ca64052b2eb0265176585ac9be23986e05d8677a32362016d83e0fc11e23b7dc3095eba46a89baad1229dd968948ab521711aff6c665d3383b10b67e75339a7b06623a2baf1272070f7a,
     'pnbits':1240,
     'rnbits':1240,
     'cost_S':163, # h = y^2 + 5*y + 1 # [1, 5, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':150,# h = y^2 + 4*y + 1 # [1, 4, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':139,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':139,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_1240"},#70
    {'k':6,
     'u':0x5d6cb2ca84c61e57d0604800d90ecb6e5d50e4bb791c09efdd1b68843c181706dbe18b5573064f1b06537ff184e18dadc154e2062e2da2a8d2926e7ec6d1a8bf52270d7df4b3f668ed3e65a204b,
     'D':29171923,
     'c':1,
     'a':-3,
     'b':0x7ec289b0e8eca22aacaf43447b60e4cb837db4499580332b075f70b956d47d0589a61cd4c2caeaf602f40d796a383e133d05e8e0e10e3ac4d735465623858a24e9ed7e8469e99f8940eb2012a6d402be3596cc21e09ae3883d7a04408efc10280349cc7758a803d8c6ad1ba78eafa122b83ea998e4b2a6390208ba4d4786ad60427869e9e11e6a7a160b00791107465f7dfe198e05911a40cba2a2,
     'pnbits':1240,
     'rnbits':1240,
     'cost_S':175,
     'deg_h_S':2,
     'cost_SS':177,# h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':179,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':179,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':186,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':177,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_1240"},#71
    #STNFS (dh=2,cost=175(try174),Px=4*x^2+1) SNFS (d=2,cost=179) GJL (df=7,d=2,cost=186) (df=8,d=2,cost=188) SSingh (df=9,dg=6,sieving_dim=2,cost=177)
    {'k':4,
     'u':-0x4b8cd85876cc976eae11ec31cdcc54dc706afdb6a464630c866a3f319763c6a25d3c3377aba70c9898f57ae0fff38b7152d499b0289995711528f79cceff5ac6be8338edb4870b4424cc0ea48307f7e328a943dd4a,
     'D':972483,
     'c':1,
     'a':-3,
     'b':0x136e44805398f6331866fe2bcf51c758e7f154373e3bb98f9961ba6725fb2ac27f515c95a86d00cc6e5c1a9342e4210edfc6017e8ac1929cceca2d4189d9d263a9e3bffeff4f267c53ed757b8819f20ebc2196c4d963223607191b818ebb188a0864d733c6a4ad12a850b359657bbdae5c2d1f87671569f5e200ce027b46879d666854d48c3014fad3938ab391bb76e9ea05883f48f27bd13ff1d1b2cf296fd54c1f2d4d2759e239ced6,
     'pnbits':1357,
     'rnbits':1357,
     'cost_S':171, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':156,# h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':144,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':144,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':None,
     'SS_deg_f':None,
     'SS_deg_g':None,
     'sieving_dim_NFS_SS':None,
     'label':"MNT4_1357"},#72
    {'k':6,
     'u':0x25c66c2c3b664bb75708f618e6e62a6e38357edb5232318643351f98cbb1e3512e9e19bbd5d3864c4c7abd707ff9c5b8a96a4cd8144ccab88a947bce677fad635f419c76da4385a2126607524183fbf19454a1eea5,
     'D':972483,
     'c':1,
     'a':-3,
     'b':0x3c45c16c3655ae5c37bf01f47ee3c6ee0439d4315d95fdac6d844781d7190aa52ea22355085c5c95987bcb00071ab0217cedeaa4d5a8e02664c055cf84c24713bf777bd2ebe208bdfe4639832988d72c2fed917388eb8689c5ef0eeff20dcd30128812c914a97a8144d6a45a736f30d9299edecc2dea4c06dcbbbd728faceb7050839fcf5f2a48e7f3b0401f6a5381660ba2cb023aded085cad36922043b4818926ef25485aef50eec4,
     'pnbits':1357,
     'rnbits':1357,
     'cost_S':182,
     'deg_h_S':2,
     'cost_SS':183, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':184,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':184,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':194,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':183,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_1357"},#73
    #STNFS (dh=2,cost=182,Px=4*x^2+1) SNFS (d=2,cost=184) GJL (df=7,cost=194) (df=8,cost=195) SSingh (df=9,dg=6,sieving_dim=2,cost=183)
    {'k':4,
     'u':0x8b391da231f3bfd8a3887bd844a8018848f9ea5654338c9cdf35e5c6752419293b9b14b891dfc84fffc84aada0241bd2a1d18a704293cb3e624367e3d7cb89382f5abc2d98479128158ea1162d2af11355ff210ed0d8e94e43b1c98581a710ed4e7b1d956a773c821e82de451ad7e0f08301c626a908d2edc4db4fc5c0f279d0,
     'D':18271507,
     'c':1,
     'a':-3,
     'b':0x135fa94331c6ed89bc4e625bcf07fbe285200785536d30c3b238a4ba9134625509e095b861522207563c1b0f509cf22a864e2b58a80ad5168678a1f9d86c5419a4b05644feff07257ed7127a6f9ce80a2a9c9ba1e64fefd0fae2e590981835191f9c7cd9fee5c7228689c267b32d62d82e34b13235a4de009db1ca0b7b9438fac97ee5d3e50b41a11728f9e1ee135d511d1fb872e31bf78e974157d79585b297b0b370cd4b1f7e7383ef63b59a06559d26fe08ec964f91aa132a632f2edb2cec28a0534c72ae701549f56c4c14e3a05bd8c602239a514ee26674990da65713e2b9176df731cd960a317ab983391e2da818e143f9f546338ed1b1ac445b7be526,
     'pnbits':2047,
     'rnbits':2047,
     'cost_S':213, # h = y^2 + 3*y + 1 # [1, 3, 1] class number Kh = 1
     'deg_h_S':2,
     'cost_SS':188,
     'deg_h_SS':2,
     'deg_aux_SS':4,#if deg_aux_SS==3: cost_SS==190, h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'cost_SNFS':174,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':174,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':None,
     'GJL_deg_f':None,
     'sieving_dim_GJL':None,
     'cost_NFS_SS':199, # 199 with (f,g)=(10,8)
     'SS_deg_f':12,
     'SS_deg_g':8,
     'sieving_dim_NFS_SS':2,
     'label':"MNT4_2047"},#74
    {'k':6,
     'u':-0x459c8ed118f9dfec51c43dec225400c4247cf52b2a19c64e6f9af2e33a920c949dcd8a5c48efe427ffe42556d0120de950e8c5382149e59f3121b3f1ebe5c49c17ad5e16cc23c8940ac7508b16957889aaff9087686c74a721d8e4c2c0d38876a73d8ecab53b9e410f416f228d6bf0784180e31354846976e26da7e2e0793ce8,
     'D':18271507,
     'c':1,
     'a':-3,
     'b':0x41d434eeaebea54cfd369e20392cd5bc4262c024a87966498789efc820bfd34a38ce45b22a64852ab471720fa83ade841cdc46e4d08bfcb630875cd224bda9b8f063fb940ccc3c01f365df79a51f6f6a90594ec6e139e3bb0163be0298d1507caa30d0f84ca9742f7fe24bb99170ea6b631856f80c069ec4c3fa1987b42b4491b0bb8fe10419ce31563c1f5f97449be877715df97861c1b0abdc2384dd5748847d0464d2af757e68e3e9baa2cf0dbc5ddcea5bf1f955d9f1b64052cb95a46e67bbf7ba0498ce6ff07ac2701b814895adb754d673b81b08cd3e3d495ec2e3ee93f3094abd2638dc7f6685311c62a2ce8069b195442c99396efa6b6179e32f12d0,
     'pnbits':2047,
     'rnbits':2047,
     'cost_S':225,
     'deg_h_S':2,
     'cost_SS':217, # h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':212,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':212,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':232,
     'GJL_deg_f':8,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':217,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_2047"},#75
    #SNFS (d=2,cost=212) GJL (df=7,sdim=2,cost=233) (df=8,sdim=2,cost=232) (df=9,sdim=2,cost=234) (df=10,sdim=2,cost=238) SSingh (df=8,dg=6,sdim=2,cost>=223) (df=9,dg=6,sdim=2,cost=217) (df=10,dg=8,sdim=2,cost=), (df=12,dg=9,sdim=2,cost=)
]

"""
Computation of curve parameters a, b with Sutherland code

classpoly D 0 p
where 0 means computes the j-invariant

classpoly 21048712401 0 6210044120409721004947206240885978274523751269793792001

alternatively, the Edwards coordinates are given at
https://github.com/scipr-lab/libff/blob/develop/libff/algebra/curves/edwards/edwards_init.cpp
a=1
d=600581931845324488256649384912508268813600056237543024
That gives the Montgomery form y^2 = x^3 + 2*(a+d)*x^2 + (a-d)^2*x
a2 = 2*(a+d)
a4 = (a-d)**2
a6 = 0
then translating by f(x-a2/3) one gets
a = 5146191432578561295273296286096700125072183936719675670
b = 1608172008978859923628383951114075654993087196492545460
"""

test_vector_MNTG = [
    #t = u+1
    #r = Phi_k(u)/d
    #r = (px+1-tx)/c
    {'k':6,
     'u0':0x8eed757d90615e40000000, # t-1 = -26*u0-1
     'u':-0xe841deec0a9e39280000003,
     'D':21048712401,
     'y':0x3266401cfffea660000,
     'c':4,
     'd':13,
     'a':5146191432578561295273296286096700125072183936719675670,
     'b':1608172008978859923628383951114075654993087196492545460,
     'pnbits':183,
     'rnbits':181,
     # estimates should be close to MNT6 of the same size, above there is MNT6-177
     'cost_S':71,
     'deg_h_S':2,
     'cost_Conj':71,# TNFS-Conj
     'deg_h_Conj':2,
     'cost_SS':96,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':83,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':83,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':89,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':83,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNTG_183 SNARK"},
]

test_vector_MNT6_SB06 = [
    #r = Phi_6(u)/3
    #c = 3 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = x^2 - x + 1
    #px = x^2 + 1
    {'k':6,
     'u':-0xa7861529d98d1433ef32,
     'D':62003,
     'c':3,
     'd':3,
     'a':-3,
     'b':0x4a43bdd92bd614b40ed8ea5d5dda0db4876dab1e,
     'pnbits':159,
     'rnbits':158,
     'label':"MNT6_159 Scott-Barreto DCC 06"},# 0
    #r = Phi_6(u)/3
    #c = 2 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (2*x^2 - 2*x + 2)/3
    #px = (2*x^2 + x + 2)/3
    {'k':6,
     'u':-0xdd128bdd2db67e8a9579,
     'D':7847065,
     'c':2,
     'd':3,
     'a':-3,
     'b':0x3c152bb439ab18dcf125ae76d03ea9aa9d4205c9,
     'pnbits':159,
     'rnbits':158,
     'label':"MNT6_159 Scott-Barreto DCC 06"},# 1
    #r = Phi_6(u)/13
    #c = 4 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (4*x^2 - 4*x + 4)/13
    #px = (4*x^2 + 9*x + 4)/13
    {'k':6,
     'u':0x1a62bf56cfa9e1e5f7b8f,
     'D':717595,
     'c':4,
     'd':13,
     'a':-3,
     'b':0x143791687b84b21278ca3fedcf598cdc503c2ba9,
     'pnbits':160,
     'rnbits':158,
     'label':"MNT6_160 Scott-Barreto DCC 06"},# 2
    #r = Phi_6(u)/13
    #c = 4 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (4*x^2 - 4*x + 4)/13
    #px = (4*x^2 + 9*x + 4)/13
    {'k':6,
     'u':-0x1ca4343626d14109ab96b,
     'D':1397298,
     'c':4,
     'd':13,
     'a':-3,
     'b':0x983e982cd9f263852e293b8374f1bfe0acd4fdc9,
     'pnbits':160,
     'rnbits':158,
     'label':"MNT6_160 Scott-Barreto DCC 06"},# 3
    #r = Phi_6(u)/13
    #c = 4 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (4*x^2 - 4*x + 4)/13
    #px = (4*x^2 + 9*x + 4)/13
    {'k':6,
     'u':-0x19283002199000dd55b83,
     'D':1523371,
     'c':4,
     'd':13,
     'a':-3,
     'b':0x29316f28dd69b24996be73bbeff61871cefd95ac,
     'pnbits':160,
     'rnbits':158,
     'label':"MNT6_160 Scott-Barreto DCC 06"},# 4
    #r = Phi_6(u)/13
    #c = 4 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (4*x^2 - 4*x + 4)/13
    #px = (4*x^2 + 9*x + 4)/13
    {'k':6,
     'u':-0x15d7917b7999921213331,
     'D':1983787,
     'c':4,
     'd':13,
     'a':-3,
     'b':0x7ba12469ab49a608279ab58e04c1d4e4546aa53,
     'pnbits':160,
     'rnbits':158,
     'label':"MNT6_160 Scott-Barreto DCC 06"},# 5
    #r = Phi_6(u)/13
    #c = 4 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (4*x^2 - 4*x + 4)/13
    #px = (4*x^2 + 9*x + 4)/13
    {'k':6,
     'u':-0x171711f606609fc1c2a23,
     'D':8807457,
     'c':4,
     'd':13,
     'a':-3,
     'b':0x3895102b475482f97e904fd72fa1d101697aed4f,
     'pnbits':160,
     'rnbits':158,
     'label':"MNT6_160 Scott-Barreto DCC 06"},# 6
    #r = Phi_6(u)/13
    #c = 4 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (4*x^2 - 4*x + 4)/13
    #px = (4*x^2 + 9*x + 4)/13
    {'k':6,
     'u':-0x19b28d88986138151eea7,
     'D':9154385,
     'c':4,
     'd':13,
     'a':-3,
     'b':0xbb09072fe88ce4512d9147a3a633f1277fd958ba,
     'pnbits':160,
     'rnbits':158,
     'label':"MNT6_160 Scott-Barreto DCC 06"},# 7
    #r = Phi_6(u)/3
    #c = 2 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (2*x^2 - 2*x + 2)/3
    #px = (2*x^2 + x + 2)/3
    {'k':6,
     'u':-0xf18622957e05c6336a5b,
     'D':85700746,
     'c':2,
     'd':3,
     'a':-3,
     'b':0x2221ae162e97d174aee464c97290f882ff89024f,
     'pnbits':160,
     'rnbits':159,
     'label':"MNT6_160 Scott-Barreto DCC 06"},# 8
    #r = Phi_6(u)/1
    #c = 1 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = x^2 - x + 1
    #px = x^2 + 1
    {'k':6,
     'u':-0xb50cbcfdca050ce874be,
     'D':1173931627,
     'c':1,
     'd':1,
     'a':-3,
     'b':0x2d851351189ba68859f3a34893d0512b2f025bf6,
     'pnbits':160,
     'rnbits':160,
     'label':"MNT6_160 Scott-Barreto DCC 06"},# 9
    #r = Phi_6(u)/1
    #c = 1 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = x^2 - x + 1
    #px = x^2 + 1
    {'k':6,
     'u':-0xbd9e196498f480b4e92a,
     'D':1175123707,
     'c':1,
     'd':1,
     'a':-3,
     'b':0x299ce219b7b01348fc2b5007b6ab1ee1005676f7,
     'pnbits':160,
     'rnbits':160,
     'label':"MNT6_160 Scott-Barreto DCC 06"},#10
    #r = Phi_6(u)/3
    #c = 2 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = (2*x^2 - 2*x + 2)/3
    #px = (2*x^2 + x + 2)/3
    {'k':6,
     'u':-0x10f0cefd25e66bcdab6b85a05,
     'D':3371809,
     'c':2,
     'd':3,
     'a':-3,
     'b':0x7eeafaf4178e7349192e71fa4eb40c681a11a9b5b4f2c0c9,
     'pnbits':192,
     'rnbits':191,
     'label':"MNT6_192 Scott-Barreto DCC 06"},#11
    #r = Phi_6(u)/1
    #c = 1 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = x^2 - x + 1
    #px = x^2 + 1
    {'k':6,
     'u':0xbf22d3eb2dc505cb471d89b1ed04,
     'D':496659,
     'c':1,
     'd':1,
     'a':-3,
     'b':0x3482e3fe1715cc7bc4e3ee50bcb16f371ba7eedf122a71c29e80ea2,
     'pnbits':224,
     'rnbits':224,
     'label':"MNT6_224 Scott-Barreto DCC 06"},#12
    #r = Phi_6(u)/1
    #c = 1 # s.t. n=r*c = p+1-tr
    #tx = x + 1
    #nx = x^2 - x + 1
    #px = x^2 + 1
    {'k':6,
     'u':-0xfb1d5f420f3c9a934b42b272ed578dbe,
     'D':56415963,
     'c':1,
     'd':1,
     'a':-3,
     'b':0x6e974d68ef44f266ae3dd5d1f97c497c1d5452d1b074a6c06a25d4e5819ccd1c,
     'pnbits':256,
     'rnbits':256,
     'label':"MNT6_256 Scott-Barreto DCC 06"},#13
]

# u = (p+1-r-1)//2

# c such that r*c = p+1-t and r is prime
# snfs is roughly the same as NFS-Conjugation because p is given by a quadratic polynomial
test_vector_MNT6_with_c = [
    {'k':6,
     'u':-0x53c30a94ecc68a19f799,
     'D':62003,
     'c':3,
     'a':-3,
     'b':0x6ac47cceddc1b6d43bca9f51762c9977b746e3a2,
     'pnbits':159,
     'rnbits':158,
     'cost_S':65,
     'deg_h_S':2,
     'cost_SS':93,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':79,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':78,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':84,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':80,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_159"},# 0
    #STNFS (dh=2,cost=65,Px=4*x^2+1) SNFS (d=4,cost=79) GJL (df=7,d=2,cost=86) (df=7,d=3,cost=84) SSingh (df=9,dg=6,dim=3,cost=80)
    #inv_zeta_Kh, w = 0.860829, 2
    #h = y^2 + y - 1 # [-1, 1, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 + 395554295667128312395673*x^2 + 395554295667128312395670*x - 1
    #max_fi =     37, log_2 max_gi =  78.39
    #aut = 6
    #alpha_f = 0.3924 alpha_g = 2.3023 sum = 2.6947
    {'k':6,
     'u':0x19c8dc2b799b9aff7d9b11,
     'D':3447443,
     'c':397,
     'a':-3,
     'b':0x32a6bf082e582c2d246f62f3fb90e93b0069b5701f8,
     'pnbits':172,
     'rnbits':163,
     'cost_S':68,
     'deg_h_S':2,
     'cost_SS':96,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':81,
     'sieving_dim_SNFS':4,
     'cost_NFS_Conj':80,
     'sieving_dim_Conj':4,
     'cost_NFS_GJL':87,
     'GJL_deg_f':7,
     'sieving_dim_GJL':3,
     'cost_NFS_SS':82,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_172"},# 1
    #STNFS (dh=2,cost=68,Px=4*x^2+1) SNFS (d=4,cost=81) GJL (df=7,d=2,cost=88) (df=7,d=3,cost=87) SSingh (df=9,dg=6,dim=3,cost=82)
    #inv_zeta_Kh, w = 0.696878, 2
    #h = y^2 - 2 # [-2, 0, 1] class number Kh = 1
    #f = 4*x^6 - 23*x^4 - 6*x^3 + 37*x^2 + 24*x + 4
    #g = x^3 - 31171680203341980338330385*x^2 - 31171680203341980338330388*x - 1
    #max_fi =     37, log_2 max_gi =  84.69
    #aut = 6
    #alpha_f = -0.5586 alpha_g = 2.2617 sum = 1.7032    
    {'k':6,
     'u':-0x1943aab18399792ff771e8800f26d832,
     'D':205483,
     'c':4759,
     'a':-3,
     'b':0x4dd4046252d362f8bd116c51a5724fa01eca6d2a52488bd087d5c140dde5aec,
     'pnbits':252,
     'rnbits':240,
     'cost_S':81,
     'deg_h_S':2,
     'cost_SS':106,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':96,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':95,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':100,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':97,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':3,
     'label':"MNT6_252"},# 2
    #STNFS (dh=2,cost=81,Px=4*x^2+1) SNFS (d=3,cost=96) (d=4,cost=97) (more rels with d=3) GJL (df=7,d=2,cost=100) SSingh (df=9,dg=6,dim=3,cost=97)    
    {'k':6,
     'u':-0x4fe91b90242b3d996f19c6ad3bbc48e55fd8d6b3eb8ab728754,
     'D':238859,
     'c':14741701,
     'a':-3,
     'b':0x26bfa208a02f469e8951497ef275d210ac2120046d4cfe6e12e31264015e378cad80e58fb8920a5a0747bc334306c74ee47808,
     'pnbits':407,
     'rnbits':383,
     'cost_S':100,
     'deg_h_S':2,
     'cost_SS':118,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':116,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':115,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':118,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':119,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_407"},# 3
    #STNFS (dh=2,cost=100,Px=4*x^2+1) SNFS (d=3,cost=116) GJL (df=8,d=2,cost=126) (df=8,d=3,cost=>130) (df=7,d=2,cost=118) (df=7,d=3,cost>=130) SSingh (df=9,dg=6,dim=2,cost=119) (df=9,dg=6,dim=3,cost=121)
    {'k':6,
     'u':-0x26287104305513a2e0692b2c31dde349ee9e0589bbead38e41056d7,
     'D':140963,
     'c':152841,
     'a':-3,
     'b':0x1492f502f5082d4ba6139caaac83f874a62313a53fe9771d0748e32458b96f83c5861265b1d3ae224c6e52cf0a6be0a8d28e7b87067156,
     'pnbits':437,
     'rnbits':420,
     'cost_S':104,
     'deg_h_S':2,
     'cost_SS':122,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':120,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':119,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':122,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':122,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_437"},# 4
    #STNFS (dh=2,cost=104,Px=4*x^2+1) SNFS (d=3,cost=120) GJL (df=8,d=2,cost=129) (df=7,d=2,cost=122) SSingh (df=9,dg=6,dim=2,cost=122) (df=9,dg=6,dim=3,cost=125)    
    {'k':6,
     'u':-0xad6704de3eb4ee99ad4ec3d7e74dced8f628f5d11be3321af96cdfc045a,
     'D':2212219,
     'c':392551,
     'a':-3,
     'b':0x14280481f28d70c450aa814670a041e3f11161089f9036959e435d0d1f8d1ad3364ad179434baa56fa4598b63e0195aab4f89fbcc1f3290372d03f4,
     'pnbits':473,
     'rnbits':455,
     'cost_S':108,
     'deg_h_S':2,
     'cost_SS':124,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':124,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':123,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':125,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':125,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_473"},# 5
    #STNFS (dh=2,cost=108,Px=4*x^2+1) SNFS (d=3,cost=124) GJL (df=8,d=2,cost=132) (df=7,d=2,cost=125) SSingh (df=9,dg=6,dim=2,cost=125)
    {'k':6,
     'u':-0x361a6a3f4d378270ea50daf96acb3eb05bef23fd6a477b56a62b138c99ceb,
     'D':873867,
     'c':16462302129913,
     'a':-3,
     'b':0x291143ac51d2313e798f87b2fb20df39e23d56fc36adda36862d808ac2da90ae0c44daae8b2fe4614338333ccfd316c3341232283ceb652595f519d15a,
     'pnbits':486,
     'rnbits':442,
     'cost_S':109,
     'deg_h_S':2,
     'cost_SS':126,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':125,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':124,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':126,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':126,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_486"},# 6
    #STNFS (dh=2,cost=109,Px=4*x^2+1) SNFS (d=3,cost=125) GJL (df=8,d=2,cost=134) (df=7,d=2,cost=126) SSingh (df=9,dg=6,dim=2,cost=126)
    {'k':6,
     'u':-0xc4df9bbc061c04db96ba9b8e881a4043247086574b86abb5774175dc3931763d2,
     'D':311387,
     'c':201,
     'a':-3,
     'b':0xdab3542f99b51fc01804e6217969eb06f7986a46e755d53b3bc500bc0987234c7a3911ed07d3cea444b8a2c2ba2bd5bd166956e45f4f215a41157090aef02d117c,
     'pnbits':522,
     'rnbits':514,
     'cost_S':113,
     'deg_h_S':2,
     'cost_SS':128,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':129,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':128,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':130,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':129,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_522"},# 7
    #STNFS (dh=2,cost=113,Px=4*x^2+1) SNFS (d=3,cost=129) GJL (df=7,d=2,cost=) SSingh (df=9,dg=6,sieving_dim=2,cost=129)    
    {'k':6,
     'u':-0x28d22dbdcdae41854574b485801eabeefbca4c0fa267752ce63e681e5c18047c1ec62f7988bb004359fba,
     'D':594739,
     'c':135683527,
     'a':-3,
     'b':0x14d8dbec3e3376781071b641465fcc5a5032df103d2848488d27d67e4bc6d0c04415b23c3f61896b30c49e59ab68eb9907b7880396604ab0dfc21b98f51c52b2a81a4c20e99984b213c1277f46381c5b4b66c87cc,
     'pnbits':677,
     'rnbits':650,
     'cost_S':129,
     'deg_h_S':2,
     'cost_SS':141,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':144,
     'sieving_dim_SNFS':3,
     'cost_NFS_Conj':143,
     'sieving_dim_Conj':3,
     'cost_NFS_GJL':144,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':141,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_677"},# 8
    #STNFS (dh=2,cost=129,Px=4*x^2+1) SNFS (d=2,cost=151), (d=3,cost=144) GJL (df=7,d=2,cost=144) SSingh (df=9,dg=6,sieving_dim=2,cost=141)
    {'k':6,
     'u':0x292d3f531712830034c874bd782a42592ed5c3ebfb65542be465891a398b0b66b8aa1db020f3acf8b865e1e7d7ee6c32219002a6589d96d2496435f51390bbe9,
     'D':4325483,
     'c':3229,
     'a':-3,
     'b':0xdb782ddfa2de602e8f05262c487f2a1d3d2bb5c96611d502d3b0caf5c62f3df0f6cfcdaf621c3c2e0dfd1ca60157c379802d4980bd19afd31842252b8dbd876280f06ff03ed7ddb6baae1accab542f1ae8b98fb95d2323503cc3df7a9ece2dd9ef0e64106023b7e8ae43a90b461bed73460fa4d13a46ab8bcda9df450e8f7a,
     'pnbits':1021,
     'rnbits':1010,
     'cost_S':158,
     'deg_h_S':2,
     'cost_SS':164,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':169,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':168,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':171,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':164,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_1021"},# 9
    #STNFS (dh=2,cost=158(try157),Px=4*x^2+1) SNFS (d=2,cost=169) GJL (df=8,d=2,cost=174) (df=7,d=2,cost=171) SSingh (df=8,dg=6,sieving_dim=2,cost=167) (df=10,dg=8,sieving_dim=2,cost>=178) (df=9,dg=6,sdim=2,cost=164), (df=12,dg=9,cost>=189)
    {'k':6,
     'u':-0x6d0b786d21f34fe63357057952aa560999b8b89b5510b06bee68b7640c45ddf7f09621cd10d85deba084d74b87721f35be99607d42078d558ec78d5a693887275a2dfceeb,
     'D':26933363,
     'c':43742199,
     'a':-3,
     'b':0xb58353b0674a7efe3551d519e5d9cb9b66dba486708912220690decb303f54b121bc4cbd4f27c565a7f6dda011e65190ac733eb65eb1142ea66742c2a7b9cd2cd89b0faf24205cff633ceec0d6fbfd41d8e0230fc81bad42eb09e2a0702cbc30eca2bb86b56304231d4f4f07e793f9b32f1ac6e5c2a6834eaa20fd90ad596429257c704f050dd20d54,
     'pnbits':1096,
     'rnbits':1071,
     'cost_S':164,
     'deg_h_S':2,
     'cost_SS':168,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':172,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':171,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':176,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':168,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_1096"},#10
    #STNFS (dh=2,cost=164(try163),Px=4*x^2+1) SNFS (d=2,cost=172) GJL (df=7,d=2,cost=176) (df=8,d=2,cost=179) SSingh (df=9,dg=6,sieving_dim=2,cost=168)
    {'k':6,
     'u':-0x7f8ba6a3b263ed73fc5ed8d34c99620e5310e2e8416e67a7646a9a2023f3342b821f63dac53a531265b8a46444e46d508d16764b69251f3bf9236c0ea215c9b3a88643de5577681f28717f20,
     'D':3211771,
     'c':2191,
     'a':-3,
     'b':0xd21d8f00625efbb48b378c9ccb5a4d7a195fbc026e6c06c39180cec50304cb7edcb00ee29c8c33df7be4db4d4da834e43a6e48e64f7c84a533d7311e0e61d2b5478aaf5b1d5eb283562fb2788b839e68fd57b40d714f5dd9446626f72058a16392fad4b1d0f4d5729c3a3b755a8871544040985a3e27a435cec61233367d27b0cf463630e8cfcf6a2f7cd106a709a15c590caeb9ebea8f5c,
     'pnbits':1216,
     'rnbits':1205,
     'cost_S':172,
     'deg_h_S':2,
     'cost_SS':175,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':178,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':177,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':185,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':176,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_1216"},#11
    #STNFS (dh=2,cost=172,Px=4*x^2+1) SNFS (d=2,cost=178) GJL (df=7,d=2,cost=185) (df=8,d=2,cost=187) SSingh (df=9,dg=6,sieving_dim=2,cost=176)
    {'k':6,
     'u':0xbbcd8209c9a43befabe98c73d2b120c0aefd993ec4b0e9919376490cbc4d538e1d1f1b3e6d6ed20ac42f0bb2ac2ea08c6c2dea6d6e9ea84c9b048ddb9f9ba647b64d7045f5a2290935d379d5229993e3aa183838542be9f4fbd1e4f03a853,
     'D':20685347,
     'c':27453,
     'a':-3,
     'b':0x1525738d13a481fb37644a63922ae3eb1ab6b87030271ff5b0dd7f8e12b9ed5186fe4c5d2aca550572e1bc7e0af686adadb465745c0ea2de5c3fbcfffabde049aff3be1390d7e08d90d6a70e804539723eaf9b43c2ea38823825058f365cf0da7b03349a31ae80ac959902640618d5aeb879c052303643b84108969a47ffde9c7641bc4d628c649bf78caea6dbf15d8630dab24838d0125663812ae3c9e898831c0bdeb86c50a469bf4f0252787ead4043d0e7e55ff1042a7ee0dcaa100,
     'pnbits':1514,
     'rnbits':1499,
     'cost_S':193,
     'deg_h_S':2,
     'cost_SS':192,
     'deg_h_SS':2,
     'deg_aux_SS':3,
     'cost_SNFS':191,
     'sieving_dim_SNFS':2,
     'cost_NFS_Conj':190,
     'sieving_dim_Conj':2,
     'cost_NFS_GJL':203,
     'GJL_deg_f':7,
     'sieving_dim_GJL':2,
     'cost_NFS_SS':191,
     'SS_deg_f':9,
     'SS_deg_g':6,
     'sieving_dim_NFS_SS':2,
     'label':"MNT6_1514"},#12
    #STNFS (dh=2,cost=193,Px=4*x^2+1) SNFS (d=2,cost=191) GJL(df=8,sdim=2,cost=204) (df=7,cost=203)
    # SSingh: (df=9,dg=6,sdim=2,cost=191), (df=3*4=12,dg=3*3=9,cost>=210), (df=2*4=8,dg=2*3=6,cost>=197) (df=2*5=10,dg=2*4=8,cost>=203)
]
