"""
A list of curves and parameters for STNFS

1. Curve parameters from eprint 2019/485, revision of September 2019
2. STNFS parameters from eprint 2019/1371, revision of December 2019

"""
import tnfs
from tnfs.alpha.alpha_tnfs_2d import alpha_TNFS_2d
from tnfs.curve.cyclo_kDe import Cyclo_kDe, polynomial_params
from tnfs.simul.simulation_tnfs import Simulation_TNFS
from tnfs.simul.polyselect_utils import automorphism_factor, divide_degree_by_d, divide_degree_by_two_palindrome, divide_degree_by_two_even
from tnfs.simul.polyselect import Polyselect

ZZy.<y> = ZZ[]
ZZxy.<x> = ZZy[]

samples = 100000

curve_test_vector = [
    # k=9
    # BLS-9
    {'k':9,'D':3,'e0':1,'u':2**43+2**37+2**7+1       ,'a':0,'b':1,'pnbits':343,'rnbits':257,'cofact_r':1,'deg_h_S':9,'cost_S': 99,'div_deg_P':1,'h_coeff_S':[-1,0,1,0,1,-1,0,0,0,1],'inv_zeta_KhS':0.96470018, 'label':"eprint 2016/1187 Sec 4.2 u=2^43+2^37+2^7+1"},
    {'k':9,'D':3,'e0':1,'u':2**70+2**59+2**46+2**41+1,'a':0,'b':1,'pnbits':559,'rnbits':419,'cofact_r':1,'deg_h_S':9,'cost_S':125,'div_deg_P':1,'h_coeff_S':[-1,-1,-1,-1,0,0,1,1,0,1],'inv_zeta_KhS':0.98048671, 'label':"eprint 2016/1187 Sec 4.2 u=2^70+2^59+2^46+2^41+1"},
    {'k':9,'D':3,'e0':1,'u':2**74+2**35-2**22+2      ,'a':0,'b':9,'pnbits':591,'rnbits':443,'cofact_r':1,'deg_h_S':9,'cost_S':128,'div_deg_P':1,'h_coeff_S':[-1,1,-1,0,1,-1,1,0,0,1],'inv_zeta_KhS':0.97538522, 'label':"eprint 2019/485 Tab 23 u=2^74+2^35-2^22+2"},
    # FST 6.6(-x)
    {'k':9,'D':3,'e0':7,'u':0x92d789aaa947259ee,'a':0,'b':9,'pnbits':536,'rnbits':402,'cofact_r':1,'deg_h_S':9,'cost_S':122,'div_deg_P':1,'h_coeff_S':[-1,-1,0,-1,0,0,0,1,0,1],'inv_zeta_KhS':0.98075428, 'label':"Cyclo_k9_D3_e7_FST66_p536"},
    {'k':9,'D':3,'e0':7,'u':0x7b7a96df6a30d1b5f,'a':0,'b':-3,'pnbits':535,'rnbits':401,'cofact_r':1,'deg_h_S':9,'cost_S':122,'div_deg_P':1,'h_coeff_S':[-1,-1,0,0,-1,0,1,0,0,1],'inv_zeta_KhS':0.95365230, 'label':"Cyclo_k9_D3_e7_p535"},
    {'k':9,'D':3,'e0':7,'u':0x86a79794f28a0dce8,'a':0,'b':-3,'pnbits':535,'rnbits':401,'cofact_r':1,'deg_h_S':9,'cost_S':122,'div_deg_P':1,'h_coeff_S':[-1,0,-1,-1,0,-1,0,0,0,1],'inv_zeta_KhS':0.86783691, 'label':"Cyclo_k9_D3_e7_p535"},
    {'k':9,'D':3,'e0':7,'u':-0x7b7a96df6a30d4a8b,'a':0,'b':9,'pnbits':535,'rnbits':401,'cofact_r':1,'deg_h_S':9,'cost_S':122,'div_deg_P':1,'h_coeff_S':[-1,0,0,1,-1,0,-1,0,0,1],'inv_zeta_KhS':0.94271892, 'label':"Cyclo_k9_D3_e7_p535"},
    {'k':9,'D':3,'e0':7,'u':-0x86a79794f28a0d296,'a':0,'b':7,'pnbits':535,'rnbits':401,'cofact_r':1,'deg_h_S':9,'cost_S':122,'div_deg_P':1,'h_coeff_S':[-1,0,1,1,1,-1,-1,0,0,1],'inv_zeta_KhS':0.98155567, 'label':"Cyclo_k9_D3_e7_p535"},
    # FST 6.2
    {'k':9,'D':1,'e0':5,'u':0x400637,'a': 7,'b':0,'pnbits':483,'rnbits':265,'cofact_r':1,'deg_h_S':9,'cost_S':116,'div_deg_P':2,'h_coeff_S':[-1,0,1,1,1,-1,-1,0,0,1],'inv_zeta_KhS':0.98155567, 'label':"eprint 2019/485 Table 6 u=-1+2^3+2^4+2^5+2^9+2^10+2^22"},# Special P(x^2)=4*p(x)
    {'k':9,'D':1,'e0':5,'u':0x42152d,'a': 7,'b':0,'pnbits':484,'rnbits':265,'cofact_r':1,'deg_h_S':9,'cost_S':116,'div_deg_P':2,'h_coeff_S':[-1,-1,-1,-1,0,0,1,1,0,1],'inv_zeta_KhS':0.98048671, 'label':"Cyclo_k9_D1_e5_p484"},
    {'k':9,'D':1,'e0':5,'u':0x4419b1,'a':17,'b':0,'pnbits':484,'rnbits':266,'cofact_r':1,'deg_h_S':9,'cost_S':116,'div_deg_P':2,'h_coeff_S':[-1,-1,0,-1,0,0,0,1,0,1],'inv_zeta_KhS':0.98075428, 'label':"Cyclo_k9_D1_e5_p484"},
    # FST 6.7
    {'k':9,'D':2,'e0':1,'u':-0xa2f,'a':-30,'b':56,'pnbits':520,'rnbits':273,'cofact_r':1,'deg_h_S':9,'cost_S':140,'div_deg_P':4,'h_coeff_S':[-1,0,1,0,1,-1,0,0,0,1],'inv_zeta_KhS':0.96470018, 'label':"eprint 2019/485 Table 19 u=-2^11-2^9-2^6+2^4+1"},# Special P(u^4) = 8*p(u)
    # k=10
    # FST 6.3
    {'k':10,'D':1,'e0':1,'u': 0x800023e9,'a':1,'b': 0,'pnbits':433,'rnbits':249,'cofact_r':1,'deg_h_S':10,'cost_S':120,'div_deg_P':2,'h_coeff_S':[-1,0,-1,-1,0,0,0,1,0,0,1],'inv_zeta_KhS':0.97128016, 'label':"eprint 2019/485 Table 7 u=1+2^3-2^5+2^10+2^13+2^31"},
    {'k':10,'D':1,'e0':1,'u': 0xffffc92b,'a':1,'b': 0,'pnbits':446,'rnbits':256,'cofact_r':1,'deg_h_S':10,'cost_S':121,'div_deg_P':2,'h_coeff_S':[-1,0,0,0,-1,1,1,0,0,0,1],'inv_zeta_KhS':0.96293584, 'label':"FST 6.3 p446"},
    # FST 6.6
    {'k':10,'D':3,'e0':1,'u':-0xfb705096,'a':0,'b':-2,'pnbits':511,'rnbits':256,'cofact_r':1,'deg_h_S':10,'cost_S':150,'div_deg_P':2,'h_coeff_S':[1,0,-1,0,-1,0,0,0,-1,0,1],'inv_zeta_KhS':0.94900778, 'label':"FST 6.6 p511"},
    # k=11
    # FST 6.2
    {'k':11,'D':1,'e0':6,'u':0x1f63,'a':19,'b':0,'pnbits':336,'rnbits':249,'cofact_r':1409,'deg_h_S':11,'cost_S':118,'div_deg_P':2,'h_coeff_S':[-1,0,-1,-1,0,0,0,0,1,0,0,1],'inv_zeta_KhS':0.96968581, 'label':"FST 6.2 p336"},
    {'k':11,'D':1,'e0':6,'u':0x2017,'a':11,'b':0,'pnbits':337,'rnbits':189,'cofact_r':3986392263936102701341,'deg_h_S':11,'cost_S':118,'div_deg_P':2,'h_coeff_S':[-1,0,-1,-1,0,0,0,0,1,0,0,1],'inv_zeta_KhS':0.96968581,'label':"FST 6.2 u=2^13+2^5-2^3-1"},
    {'k':11,'D':1,'e0':6,'u':0x40ff,'a':19,'b':0,'pnbits':363,'rnbits':281,'cofact_r':1,   'deg_h_S':11,'cost_S':122,'div_deg_P':2,'h_coeff_S':[1,1,1,0,1,0,0,0,0,0,0,1],'inv_zeta_KhS':0.94992919, 'label':"eprint 2019/485 Table 6 u=2^14+2^8-1 FST62 p363"},
    # FST 6.6
    {'k':11,'D':3,'e0':4,'u': 0x17c2,'a':0,'b':4,'pnbits':301,'rnbits':252,'cofact_r':  1,'deg_h_S':11,'cost_S':115,'div_deg_P':2,'h_coeff_S':[1,0,0,-1,0,-1,0,1,0,0,0,1],'inv_zeta_KhS':0.92373456, 'label':"u=2^13-2^11-2^6+2"},
    {'k':11,'D':3,'e0':4,'u':-0x208a,'a':0,'b':4,'pnbits':311,'rnbits': 78,'cofact_r':12570299770722199504710162467421029765488569491384377561,'deg_h_S':11,'cost_S':114,'div_deg_P':2,'h_coeff_S':[-1,-1,0,0,1,0,0,0,0,1,0,1],'inv_zeta_KhS':0.98231804,'label':"FST 6.6 p311 u=-2^13-2^7-2^3-2"},
    {'k':11,'D':3,'e0':4,'u':-0x2282,'a':0,'b':4,'pnbits':314,'rnbits':257,'cofact_r': 67,'deg_h_S':11,'cost_S':115,'div_deg_P':2,'h_coeff_S':[-1,0,-1,-1,0,0,0,0,1,0,0,1],'inv_zeta_KhS':0.96968581, 'label':"FST 6.6 p314 u=-2^13-2^9-2^7-2"},
    {'k':11,'D':3,'e0':4,'u': 0x2470,'a':0,'b':4,'pnbits':315,'rnbits':255,'cofact_r':463,'deg_h_S':11,'cost_S':115,'div_deg_P':2,'h_coeff_S':[1,-1,0,0,0,1,0,-1,0,0,0,1],'inv_zeta_KhS':0.97056943,'label':"u=2^13+2^10+2^7-2^4"},
    {'k':11,'D':3,'e0':4,'u':-0x46d0,'a':0,'b':4,'pnbits':338,'rnbits':283,'cofact_r':  1,'deg_h_S':11,'cost_S':118,'div_deg_P':2,'h_coeff_S':[-1,0,0,0,-1,-1,1,0,0,0,0,1],'inv_zeta_KhS':0.96693036, 'label':"eprint 2019/485 Table 16 u=-(2^4+2^6+2^7+2^9+2^10+2^14)=-2^14-2^11+2^8+2^6-2^4 FST 6.6 p338"},
    # k=12
    # BN
    {'k':12,'D':3,'e0':0,'u':2**114+2**101-2**14-1,'a':0,'b':-4,'pnbits':462,'rnbits':462,'cofact_r':1,'deg_h_S':6,'cost_S':135,'div_deg_P':1,'h_coeff_S':[1,1,0,1,0,1,1],'inv_zeta_KhS':0.88441183,'label':"Barbulescu-Duquesne JoC 18"},
    # BLS12
    {'k':12,'D':3,'e0':1,'u':-2**77-2**59+2**9, 'a':0,'b': 4,'pnbits':461,'rnbits':309,'cofact_r':1,'cost_S':135,'deg_h_S':6,'div_deg_P':1,'h_coeff_S':[1,0,-2,0,-1,0,1],'inv_zeta_KhS':0.96356914, 'label':"Barbulescu-Duquesne JoC 18"},
    {'k':12,'D':3,'e0':1,'u':-2**77+2**50+2**33,'a':0,'b': 4,'pnbits':461,'rnbits':308,'cofact_r':1,'cost_S':134,'deg_h_S':6,'div_deg_P':1,'h_coeff_S':[-1,1,0,0,0,0,1],'inv_zeta_KhS':0.93897437, 'label':"Barbulescu-Duquesne JoC 18"},
    # FST 6.7
    {'k':12,'D':2,'e0':1,'u':0xfffff9d7, 'a':-30,'b':56,'pnbits':445,'rnbits':256,'cofact_r':1,'deg_h_S':12,'cost_S':134,'div_deg_P':2,'h_coeff_S':[-1,-1,0,0,-1,0,0,0,0,0,0,1,1],'inv_zeta_KhS':0.91612617, 'label':"Cyclo_k12_D2_e1_p445"},
    {'k':12,'D':2,'e0':1,'u':0x100024001,'a':-30,'b':56,'pnbits':446,'rnbits':257,'cofact_r':1,'deg_h_S':12,'cost_S':134,'div_deg_P':2,'h_coeff_S':[1,-1,0,0,-1,0,0,0,0,1,0,0,1],'inv_zeta_KhS':0.98717286, 'label':"eprint 2019/485 u=1+2^14+2^17+2^32"},
    # FST 6.4
    {'k':12,'D':1,'e0':1,'u':0xffffffffffffff03, 'a': 2,'b':0,'pnbits':510,'rnbits':256,'cofact_r':1,'deg_h_S':12,'cost_S':138, 'div_deg_P':1,'h_coeff_S':[-1,0,0,-1,0,1,0,0,0,1,0,0,1],'inv_zeta_KhS':0.95100795,'label':"Cyclo_k12_D4_e1_p510"},
    {'k':12,'D':1,'e0':1,'u':0x10000000000000b0b,'a': 6,'b':0,'pnbits':511,'rnbits':257,'cofact_r':1,'deg_h_S':12,'cost_S':138, 'div_deg_P':1,'h_coeff_S':[-1,1,0,0,-1,1,0,0,0,0,0,0,1],'inv_zeta_KhS':0.83929621, 'label':" eprint 2019/485 Table 8u=1+2+2^3+2^8+2^9+2^11+2^64"},
    # k=13
    # FST 6.2
    {'k':13,'D':1,'e0':7,'u':0x801,'a':51,'b':0,'pnbits':329,'rnbits':218,'cofact_r':102564562798321,'deg_h_S':13,'cost_S':140,'div_deg_P':4,'h_coeff_S':[-1,0,0,0,0,-1,1,0,1,0,0,0,0,1],'inv_zeta_KhS':0.94074586, 'label':"u=-2^11-1"},
    {'k':13,'D':1,'e0':7,'u':0x10451b,'a':17,'b':0,'pnbits':599,'rnbits':481,'cofact_r':1, 'deg_h_S':13,'cost_S':162,'div_deg_P':2,'h_coeff_S':[-1,0,0,-1,0,-1,0,0,0,0,1,0,0,1],'inv_zeta_KhS':0.94377139,'label':"eprint 2019/485 u=1+2+2^3+2^4+2^8+2^10+2^14+2^20=2^20+2^14+2^10+2^8+2^5-2^2-1"},
    # FST 6.6
    {'k':13,'D':3,'e0':9,'u': -0x863,'a':0,'b':9,'pnbits':309,'rnbits':218,'cofact_r':407609994343423,'deg_h_S':13,'cost_S':139,'div_deg_P':3,'h_coeff_S':[-1,0,0,0,1,-1,0,0,1,0,0,0,0,1],'inv_zeta_KhS':0.94692825, 'label':"u=-2^11-2^7+2^5-2^2+1"},
    {'k':13,'D':3,'e0':9,'u':-0x2c90,'a':0,'b':9,'pnbits':376,'rnbits':324,'cofact_r':1, 'deg_h_S':13,'cost_S':152,'div_deg_P':3,'h_coeff_S':[-1,0,-1,1,0,0,0,0,1,0,0,0,0,1],'inv_zeta_KhS':0.87011688, 'label':"u=-2^14+2^12+2^10-2^7-2^4"},
    # k=14
    # FST 6.3
    {'k':14,'D':1,'e0':1,'u':0x37723d,'a':1,'b':0,'pnbits':391,'rnbits':262,'cofact_r':1, 'deg_h_S':14,'cost_S':131,'div_deg_P':2,'h_coeff_S':[-1,0,0,0,1,-1,0,0,0,1,0,0,0,0,1],'inv_zeta_KhS':0.95645440, 'label':"eprint 2019/485 u=1-2^2+2^6+2^9-2^12-2^15-2^19+2^22"},
    {'k':14,'D':1,'e0':1,'u': 0x3ff94d,'a':1,'b':0,'pnbits':394,'rnbits':264,'cofact_r':1, 'deg_h_S':14,'cost_S':132,'div_deg_P':2,'h_coeff_S':[1,0,0,1,0,0,-1,0,-1,0,0,0,0,0,1],'inv_zeta_KhS':0.97739583, 'label':"Cyclo_k14_D1_e1_p394"},
    # FST 6.6
    {'k':14,'D':3,'e0':5,'u':0x4226bf,'a':0,'b':2,'pnbits':352,'rnbits':265,'cofact_r':1, 'deg_h_S':14,'cost_S':150,'div_deg_P':1,'h_coeff_S':[1,0,0,0,0,0,-1,0,0,0,0,0,0,0,1],'inv_zeta_KhS':0.96524213, 'label':"eprint 2019/485 -1+2^6+2^7+2^9+2^10+2^13+2^17+2^22 (tab.15) = 2^22+2^17+2^13+2^11-2^8-2^6-1 (HW-2-NAF)"},
    {'k':14,'D':3,'e0':5,'u':-0x3f1e6e,'a':0,'b': 2,'pnbits':351,'rnbits':264,'cofact_r':1, 'deg_h_S':14,'cost_S':150,'div_deg_P':1,'h_coeff_S':[1,0,-1,0,0,0,1,0,0,0,-1,0,0,0,1],'inv_zeta_KhS':0.93175030, 'label':"Cyclo_k14_D3_e5_p351"},
    {'k':14,'D':3,'e0':5,'u': 0x3ee129,'a':0,'b':11,'pnbits':351,'rnbits':264,'cofact_r':1, 'deg_h_S':14,'cost_S':151,'div_deg_P':1,'h_coeff_S':[-1,-1,0,1,1,0,0,0,0,0,0,0,0,0,1],'inv_zeta_KhS':0.91414686, 'label':"Cyclo_k14_D3_e5_p351"},
    {'k':14,'D':3,'e0':5,'u': 0x4190d7,'a':0,'b': 2,'pnbits':351,'rnbits':265,'cofact_r':1, 'deg_h_S':14,'cost_S':150,'div_deg_P':1,'h_coeff_S':[1,0,-1,0,0,0,1,0,0,0,0,0,-1,0,1],'inv_zeta_KhS':0.98515689, 'label':"Cyclo_k14_D3_e5_p351"},
    {'k':14,'D':3,'e0':5,'u':-0x41a12e,'a':0,'b':-4,'pnbits':351,'rnbits':265,'cofact_r':1, 'deg_h_S':14,'cost_S':151,'div_deg_P':1,'h_coeff_S':[1,0,0,0,0,0,-1,0,0,0,0,0,0,0,1],'inv_zeta_KhS':0.96524213, 'label':"Cyclo_k14_D3_e5_p351"},
    # k=15
    # BLS
    {'k':15,'D':3,'e0':1,'u': 0x80080024,'a':0,'b':2,'pnbits':371,'rnbits':249,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[1,0,1,0,0,1,0,0,0,0,1,0,0,0,0,1],'inv_zeta_KhS':0.9455, 'label':"eprint 2016/1187 Sec 8.1 2^2+2^5+2^19+2^31"},
    {'k':15,'D':3,'e0':1,'u':0x100090402,'a':0,'b':3,'pnbits':383,'rnbits':257,'cofact_r':1, 'deg_h_S':15,'cost_S':138,'div_deg_P':1,'h_coeff_S':[-1,-1,0,0,0,0,0,1,0,0,1,0,0,0,0,1],'inv_zeta_KhS':0.9553, 'label':"eprint 2019/485 2+2^10+2^16+2^19+2^32"},
    {'k':15,'D':3,'e0':1,'u': 0xeac0d34e,'a':0,'b':2,'pnbits':381,'rnbits':256,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[-1,1,0,0,0,0,0,-1,1,0,0,0,0,0,0,1],'inv_zeta_KhS':0.9476, 'label':"Cyclo_k15_D3_e1_p381"},
    {'k':15,'D':3,'e0':1,'u':-0xeac0f212,'a':0,'b':2,'pnbits':381,'rnbits':256,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[-1,1,0,0,0,-1,0,0,1,0,0,0,0,0,0,1],'inv_zeta_KhS':0.9488, 'label':"Cyclo_k15_D3_e1_p381"},
    {'k':15,'D':3,'e0':1,'u':-0xfffff9ec,'a':0,'b':2,'pnbits':383,'rnbits':256,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[1,-1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1],'inv_zeta_KhS':0.9676, 'label':"Cyclo_k15_D3_e1_p383"},
    {'k':15,'D':3,'e0':1,'u': 0xffffee69,'a':0,'b':1,'pnbits':383,'rnbits':256,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1],'inv_zeta_KhS':0.8732, 'label':"Cyclo_k15_D3_e1_p383"},
    {'k':15,'D':3,'e0':1,'u':0x108cc1b8b,'a':0,'b':1,'pnbits':383,'rnbits':257,'cofact_r':1, 'deg_h_S':15,'cost_S':138,'div_deg_P':1,'h_coeff_S':[1,-1,1,0,0,0,0,1,0,0,0,0,0,0,0,1],'inv_zeta_KhS':0.8655, 'label':"Cyclo_k15_D3_e1_p383"},
    # FST 6.6
    {'k':15,'D':3,'e0':11,'u':-0x100011005,'a':0,'b':15,'pnbits':383,'rnbits':257,'cofact_r':1, 'deg_h_S':15,'cost_S':138,'div_deg_P':1,'h_coeff_S':[1,0,0,-1,-1,0,0,0,0,0,1,0,0,0,0,1],'inv_zeta_KhS':0.9808, 'label':"eprint 2019/485 u=1+2^2+2^12+2^16+2^32"},
    {'k':15,'D':3,'e0':11,'u': 0xeac0dd8f,'a':0,'b':67,'pnbits':381,'rnbits':256,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[1,0,0,1,0,0,0,-1,0,-1,0,0,0,0,0,1],'inv_zeta_KhS':0.929, 'label':"Cyclo_k15_D3_e11_p381"},
    {'k':15,'D':3,'e0':11,'u':-0xeac0cea8,'a':0,'b':13,'pnbits':381,'rnbits':256,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[-1,0,-1,0,0,0,1,0,0,0,1,0,0,0,0,1],'inv_zeta_KhS':0.864, 'label':"Cyclo_k15_D3_e11_p381"},
    {'k':15,'D':3,'e0':11,'u':-0xffffc7fa,'a':0,'b':-7,'pnbits':383,'rnbits':256,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[-1,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,1],'inv_zeta_KhS':0.946, 'label':"Cyclo_k15_D3_e11_p383"},
    {'k':15,'D':3,'e0':11,'u': 0xffffed6d,'a':0,'b': 7,'pnbits':383,'rnbits':256,'cofact_r':1, 'deg_h_S':15,'cost_S':137,'div_deg_P':1,'h_coeff_S':[-1,1,0,0,0,0,0,0,0,0,1,0,0,-1,0,1],'inv_zeta_KhS':0.924, 'label':"Cyclo_k15_D3_e11_p383"},
    # k=16
    # KSS16
    {'k':16,'D':1,'e0':0,'u':-2**34+2**27-2**23+2**20-2**11+1,'a':1,'b':0,'pnbits':330,'rnbits':257,'cofact_r':1,'cost_S':140,'deg_h_S':16,'div_deg_P':1,'h_coeff_S':[1, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 1],'inv_zeta_KhS':0.86627556, 'label':"Barbulescu-Duquesne JoC 18"},
    # k=17
    # FST 6.2
    {'k':17,'D':1,'e0':9,'u': 0x105,'a':11,'b':0, 'pnbits':304,'rnbits':135, 'cofact_r':9319010681331192024888780409936788353,'deg_h_S':17,'cost_S':153,'div_deg_P':4,'h_coeff_S':[1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,0,0,0,1],'inv_zeta_KhS':0.9725, 'label':"u=+2^8+2^2+1"},
    {'k':17,'D':1,'e0':9,'u': 0x24b,'a':47,'b':0, 'pnbits':348,'rnbits':243, 'cofact_r':4798889260244533,   'deg_h_S':17,'cost_S':162,'div_deg_P':4,'h_coeff_S':[-1,-1,0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0,1],'inv_zeta_KhS':0.972, 'label':"u=+2^9+2^6+2^4-2^2-1"},
    {'k':17,'D':1,'e0':9,'u': 0x2c7,'a':43,'b':0, 'pnbits':359,'rnbits':254, 'cofact_r':1130990220889301,   'deg_h_S':17,'cost_S':163,'div_deg_P':4,'h_coeff_S':[1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,1],'inv_zeta_KhS':0.9128, 'label':"u=+2^10-2^8-2^6+2^3-1"},
    {'k':17,'D':1,'e0':9,'u': 0x44d,'a':39,'b':0, 'pnbits':382,'rnbits':262, 'cofact_r':4340330150541605573,'deg_h_S':17,'cost_S':167,'div_deg_P':4,'h_coeff_S':[-1,-1,0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0,1],'inv_zeta_KhS':0.973, 'label':"u=+2^10+2^6+2^4-2^2+1"},
    # FST 6.6
    {'k':17,'D':3,'e0':6,'u':-0x344,'a':0,'b':4,'pnbits':348,'rnbits':249,'cofact_r':5338718021784347281,'deg_h_S':17,'cost_S':168, 'div_deg_P':4,'h_coeff_S':[-1,0,0,0,0,0,0,1,1,-1,0,0,0,0,0,0,0,1],'inv_zeta_KhS':0.9761,'label':"u=-2^10+2^8-2^6-2^2"},
]

for item in curve_test_vector:
    k,D,e0,u,a,b,pnbits,rnbits,cofact_r,deg_h_S,cost_S,div_deg_P,h_coeff_S,inv_zeta_KhS,label=item['k'],item['D'],item['e0'],item['u'],item['a'],item['b'],item['pnbits'],item['rnbits'],item['cofact_r'],item['deg_h_S'],item['cost_S'],item['div_deg_P'],item['h_coeff_S'],item['inv_zeta_KhS'],item['label']
    inv_zeta_KhS = float(inv_zeta_KhS)
    print("")
    print("k={} D={} e0={}".format(k,D,e0))
    print(label)
    E = Cyclo_kDe(k,D,e0, u=u, a=a,b=b, cofactor_r=cofact_r)
    px, rx, tx, cx, yx, betax, lambx, D = polynomial_params(k,D,e0)
    E.print_parameters()
    p = E.p()
    r = E.r()
    Fp = E.Fp()
    Fpz,z = E.Fpz()
    print("")
    px_denom = lcm([ci.denom() for ci in px.list()])
    rx_denom = lcm([ci.denom() for ci in rx.list()])
    tx_denom = lcm([ci.denom() for ci in tx.list()])
    if px_denom > 1:
        print("px = ({})/{}".format(px*px_denom, px_denom))
    else:
        print("px = {}".format(px))
    if rx_denom > 1:
        print("rx = ({})/{}".format(rx*rx_denom, rx_denom))
    else:
        print("rx = {}".format(rx))
    if tx_denom > 1:
        print("tx = ({})/{}".format(tx*tx_denom, tx_denom))
    else:
        print("tx = {}".format(tx))
    poly_p = ZZy(px_denom*px)
    poly_P = poly_p
    poly_u = [-u, 1]
    # special modification of polynomial p(x) --> reduce its degree while increasing its coefficients
    if div_deg_P > 1:
        cp = poly_p.list()
        if div_deg_P == 2 and tnfs.simul.polyselect_utils.is_even(cp):
            cP, poly_u = divide_degree_by_two_even(cp, u)
            poly_P = ZZy(cP)            
            print("poly P = {} where y = x^2".format(poly_P))
        elif div_deg_P == 2 and tnfs.simul.polyselect_utils.is_palindrome(poly_p.list()):
            cP, poly_u = divide_degree_by_two_palindrome(poly_p.list(), u)
            poly_P = ZZy(cP)
            print("poly P = {} where y = x+1/x".format(poly_P))
        elif div_deg_P == 4 and tnfs.simul.polyselect_utils.is_palindrome(poly_p.list()):
            cP, poly_u = divide_degree_by_two_palindrome(poly_p.list(), u)
            poly_P = ZZy(cP)
            # now we have P(x+1/x)*x^deg(P) = p(x), we need P((u+1/u)^2)*u^i = p(u)
            print("poly P = {} where y = (x+1/x)".format(poly_P))
            print("poly u = {} = u*(x-(u+1/u)) = u*x - u^2-1".format(poly_u))
            ccP, poly_u = divide_degree_by_d(cP, 2, (u+1/u))
            ccP = [u*ci for ci in ccP]
            poly_P = ZZy(ccP)
            print("poly P = {} where y = (u+1/u)^2".format(poly_P))
            poly_u = [-(u**2+1)**2, u**2]
            print("poly u = {} = u^2*(x-(u+1/u)^2) = u^2*x - u^4-2*u^2-1".format(poly_u))
        else:
            cP, poly_u = divide_degree_by_d(poly_p.list(), div_deg_P, u)
            poly_P = ZZy(cP)
            print("poly P = {} where y = u^{}".format(poly_P, div_deg_P))
        assert (poly_P.resultant(ZZy(poly_u)) % p) == 0

    if (len(h_coeff_S) == 0) or (inv_zeta_KhS == 0):
        print("no valid poly h given")
        continue
    h = ZZy(h_coeff_S)
    Kh = NumberField(h, 'ah')
    number_roots_unity = len(Kh.roots_of_unity())
    w = number_roots_unity
    deg_h = h.degree()
    deg_g=(k//h.degree())
    ps = Polyselect(E)
    res = ps.Resultant(deg_g, poly_P, poly_u, h=h, with_y=(gcd(deg_g,deg_h) > 1))
    f, g, max_fi, max_gi, aut_fg = res[:5]

    clKh = Kh.class_number()
    if clKh == 1:
        alpha_f = float(alpha_TNFS_2d(f,h,1000,test_principal=True))
        alpha_g = float(alpha_TNFS_2d(g,h,1000,test_principal=True))
    else:
        alpha_f = float(alpha_TNFS_2d(f,h,1000,test_principal=False))
        alpha_g = float(alpha_TNFS_2d(g,h,1000,test_principal=False))
    sum_alpha = alpha_f + alpha_g
    aut = g.degree()
    if gcd(h.degree(), k//h.degree()) == 1:
        aut *= automorphism_factor(h.list())
    print("f = {}".format(f))
    print("g = {}".format(g))
    print("max_fi = {:6d}, log_2 max_gi = {:6.2f}".format(max_fi, float(log(max_gi,2))))
    print("aut = {}".format(aut))
    print("alpha_f = {:.4f} alpha_g = {:.4f}".format(alpha_f,alpha_g))
    print("    ({:.8f},{:2d}, {:40s}, {:.4f}, {:.4f}, {:.4f}),".format(float(inv_zeta_KhS),int(w),str(h),float(alpha_f),float(alpha_g),float(sum_alpha)))
    print("{} random samples".format(samples))
    simul = Simulation_TNFS(p,r,Fp,Fpz,h,f,g,ZZxy,cost_S,aut, inv_zeta_Kh=inv_zeta_KhS, alpha_f=alpha_f, alpha_g=alpha_g, weight_per_row=200, cst_filtering=20, count_sieving=True)
    simul.print_params()
    simul.simulation(samples=samples)
    simul.print_results()

    print("\n# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")
