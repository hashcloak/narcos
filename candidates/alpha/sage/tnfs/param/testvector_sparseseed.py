"""
Seeds of pairing-friendly curves, BN, BLS12, BLS24, BLS48, Fotiadis-Martindale, KSS16, KSS18, KSS32, KSS36, KSS40 and KSS54
"""

# parameters u from Table 1 in https://eprint.iacr.org/2010/429
test_vector_bn =[
    {'u':              -(2**38+2**28+1), 'b':  17, 'pnbits': 158,  'rnbits':158, 'cost_S': 86, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':        -(2**46+2**23+2**22+1), 'b':4097, 'pnbits': 190,  'rnbits':190, 'cost_S': 92, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':                 2**54-2**44+1, 'b': 257, 'pnbits': 222,  'rnbits':222, 'cost_S': 98, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':             2**62-2**54+2**44, 'b':   5, 'pnbits': 254,  'rnbits':254, 'cost_S':102, 'deg_h_S':6, 'label':"Nogami et al PAIRING'2008, TEPLA"},
    {'u':              -(2**62+2**55+1), 'b':   2, 'pnbits': 254,  'rnbits':254, 'cost_S':102, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':            0x44e992b44a6909f1, 'b':   3, 'pnbits': 254,  'rnbits':254, 'cost_S':103, 'deg_h_S':6, 'label':"Ethereum BN254"},
    {'u':2**62+2**59+2**55+2**15+2**10-1,'b':   5, 'pnbits': 254,  'rnbits':254, 'cost_S':103, 'deg_h_S':6, 'label':"CBMNPZ15 eprint 2015/247"},
    {'u':                    1868033**3, 'b':   3, 'pnbits': 256,  'rnbits':256, 'cost_S':103, 'deg_h_S':6, 'label':"Naehrig-Niederhagen-Schwabe LATINCRYPT 2010 e2010/186"},
    {'u':        -(2**70+2**58+2**38+1), 'b':   2, 'pnbits': 286,  'rnbits':286, 'cost_S':109, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':               2**78+2**62+2+1, 'b':   2, 'pnbits': 318,  'rnbits':318, 'cost_S':114, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':        -(2**86-2**69+2**28+1), 'b':   2, 'pnbits': 350,  'rnbits':350, 'cost_S':119, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':        -(2**94+2**76+2**72+1), 'b':   2, 'pnbits': 382,  'rnbits':382, 'cost_S':123, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':    0x49e69d16fdc80216226909f1, 'b':  -4, 'pnbits': 383,  'rnbits':383, 'cost_S':123, 'deg_h_S':6, 'label':"El Housni Feb 2022"},
    {'u':       -(2**102+2**84-2**55+1), 'b':   2, 'pnbits': 414,  'rnbits':414, 'cost_S':129, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':                2**110+2**36+1, 'b': 257, 'pnbits': 446,  'rnbits':446, 'cost_S':132, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':-0x4000000000001000008780000000,'b':  57, 'pnbits': 446,  'rnbits':446, 'cost_S':132, 'deg_h_S':6, 'label':"Pluto curve https://github.com/daira/pluto-eris/"},
    {'u':         2**114+2**101-2**14-1, 'b':  -4, 'pnbits': 462,  'rnbits':462, 'cost_S':135, 'deg_h_S':6, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':       -(2**118-2**55-2**19+1), 'b':   2, 'pnbits': 478,  'rnbits':478, 'cost_S':138, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':       -(2**126+2**53-2**50+1), 'b': 257, 'pnbits': 510,  'rnbits':510, 'cost_S':142, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':      -(2**134+2**114+2**30+1), 'b':   2, 'pnbits': 542,  'rnbits':542, 'cost_S':146, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':      -(2**142+2**120-2**99+1), 'b':   2, 'pnbits': 574,  'rnbits':574, 'cost_S':150, 'deg_h_S':6, 'label':"Pereira et al eprint 2010/429"},
    {'u':        -(2**150-2**95+2**8+1), 'b':   2, 'pnbits': 606,  'rnbits':606, 'cost_S':153, 'deg_h_S':4, 'label':"Pereira et al eprint 2010/429"},
    {'u':         2**158-2**128-2**68+1, 'b': 257, 'pnbits': 638,  'rnbits':638, 'cost_S':157, 'deg_h_S':4, 'label':"Pereira et al eprint 2010/429"},
    {'u':         2**158-2**128-2**68+1, 'b':   5, 'pnbits': 638,  'rnbits':638, 'cost_S':157, 'deg_h_S':4, 'label':"Pereira et al eprint 2010/429"},
    # Hw 3
    {'u':         -2**254+2**33+2**6   , 'b':  10, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^33+2^6"},
    {'u':         -2**254+2**178-2**101, 'b':   7, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^178-2^101"},
    {'u':         -2**254+2**184-2**96 , 'b':  15, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^184-2^96"},
    {'u':         +2**254-2**231+2**193, 'b':  10, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^231+2^193"},
    {'u':         +2**254-2**242-2**220, 'b':  26, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^242-2^220"},
    # Hw 4
    {'u':    +2**254-2**248+2**72+2**44, 'b':  26, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^248+2^72+2^44"},
    {'u':   -2**254+2**248-2**101+2**78, 'b':  19, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^248-2^101+2^78"},
    {'u':   -2**254+2**248+2**128-2**59, 'b':  13, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^248+2^128-2^59"},
    {'u':  -2**254+2**248+2**145-2**133, 'b':   5, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^248+2^145-2^133"},
    {'u':   -2**254+2**248-2**190+2**57, 'b':   7, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^248-2^190+2^57"},
    {'u':  -2**254+2**248-2**205+2**132, 'b':  13, 'pnbits':1022, 'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^248-2^205+2^132"},    
    # Hw 5
    {'u':+2**254-2**249+2**246+2**73-2**16,  'b':10,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^246+2^73-2^16"},
    {'u':+2**254-2**249+2**246+2**136+2**87, 'b':19,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^246+2^136+2^87"},
    {'u':+2**254-2**249+2**246-2**168-2**5,  'b': 5,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^246-2^168-2^5"},
    {'u':-2**254+2**249-2**246-2**179+2**42, 'b':13,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^249-2^246-2^179+2^42"},
    {'u':+2**254-2**249+2**246-2**180-2**5,  'b': 5,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^246-2^180-2^5"},
    {'u':+2**254-2**249+2**246-2**185+2**63, 'b': 5,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^246-2^185+2^63"},
    {'u':+2**254-2**249+2**246-2**196+2**85, 'b':10,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^246-2^196+2^85"},
    {'u':+2**254-2**249+2**247-2**57+2**17,  'b': 7,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^247-2^57+2^17"},
    {'u':+2**254-2**249+2**247+2**59-2**6,   'b':14,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^247+2^59-2^6"},
    {'u':-2**254+2**249-2**247-2**77-2**44,  'b':10,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^249-2^247-2^77-2^44"},
    {'u':-2**254+2**249-2**247+2**137+2**111,'b':17,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^249-2^247+2^137+2^111"},
    {'u':+2**254-2**249+2**247+2**145-2**102,'b':10,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^247+2^145-2^102"},
    {'u':-2**254+2**249-2**247-2**171-2**133,'b': 7,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^249-2^247-2^171-2^133"},
    {'u':-2**254+2**249-2**247+2**180+2**52, 'b':10,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, -2^254+2^249-2^247+2^180+2^52"},
    {'u':+2**254-2**249+2**247+2**181+2**65, 'b':22,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^247+2^181+2^65"},
    {'u':+2**254-2**249+2**247+2**207+2**69, 'b':23,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^247+2^207+2^69"},
    {'u':+2**254-2**249+2**247+2**216-2**115,'b': 5,'pnbits':1022,'rnbits':1022, 'cost_S':191, 'deg_h_S':4, 'label':"July 23, 2019, +2^254-2^249+2^247+2^216-2^115"},
]

test_vector_bls12 = [
    {'u':   2**63+2**58+2**56+2**51+2**47+2**46+1, 'b': 1,'pnbits': 377, 'rnbits':253, 'cost_S':126, 'deg_h_S':6, 'label':"Zexe eprint 2018/962 Table 16"},
    {'u':  -(2**63+2**62+2**60+2**57+2**48+2**16), 'b': 4,'pnbits': 381, 'rnbits':255, 'cost_S':126, 'deg_h_S':6, 'label':"Zcash https://blog.z.cash/new-snark-curve/"},
    #
    {'u':-2**63+2**54+2**51-2**44+1,               'b': 1, 'pnbits':377, 'rnbits':252, 'cost_S':126, 'deg_h_S':6, 'label':"Sept 29, 2020, -2^63+2^54+2^51-2^44+1"},
    {'u': 2**63+2**61-2**58-2**56+2**50+1,         'b': 1, 'pnbits':379, 'rnbits':254, 'cost_S':126, 'deg_h_S':6, 'label':"Oct 5, 2020, 2^63+2^61-2^58-2^56+2^50+1"},
    {'u':-2**64+2**50+2**46-2**42+1,               'b': 1, 'pnbits':383, 'rnbits':256, 'cost_S':126, 'deg_h_S':6, 'label':"Sept 29, 2020, -2^64+2^50+2^46-2^42+1"},
    {'u':-2**64+2**51+2**46-2**41+1,               'b': 1, 'pnbits':383, 'rnbits':256, 'cost_S':126, 'deg_h_S':6, 'label':"Sept 29, 2020, -2^64+2^51+2^46-2^41+1"},
    {'u':-2**64+2**54-2**50+2**46+1,               'b': 1, 'pnbits':383, 'rnbits':256, 'cost_S':126, 'deg_h_S':6, 'label':"Sept 29, 2020, -2^64+2^54-2^50+2^46+1"},
    {'u': 2**64+2**59-2**57-2**55+2**53+2**51+1,   'b': 1, 'pnbits':383, 'rnbits':257, 'cost_S':126, 'deg_h_S':6, 'label':"Oct 5, 2020, 2^64+2^59-2^57-2^55+2^53+2^51+1"},
    #
    {'u':     2**64+2**16-2**13+2**11+2**9+2**5+1, 'b': 1,'pnbits': 383, 'rnbits':257, 'cost_S':126, 'deg_h_S':6, 'label':"eprint 2018/969 p.26 and http://www.martindale.info/research/TNFS-secure-pairing-implementations/"},
    {'u':          2**64+2**18+2**12+2**10+2**3+1, 'b': 1,'pnbits': 383, 'rnbits':257, 'cost_S':126, 'deg_h_S':6, 'label':"eprint 2018/969 p.26"},
    {'u':            2**64+2**51+2**24+2**12+2**9, 'b': 4,'pnbits': 383, 'rnbits':257, 'cost_S':126, 'deg_h_S':6, 'label':"eprint 2019/077 Miracl"},
    {'u':              -(2**73+2**72+2**50+2**24), 'b': 9,'pnbits': 440, 'rnbits':295, 'cost_S':132, 'deg_h_S':6, 'label':"Zhaohui-Cheng in Barbulescu-Duquesne JoC 18"},
    {'u':              -(2**12-2**48-2**71+2**74), 'b': 7,'pnbits': 442, 'rnbits':296, 'cost_S':132, 'deg_h_S':6, 'label':"Zhaohui-Cheng in Barbulescu-Duquesne JoC 18"},
    #
    {'u': 2**74+2**62-2**53+1,                     'b': 1, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Sept 29, 2020, 2^74+2^62-2^53+1"},
    {'u':-2**74-2**54-2**49-2**41+1,               'b': 1, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Sept 29, 2020, -2^74-2^54-2^49-2^41+1"},
    {'u': 2**74+2**61-2**57-2**40+1,               'b': 1, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Sept 29, 2020, +2^74+2^61-2^57-2^40+1"},
    {'u':-2**74-2**63-2**60-2**43+1,               'b': 1, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Sept 29, 2020, -2^74-2^63-2^60-2^43+1"},
    {'u': 2**74+2**67+2**49-2**41+1,               'b': 1, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Sept 29, 2020, +2^74+2^67+2^49-2^41+1"},
    #
    {'u':               2**74-2**66+2**50+2**24-1, 'b': 1, 'pnbits':443, 'rnbits':296, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74-2^66+2^50+2^24-1"},
    {'u':               2**74-2**39+2**23-2**19-1, 'b': 1, 'pnbits':443, 'rnbits':296, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74-2^39+2^23-2^19-1"},
    {'u':            2**74+2**40-2**36+2**18+2**9, 'b':-2, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^40-2^36+2^18+2^9"},
    {'u':           2**74+2**29+2**21+2**18+2**16, 'b': 4, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^29+2^21+2^18+2^16"},
    {'u':            2**74+2**44-2**38+2**28-2**4, 'b': 5, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^44-2^38+2^28-2^4"},
    {'u':            2**74+2**39-2**31-2**11-2**2, 'b': 2, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^39-2^31-2^11-2^2"},
    {'u':            2**74+2**45+2**39+2**34-2**5, 'b':-3, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^45+2^39+2^34-2^5"},
    {'u':                2**74-2**46+2**31-2**5-2, 'b': 3, 'pnbits':443, 'rnbits':296, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74-2^46+2^31-2^5-2"},
    {'u':            2**74+2**48+2**35-2**23+2**3, 'b':-3, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^48+2^35-2^23+2^3"},
    {'u':           2**74+2**50+2**27+2**18+2**15, 'b': 4, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^50+2^27+2^18+2^15"},
    {'u':           2**74+2**56+2**38-2**36+2**33, 'b':-2, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^56+2^38-2^36+2^33"},
    {'u':           2**74+2**62+2**43-2**31-2**26, 'b':-3, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^62+2^43-2^31-2^26"},
    {'u':            2**74+2**62-2**60+2**28+2**7, 'b':-3, 'pnbits':443, 'rnbits':297, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74+2^62-2^60+2^28+2^7"},
    {'u':            2**74-2**70-2**63-2**15-2**4, 'b':-2, 'pnbits':442, 'rnbits':296, 'cost_S':132, 'deg_h_S':6, 'label':"Oct 18, 2019, 2^74-2^70-2^63-2^15-2^4"},
    {'u':-(2**74+2**73+2**63+2**57+2**50+2**17+1), 'b': 1,'pnbits': 446, 'rnbits':299, 'cost_S':132, 'deg_h_S':6, 'label':"NancyMay24"},
    {'u':           2**76 + 2**53 + 2**31 + 2**11, 'b':10,'pnbits': 455, 'rnbits':305, 'cost_S':133, 'deg_h_S':6, 'label':"relic/src/epx/relic_ep2_curve.c"},
    {'u':                   -2**77 - 2**59 + 2**9, 'b': 4,'pnbits': 461, 'rnbits':309, 'cost_S':135, 'deg_h_S':6, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':                  -2**77 + 2**50 + 2**33, 'b': 4,'pnbits': 461, 'rnbits':308, 'cost_S':135, 'deg_h_S':6, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':-2**77-2**71-2**64+2**37+2**35+2**22-2**5,'b':-2,'pnbits': 461, 'rnbits':308, 'cost_S':135, 'deg_h_S':6, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':                    2**106-2**72+2**69-1, 'b': 1,'pnbits': 635, 'rnbits':424, 'cost_S':152, 'deg_h_S':6, 'label':"BosCostelloNaehrig SAC2013 eprint 2013/458"},
    {'u':    -2**106-2**92-2**60-2**34+2**12-2**9, 'b':-2,'pnbits': 635, 'rnbits':425, 'cost_S':152, 'deg_h_S':6, 'label':"CBMNPZ15 eprint 2015/247"},
    {'u':               -2**107+2**105+2**93+2**5, 'b': 4,'pnbits': 638, 'rnbits':427, 'cost_S':152, 'deg_h_S':6, 'label':"Aranha et al Pairing 2012 eprint 2012/232"},
    {'u':    -2**192+2**188-2**115-2**110-2**44-1, 'b': 1,'pnbits':1150, 'rnbits':768, 'cost_S':193, 'deg_h_S':6, 'label':"July 23, 2019, -2^192+2^188-2^115-2^110-2^44-1"},
]

# Fotiados-Martindale curves (family of Fotiadis-Konstantinou)
test_vector_fm12_17 = [
    {'u':0x41d6919153c8023f, 'b': 5, 'pnbits':384,'rnbits':254, 'deg_h_S':6,'cost_S':129, 'label':"FM12_p384"},#
    {'u':0x49e69d1640cc2ee0, 'b':-2, 'pnbits':384,'rnbits':254, 'deg_h_S':6,'cost_S':129, 'label':"FM12_p384"},#
    {'u':0x57e2266168ce8f7e, 'b': 2, 'pnbits':386,'rnbits':256, 'deg_h_S':6,'cost_S':129, 'label':"FM12_r256_p386"},#
    {'u':0x6882f5c030b0e454, 'b':-2, 'pnbits':387,'rnbits':256, 'deg_h_S':6,'cost_S':129, 'label':"FM12_r256_p387"},#
    {'u':-2**64-2**63-2**11-2**10, 'b':-2, 'pnbits': 399, 'rnbits':264, 'cost_S':130, 'deg_h_S':6, 'label':"eprint 2019/555"},
    {'u':-2**72-2**71-2**36, 'b':-2, 'pnbits': 447, 'rnbits':296, 'cost_S':135, 'deg_h_S':6, 'label':"eprint 2019/555"},
]

test_vector_kss16 = [
    {'u': 0x38fab7583, 'a':-2, 'pnbits':329, 'rnbits':255, 'cost_S':138, 'deg_h_S':16, 'label':"El Housni Feb 2022"},
    {'u':-0x3eff7803f, 'a': 1, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"-2^34+2^28+2^19+2^15-2^6+1 Hw2naf 6"},
    {'u': 0x3eff008ff, 'a':10, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"+2^34-2^28-2^20+2^11+2^8-1 Hw2naf 6"},
    {'u': 0x3ef9ffff9, 'a': 1, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"+2^34-2^28-2^23+2^21-2^3+1 Hw2naf 6"},
    {'u': 0x3ebff0011, 'a': 1, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"+2^34-2^28-2^26-2^16+2^4+1 Hw2naf 6"},
    {'u':-0x3dffb7001, 'a':14, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"-2^34+2^29+2^18+2^15+2^12-1 Hw2naf 6"},
    {'u':-0x3e8080011, 'a': 5, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"-2^34+2^29-2^27-2^19-2^4-1 Hw2naf 6"},
    {'u': 0x3d81fe001, 'a': 1, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"+2^34-2^29-2^27+2^21-2^13+1 Hw2naf 6"},
    {'u':-0x3c01d0001, 'a': 5, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"-2^34+2^30-2^21+2^18-2^16-1 Hw2naf 6"},
    {'u':-0x3c7804001, 'a':13, 'pnbits':330, 'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"-2^34+2^30-2^27+2^23-2^14-1 Hw2naf 6"},
    #
    {'u':-2**34+2**30-2**20-2**17-2**14+2**12+1,        'a': 1,'pnbits':330,'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"July 22, 2019, -2^34+2^30-2^20-2^17-2^14+2^12+1"},
    {'u':-2**34+2**30+2**23-2**20+2**18+2**7+1,         'a': 1,'pnbits':330,'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"July 22, 2019, -2^34+2^30+2^23-2^20+2^18+2^7+1"},
    {'u':+2**34-2**30+2**26+2**23+2**14-2**5+1,         'a': 1,'pnbits':330,'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"July 22, 2019, +2^34-2^30+2^26+2^23+2^14-2^5+1"},
    {'u':-2**34+2**30-2**26+2**24-2**12+2**3+1,         'a': 1,'pnbits':330,'rnbits':256, 'cost_S':141, 'deg_h_S':16, 'label':"July 22, 2019, -2^34+2^30-2^26+2^24-2^12+2^3+1"},
    #
    {'u':-2**34+2**27-2**23+2**20-2**11+1,              'a': 1,'pnbits':330,'rnbits':257, 'cost_S':141, 'deg_h_S':16, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':        2**35-2**32-2**18+2**8+1,              'a': 1,'pnbits':339,'rnbits':263, 'cost_S':141, 'deg_h_S':16, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':        2**49+2**26+2**15-2**7-1,              'a': 3,'pnbits':481,'rnbits':377, 'cost_S':163, 'deg_h_S':16, 'label':"Zhang Lin INDOCRYPT 2012, a=3,-3"},
    {'u':2**50+2**47-2**38+2**32+2**25-2**15-2**5-1,    'a':17,'pnbits':492,'rnbits':386, 'cost_S':163, 'deg_h_S':16, 'label':"Ghammam-Fouotsa eprint 2016/472"},
    #
    {'u':-2**77-2**33-2**28-2**10-1, 'a': 3, 'pnbits':761, 'rnbits':601, 'cost_S':194, 'deg_h_S':16, 'label':"-2^77-2^33-2^28-2^10-1"},
    {'u':-2**77+2**57-2**35+2**25+1, 'a': 1, 'pnbits':761, 'rnbits':601, 'cost_S':194, 'deg_h_S':16, 'label':"-2^77+2^57-2^35+2^25+1"},
    {'u':-2**77+2**63-2**36-2**24+1, 'a': 1, 'pnbits':761, 'rnbits':601, 'cost_S':194, 'deg_h_S':16, 'label':"-2^77+2^63-2^36-2^24+1"},
    {'u': 2**77-2**63-2**60+2**43-1, 'a': 5, 'pnbits':761, 'rnbits':601, 'cost_S':194, 'deg_h_S':16, 'label':"+2^77-2^63-2^60+2^43-1"},
    {'u':-2**77+2**67-2**27-2**18+1, 'a': 1, 'pnbits':761, 'rnbits':601, 'cost_S':194, 'deg_h_S':16, 'label':"-2^77+2^67-2^27-2^18+1"},
    {'u':-2**77+2**71+2**56-2**31+1, 'a': 1, 'pnbits':760, 'rnbits':600, 'cost_S':194, 'deg_h_S':16, 'label':"-2^77+2^71+2^56-2^31+1"},
    {'u':-2**77-2**75-2**31-2**11+1, 'a': 1, 'pnbits':764, 'rnbits':603, 'cost_S':194, 'deg_h_S':16, 'label':"-2^77-2^75-2^31-2^11+1"},
    #
    {'u':2**78-2**76-2**28+2**14+2**7+1,                'a': 1,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, 2^78-2^76-2^28+2^14+2^7+1, a=1,-3, in paper"},
    {'u':2**78-2**76+2**34-2**29-2**4+1,                'a': 1,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, 2^78-2^76+2^34-2^29-2^4+1"},
    {'u':-(2**78-2**76-2**72+2**49-2**14+2**4-1),       'a': 1,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, -(2^78-2^76-2^72+2^49-2^14+2^4-1)"},
    {'u':-(2**78-2**76-2**72-2**23-2**13+2**5+1),       'a': 3,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, -(2^78-2^76-2^72-2^23-2^13+2^5+1)"},
    {'u':2**78-2**76-2**73+2**49+2**16+2**5-1,          'a': 3,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, 2^78-2^76-2^73+2^49+2^16+2^5-1"},
    {'u':-(2**78-2**76+2**70+2**54+2**37+2**35+2**3+1), 'a':11,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, -(2^78-2^76+2^70+2^54+2^37+2^35+2^3+1)"},
    {'u':-(2**78-2**76-2**72+2**41+2**25+2**21+2**8+1), 'a': 3,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, -(2^78-2^76-2^72+2^41+2^25+2^21+2^8+1)"},
    {'u':-(2**78-2**76-2**72+2**51-2**49+2**47+2**3+1), 'a': 3,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, -(2^78-2^76-2^72+2^51-2^49+2^47+2^3+1)"},
    {'u':-(2**78-2**76-2**72+2**67+2**28+2**20+2**5+1), 'a': 5,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, -(2^78-2^76-2^72+2^67+2^28+2^20+2^5+1)"},
    {'u':-(2**78-2**76-2**73+2**48+2**15+2**12+2**7+1), 'a': 7,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, -(2^78-2^76-2^73+2^48+2^15+2^12+2^7+1)"},
    {'u':-(2**78-2**76-2**73+2**60+2**55+2**28+2**15+1),'a':17,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, -(2^78-2^76-2^73+2^60+2^55+2^28+2^15+1)"},
    {'u':2**78-2**76-2**73+2**68+2**54+2**44+2**29+1,   'a': 1,'pnbits':766,'rnbits':605, 'cost_S':194, 'deg_h_S':16, 'label':"July 22, 2019, 2^78-2^76-2^73+2^68+2^54+2^44+2^29+1"},
]

test_vector_kss18 = [
    {'u':0xc0c44000000, 'b': 2,'pnbits': 345, 'rnbits': 256, 'cost_S':150, 'deg_h_S':18, 'label':"El Housni Feb 2022"},
    {'u':2**44 + 2**22 - 2**9 + 2, 'b': 3,'pnbits': 348, 'rnbits': 256, 'cost_S':151, 'deg_h_S':18, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':-2**44-2**36+2**33+2**17, 'b':28,'pnbits': 348, 'rnbits': 256, 'cost_S':151, 'deg_h_S':18, 'label':"-2^44-2^36+2^33+2^17"},
    {'u':2**44-2**36-2**29-2**24+2**19, 'b': 6, 'pnbits':348, 'rnbits':256, 'cost_S':151, 'deg_h_S': 18, 'label':"2^44-2^36-2^29-2^24+2^19"},
    {'u':2**44-2**37-2**24+2**19+2**17, 'b':14, 'pnbits':348, 'rnbits':256, 'cost_S':151, 'deg_h_S': 18, 'label':"2^44-2^37-2^24+2^19+2^17"},
    {'u':-2**44-2**37-2**34-2**28-2**23,'b': 2, 'pnbits':348, 'rnbits':256, 'cost_S':151, 'deg_h_S': 18, 'label':"-2^44-2^37-2^34-2^28-2^23"},
    {'u':-2**44+2**37+2**35+2**22-2**17,'b': 6, 'pnbits':348, 'rnbits':256, 'cost_S':151, 'deg_h_S': 18, 'label':"-2^44+2^37+2^35+2^22-2^17"},
    {'u':2**64-2**51+2**47+2**28,  'b': 2,'pnbits': 508, 'rnbits': 376, 'cost_S':179, 'deg_h_S':18, 'label':"BosCostelloNaehrig SAC2013 eprint 2013/458"},
    {'u':2**64+2**47+2**43+2**37+2**26+2**25+2**19-2**13-2**7,
                                   'b': 2,'pnbits': 508, 'rnbits': 376, 'cost_S':179, 'deg_h_S':18, 'label':"CBMNPZ15 eprint 2015/247"},
    {'u':-2**64-2**51+2**46+2**12, 'b': 2,'pnbits': 508, 'rnbits': 376, 'cost_S':179, 'deg_h_S':18, 'label':"Aranha et al Pairing 2012 eprint 2012/232"},
    #
    {'u':-2**80+2**37+2**34      , 'b': 7, 'pnbits':636, 'rnbits': 472, 'cost_S':193, 'deg_h_S': 9, 'label':"-2^80+2^37+2^34"},
    {'u':-2**80-2**26-2**24+2**7 , 'b': 2, 'pnbits':636, 'rnbits': 472, 'cost_S':193, 'deg_h_S': 9, 'label':"-2^80-2^26-2^24+2^7"},
    {'u': 2**80-2**33+2**30-2**2 , 'b':17, 'pnbits':636, 'rnbits': 472, 'cost_S':193, 'deg_h_S': 9, 'label':"2^80-2^33+2^30-2^2"},
    {'u': 2**80-2**43+2**31-2**5 , 'b':14, 'pnbits':636, 'rnbits': 472, 'cost_S':193, 'deg_h_S': 9, 'label':"2^80-2^43+2^31-2^5"},
    {'u': 2**80-2**49-2**18-2**3 , 'b': 2, 'pnbits':636, 'rnbits': 472, 'cost_S':193, 'deg_h_S': 9, 'label':"2^80-2^49-2^18-2^3"},
    {'u': 2**80+2**51-2**27-2**17, 'b':14, 'pnbits':636, 'rnbits': 472, 'cost_S':193, 'deg_h_S': 9, 'label':"2^80+2^51-2^27-2^17"},
    {'u': 2**80+2**56-2**40+2**24, 'b': 2, 'pnbits':636, 'rnbits': 472, 'cost_S':193, 'deg_h_S': 9, 'label':"2^80+2^56-2^40+2^24"},
    {'u':-2**80-2**73-2**13+2**6 , 'b': 2, 'pnbits':636, 'rnbits': 472, 'cost_S':193, 'deg_h_S': 9, 'label':"-2^80-2^73-2^13+2^6"},
    #
    {'u':0x12fffdfdfffffffffc000,  'b': 6,'pnbits': 638, 'rnbits': 474, 'cost_S':193, 'deg_h_S': 9, 'label':"July21,2019, 2^80+2^77+2^76-2^61-2^53-2^14"},
    {'u':0x12ffffffffe9ff0000000,  'b': 7,'pnbits': 638, 'rnbits': 474, 'cost_S':193, 'deg_h_S': 9, 'label':"July21,2019, 2^80+2^77+2^76-2^40-2^38-2^37-2^28"},
    {'u':-0x130000000040100020080, 'b': 6,'pnbits': 638, 'rnbits': 474, 'cost_S':193, 'deg_h_S': 9, 'label':"July21,2019, -(2^80+2^77+2^76+2^42+2^32+2^17+2^7)"},
    {'u':-2**85-2**31-2**26+2**6,  'b': 2,'pnbits': 676, 'rnbits': 502, 'cost_S':196, 'deg_h_S': 9, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':2**186-2**75-2**22+2**4,  'b':-2,'pnbits':1484, 'rnbits':1108, 'cost_S':256, 'deg_h_S': 6, 'label':"Barbulescu-Duquesne JoC 18"}]

test_vector_bls24 = [
    # with 2^n | (u-1) and n is maximal
    {'u': 2**31-2**29+2**22-2**20+1            ,'b':1,'pnbits':305, 'rnbits':245, 'cost_S':158, 'deg_h_S':24, 'label':"Sept 29, 2020,  2^31-2^29+2^22-2^20+1"},
    {'u':-2**31-2**28-2**26-2**24-2**20+1      ,'b':1,'pnbits':311, 'rnbits':250, 'cost_S':159, 'deg_h_S':24, 'label':"Sept 29, 2020, -2^31-2^28-2^26-2^24-2^20+1"},
    {'u': 2**31+2**29-2**23+2**21-2**18+1      ,'b':1,'pnbits':312, 'rnbits':251, 'cost_S':159, 'deg_h_S':24, 'label':"Sept 29, 2020, 2^31+2^29-2^23+2^21-2^18+1"},
    {'u':-2**32+2**30+2**22-2**20+1            ,'b':1,'pnbits':315, 'rnbits':253, 'cost_S':160, 'deg_h_S':24, 'label':"Sept 29, 2020, -2^32+2^30+2^22-2^20+1"},
    {'u':-2**32+2**30-2**27-2**24-2**20+2**18+1,'b':1,'pnbits':315, 'rnbits':254, 'cost_S':160, 'deg_h_S':24, 'label':"Sept 29, 2020, -2^32+2^30-2^27-2^24-2^20+2^18+1"},

    {'u': 2**32-2**29-2**27+2**24+2**16+2**15  ,'b':4,'pnbits':317, 'rnbits':255, 'cost_S':160, 'deg_h_S':24, 'label':"El Housni Feb 2022, 2^32-2^29-2^27+2^24+2^16+2^15"},
    {'u': 2**32-2**29+2**25-2**23+2**21-2**18+1,'b':1,'pnbits':317, 'rnbits':255, 'cost_S':160, 'deg_h_S':24, 'label':"Sept 29, 2020, 2^32-2^29+2^25-2^23+2^21-2^18+1"},
    {'u':-2**32-2**26-2**23-2**19+1            ,'b':1,'pnbits':319, 'rnbits':257, 'cost_S':161, 'deg_h_S':24, 'label':"Sept 29, 2020, -2^32-2^26-2^23-2^19+1"},
    #
    {'u': 2**32-2**28+2**16-2**10+2**5+2**2+1  ,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, 2^32-2^28+2^16-2^10+2^5+2^2+1"},
    {'u': 2**32-2**28-2**21-2**18+2**16+2**4-1 ,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, 2^32-2^28-2^21-2^18+2^16+2^4-1"},
    {'u':-2**32+2**28+2**21+2**19-2**16+2**9-1 ,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, -2^32+2^28+2^21+2^19-2^16+2^9-1"},
    {'u':-2**32+2**28+2**22+2**20+2**17-2**11-1,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, -2^32+2^28+2^22+2^20+2^17-2^11-1"},
    {'u':-2**32+2**28-2**23+2**21+2**14+2**12-1,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, -2^32+2^28-2^23+2^21+2^14+2^12-1"},
    {'u':-2**32+2**28-2**23+2**21+2**18+2**12-1,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, in paper, -2^32+2^28-2^23+2^21+2^18+2^12-1"},
    {'u':-2**32+2**28+2**24-2**13-2**11+2**7-1 ,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, -2^32+2^28+2^24-2^13-2^11+2^7-1"},
    {'u':-2**32+2**28+2**24-2**18+2**9-2**3+1  ,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, -2^32+2^28+2^24-2^18+2^9-2^3+1"},
    {'u': 2**32-2**28+2**24+2**19+2**13-2**3+1 ,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, 2^32-2^28+2^24+2^19+2^13-2^3+1"},
    {'u': 2**32-2**28-2**25+2**14+2**10-2**4-1 ,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, 2^32-2^28-2^25+2^14+2^10-2^4-1"},
    {'u': 2**32-2**28+2**26-2**19-2**8-2**6+1  ,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, 2^32-2^28+2^26-2^19-2^8-2^6+1"},
    {'u': 2**32-2**28+2**26+2**21-2**15-2**13-1,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, 2^32-2^28+2^26+2^21-2^15-2^13-1"},
    {'u': 2**32-2**28+2**26-2**21+2**16-2**10-1,'b':1,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"July 23, 2019, 2^32-2^28+2^26-2^21+2^16-2^10-1"},#12
    #
    {'u': -2**32+2**28+2**12,'b':5,'pnbits':318,'rnbits':256,'cost_S':161,'deg_h_S':24,'label':"-2^32+2^28+2^12 Barbulescu-El Mrabet-Ghammam 2019/485 Tab.8"},
    #
    {'u':-2**48+2**45+2**31-2**7,    'b': 4, 'pnbits': 477, 'rnbits': 383, 'cost_S':188.6, 'deg_h_S':24, 'label':"Aranha et al Pairing 2012 eprint 2012/232"},
    {'u':-1+2**11-2**28-2**51,       'b': 1, 'pnbits': 509, 'rnbits': 409, 'cost_S':193.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-2**16+2**46-2**48+2**51,   'b':-3, 'pnbits': 507, 'rnbits': 407, 'cost_S':193.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':2**15-2**22+2**47-2**51,    'b': 4, 'pnbits': 508, 'rnbits': 408, 'cost_S':193.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-2**4-2**6-2**8-2**51,      'b': 4, 'pnbits': 509, 'rnbits': 409, 'cost_S':193.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':2**39+2**45-2**48-2**51,    'b': 4, 'pnbits': 510, 'rnbits': 410, 'cost_S':193.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**8-2**15+2**51,        'b': 1, 'pnbits': 509, 'rnbits': 408, 'cost_S':193.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-2**56-2**43+2**9-2**6,     'b':-2, 'pnbits': 559, 'rnbits': 449, 'cost_S':201.0, 'deg_h_S':24, 'label':"Barbulescu-Duquesne JoC 18"},
    {'u':2**63-2**47+2**38,          'b': 4, 'pnbits': 629, 'rnbits': 504, 'cost_S':212.0, 'deg_h_S':24, 'label':"BosCostelloNaehrig SAC2013 eprint 2013/458"},
    {'u':-(2**63-2**47-2**31-2**26-2**24+2**8-2**5+1),
                                     'b': 1, 'pnbits': 629, 'rnbits': 504, 'cost_S':212.0, 'deg_h_S':24, 'label':"CBMNPZ15 eprint 2015/247"},
    {'u':-2**103-2**101+2**68+2**50, 'b':-2, 'pnbits':1032, 'rnbits': 827, 'cost_S':265.0, 'deg_h_S':24, 'label':"Barbulescu-Duquesne JoC 18"},
    #
    #u = 1 mod 3 = 4 mod 6 = 40 mod 72 u =  0x4d000000000001000000008200 # Hw_2NAF: 7 Hw: 7 c = 3*(1*2033522837866118092059286973269)^2 c1 is prime: True
    {'u':2**102+2**100-2**98+2**96+2**48+2**15+2**9,  'b':22,'pnbits':1022,'rnbits':819,'cost_S':264,'deg_h_S':24,'label':"July 24, 2019, 2^102+2^100-2^98+2^96+2^48+2^15+2^9"},
    #u = 1 mod 3 = 4 mod 6 = 16 mod 72 u =  0x4cbfffffffffffffbfffffffe0 # Hw_2NAF: 7 Hw: 7 c = 3*(1*2026920490989929303343206760437)^2 c1 is prime: True
    {'u':2**102+2**100-2**98+2**96-2**94-2**38-2**5,  'b': 4,'pnbits':1022,'rnbits':819,'cost_S':264,'deg_h_S':24,'label':"July 24, 2019, 2^102+2^100-2^98+2^96-2^94-2^38-2^5"},
    #u = 1 mod 3 = 4 mod 6 = 52 mod 72 u =  0x4ffffffffffffffffffffffc84 # Hw_2NAF: 5 Hw: 5 tw2 is prime: True
    {'u':2**102+2**100-2**10+2**7+2**2,               'b': 4,'pnbits':1022,'rnbits':819,'cost_S':264,'deg_h_S':24,'label':"July 24, 2019, 2^102+2^100-2^10+2^7+2^2"},
    #u = 1 mod 3 = 4 mod 6 = 64 mod 72 u = -0x4cffffffffdffffffffdf80000 # Hw_2NAF: 7 Hw: 7 tw2 is prime: True
    {'u':-2**102-2**100+2**98-2**96+2**61+2**25+2**19,'b':-2,'pnbits':1022,'rnbits':819,'cost_S':264,'deg_h_S':24,'label':"July 24, 2019, -2^102-2^100+2^98-2^96+2^61+2^25+2^19"},
    #u = 1 mod 3 = 4 mod 6 = 16 mod 72 u = -0x4cffffff00000007c000000000 # Hw_2NAF: 7 Hw: 7 tw2 is prime: True
    {'u':-2**102-2**100+2**98-2**96+2**72-2**43+2**38,'b': 4,'pnbits':1022,'rnbits':819,'cost_S':264,'deg_h_S':24,'label':"July 24, 2019, -2^102-2^100+2^98-2^96+2^72-2^43+2^38"},
]

# for u= 7 mod 72, b=1
# for u=16 mod 72, b=4
# for u=31 mod 72, b=1
# for u=64 mod 72, b=-2
test_vector_bls24_cln = [ #https://eprint.iacr.org/2011/465
    {'u':-1-2**8+2**38+2**45,        'b':1, 'pnbits':449, 'rnbits':361, 'cost_S':186.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**3-2**5-2**19+2**46,   'b':1, 'pnbits':459, 'rnbits':368, 'cost_S':187.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**11-2**26-2**35-2**47, 'b':1, 'pnbits':469, 'rnbits':377, 'cost_S':188.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**19-2**24+2**27-2**48, 'b':1, 'pnbits':479, 'rnbits':384, 'cost_S':189.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**11-2**28-2**35-2**49, 'b':1, 'pnbits':489, 'rnbits':393, 'cost_S':190.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**4-2**21-2**50,        'b':1, 'pnbits':499, 'rnbits':401, 'cost_S':191.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**11-2**28-2**51,       'b':1, 'pnbits':509, 'rnbits':409, 'cost_S':192.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**22-2**26-2**36-2**52, 'b':1, 'pnbits':519, 'rnbits':417, 'cost_S':193.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**44+2**51+2**53,       'b':1, 'pnbits':532, 'rnbits':427, 'cost_S':194.0, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    #
    {'u':-1-2**3-2**29-2**38-2**55,  'b':1, 'pnbits':549, 'rnbits':441, 'cost_S':197.5, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**3-2**11-2**51+2**56,  'b':1, 'pnbits':558, 'rnbits':448, 'cost_S':199.2, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**15-2**22-2**56,       'b':1, 'pnbits':559, 'rnbits':449, 'cost_S':199.2, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**16+2**23-2**28+2**57, 'b':1, 'pnbits':569, 'rnbits':456, 'cost_S':199.2, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**28+2**51+2**58,       'b':1, 'pnbits':579, 'rnbits':465, 'cost_S':201.5, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**12-2**28-2**50-2**58, 'b':1, 'pnbits':579, 'rnbits':465, 'cost_S':201.5, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    #
    # DL cost is not precise for the data below
    {'u':-1+2**15+2**43-2**61,       'b':1, 'pnbits':609, 'rnbits':488, 'cost_S':207, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**49+2**55+2**62,       'b':1, 'pnbits':619, 'rnbits':497, 'cost_S':208, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**19-2**23-2**26-2**63, 'b':1, 'pnbits':629, 'rnbits':505, 'cost_S':209, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**8-2**35-2**61-2**63,  'b':1, 'pnbits':632, 'rnbits':507, 'cost_S':209, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**17-2**54+2**61-2**64, 'b':1, 'pnbits':637, 'rnbits':511, 'cost_S':210, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**35+2**60-2**64,       'b':1, 'pnbits':638, 'rnbits':512, 'cost_S':210, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**10+2**14-2**18+2**64, 'b':1, 'pnbits':639, 'rnbits':512, 'cost_S':210, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**3-2**37-2**50-2**65,  'b':1, 'pnbits':649, 'rnbits':521, 'cost_S':211, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    #
    {'u':-1-2**52+2**59-2**81,       'b':1, 'pnbits':809, 'rnbits':648, 'cost_S':236, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**26-2**74+2**82,       'b':1, 'pnbits':819, 'rnbits':656, 'cost_S':238, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**21-2**73-2**82,       'b':1, 'pnbits':819, 'rnbits':657, 'cost_S':238, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**7-2**13-2**27-2**82,  'b':1, 'pnbits':819, 'rnbits':657, 'cost_S':238, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**11-2**23-2**32-2**83, 'b':1, 'pnbits':829, 'rnbits':665, 'cost_S':240, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**48-2**52-2**72-2**84, 'b':1, 'pnbits':839, 'rnbits':673, 'cost_S':241, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**8-2**12+2**16-2**85,  'b':1, 'pnbits':849, 'rnbits':680, 'cost_S':242, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**3+2**31-2**86,        'b':1, 'pnbits':859, 'rnbits':688, 'cost_S':243, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**16-2**20-2**87,       'b':1, 'pnbits':869, 'rnbits':697, 'cost_S':244, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**53-2**56+2**88,       'b':1, 'pnbits':879, 'rnbits':704, 'cost_S':246, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**23+2**67+2**90,       'b':1, 'pnbits':899, 'rnbits':721, 'cost_S':250, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    #
    {'u':-1-2**12-2**93-2**95-2**107,'b':1, 'pnbits':1069, 'rnbits':857, 'cost_S':270, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**65-2**75+2**109,      'b':1, 'pnbits':1089, 'rnbits':872, 'cost_S':272, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**64-2**100+2**110,     'b':1, 'pnbits':1099, 'rnbits':880, 'cost_S':273, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1+2**15+2**93-2**111,      'b':1, 'pnbits':1109, 'rnbits':888, 'cost_S':274, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
    {'u':-1-2**13+2**57-2**112,      'b':1, 'pnbits':1119, 'rnbits':896, 'cost_S':275, 'deg_h_S':24, 'label':"CostelloLauterNaehrig Indocrypt 2011 eprint 2011/465"},
]

test_vector_bls48 = [
    {'u': 0xf5ef, 'b': 1, 'pnbits':286, 'rnbits':256, 'cost_S':230, 'deg_h_S':16, 'cost_C':222, 'deg_h_C':8, 'label':"MIRACL Mike Scott 2^16-2^11-2^9-2^4-1 Hw2naf 5"},
    {'u':-0xf5e1, 'b': 1, 'pnbits':286, 'rnbits':256, 'cost_S':230, 'deg_h_S':16, 'cost_C':222, 'deg_h_C':8, 'label':"-2^16+2^11+2^9+2^5-1 Hw2naf 5"},
    {'u': 0xf0c4, 'b': 5, 'pnbits':285, 'rnbits':255, 'cost_S':230, 'deg_h_S':16, 'cost_C':222, 'deg_h_C':8, 'label':"2^16-2^12+2^8-2^6+2^2 Hw2naf 5"},
]

test_vector_fk12d3 = [
    {'u':0x6e94899bd5a22c8d,       'b':2, 'pnbits':388, 'rnbits': 257, 'cost_S':127, 'deg_h_S':6, 'label': "Fotiadis-Konstantinou eprint 2018/1017 Table 8"},
    {'u':-2**64-2**63-2**11-2**10, 'b':7, 'pnbits':399, 'rnbits': 264, 'cost_S':130, 'deg_h_S':6, 'label': "Fotiadis-Martindale eprint 2019/555 Family 17 (a)" },
    {'u':-2**66+2**25+2**17,       'b':5, 'pnbits':407, 'rnbits': 270, 'cost_S':131, 'deg_h_S':6, 'label': "Fotiadis-Martindale eprint 2019/555 Family 17 (b)" },
]
test_vector_fm12d3 = [
    {'u':0x2977b39a701cd0b50,                    'b':2, 'pnbits':388, 'rnbits': 257, 'cost_S':127, 'deg_h_S':6, 'label': "Fotiadis-Konstantinou eprint 2018/1017 Table 8 u=6*u0+2"},
    {'u':-2**67-2**64-2**14-2**11+2,             'b':7, 'pnbits':399, 'rnbits': 264, 'cost_S':130, 'deg_h_S':6, 'label': "Fotiadis-Martindale eprint 2019/555 Family 17 (a), u=6*u0+2"},
    {'u':-2**68-2**67+2**27+2**26+2**19+2**18+2, 'b':5, 'pnbits':407, 'rnbits': 270, 'cost_S':131, 'deg_h_S':6, 'label': "Fotiadis-Martindale eprint 2019/555 Family 17 (b), u=6*u0+2"},
]



# KSS32 curves with seed u = [325, 5889] mod 6214
test_vector_kss32 = [
    {'u':-0xce3eb,  'a':-2, 'pnbits':333,'rnbits':251, 'cofact_r': 308993, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^20+2^18-2^16+2^13-2^10+2^4+2^2+1 Hw2NAF=8"},#
    {'u':-0x1e9f0d, 'a':-2, 'pnbits':356,'rnbits':276, 'cofact_r':   9409, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^21+2^17-2^15-2^13+2^8-2^4+2^2-1 Hw2NAF=8"},#
    {'u':-0x257547, 'a':-3, 'pnbits':361,'rnbits':282, 'cofact_r':   2657, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^21-2^19+2^17+2^15+2^12-2^10-2^8-2^6-2^3+1 Hw2NAF=10"},#
    {'u':-0x41e74f, 'a': 1, 'pnbits':376,'rnbits':307, 'cofact_r':      1, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^22-2^17+2^13-2^11+2^8-2^6-2^4+1 Hw2NAF=8"},#
    {'u': 0x4599f9, 'a':-3, 'pnbits':377,'rnbits':301, 'cofact_r':     97, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^22+2^19-2^17-2^15+2^13-2^11+2^9-2^3+1 Hw2NAF=9"},#
    {'u':-0x462b9d, 'a':-2, 'pnbits':377,'rnbits':308, 'cofact_r':      1, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^22-2^19+2^17-2^14+2^12+2^10+2^7-2^5+2^2-1 Hw2NAF=10"},#
    {'u':-0x4c24d7, 'a': 1, 'pnbits':380,'rnbits':310, 'cofact_r':      1, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^22-2^20+2^18-2^13-2^10-2^8+2^5+2^3+1 Hw2NAF=9"},#
    #
    {'u':0x2a1953b, 'a': 2, 'pnbits':436,'rnbits':360, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^25+2^23+2^21+2^17-2^15+2^12+2^10+2^8+2^6-2^2-1 Hw2NAF=11"},#
    {'u':0x34617a5, 'a': 2, 'pnbits':442,'rnbits':365, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^26-2^24+2^22+2^19-2^17+2^13-2^11-2^7+2^5+2^2+1 Hw2NAF=11"},#
    {'u':0x59c067b, 'a': 2, 'pnbits':456,'rnbits':378, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27-2^25-2^23+2^21-2^18+2^11-2^9+2^7-2^2-1 Hw2NAF=10"},#
    {'u':0x5e5f939, 'a': 1, 'pnbits':457,'rnbits':379, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27-2^25-2^21+2^19-2^17-2^11+2^8+2^6-2^3+1 Hw2NAF=10"},#
    {'u':0x5fa8f19, 'a': 1, 'pnbits':457,'rnbits':379, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27-2^25-2^19+2^17+2^15+2^12-2^8+2^5-2^3+1 Hw2NAF=10"},#
    {'u':0x79611d1, 'a': 1, 'pnbits':464,'rnbits':385, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27-2^23+2^21-2^19-2^17+2^12+2^9-2^6+2^4+1 Hw2NAF=10"},#
    {'u':0x7d3f9d5, 'a': 2, 'pnbits':464,'rnbits':386, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27-2^22+2^20+2^18-2^11+2^9-2^6+2^4+2^2+1 Hw2NAF=10"},#
    {'u':0x807a965, 'a': 2, 'pnbits':465,'rnbits':386, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27+2^19-2^15+2^13+2^11+2^9-2^7-2^5+2^2+1 Hw2NAF=10"},#
    {'u':0x9aaa9b7, 'a': 5, 'pnbits':470,'rnbits':390, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27+2^25-2^23+2^21+2^19+2^17+2^15+2^13+2^11+2^9-2^6-2^3-1 Hw2NAF=13"},#
    {'u':0x9b58e9f, 'a': 3, 'pnbits':470,'rnbits':391, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27+2^25-2^22-2^19-2^17-2^15+2^12-2^9+2^7+2^5-1 Hw2NAF=11"},#
    {'u':0xa822cb3, 'a':-2, 'pnbits':472,'rnbits':392, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^27+2^25+2^23+2^17+2^14-2^12-2^10+2^8-2^6-2^4+2^2-1 Hw2NAF=12"},#
    {'u':0xb6ffd05, 'a': 2, 'pnbits':474,'rnbits':394, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^28-2^26-2^23-2^20-2^10+2^8+2^2+1 Hw2NAF=8"},#
    {'u':0xc7ace79, 'a': 1, 'pnbits':477,'rnbits':396, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^28-2^26+2^23-2^18-2^16-2^14+2^12-2^9+2^7-2^3+1 Hw2NAF=11"},#
    #
    {'u':-0x7037b15, 'a': 2, 'pnbits':462, 'rnbits':383, 'cost_S':None, 'deg_h_S':None, 'label':"-2^27+2^24-2^18+2^15+2^10+2^8-2^4-2^2-1 Hw 9"},
    {'u':-0x7afc161, 'a': 3, 'pnbits':464, 'rnbits':385, 'cost_S':None, 'deg_h_S':None, 'label':"-2^27+2^22+2^20+2^14-2^9+2^7+2^5-1 Hw 8"},
    {'u':-0x994005d, 'a':-2, 'pnbits':470, 'rnbits':390, 'cost_S':None, 'deg_h_S':None, 'label':"-2^27-2^25+2^23-2^20-2^18-2^7+2^5+2^2-1 Hw 9"},
    #
    {'u':-0x676311f, 'a': 1, 'pnbits':459,'rnbits':381, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^27+2^25-2^23+2^19+2^17-2^14+2^12-2^8-2^5+1 Hw2NAF=10"},#
    {'u':-0x6c95b77, 'a': 1, 'pnbits':461,'rnbits':382, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^27+2^24+2^22-2^19-2^17+2^15+2^13+2^10+2^7+2^3+1 Hw2NAF=11"},#
    {'u':-0x6db5f6b, 'a': 2, 'pnbits':461,'rnbits':383, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^27+2^24+2^21+2^18+2^15+2^13+2^7+2^4+2^2+1 Hw2NAF=10"},#
    {'u':-0x7521339, 'a':17, 'pnbits':463,'rnbits':384, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^27+2^24-2^22-2^20-2^17-2^12-2^10+2^8-2^6+2^3-1 Hw2NAF=11"},#
    #
    {'u':-0x73c6c0203, 'a':19, 'pnbits':606,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^35+2^32-2^30+2^26-2^23+2^20+2^18-2^9-2^2+1 Hw2NAF=10"},#
    {'u':-0x73c6c0203, 'a':19, 'pnbits':606,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^35+2^32-2^30+2^26-2^23+2^20+2^18-2^9-2^2+1 Hw2NAF=10"},#
    {'u':-0x750051bdf, 'a': 1, 'pnbits':607,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^35+2^32-2^30-2^28-2^18-2^16-2^13+2^10+2^5+1 Hw2NAF=10"},#
    {'u':-0x75f11dcff, 'a': 1, 'pnbits':607,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^35+2^31+2^29+2^24-2^20-2^17+2^13+2^10-2^8+1 Hw2NAF=10"},#
    {'u': 0x7341c0141, 'a': 1, 'pnbits':606,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^35-2^32+2^30-2^28+2^26+2^21-2^18+2^8+2^6+1 Hw2NAF=10"},#
    {'u': 0x744fdcdff, 'a': 3, 'pnbits':607,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^35-2^32+2^30+2^26+2^24-2^17-2^14+2^12-2^9-1 Hw2NAF=10"},#
    {'u': 0x7341c0141, 'a': 1, 'pnbits':606,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^35-2^32+2^30-2^28+2^26+2^21-2^18+2^8+2^6+1 Hw2NAF=10"},#
    {'u': 0x744fdcdff, 'a': 3, 'pnbits':607,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^35-2^32+2^30+2^26+2^24-2^17-2^14+2^12-2^9-1 Hw2NAF=10"},#
    {'u': 0x74fa403df, 'a': 3, 'pnbits':607,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^35-2^32+2^30+2^28-2^23+2^21+2^18+2^10-2^5-1 Hw2NAF=10"},#
    {'u': 0x757eb7c05, 'a': 2, 'pnbits':607,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^35-2^31-2^29-2^27-2^20-2^18-2^15-2^10+2^2+1 Hw2NAF=10"},#
]
# KSS36 curves with seed u = [287, 308, 497, 539, 728, 749] mod 777
test_vector_kss36 = [
    {'u':-0xbbcdf8, 'b': 3, 'pnbits':315,'rnbits':256, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^24+2^22+2^18+2^14-2^12+2^9+2^3  Hw2NAF=7"},#
    {'u':-0xbdf412, 'b': 5, 'pnbits':316,'rnbits':256, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^24+2^22+2^17+2^12-2^10-2^4-2    Hw2NAF=7"},#
    {'u': 0xbfef0d, 'b': 2, 'pnbits':316,'rnbits':256, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^24-2^22-2^12-2^8+2^4-2^2+1      Hw2NAF=7"},#
    #
    {'u':0x49fbfd004, 'b':-4, 'pnbits':465, 'rnbits':384, 'cost_S':None, 'deg_h_S':None, 'label':"+2^34+2^31+2^29-2^22-2^14+2^12+2^2 Hw 7"},
    {'u':0x4beffbc10, 'b':-7, 'pnbits':465, 'rnbits':384, 'cost_S':None, 'deg_h_S':None, 'label':"+2^34+2^32-2^30-2^24-2^14-2^10+2^4 Hw 7"},
    {'u':0x4affff7ec, 'b':-5, 'pnbits':465, 'rnbits':384, 'cost_S':None, 'deg_h_S':None, 'label':"+2^34+2^32-2^30-2^28-2^11-2^4-2^2 Hw 7"},
    #
    {'u':-0x1dff7f7f7fe0, 'b': 3, 'pnbits':614, 'rnbits':512, 'cost_S':None, 'deg_h_S':None, 'label':"-2^45+2^41+2^31+2^23+2^15+2^5 Hw 6"},
    {'u':-0x1e0102002080, 'b': 3, 'pnbits':614, 'rnbits':512, 'cost_S':None, 'deg_h_S':None, 'label':"-2^45+2^41-2^32-2^25-2^13-2^7 Hw 6"},
    {'u': 0x1e07fff08080, 'b': 3, 'pnbits':614, 'rnbits':512, 'cost_S':None, 'deg_h_S':None, 'label':"+2^45-2^41+2^35-2^20+2^15+2^7 Hw 6"},
    {'u':-0x1e077ffc2000, 'b': 2, 'pnbits':614, 'rnbits':512, 'cost_S':None, 'deg_h_S':None, 'label':"-2^45+2^41-2^35+2^31+2^18-2^13 Hw 6"},
    {'u':-0x1d8008200002, 'b':17, 'pnbits':614, 'rnbits':512, 'cost_S':None, 'deg_h_S':None, 'label':"-2^45+2^41+2^39-2^27-2^21-2 Hw 6"},
    {'u': 0x1cfffffefc02, 'b': 2, 'pnbits':614, 'rnbits':512, 'cost_S':None, 'deg_h_S':None, 'label':"+2^45-2^42+2^40-2^16-2^10+2 Hw 6"},
]
# KSS40 curves with seed u = [415, 1165, 1205, 1955] mod 2370
test_vector_kss40 = [
    {'u': 0x423f1,  'a': 3, 'pnbits':377,'rnbits':258, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^18+2^13+2^10-2^4+1              Hw2NAF=5"},#
    {'u':-0x6d12d,  'a': 2, 'pnbits':393,'rnbits':270, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^19+2^16+2^14-2^12-2^8-2^6+2^4+2^2-1 Hw2NAF=9"},#
    #
    {'u':0x3c71c0f, 'a': 1, 'pnbits':551,'rnbits':384, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^26-2^22+2^19-2^16+2^13-2^10+2^4-1 Hw2NAF=8"},#
    #
    {'u':-0x3baff7fc7, 'a':15, 'pnbits':726,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^34+2^30+2^26+2^24+2^15+2^6-2^3+1 Hw2NAF=8"},#
    {'u':-0x3c83ec03f, 'a': 6, 'pnbits':727,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^34+2^30-2^27-2^22+2^16+2^14-2^6+1 Hw2NAF=8"},#
    {'u': 0x3c9402041, 'a': 6, 'pnbits':727,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^34-2^30+2^27+2^24+2^22+2^13+2^6+1 Hw2NAF=8"},#
    {'u':-0x3b4180409, 'a': 1, 'pnbits':726,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^34+2^30+2^28-2^26-2^21+2^19-2^10-2^3-1 Hw2NAF=9"},#
    {'u':-0x3c0e03819, 'a': 1, 'pnbits':726,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^34+2^30-2^24+2^21-2^14+2^11-2^5+2^3-1 Hw2NAF=9"},#
    {'u':-0x3ceff89f9, 'a': 1, 'pnbits':727,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^34+2^30-2^28+2^24+2^15-2^11-2^9+2^3-1 Hw2NAF=9"},#
    {'u':-0x3d010fb61, 'a': 1, 'pnbits':727,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^34+2^30-2^28-2^20-2^16+2^10+2^7+2^5-1 Hw2NAF=9"},#
    {'u': 0x3ba0fceff, 'a': 1, 'pnbits':726,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^34-2^30-2^27+2^25+2^20-2^14+2^12-2^8-1 Hw2NAF=9"},#
    {'u': 0x3bfc374ff, 'a': 1, 'pnbits':726,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^34-2^30-2^22+2^18-2^15-2^12+2^10+2^8-1 Hw2NAF=9"},#
    {'u': 0x3c4030247, 'a': 1, 'pnbits':726,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^34-2^30+2^26+2^18-2^16+2^9+2^6+2^3-1 Hw2NAF=9"},#
    {'u': 0x3c7811017, 'a': 1, 'pnbits':727,'rnbits':512, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^34-2^30+2^27-2^23+2^16+2^12+2^5-2^3-1 Hw2NAF=9"},#
]
test_vector_kss54 = [
    {'u':-0x2886,   'b':-9, 'pnbits':283,'rnbits':255, 'deg_h_S':None,'cost_S':None, 'deg_h_C':9, 'cost_C':232, 'label':"Mike Scott, u=-2^13-2^11-2^7-2^3+2 Hw2NAF=5"},#
    {'u': 0x31be,   'b':-4, 'pnbits':289,'rnbits':260, 'deg_h_S':None,'cost_S':None, 'deg_h_C':9, 'cost_C':232, 'label':"u=+2^14-2^12+2^9-2^6-2               Hw2NAF=5"},#
    #
    {'u':-0x16e7c2, 'b': 3, 'pnbits':427,'rnbits':384, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^21+2^19+2^16+2^13-2^11+2^6-2    Hw2NAF=7"},#
    {'u':-0x170c04, 'b':15, 'pnbits':427,'rnbits':384, 'deg_h_S':None,'cost_S':None, 'label':"u=-2^21+2^19+2^16-2^12+2^10-2^2      Hw2NAF=6"},#
    {'u':0x170106, 'b':-17, 'pnbits':427,'rnbits':384, 'deg_h_S':None,'cost_S':None, 'label':"u=+2^21-2^19-2^16+2^8+2^3-2          Hw2NAF=6"},#
    #
    {'u':-0xc400240, 'b':28, 'pnbits':569, 'rnbits':512, 'cost_S':None, 'deg_h_S':None, 'label':"-(2^27+2^26+2^22+2^9+2^6) Hw 5"},
    {'u': 0xc404042, 'b':-4, 'pnbits':569, 'rnbits':512, 'cost_S':None, 'deg_h_S':None, 'label':"eprint 2018/193 2^27+2^26+2^22+2^14+2^6+2 Hw 6"},
]
