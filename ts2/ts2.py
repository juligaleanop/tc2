#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 17:18:40 2026

@author: Julian 
"""
import numpy as np
import scipy.signal as sp
import sympy as syp
from pytc2.sistemas_lineales import bodePlot, pzmap, sos2tf_analog

z,p,k = sp.buttap(6)
num,den = sp.zpk2tf(z,p,k)
T = sp.TransferFunction(num,den)
pzmap(T)
bodePlot(T)

VI,VO,VA,V1 = syp.symbols("VI,VO,VA,V1")
G,C,K = syp.symbols("G,C,K")
s = syp.symbols("s")
so_1 = syp.solve([
    3*G*V1-G*VI-G*VO-G*VA,
    ((s*C)**2/G**2)*VO-VA,
    ((s*C/G)/K)*VO+V1
    ],
    [VO,VI,VA,V1])

T_KHN = so_1[VO]/so_1[VI]
print(T_KHN) #verifica!!!

VI,VO,VA = syp.symbols("VI,VO,VA")
G1,G2,G3,C = syp.symbols("G1,G2,G3,C")
s = syp.symbols("s")
so_2 = syp.solve([
    (G1+G2+G3+s*C)*VA-G3*VI-G1*VO,
    (s*C)/G2*VO+VA,
    ],
    [VO,VI,VA])

T_MFB = so_2[VO]/so_2[VI]
print(T_MFB) #verifica tambien!!!

#Variables de la SOS1 - Sallen-Key
G_SK = 1
G_4_SK = 0.07
K_SK = 1+G_4_SK/G_SK
C_SK = 1

#Variables de la SOS2 - KHN
G_KHN = 1
G_7_KHN = 1.1213
M_KHN = 1+G_7_KHN/G_KHN
C_KHN = 1

#Variables de la SOS3 - Multiple Feedback
G_MFB = 1
G_3_MFB = 2.956
K_MFB = G_3_MFB/G_MFB
C_1_MFB = 9.574
C_2_MFB = 0.10455

SOS = np.array([[0,0,K_SK*(G_SK**2)/(C_SK**2),1,G_SK*(3-K_SK)/C_SK,(G_SK**2)/(C_SK**2)],
                [0,0,-(G_KHN**2)/(C_KHN**2),1,3*G_KHN/(C_KHN*M_KHN),(G_KHN**2)/(C_KHN**2)],
                [0,0,-K_MFB*(G_MFB**2)/(C_1_MFB*C_2_MFB),1,(2*G_MFB+G_3_MFB)/C_1_MFB,(G_MFB**2)/(C_1_MFB*C_2_MFB)]])

TF = sos2tf_analog(SOS)
pzmap(TF)
w, TF_w0 = sp.freqresp(TF,w=[0.5,1])
print(f"Módulo en w=1/2: {20*np.log10(np.abs(TF_w0[0]))}")
print(f"Módulo en w=1: {20*np.log10(np.abs(TF_w0[1]))}")