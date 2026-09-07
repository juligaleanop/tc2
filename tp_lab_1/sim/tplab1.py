#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 18:02:53 2026

@author: julian
"""

import numpy as np
import scipy.signal as sp
import sympy as syp
from pytc2.sistemas_lineales import analyze_sys

# BANDPASS FILTER
# RIPPLE = 2.5 dB, BW = 9.9

amax = 2.5
amin = 15
B = 9.9

#PROTOTIPO PASABAJOS
ws=10
epsq = 10**(amax/10) - 1
print(f"E^2 = {epsq} => E={np.sqrt(epsq)}")
for n in np.arange(1,4):
    amin=10*np.log10(1+epsq*np.cosh(n*np.acosh(ws))**2)
    print(f"Atenuación para orden {n}: {amin} dB")
    
# CON orden 1 alcanza =>

VI,VO,VA = syp.symbols("VI,VO,VA")
G1,G2,G3,C1,C2 = syp.symbols("G1,G2,G3,C1,C2")
s = syp.symbols("s")
so_2 = syp.solve([
    (G1+G2+s*C1+s*C2)*VA-G1*VI-s*C1*VO,
    -s*C2*VA-G3*VO,
    ],
    [VO,VI,VA])

T_MFB = so_2[VO]/so_2[VI]
print(T_MFB) #verifica

