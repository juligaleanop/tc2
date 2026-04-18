#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 17:18:40 2026

@author: Julian 
"""
import numpy as np
from scipy.signal import TransferFunction, bode
import sympy as syp
from pytc2.sistemas_lineales import bodePlot, pzmap
from matplotlib import pyplot as plt


V1, V2, Va = syp.symbols("V1, V2, Va")
G1, G2, G3, C1 = syp.symbols("G1, G2, G3, C1")
S = syp.symbols("S")
so = syp.solve([
        Va*(G1+G2)-V1*G1-V2*G2,
        Va*(S*C1+G3)-V1*S*C1
        ],
        [V1,V2,Va])

T_si = so[V2]/so[V1]
print(T_si)

w0 = 1
p = np.array([1, -w0])
q = np.array([1, w0])

T_nu = TransferFunction(p, q)

w=np.arange(0,10**3,0.01)
w, mag, phase = bode(T_nu, w=w)
phase_rad = phase*2*np.pi/360

plt.figure()
plt.subplot(211)
plt.semilogx(w, mag)
plt.ylim(-1,1)
plt.xlim(10**-2,10**2)
plt.grid()
plt.title('Diagrama de Bode')
plt.xlabel('Frecuencia [w]')
plt.ylabel('Magnitud [dB]')

plt.subplot(212)
plt.semilogx(w, phase_rad)
plt.ylim(0,np.pi)
plt.xlim(10**-2,10**2)
plt.grid()
plt.title('Diagrama de Bode')
plt.xlabel('Frecuencia [w]')
plt.ylabel('Fase [grados]')

plt.tight_layout()
plt.show()