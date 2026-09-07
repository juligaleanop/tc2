#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 20:56:53 2026

@author: julian
"""
import pandas as pd
import numpy as np

f = np.array([100,600,2800,3900,5000,5250,5500,5750,6000,6250,6500,6750,7000,21000,35000,50000,60000,100000])
Vi = np.array([1.55,1.58,1.63,1.15,1.62,1.53,1.44,1.34,1.25,1.17,1.13,1.12,1.1,1.24,1.25,1.26,1.26,1.26])
Vo = np.array([0.084,0.416,1.85,2.75,6.78,7.11,7.23,7.15,6.88,6.48,6.05,5.65,5.3,1.03,0.61,0.430,0.362,0.224])
Dtu = np.array([2600,440,108,84,78,80,80,82,82,82,84,84,84,34.8,21,15,12.4,7.5])
Dt = Dtu*10**(-6)

A = 20*np.log10(Vo/Vi)
An = A-20*np.log10(Vo[8]/Vi[8])
Ph = f*Dt*360

print(f)
print('Ganancia:')
print(f'{An} dB')
print('\nFase:')
print(f'{Ph}Pi')


datos = {
    'f[Hz]' : f,
    'Vi[V]' : Vi,
    'Vo[V]' : Vo,
    'Dt[us]': Dt,
    'A[dB]' : A,
    'An[dB]': An,
    'Ph[°]' : Ph,
}

df=pd.DataFrame(datos)

df.to_csv('../mediciones/osc.csv', index=False)