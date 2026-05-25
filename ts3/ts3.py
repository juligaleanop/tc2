#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 12:29:52 2026

@author: julian
"""

import numpy as np
import scipy.signal as sp
from pytc2.sistemas_lineales import bodePlot, pzmap, sos2tf_analog, tf2sos_analog, pretty_print_SOS

ripple = 0.5

z,p,k = sp.cheb1ap(3, ripple)
num,den = sp.zpk2tf(z,p,k)
T_nu = sp.TransferFunction(num,den)
pzmap(T_nu)

num_hp,den_hp = sp.lp2hp(num,den)
SOS_hp = tf2sos_analog(num_hp, den_hp)

print(num_hp, den_hp)
print(SOS_hp)
pzmap(sp.TransferFunction(num_hp,den_hp))