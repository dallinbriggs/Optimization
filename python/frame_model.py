#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

            #name :  [rotor_size(in), frame_mass(kg)]
Frames = {'Quanum AXE FPV 180' : [4, 0.074],
           'RotorBits Quad' : [10, 0.304],
           'Scarab 360' : [7, 0.240],
           'SK Lightning 480' : [11, 0.510],
           'tarot mini 300' : [6, 0.135],
           'H-King Color 250 FPV' : [5, 0.130],
           'S500 Glass Fiber Quad' : [10, .405],
#           'Diaton ET 180' : [4.5, 0.036],
#           'Diaton ET 200' : [5, 0.039],
#           'Diaton Blade 150' : [3, 0.0215],
#           'Diaton Blade 180' : [4, 0.024],
#           'Diaton Blade 200' : [6, 0.0275],
#           'Totem Q380' : [9, 0.160],
#           'Totem Q330' : [6, 0.150]
           }


data = np.array(Frames.values())
data[:,0] = 0.0254*data[:,0]/2

f = data[:,1]
A = np.array([data[:,0], np.ones(f.size)])
a_w, b_w  = np.linalg.lstsq(A.T,f)[0]
print a_w
plt.figure()
plt.plot(data[:,0], data[:,1], 'ko')
plt.plot([0.05,0.14], [a_w*0.05+b_w,a_w*0.14+b_w], '-k')
plt.show()
