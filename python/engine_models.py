#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

            #name :  [displacment(cc), mass(kg), power(hp)]
Engines = {'DA-35' : [35, 0.935, -1],
          'DA-50' : [50, 1.33, -1],
          'DA-60' : [60.5, 1.41, -1],
          'DA-70' : [70, 1.61, -1],
          'DA-85' : [85.9, 1.95, -1],
          'DA-100L' : [100, 2.53, 9.8],
          'DA-100inline' : [100, 3.12, -1],
          'DA-120' : [121, 2.25, -1],
          'DA-150L' : [150, 3.3, -1],
          'DA-170' : [171.5, 3.56, -1],
          'DA-200' : [200, 4.95, 19],
          'Evilution NX' : [8.52, 0.435, -1],
          'OS 50SX-H' : [8.17, 0.406, 1.87],
          'OS 105HZ-R' : [17.17, 0.608, 3.75],
          'HK MLD-70' : [70, 2.18, 6.6]}


data = np.array(Engines.values())
f = data[:,1]
A = np.array([data[:,0], np.ones(f.size)])

a_w, b_w  = np.linalg.lstsq(A.T,f)[0]
print a_w,b_w
plt.figure()
plt.xlabel('Engine Size (cc)')  # labels for axes
plt.ylabel('Engine Mass (kg)')
plt.title('Mass Model')
plt.plot(data[:,0], data[:,1], 'ko', label='Data')
plt.plot([0,200], [b_w,a_w*200+b_w], '-k', label='Best Fit')
plt.legend(loc='upper left')
#plt.show()

powerdata = np.array([0])
for row in data[:]:
    if row[2] > 0 and powerdata.size < 3:
        powerdata = np.array(row)
    elif row[2] > 0:
        powerdata = np.append(powerdata, row)

powerdata = powerdata.reshape(5,3)
f = powerdata[:,2]*745.699872 #put the powerdata in Watts

plt.figure(2)
plt.xlabel('Engine Size (cc)')  # labels for axes
plt.ylabel('Engine Power (W)')
plt.title('Power Model')
# data points
plt.plot(powerdata[:,0], f, 'ko', label='Data')

# linear model
A = np.array([powerdata[:,0], np.ones(f.size)])
a, b  = np.linalg.lstsq(A.T,f)[0]
#print a,b
plt.plot([0,200], [b,a*200+b], '-k', label='Linear Fit')

# data from http://www.barnardmicrosystems.com/UAV/engines/2_stroke.html
a_ref = 76.0
b_ref = 502.0
plt.plot([0,200], [b_ref,a_ref*200+b_ref], '--k', label='Barnard Microsystems Fit')

# quadratic model
A = np.array([np.multiply(powerdata[:,0],powerdata[:,0]), powerdata[:,0]])
a, b  = np.linalg.lstsq(A.T,f)[0]
# print a,b
x = np.linspace(0,200, num=100)
y = a*np.multiply(x,x) + b*x
plt.plot(x, y, '-.k', label='Quadratic Fit thru Origin')
plt.legend(loc='upper left')
plt.show()

# print a_ref/a_w  #Watts/kg
