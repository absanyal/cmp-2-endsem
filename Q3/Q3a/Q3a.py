# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 20:06:21 2017

@author: amit
"""

import numpy as np

A = 1
a = 1

I = complex(0, 1)

atoms = 10000

sample = 50

def V(x):
    return -A * (1 + np.cos(2 * np.pi * x / a ))

x = list(np.linspace( 0, 1, num = sample + 1 ))
x.pop(-1)


kl = list(np.linspace( -np.pi/a, np.pi/a, atoms * a ))
kl.pop(-1)
#kl = [-np.pi/a]

h = x[1] - x[0]

H = np.zeros((sample, sample), dtype = np.complex)

f = open('spectrum.dat', 'w')

for k in kl:
    p = ''
    for i in range( len(x) ):
        H[i][i] = 1/h**2 + (k**2) / 2 + V( x[i] )
    for i in range( 0, len(x) - 1 ):
        H[i][i+1] = -1/(2 * h **2) - I * k / (2 * h)
    for i in range( 1, len(x) ):
        H[i][i-1] = -1/(2 * h **2) + I * k / (2 * h)
    H[0][sample - 1] = -1/(2 * h **2) + I * k / (2 * h)
    H[sample - 1][0] = -1/(2 * h **2) - I * k / (2 * h)

    eH = list(np.real(np.linalg.eigvals(H)))
    eH.sort()

    p += str(k) + '\t'

    for ev in eH:
        p += str(ev) + '\t'

    p += '\n'

    f.write(p)

f.close()
