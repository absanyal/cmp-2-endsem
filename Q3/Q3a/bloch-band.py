# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 20:18:17 2017

@author: amit
"""

import numpy as np

A = 1
a = 1
h = 1
m = 1
c = h**2 / (2 * m)

G = np.linspace(-2, 2, 5)

f = open('spectrum-bloch.dat', 'w')

ql = list(np.linspace(-np.pi/a, np.pi/a, 10000))
ql.pop(-1)

for q in ql:
    p = ''
    H = np.zeros((5,5))
    for i in range(5):
        H[i][i] = c * (q - (2 * np.pi / a) * G[i])**2 - A
    for i in range(0, 4):
        H[i][i+1] = -A/2
    for i in range(1, 5):
        H[i][i-1] = -A/2

    eH = list(np.real(np.linalg.eigvals(H)))
    eH.sort()

    p += str(q) + '\t'

    for ev in eH:
        p += str(ev) + '\t'

    p += '\n'

    f.write(p)

f.close()
