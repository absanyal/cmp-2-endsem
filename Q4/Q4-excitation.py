# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 17:08:04 2017

@author: amit
"""

import numpy as np
import matplotlib.pyplot as plt

A = 1
a = 1

I = complex(0, 1)

atoms = 20 + 1

sample = 20

def V(x):
    return -A * (1 + np.cos(2 * np.pi * x / a ))

x = list(np.linspace( 0, 1, num = sample + 1 ))
x.pop(-1)

def kdel(x, y):
    if (x == y):
        return 1
    else:
        return 0

kl = list(np.linspace( -np.pi/a, np.pi/a, atoms * a ))
kl.pop(-1)
#kl = [-np.pi/a]

h = x[1] - x[0]

k_spectrum_b0 = []
k_spectrum_b1 = []

H = np.zeros((sample, sample), dtype = np.complex)

evl = []

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

    eVal, eVec =  np.linalg.eig(H)
    eVec = np.transpose(eVec)
    ev_list = list(zip( eVal, eVec ))
    ev_list.sort(key=lambda tup:tup[0], reverse=False)
    eVal, eVec = zip(*ev_list)

    k_spectrum_b0.append([k, np.real(eVal[0]), eVec[0]])
    k_spectrum_b1.append([k, np.real(eVal[1]), eVec[1]])

def u0(k_i, x_i):
    return k_spectrum_b0[k_i][2][x_i % sample]

def u1(k_i, x_i):
    return k_spectrum_b1[k_i][2][x_i % sample]

#xn = [h * i for i in range(atoms*sample)]
#xn.pop(-1)

def fd(k1_i, k2_i, x1_i, x2_i):
    return np.real(np.conjugate(u1(k1_i, x1_i)) \
    * np.conjugate(u0(k2_i, x2_i)) \
    * u1(k1_i, x1_i) * u0(k2_i, x2_i) \
    * 1/( abs(x[x1_i] - x[x2_i] ) + h/2 )) / (4 * np.pi)

def fe(k1_i, k2_i, x1_i, x2_i):
    return np.real(np.conjugate(u1(k1_i, x1_i)) \
    * np.conjugate(u0(k2_i, x2_i)) \
    * u1(k1_i, x2_i) * u0(k2_i, x1_i) \
    * 1/( abs(x[x1_i] - x[x2_i] ) + h/2 )) \
    * np.exp( I * ( kl[k1_i] - kl[k2_i] ) * ( x[x1_i] - x[x2_i] ) )\
    / (4 * np.pi)

J = 0
K = 0
k1_i = 0
k2_i = 0
for x1_i in range(len(x)):
    for x2_i in range(len(x)):
        J += fd(k1_i, k2_i, x1_i, x2_i) #* (1 - kdel(k1_i, k2_i))
        K += fe(k1_i, k2_i, x1_i, x2_i) #* (1 - kdel(k1_i, k2_i))

K = np.real(K)

correction = -J + 2 *K

print(correction)