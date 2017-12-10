# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 16:03:40 2017

@author: amit
"""

import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
import scipy.constants as con
import scipy.integrate as gra
import matplotlib.pyplot as plt
from numpy.linalg import inv
from numpy import *
from scipy.integrate import nquad
from numpy import linalg as LA
from numpy.linalg import *

no_of_electrons_taken=19

A = 1
a = 1

energyold=np.zeros(no_of_electrons_taken,dtype=complex)
energynew=np.zeros(no_of_electrons_taken,dtype=complex)
energyini=np.zeros(no_of_electrons_taken,dtype=complex)

A = 1
a = 1

band = 0

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

    k_spectrum_b0.append([k, np.real(eVal[band]), eVec[band]])
    #evl.append(eVal[band])

def u(k_i, x_i):
    return k_spectrum_b0[k_i][2][x_i]

def rho_H(x_i):
    rh = 0
    for k_i in range(len(kl)):
        rh += -np.conjugate(u(k_i, x_i))
    return rh

def rho_HF(x1_i, x2_i):
    rhf1 = 0
    rhf2 = 0
    for k1_i in range(len(kl)):
        for k2_i in range(len(kl)):
            rhf1 += -np.conjugate(u(x1_i, k2_i)) \
            * np.conjugate(u(x2_i, k1_i)) \
            * u(x1_i, k1_i) * u(x2_i, k2_i)
    for k_i in range(len(kl)):
        rhf2 += np.conjugate(u(k_i, x1_i))*u(k_i, x1_i)
    if (rhf2 == 0):
        return 0
    else:
        return rhf1/rhf2

def deno(x1_i, x2_i):
    return 1/( abs( x[x1_i] - x[x2_i] ) + h/2 )

H_i = np.zeros((sample, sample), dtype = np.complex)
for loops in range(10):
    for i in range(sample):
        Sum=0
        k=0
        Sum = Sum + ((rho_H(0)-rho_HF(i,0))/deno(i,0))
        Sum = Sum + ((rho_H(sample-1)-rho_HF(i,sample-1))/deno(i,sample-1))
        for j in range(sample):
            if (k == 0):
                Sum += 4 * ( rho_H(i) - rho_HF(i, j) / deno(i, j) )
                k = 1
            elif (k == 1):
                Sum += 2 * ( rho_H(i) - rho_HF(i, j) / deno(i, j) )
                k = 0
        Sum *= 0.01/3
        H_i[i][i] -= np.real(Sum)

    H_new = H + H_i
    k_spectrum_b0 = []
    for k in kl:
        eVal, eVec =  np.linalg.eig(H_new)
        eVec = np.transpose(eVec)
        ev_list = list(zip( eVal, eVec ))
        ev_list.sort(key=lambda tup:tup[0], reverse=False)
        eVal, eVec = zip(*ev_list)
        k_spectrum_b0.append([k, np.real(eVal[band]), eVec[band]])

    print(eVal[0])