# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 18:57:41 2017

@author: amit
"""

import numpy as np
import matplotlib.pyplot as plt

A = 1
a = 1

band = 4

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
    evl.append(np.real(eVal[band]))

n = 2/a
E_f = np.real(max(evl))
l = np.sqrt( 3 * n / ( 2 * E_f) )

def u(k_i, x_i):
    return k_spectrum_b0[k_i][2][x_i]

#xn = [h * i for i in range(atoms*sample)]
#xn.pop(-1)

def fd(k1_i, k2_i, x1_i, x2_i):
    return np.real(np.conjugate(u(k1_i, x1_i)) \
    * np.conjugate(u(k2_i, x2_i)) \
    * u(k1_i, x1_i) * u(k2_i, x2_i) \
    * 1/( abs(x[x1_i] - x[x2_i] ) + h/2 )) \
    * np.exp(- l * abs(x[x1_i] - x[x2_i] ) ) \
    / (4 * np.pi)

def fe(k1_i, k2_i, x1_i, x2_i):
    return np.real(np.conjugate(u(k1_i, x1_i)) \
    * np.conjugate(u(k2_i, x2_i)) \
    * u(k1_i, x2_i) * u(k2_i, x1_i) \
    * 1/( abs(x[x1_i] - x[x2_i] ) + h/2 )) \
    * np.exp( I * ( kl[k1_i] - kl[k2_i] ) * ( x[x1_i] - x[x2_i] ) ) \
    * np.exp(- l * abs(x[x1_i] - x[x2_i] ) ) \
    / (4 * np.pi)

J_k = []
K_k = []
for k1_i in range(len(kl)):
    J = 0
    K = 0
    print(kl[k1_i])
    for k2_i in range(len(kl)):
        for x1_i in range(len(x)):
            for x2_i in range(len(x)):
                J += fd(k1_i, k2_i, x1_i, x2_i) * (1 - kdel(k1_i, k2_i))
                K += fe(k1_i, k2_i, x1_i, x2_i) * (1 - kdel(k1_i, k2_i))
    J_k.append(J)
    K_k.append(K)


J_k = np.real(np.array(J_k))
K_k = np.real(np.array(K_k))
correction = np.real(2 * J_k - K_k)


print(np.real(J_k))
print(np.real(K_k))
print(correction)

evl = np.real(np.array(evl))
kl = np.array(kl)

plt.plot(kl, evl)
plt.plot(kl, evl + correction)
plt.show()

filename = 'spectrum-screen-' + str(band) + '.dat'
f = open(filename, 'w')

ec = evl + correction

for i in range(len(kl)):
    p = str(kl[i]) + '\t' + str(ec[i]) + '\n'
    f.write(p)

f.close()

#x = np.array(x)
#
#u0 = np.array([ np.absolute(u( 0, xi )) \
#for xi in range(len(x)) ])
#
#plt.plot(x, u0)
#
#plt.show()


#J = 0
#K = 0
#k1_i = 1
#k2_i = 2
#for x1_i in range(len(x)):
#            for x2_i in range(len(x)):
#                J += fd(k1_i, k2_i, x1_i, x2_i) * (1 - kdel(k1_i, k2_i))
#                K += fe(k1_i, k2_i, x1_i, x2_i) * (1 - kdel(k1_i, k2_i))
#            print(x[x1_i], np.real(J), sep = '\t')
#
#print(J, K)
