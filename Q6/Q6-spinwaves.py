# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 13:15:21 2017

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


psi_x_1 = np.zeros(len(x))

xn = [x[0] + i * h for i in range((atoms - 1) * sample)]

def u_0(k_i, x_i):
    return k_spectrum_b0[k_i][2][x_i % sample]

def W_0_n(n):
    W_x_0 = np.zeros(len(xn), dtype = np.complex)
    for i in range(len(xn)):
          for j in range(len(kl)):
              W_x_0[i] += u_0(j, i) * np.exp(I * kl[j] * xn[i]) \
              * np.exp( -I * kl[j] * n * a ) / \
              np.sqrt(len(kl))
    return W_x_0

def u_1(k_i, x_i):
    return k_spectrum_b1[k_i][2][x_i % sample]

def W_1_n(n):
    W_x_1 = np.zeros(len(xn), dtype = np.complex)
    for i in range(len(xn)):
          for j in range(len(kl)):
              W_x_1[i] += u_1(j, i) * np.exp(I * kl[j] * xn[i]) \
              * np.exp( -I * kl[j] * n * a ) / \
              np.sqrt(len(kl))
    return W_x_1

n = 10
p0 = np.real(np.conjugate(W_0_n(n)) * W_0_n(n))
p1 = np.real(np.conjugate(W_1_n(n)) * W_1_n(n))
plt.plot(xn, p0, label = 'Wannier for ground state' )
plt.plot(xn, p1, label = 'Wannier for first excited state')
plt.legend()
plt.xlabel('x')
plt.ylabel('Squared norm of Wannier function')
#plt.savefig('wannier.pdf')

"""

##################################################

THIS INTEGRATION WORKS
It runs very slowly. So the code has been run once
and the value has been saved to let the code do
the remaining calculations quickly.

##################################################
def fd_0(x1_i, x2_i):
    return \
    np.real(\
    np.conjugate(W_0_n(4)[x1_i] \
    * np.conjugate(W_0_n(5)[x2_i]) \
    * W_0_n(4)[x1_i] * W_0_n(5)[x2_i] \
    * 1/( abs(xn[x1_i] - xn[x2_i] ) + h/2 )) ) \
    / (4 * np.pi)

def fd_1(x1_i, x2_i):
    return \
    np.real(\
    np.conjugate(W_1_n(4)[x1_i] \
    * np.conjugate(W_1_n(5)[x2_i]) \
    * W_1_n(4)[x1_i] * W_1_n(5)[x2_i] \
    * 1/( abs(xn[x1_i] - xn[x2_i] ) + h/2 )) ) \
    / (4 * np.pi)

J0 = 0
J1 = 0
for x1_i in range(2*len(xn)):
    for x2_i in range(2*len(xn)):
        print(x1_i,x2_i)
        J0 += fd_0(x1_i, x2_i)
        J1 += fd_1(x1_i, x2_i)

print(J0, J1)
"""
J0= 15.49
J1= 13.42

nu = 2
s = 0.5

def gamma_k(k):
    return ( np.exp(I * k * (a) ) + np.exp(I * k * (-a) ) ) / nu

#For ground state
E0_0 = -J0 * s**2 * nu * atoms

E_k_0 = np.array([ \
E0_0 + 2 * J0 * s * (1 - gamma_k(k)) for k in kl \
 ])

#For exc state
E0_1 = -J1 * s**2 * nu * atoms

E_k_1 = np.array([ \
E0_1 + 2 * J1 * s * (1 - gamma_k(k)) for k in kl \
 ])

#plt.plot(kl, E_k_0, label = 'Dispersion for ground state')
#plt.plot(kl, E_k_1, label = 'Dispersion for first excited state')
#plt.xlabel('k')
#plt.ylabel('Energy')
#plt.legend(framealpha = 0.2)
#plt.savefig('spin-wave-disp.pdf')

