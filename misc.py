# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 17:56:40 2021

@author: wes_c
"""

from sympy import Symbol
from sympy.matrices import Matrix
from sympy.integrals import integrate

x = Symbol('x')
n = Symbol('n')

N1 = ((1-x)/2)*((1-n)/2)
N2 = ((1+x)/2)*((1-n)/2)
N3 = ((1+x)/2)*((1+n)/2)
N4 = ((1+x)/2)*((1+n)/2)
dN1dx = (n-1)/4
dN1dn = (x-1)/4
dN2dx = (1-n)/4
dN2dn = (-x-1)/4
dN3dx = (n+1)/4
dN3dn = (x+1)/4
dN4dx = (-n-1)/4
dN4dn = (1-x)/4

SHAPE_FNS = [N1, N2, N3, N4]
DX_FNS = [dN1dx, dN2dx, dN3dx, dN4dx]
DN_FNS = [dN1dn, dN2dn, dN3dn, dN4dn]

def jacobian(x_vec, y_vec):
    J11 = sum(x_vec[i]*DX_FNS[i] for i in range(4))
    J12 = sum(y_vec[i]*DX_FNS[i] for i in range(4))
    J21 = sum(x_vec[i]*DN_FNS[i] for i in range(4))
    J22 = sum(y_vec[i]*DN_FNS[i] for i in range(4))
    return Matrix([[J11, J12], [J21, J22]])

def jacobian_inv(jac):
    return jac.inv()

def integral_prod(jac, j_inv, i, j):
    p1 = j_inv[0, 0]*DX_FNS[i] + j_inv[0, 1]*DN_FNS[i]
    p2 = j_inv[0, 0]*DX_FNS[j] + j_inv[0, 1]*DN_FNS[j]
    p3 = j_inv[1, 0]*DX_FNS[i] + j_inv[1, 1]*DN_FNS[i]
    p4 = j_inv[1, 0]*DX_FNS[j] + j_inv[1, 1]*DN_FNS[j]
    return p1 * p2 * p3 * p4 * jac.det()

def double_integral(int_prod):
    out = integrate(int_prod, (x, 0, 1))
    out = integrate(out, (n, 0, 1))
    return float(out)


