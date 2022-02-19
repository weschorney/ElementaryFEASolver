# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 17:47:12 2021

@author: wes_c
"""

import math

def make_elements(Lx, Ly, Nx, Ny):
    x_cuts = [i * (Lx/Nx) for i in range(Nx + 1)]
    y_cuts = [i * (Ly/Ny) for i in range(Ny + 1)]
    return x_cuts, y_cuts

def get_global_nodes(ele_number, Nx):
    level = math.ceil(ele_number/Nx)
    x_num = ele_number - Nx*((ele_number - 1)//Nx)
    xbase_ln = x_num
    xbase_rn = x_num + 1
    global_nodes = [(level-1)*(Nx+1)+xbase_ln, (level-1)*(Nx+1)+xbase_rn,
                    level*(Nx+1)+xbase_rn, level*(Nx+1)+xbase_ln] #nondim order
    return global_nodes, level, x_num

def get_local_nodes(ele_number, level, x_num, x_cuts, y_cuts):
    #returns x vec: [x1 x2 x3 x4]
    x_ln, x_rn = x_num, x_num + 1
    local_nodes = [(x_cuts[x_ln - 1], y_cuts[level - 1]), 
                   (x_cuts[x_rn - 1], y_cuts[level - 1]),
                   (x_cuts[x_rn - 1], y_cuts[level]),
                   (x_cuts[x_ln - 1], y_cuts[level])] #nondim order
    return local_nodes

if __name__ == '__main__':
    x_cuts, y_cuts = make_elements(30, 120, 3, 10)
    gn, level, x_num = get_global_nodes(9, 3)
    ln = get_local_nodes(9, level, x_num, x_cuts, y_cuts)
