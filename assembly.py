# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 16:56:42 2021

@author: wes_c
"""

import math
import numpy as np
import elements as ele
import misc as m
from itertools import product
from tqdm import tqdm

class FEMSolver:
    def __init__(self, Lx, Ly, Nx, Ny, alpha, L=5, R=7, U=11, D=math.sin):
        self.Lx = Lx
        self.Ly = Ly
        self.Nx = Nx
        self.Ny = Ny
        self.alpha = alpha
        self.num_elements = Nx * Ny
        self.num_pts = (Nx + 1)*(Ny + 1)
        self.x_cuts, self.y_cuts = ele.make_elements(Lx, Ly, Nx, Ny)
        self.global_K = np.zeros((self.num_pts, self.num_pts))
        self.L, self.R, self.U, self.D = L, R, U, D
        self.F = np.zeros((self.num_pts, 1))

    def _map_element_matrix_to_global(self, ele_matrix, global_nodes):
        glob_coord = [(i, global_nodes[i]) for i in range(len(global_nodes))]
        maps = list(product(glob_coord, glob_coord))
        #subtract one since global K starts from 0
        maps = [[(lst[0][0], lst[1][0]), (lst[0][1] - 1, lst[1][1] - 1)]\
                 for lst in maps]
        for sublist in maps:
            self.global_K[sublist[1]] += ele_matrix[sublist[0]]
        return

    def _get_element_matrix(self, element_number):
        gn, level, x_num = ele.get_global_nodes(element_number, self.Nx)
        ln = ele.get_local_nodes(element_number, level, x_num, self.x_cuts,
                                 self.y_cuts)
        x_vector = [e[0] for e in ln]
        y_vector = [e[1] for e in ln]
        element_jacobian = m.jacobian(x_vector, y_vector)
        jac_inverse = m.jacobian_inv(element_jacobian)
        k = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                integral_product = m.integral_prod(element_jacobian,
                                                   jac_inverse, i, j)
                k[i, j] = m.double_integral(integral_product)
        return k, gn

    def _get_boundary_nodes(self):
        bdry_pts = []
        #left bdry
        bdry_pts.append(list(range(1, 1 + len(self.x_cuts)*len(self.y_cuts),\
                                   len(self.x_cuts))))
        #right bdry
        bdry_pts.append(list(range(len(self.x_cuts),\
                                   1 + len(self.x_cuts)*len(self.y_cuts),\
                                   len(self.x_cuts))))
        #up bdry
        bdry_pts.append(list(range(self.num_pts - len(self.x_cuts) + 1,
                                   self.num_pts + 1)))
        #down bdry
        bdry_pts.append(list(range(1, len(self.x_cuts) + 1)))
        return bdry_pts

    def _bdry_apply(self, bdry_pts, vals):
        #need to apply to force vector
        for pt in bdry_pts:
            if callable(vals):
                value = vals(pt)
            else:
                value = vals
            #off by one so subtract
            self.F[pt - 1, 0] = value
            #replace entire vector row by 0,...,1,...,0
            new_row = np.zeros(self.num_pts)
            new_row[pt - 1] = 1
            self.global_K[pt - 1, :] = new_row
        return

    def _apply_boundary_conditions(self):
        #We apply bdry condns s.t. precedence according to d, u, r, l
        #bdry respectively
        bdry_pts = self._get_boundary_nodes()
        self._bdry_apply(bdry_pts[0], self.L)
        self._bdry_apply(bdry_pts[1], self.R)
        self._bdry_apply(bdry_pts[2], self.U)
        self._bdry_apply(bdry_pts[3], self.D)
        return

    def global_assembly(self):
        for i in tqdm(range(self.num_elements)):
            k, gn = self._get_element_matrix(i)
            k = (self.alpha**2) * k
            self._map_element_matrix_to_global(k, gn)
        return

    def solve(self):
        self.global_assembly()
        self._apply_boundary_conditions()
        K_inv = np.linalg.inv(self.global_K)
        u = np.matmul(K_inv, self.F)
        return u

if __name__ == '__main__':
    fs = FEMSolver(30, 40, 15, 10, 10)
    u = fs.solve()
