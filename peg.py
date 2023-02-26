#!/usr/bin/env python
# 
# Copyright 2013 IIT Bombay.
# Author: Manu T S
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

import numpy as np
from find_smallest import find_smallest
import copy

class peg():

    """
    Progressive edge growth algorithm for generating
    LDPC matrices. The algorithm is obtained from [1]
    """

    def __init__(self, nvar, nchk, degree_sequence):
        self.degree_sequence = degree_sequence
        self.nvar = nvar
        self.nchk = nchk
        self.H = np.zeros((nchk, nvar), dtype = np.int32)
        self.sym_degrees = np.zeros(nvar, dtype = np.int32)
        self.chk_degrees = np.zeros(nchk, dtype = np.int32)
        self.I_edge_chk2var = []
        self.I_edge_var2chk = []
        for chk in range(nchk):
            self.I_edge_chk2var.append([])
            for var in range(nvar):
                self.I_edge_chk2var[chk].append(0)
        for var in range(nvar):
            self.I_edge_var2chk.append([])
            for chk in range(nchk):
                self.I_edge_var2chk[var].append(0)
    
    def grow_edge(self, var, chk):
        self.I_edge_chk2var[chk][var] = 1
        self.I_edge_var2chk[var][chk] = 1
        self.H[chk, var] = 1
        self.sym_degrees[var] += 1
        self.chk_degrees[chk] += 1

    def bfs(self, var):
        var_list = np.zeros(self.nvar, dtype = np.int32)
        var_list[var] = 1
        cur_chk_list = []
        new_chk_list = []
        for i in range(self.nchk):
            cur_chk_list.append(0)
            new_chk_list.append(0)
        chk_Q = []
        var_Q = []
        var_Q.append(var)
        while(True):
            for _vars in var_Q:
                for i in range(self.nchk):
                    if self.H[i, _vars] == 1:
                        if cur_chk_list[i] == 0:
                            new_chk_list[i] = 1
                            chk_Q.append(i)
            var_Q = []
            for _chks in chk_Q:
                for j in range(self.nvar):
                    if self.H[_chks, j] == 1:
                        if var_list[j] == 0:
                            var_list[j] = 1
                            var_Q.append(j)
            chk_Q = []
            if new_chk_list.count(1) == self.nchk:
                new_chk = self.find_smallest_chk(cur_chk_list)
                return new_chk
            elif np.array_equal(new_chk_list, cur_chk_list):
                new_chk = self.find_smallest_chk(cur_chk_list)
                return new_chk
            else:
                cur_chk_list = copy.copy(new_chk_list)

    def find_smallest_chk(self, cur_chk_list):
        index = []
        degree = []
        for i in range(len(cur_chk_list)):
            if cur_chk_list[i] == 0:
                index.append(i)
                degree.append(self.chk_degrees[i])
        return index[find_smallest(degree)]

    def _print(self):
        print "I_edge_chk2var"
        for i in range(len(self.I_edge_chk2var)):
            print self.I_edge_chk2var[i]
        print "I_edge_var2chk"
        for i in range(len(self.I_edge_var2chk)):
            print self.I_edge_var2chk[i]

    def progressive_edge_growth(self):
        for var in range(self.nvar):
            print "edge growth at var", var
            for k in range(self.degree_sequence[var]):
                if k == 0:
                    smallest_degree_chk = find_smallest(self.chk_degrees)
                    self.grow_edge(var, smallest_degree_chk)
                else:
                    chk = self.bfs(var)
                    self.grow_edge(var, chk)

"""
References

"Regular and Irregular Progressive-Edge Growth Tanner Graphs",
Xiao-Yu Hu, Evangelos Eleftheriou and Dieter M. Arnold.
IEEE Transactions on Information Theory, January 2005.
"""
