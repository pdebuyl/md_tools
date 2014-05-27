# -*- coding: utf-8 -*-
# Copyright 2014 Pierre de Buyl
#
# This file is part of md_tools
#
# md_tools is free software and is licensed under the modified BSD license (see
# LICENSE file).

import numpy as np
cimport numpy as np
from libc.math cimport sqrt, floor, ceil
cimport cython

cdef int count_max=64
cdef int list_box=3

cdef class CellListRDF(object):
    cdef double [:,::1] data
    cdef double[::1] bounds
    cdef double[3] dl
    cdef double density
    cdef int[3] Ncells
    cdef int [:,:,:,::1] cells
    cdef int [:,:,::1] cell_count
    def __init__(self, double[:,::1] data, double[::1] bounds, double density):
        cdef int i
        self.data = np.asarray(data)
        self.bounds = np.zeros(3, dtype=np.float64)
        for i in range(3):
            self.bounds[i] = bounds[i]
        self.density = density
        l_tmp = (32./density)**(1./3.)
        for i in range(3):
            self.Ncells[i] = int(ceil(bounds[i]/l_tmp))
            self.dl[i] = bounds[i]/self.Ncells[i]
        self.cells = np.zeros([self.Ncells[0], self.Ncells[1], self.Ncells[2], count_max], np.int32)
        self.cell_count = np.zeros([self.Ncells[0], self.Ncells[1], self.Ncells[2]], np.int32)
        
    def info(self):
        print "bounds: ", self.bounds[0], self.bounds[1], self.bounds[2]
        print "dl:     ", self.dl[0], self.dl[1], self.dl[2]
        print "Ncells: ", self.Ncells[0], self.Ncells[1], self.Ncells[2]
    def fill(self):
        cdef int i, xi, yi, zi, local_count
        self.cell_count = np.zeros([self.Ncells[0], self.Ncells[1], self.Ncells[2]], np.int32)
        for i in range(self.data.shape[0]):
            xi = int(floor(self.data[i,0]/self.dl[0]))
            if xi<0:
                xi=xi+self.Ncells[0]
            elif xi>=self.Ncells[0]:
                xi=xi-self.Ncells[0]
            yi = int(floor(self.data[i,1]/self.dl[1]))
            if yi<0:
                yi=yi+self.Ncells[1]
            elif yi>=self.Ncells[1]:
                yi=yi-self.Ncells[1]
            zi = int(floor(self.data[i,2]/self.dl[2]))
            if zi<0:
                zi=zi+self.Ncells[2]
            elif zi>=self.Ncells[2]:
                zi=zi-self.Ncells[2]
            if xi>=self.Ncells[0] or yi>self.Ncells[1] or zi>self.Ncells[2]:
                print "bug"
            if xi<0 or yi<0 or zi<0:
                print "neg bug"
                break
            local_count = self.cell_count[xi,yi,zi]
            if local_count>=count_max:
                print "max in", xi, yi, zi
            self.cells[xi,yi,zi,local_count] = i
            self.cell_count[xi,yi,zi] = local_count+1

    def minmax(self):
        cdef int i,j,k
        cdef int mymin, mymax, local_count
        mymin = 10**10
        mymax = -100
        for i in range(self.Ncells[0]):
            for j in range(self.Ncells[1]):
                for k in range(self.Ncells[2]):
                    local_count = self.cell_count[i,j,k]
                    if local_count>mymax:
                        mymax=local_count
                    if local_count<mymin:
                        mymin=local_count
        return mymin, mymax
    
    cdef int get_idx(self, int i, int j, int k, int count):
        if count<self.cell_count[i,j,k]:
            return self.cells[i,j,k, count]
        else:
            return -1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef double min_dist(self, double[3] x1, int p2):
        cdef np.intp_t i
        cdef double d, b, distsqr
        distsqr=0.
        for i in range(3):
            d = x1[i]-self.data[p2,i]
            b = self.bounds[i]
            if d<-b/2.:
                d=d+b
            elif d>b/2.:
                d=d-b
            distsqr += d*d
        return distsqr**0.5
            
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def compute_rdf(self, double dr, double[::1] g):
        cdef np.intp_t i,j,k,l,p1,p2,r_i
        cdef int N=g.shape[0]
        cdef double rmax = (N-1)*dr
        cdef double[3] x1
        for i in range(3):
            x1[i]=0.
        cdef int[3] idx, idx_shift, idx_loop
        for p1 in range(self.data.shape[0]):
            for i in range(3):
                x1[i] = self.data[p1,i]
                idx[i] = int(floor(x1[i]/self.dl[i]))
            for k in range(3):
                idx_shift[k]=0
            for j in range(list_box**3):
                for k in range(3):
                    if j%list_box**k==0:
                        idx_shift[2-k] = (idx_shift[2-k]+1)%list_box
                for k in range(3):
                    idx_loop[k] = idx[k]+idx_shift[k]-1
                    if idx_loop[k]>=self.Ncells[k]:
                        idx_loop[k] -= self.Ncells[k]
                    elif idx_loop[k]<0:
                        idx_loop[k] += self.Ncells[k]
                for k in range(self.cell_count[idx_loop[0], idx_loop[1], idx_loop[2]]):
                    p2 = self.cells[idx_loop[0], idx_loop[1], idx_loop[2], k]
                    r = self.min_dist(x1, p2)
                    r_i = int(floor(r/dr))
                    if r_i<N:
                        g[r_i] += 1
        return


            
