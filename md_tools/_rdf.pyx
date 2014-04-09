# -*- coding: utf-8 -*-
# Copyright 2014 Pierre de Buyl
#
# This file is part of md_tools
#
# md_tools is free software and is licensed under the modified BSD license (see
# LICENSE file).

import numpy as np
cimport numpy as np
from libc.math cimport sqrt, floor
cimport cython

def compute_rdf(r, L, N, species=None, n_species=None, state=None, n_state=None):
    """
    Compute the radial distribution function (rdf). The number of rdf
    components, N_rdf, is 1 for a single species system, n_idx*(n_idx+1)/2 for a
    multi-species system where n_idx=n_species if state is not specified and
    n_idx=sum(n_state) if state is specified.

    Arguments
    ---------

    r: [N,3] array of positions
    L: [3] sides of a cuboid box
    N: number of elements in the rdf

    Optional arguments
    ------------------

    species: [N] array of integer in (0,n_species-1) giving the species of the particles
    n_species: integer, number of species.
    state: [N] array of integer in (0,n_state[species[i]]-1) giving the state of the particles
    n_state: [n_species] array giving the number of allowed states for 

    Returns
    -------

    dx, all_rdf, count
    dx is the radius step
    all_rdf is a [N_rdf,N]
    count is a [N_rdf] array

    """
    cdef int i, n_rdf, n_idx
    cdef double[3] cy_L
    for i in range(3):
        cy_L[i] = L[i]

    x_max = np.min(L)/2.
    dx = x_max/N

    if species is not None:
        assert n_species is not None
        if state is not None:
            assert n_state is not None
            assert len(n_state)==n_species
            n_idx = np.sum(n_state)
        else:
            n_idx = n_species
    else:
        n_idx = 1

    n_rdf = (n_idx)*(n_idx+1)/2
    if n_species is None:
        n_species = 1
    if n_state is None:
        n_state = np.ones( (n_species,), dtype=np.int32)

    result, count = _compute_rdf(r, cy_L, N, dx, n_rdf, n_species, species, state, n_state)
    return dx, np.asarray(result), np.asarray(count)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _compute_rdf(double[:, ::1] r, double[3] L, int N, double dx, int n_rdf,
                 int n_species, int[::1] species, int[::1] state, int[::1] n_state):
    cdef int i, j, coord, idx, si, sj, rdf_idx, n_idx
    cdef double dist
    cdef double dist_sqr
    cdef double inv_dx = 1./dx
    cdef double x_max_sqr = (N*dx)**2
    cdef double pi = np.pi
    cdef double k

    cdef double[:,::1] result = np.zeros( (n_rdf, N) )
    cdef int[::1] idx_root = np.zeros( (n_species,), dtype=np.int32 )

    n_idx = n_state[0]
    idx_root[0] = 0
    for i in range(1,n_species):
        idx_root[i] = idx_root[i-1] + n_state[i-1]
        n_idx += n_state[i]

    cdef int[::1] count = np.zeros( (n_idx,), dtype=np.int32 )

    for i in range(r.shape[0]):
        if species is not None:
            if state is not None:
                si = idx_root[species[i]] + state[i]
            else:
                si = species[i]
        else:
            si = 0
        count[si] += 1
        for j in range(i+1, r.shape[0]):
            if species is not None:
                if state is not None:
                    sj = idx_root[species[j]] + state[j]
                else:
                    sj = species[j]
                if si>sj:
                    rdf_idx = si*(si+1)/2 + sj
                else:
                    rdf_idx = sj*(sj+1)/2 + si
            else:
                rdf_idx=0
            dist_sqr = 0.
            for coord in range(3):
                dist = r[i,coord]-r[j,coord]
                if dist<-L[coord]/2.:
                    dist += L[coord]
                elif dist>L[coord]/2.:
                    dist -= L[coord]
                dist_sqr += dist**2
            if dist_sqr<=x_max_sqr:
                idx = int(floor(sqrt(dist_sqr)*inv_dx))
                result[rdf_idx,idx] += 1

    k = 2*pi*dx
    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            result[i,j] = result[i,j] / (k*((j+0.5)*dx)**2)

    return result, count
