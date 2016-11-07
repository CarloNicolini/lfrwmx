#   This file is part of PACO-PArtitioning Clustering Optimization a program
#   to find network partitions using modular solvers and quality functions.
#   Copyright (C) 2015 Carlo Nicolini <carlo.nicolini@iit.it>
#   PACO is free software you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation either
#   version 3 of the License, or (at your option) any later version.
#
#   Alternatively, you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation either version 2 of
#   the License, or (at your option) any later version.
#
#   PACO is distributed in the hope that it will be useful, but WITHOUT ANY
#   WARRANTY without even the implied warranty of MERCHANTABILITY or FITNESS
#   FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License and a copy of the GNU General Public License along with
#   PACO. If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
import ctypes
import numpy as np
cimport numpy as np
#np.import_array() 
from ctypes import c_double

# Cython imports
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
import cython
from cython.view cimport array as cvarray

ctypedef map[string, int] params_map

cdef extern from "benchm.h":
    cdef int benchmark_py(int excess, int defect, int num_nodes, double average_k, int max_degree, double tau, double tau2, double mixing_parameter, double mixing_parameter2, double beta, int overlapping_nodes, int overlap_membership, int nmin, int nmax, int fixed_range, double ca, vector[double] &W, vector[int] &membership)

cdef extern from "set_parameters.h":
    cdef cppclass Parameters:
        Parameters() except +
        int arrange() except +
        void print_parameters()
        int num_nodes
        double average_k
        int max_degree
        double tau_degree
        double tau_commsize
        double mixing_parameter_topological
        double mixing_parameter_weights
        double beta
        int overlapping_nodes
        int overlap_membership
        int nmin
        int nmax
        int fixed_range
        int excess
        int defect
        int randomf
        double clustering_coeff

@cython.boundscheck(False)
@cython.wraparound(False)
def pylfrw(N, avgk, maxk, mut, muw, **kwargs):    
    args = ['tau_degree', 'tau_commsize', 'minc', 'maxc', 'beta', 'on', 'om', 'C', 'verbosity']
    
    args_diff = set(kwargs.keys()) - set(args)
    if args_diff:
        raise Exception("Invalid args:" + str(tuple(args_diff)) + "as graph_rep: valid arguments are " + str(args))

    cdef Parameters pars
    pars.num_nodes = int(N)
    pars.average_k = int(avgk)
    pars.max_degree = int(maxk)
    pars.mixing_parameter_topological = float(mut)
    pars.mixing_parameter_weights = float(muw)

    pars.tau_degree  = float(kwargs.get("tau_degree", pars.tau_degree))
    pars.tau_commsize = float(kwargs.get("tau_commsize", pars.tau_commsize))
    pars.nmin = int(kwargs.get("minc", pars.nmin))
    pars.nmax = int(kwargs.get("maxc", pars.nmax))
    pars.overlapping_nodes = int(kwargs.get("on",pars.overlapping_nodes))
    pars.overlap_membership = int(kwargs.get("om",pars.overlap_membership))
    pars.clustering_coeff = float(kwargs.get("C",pars.clustering_coeff))
    try:
        pars.arrange()
    except RuntimeError:
        raise

    pars.print_parameters()

    cdef vector[double] W
    cdef vector[int] M
    
    benchmark_py(pars.excess, pars.defect, pars.num_nodes, 
        pars.average_k, pars.max_degree, pars.tau_degree,
         pars.tau_commsize, pars.mixing_parameter_topological, 
         pars.mixing_parameter_weights, pars.beta, pars.overlapping_nodes, 
         pars.overlap_membership, pars.nmin, pars.nmax,
         pars.fixed_range, pars.clustering_coeff,
         W,
         M)
    return np.reshape(W,[pars.num_nodes,pars.num_nodes]),M
    