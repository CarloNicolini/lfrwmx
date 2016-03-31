# lfrwmx
Matlab wrapper for weighted Lancichinetti-Fortunato-Radicchi benchmark graphs

# Lancichinetti-Fortunato-Radicchi benchmark functions
Lancichinetti-Fortunato-Radicchi benchmark is a network generator for benchmark of community detection algorithms.

In this release I've adapted and partially rewritten the freely-available code (https://sites.google.com/site/santofortunato/inthepress2) to generate weighted networks as described in the paper by Lancichinetti A. and Fortunato S., *Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities*, Phys. Rev. E 80, 016118.

You can compile the LFR benchmark by specifying:

    $> cmake -DMATLAB_SUPPORT=True ..
    $> make lfrw_mx

The `lfrw_mx` function is self-documented:

    LFRW: Lancichinetti-Fortunato-Radicchi network generator
    This program produces weighted adjacency matrices for graph community detection benchmark
    [WeightedGraph, membership] = lfrw('argname',argvalue);
    Input:
    'N'
        Desidered number of nodes in the network
    'k'
        Desidered average degree of the network
    'maxk'
        Desidered max degree of the network
    'minc'
        Desidered minimum size of the community
    'maxc'
        Desidered maximum size of the community
    'mut'
        Desidered community topological mixing coefficient (range [0,1])
    'muw'
        Desidered weights mixing coefficient (range [0,1])
    't1'
        Desidered exponent for the degree distribution
    't2'
        Desidered exponent for the community size distribution
    'beta'
        Desidered beta exponent
    'C'
        Desidered clustering coefficient
    Error using lfrw_mx
    Error at argument: 0: Not enough input arguments. 