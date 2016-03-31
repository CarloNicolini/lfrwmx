# lfrwmx
Matlab wrapper for weighted Lancichinetti-Fortunato-Radicchi benchmark graphs

# Lancichinetti-Fortunato-Radicchi benchmark functions
Lancichinetti-Fortunato-Radicchi benchmark is a network generator for benchmark of community detection algorithms.

In this release I've adapted and partially rewritten the freely-available code (https://sites.google.com/site/santofortunato/inthepress2) to generate weighted networks as described in the paper by Lancichinetti A. and Fortunato S., *Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities*, Phys. Rev. E 80, 016118.

## Requirements
- Linux with GCC/G++ compiler (or other suitable C++ compiler)
- CMake (www.cmake.org)
- `git`

To compile clone this dataset:

    $> git clone https:/github.com/CarloNicolini/lfrwmx.git
    $> cd lfrwmx
    $> mkdir build
    $> cd build
    $> cmake ..

For the compilation of the MEX file you need MATLAB mex headers, usually installed with a standard Matlab installation.

You can compile the LFR benchmark by specifying:

    $> cmake -DMATLAB_SUPPORT=True ..
    $> make lfrw_mx

The wrapper is also available for Octave, in case you don't have MATLAB, just install the octave-mex headers (typically on Debian based distros: `sudo apt-get install octave liboctave-dev`)

    $> cmake -DOCTAVE_SUPPORT=True ..
    $> make lfrw_mx

## FAQ

### Linking libstdc++.so.6 problem

I can compile Paco for MATLAB but after calling `paco_mx`, MATLAB prompts me with the following error message:

```
Invalid MEX-file '~/paco_mx.mexa64': /usr/local/MATLAB/R2015a/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by ...
```

This problem means that the `libstdc++.so.6` inside the Matlab library folder is pointing to a version of `libstdc++` older than the system one, usually stored in `/usr/lib/x86_64` folder.

To solve the issue you need to redirect the symbolic links in the MATLAB folder to the systemwise `libstdc++`. Hereafter we assume the MATLAB folder to be `/usr/local/MATLAB/R2015a` and the system to be some Linux variant.

Two of the symlinks for libraries need to be changed:

```
$> cd /usr/local/MATLAB/R2015a/sys/os/glnxa64
$> ls -l
```

The sym links for libstdc++.so.6 and libgfortran.so.3 should point to versions in /usr/lib, not local ones.


Before changing this libraries, first make sure `g++-4.4` and `libgfortran3`are installed :

```
$> sudo apt-get install g++-4.4 libgfortran3
```

Now, modify the symlinks:

```
$> sudo ln -fs /usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0 libgfortran.so.3
$> sudo ln -fs /usr/lib/gcc/x86_64-linux-gnu/4.4/libstdc++.so libstdc++.so.6
```

This command makes the `libstdc++.so.6` point to the `g++-4.4` `libstdc++` library.


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