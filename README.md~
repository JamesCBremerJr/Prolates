Numerical algorithms for the associated Legendre functions
==========================================================

This repository contains code for the rapid numerical evaluation of the
associated Legendre functions of the first and second kinds on the
interval (-1,1).  It runs in time independent of degree and order.
The algorithm is described in the preprint 

    James Bremer, "An algorithm for the numerical evaluation of the  associated 
       Legendre functions that runs in time independent of degree and order."  
       arXiv:1707.03287

It also contains code for rapidly calculating the roots of the associated Legendre
functions on the interval (-1,1).  The evaluation of a single root is done in time 
independent of degree and order.  A forthcoming article will explain its operation.

The file alegendre_eval.f90 contains both codes.  For the most part,
it uses precomputed expansions which are found in the files alegendre_data.bin1
and alegendre_data.bin2 to perform its computations.  The precomputed expansions
take up about 138 MB (for degrees through 1,000,000).  The data was spread
across two files because of github's limitations on binary file szes.

The file bessel_eval.f90 contains code for the rapid evaluation of Bessel
functions of orders between 0 and 1,000,000,000 on the real line.  It is included
because one of the asymptotic expansion in alegendre_eval involves Bessel
functions and relies on this code.

The file alegendre_test0.f90 contains an extremely simple test code.  It can
be compiled and executed as follows:

  gfortran -o alegendre_test -Ofast bessel_eval.f90 alegendre_eval.f90 alegendre_test.f90
  ./alegendre_test

The code this package is released under version 3 of the GPL.  See the file
LICENSE for details.

