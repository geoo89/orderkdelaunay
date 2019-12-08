# Order-k Delaunay Mosaics

Algorithm to compute higher order Delaunay mosaics in C++ and python.
The python version is dimension-agnostic but uses floating point arithmetics
and might be susceptible to floating point errors. The C++ version is
specific to the 3-dimensional setting, uses the [CGAL](https://www.cgal.org/)
library for Delaunay triangulations and exact arithmetics and includes
some unit tests using [Catch2](https://github.com/catchorg/Catch2).

## Python version

_Prerequisites:_ scipy, numpy, matplotlib

The python version also provides functionality to plot 2- and 3-
dimensional order-k Delaunay mosaics. The input point set to compute
the order-k Delaunay mosaics of can be of any dimension.
main.py computes and plots order-k Delaunay mosaics for a given
example point set.

## C++ version

_Prerequisites:_ cmake, CGAL version <= 4.9, Catch2 (included)

To work with CGAL version >= 4.10, some typedefs need to be changed,
see cpp/src/orderk_delaunay.h

To build everything and run the commandline tool:
```
cmake .
make
orderk example_data/example_input.txt
```
It will build a commandline tool that accepts an input filename,
output filename and an order up to which to compute the mosaics.
Input and output are text files. The input is one point per line,
with coordinates separated by spaces. The output is a canonical
representation, see cpp/src/orderk_delaunay.h

The unit tests use Catch2 which is included as a header and comes
with its own licence.