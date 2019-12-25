# Order-k Delaunay Mosaics

Algorithm to compute higher order Delaunay mosaics in C++ and python,
based on the paper [A Simple Algorithm for Computing Higher Order Delaunay
Mosaics](http://pub.ist.ac.at/~edels/Papers/2020-P-01-SimpleAlgorithm.pdf).
The python version is dimension-agnostic but uses floating point arithmetics
and might be susceptible to floating point errors. The C++ version is
specific to the 3-dimensional setting, uses the [CGAL](https://www.cgal.org/)
library for Delaunay triangulations and exact arithmetics and includes
some unit tests using [Catch2](https://github.com/catchorg/Catch2).

## Python version

_Prerequisites:_ scipy, numpy, (matplotlib)

Order-k Delaunay mosaics are implemented by the OrderKDelaunay class in
`python/orderk_delaunay.py`, see documentation there. The input point set to
compute the order-k Delaunay mosaics of can be of any dimension.
A plotter class is provided for plotting 2- and 3-dimensional order-k
Delaunay mosaics. As an example of the usage, `python/main.py` computes
and plots order-k Delaunay mosaics for a given example point set.

### Persistence of k-fold covers

The functionality to compute persistence of k-fold covers in 2 and 3
dimensions is implemented in the `kcover_persistence` function in
`python/kcover_persistence.py`, see documentation there. 
As an example of the usage, `python/main_kcoverp.py` computes persistence
of the k-fold cover for a given example point set and plots the resulting
persistence diagram.

_Remark:_ The only aspects that are implemented in a dimension-dependent way
(limited to 2D and 3D) are computing circumspheres of points and obtaining
lower-dimensional simplices from the top-dimensional simplices of a simplicial
complex. In principle neither should be difficult to generalize, allowing the
code to be made dimension-agnostic.

## C++ version

_Prerequisites:_ cmake, CGAL version <= 4.9, Catch2 (included);
to work with CGAL version >= 4.10, some typedefs need to be changed,
see `cpp/src/orderk_delaunay.h`

Order-k Delaunay mosaics are implemented by the OrderKDelaunay_3 class
in `cpp/src/orderk_delaunay.h`, see documentation there.

The build setup builds a commandline tool and tests. To build, run:
```
cmake .
make
```
The commandline tool accepts an input filename,
output filename and an order up to which to compute the mosaics.
Example usage:
```
./orderk example_data/example_input.txt example_data/example_output.txt 4
```
Input and output are text files. The input is one point per line,
each with 3 coordinates separated by spaces. The output file contains canonical
representations for the mosaics of each order, one per line, preceded by
an integer representing the order. Canonical representations are explained
in the documentation for the method `get_canonical_representation(int order)`
in `cpp/src/orderk_delaunay.h`.

The unit tests use Catch2 which is included as a header and comes
with its own licence.