/*
 * Copyright (c) 2019-2020 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#ifndef _ORDERK_DELAUNAY_H_
#define _ORDERK_DELAUNAY_H_

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <vector>

// CGAL <4.10
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <inttypes.h>   // for uint32_t


// point index
typedef unsigned            PIndex;
// A combinatorial vertex in the order-k Delaunay is a k-tuple (k-set)
// of points, represented as a (sorted) k-tuple of point indices.
typedef std::vector<PIndex>  CVertex;
// A subvertex is only meaningful in conjunction with the cell it belongs
// to. The points of the subvertex, together with the anchor_vertex
// from the cell, make the actual combinatorial vertex.
typedef std::vector<PIndex>  CSubVertex;


// These data types are not strictly needed in the algorithm,
// but useful for the provided output formats:
// Vertex index.
typedef uint32_t            VIndex;
// Cell with indexed vertices
typedef std::vector<VIndex>       ICell;

static const int dimension = 3;


// Data structure representing a cell in the order-k Delaunay mosaic.
// Each cell is uniquely defined as a g-th generation slice of a rhomboid,
// which in turn is uniquely defined by X_in and X_on.
// This data structure stores X_in (as anchor_vertex) and all subsets of
// X_on of generation g (as subvertices). For each subvertex, its union with
// the anchor_vertex gives the combinatorial vertex in the order-k Delaunay
// mosaic.
struct Cell {
  // Indices of the points that all vertices in this cell share
  CVertex anchor_vertex; // This is X_in
  // For each vertex of the cell, the additional points to add to the common
  // points to obtain this vertex
  std::vector<CSubVertex> subvertices;  // The union of these is X_on
  // order in which this cell appears
  int k;
};


template<class K>
class OrderKDelaunay_3 {
  public:
    typedef typename K::Point_3                                                       Point;
    typedef typename K::Vector_3                                                     Vector;

    // CGAL 4.10 breaks some backward compatibility.
    // This code is for CGAL <4.10, but replacing the typedefs below
    // should make to code compile for later CGAL versions.

    // If using CGAL <4.10:
    typedef typename CGAL::Regular_triangulation_euclidean_traits_3<K>               Traits;
    typedef typename CGAL::Triangulation_vertex_base_3<Traits>                        Vbase;
    typedef typename CGAL::Triangulation_vertex_base_with_info_3<unsigned, Traits,Vbase> Vb;
    typedef typename CGAL::Regular_triangulation_cell_base_3<Traits>                     Cb;
    typedef typename CGAL::Triangulation_data_structure_3<Vb,Cb>                        Tds;
    typedef typename CGAL::Regular_triangulation_3<Traits, Tds>                     Reg_Tri;
    typedef typename Traits::Weighted_point_3                                Weighted_point;

    // If using CGAL >=4.10:
    // typedef typename K::FT                                               Weight;
    // typedef typename K::Weighted_point_3                                 Weighted_point;
    // typedef typename CGAL::Regular_triangulation_vertex_base_3<K>        Vb0;
    // typedef typename CGAL::Triangulation_vertex_base_with_info_3<unsigned, K, Vb0> Vb;
    // typedef typename CGAL::Regular_triangulation_cell_base_3<K>          Cb;
    // typedef typename CGAL::Triangulation_data_structure_3<Vb,Cb>         Tds;
    // typedef typename CGAL::Regular_triangulation_3<K, Tds>               Reg_Tri;


    /**
     * Order-k Delaunay diagrams up to a given order.
     *
     * Upon creation, this class computes the order-k Delaunay mosaics
     * of the given point set, up to the specified order.
     * 
     * Input:
     *    bpoints: vector of input Points.
     *    order: The order up to (including) which to compute the mosaics.
     */
    OrderKDelaunay_3(const std::vector<Point>& bpoints, int order);
    /* For the specified order k, return the set of combinatorial vertices.
     *
     * Each vertex is represented as a k-tuple of point indices.
     *
     * Its geometric location can be obtained as the barycenter of the
     * k points that the combinatorial vertex refers to.
     */
    std::vector<CVertex> get_vertices(int order);
    /* For the specified order k, return the set of cells.
     *
     * Each cell is represented as a Cell struct.
     */
    std::vector<Cell> get_cells(int order);
    /* For the specified order k, return the set of simplices.
     *
     * Each simplex is represented as a 4-tuple of indices
     * into the vertex vector.
     */
    std::vector<ICell> get_triangulated_cells(int order);
    /* For the specified order k, get a canonical representation.
     *
     * Each cell is represented in as a sorted tuple of combinatorial vertices,
     * with each combinatorial vertex being a sorted tuple of k indices into
     * the original vector of input points. The returned vector is sorted.
     *
     * Thus two order-k Delaunay mosaics are combinatorially equivalent if and
     * only if their canonical representations are equivalent. This is useful
     * for testing.
     */
    std::vector<std::vector<std::vector<PIndex>>>
        get_canonical_representation(int order);
  private:

    // Make set of combinatorial vertices of the first-order Delaunay mosaic.
    // Each first-order vertex is a singleton set, containing one point index.
    std::vector<CVertex> compute_first_order_vertices(const std::vector<Point>& bpoints);

    std::vector<Point> bpoints;
    std::vector<typename K::FT> squared_lengths;
    std::vector<std::vector<ICell> > diagrams_simplices;
    std::vector<std::vector<Cell> > diagrams_cells;
    std::vector<std::vector<CVertex> > diagrams_vertices;
};

#include "orderk_delaunay_impl.h"

#endif // _ORDERK_DELAUNAY_H_
