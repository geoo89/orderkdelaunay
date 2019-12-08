/*
 * Copyright (c) 2019 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#include "orderk_delaunay.h"


// Compute common intersection of a set of combinatorial vertices, which
// are sets of point indices.
CSubVertex intersection(const std::vector<CVertex> &vertices) {
    CSubVertex previous_intersection = vertices[0];
    CSubVertex new_intersection;

    for (unsigned i = 1; i < vertices.size(); ++i) {
        std::set_intersection(
            previous_intersection.begin(), previous_intersection.end(),
            vertices[i].begin(), vertices[i].end(),
            std::back_inserter(new_intersection));
        std::swap(previous_intersection, new_intersection);
        new_intersection.clear();
    }
    return previous_intersection;
}


// All subsets of {0, 1, 2, 3} of size 2 and 3 respectively.
// Used to get second and third generation cells from a first generation cell.
const std::vector<std::vector<int> >
    combinatorial_pairs {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
const std::vector<std::vector<int> >
    combinatorial_triplets {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};


// Make set of combinatorial vertices of the first-order Delaunay mosaic.
// Each first-order vertex is a singleton set, containing one point index.
std::vector<CVertex> compute_first_order_vertices(
    const std::vector<Point>& bpoints) {
  std::vector<CVertex> cvertices1;
  for (int i = 0; i < bpoints.size(); ++i) {
    CVertex cv;
    cv.push_back(i);
    cvertices1.push_back(cv);
  }

  return cvertices1;
}


/**
 * Constructor for Order-k Delaunay mosaics up to a given order.
 * 
 * Input:
 *    bpoints: vector of input Points.
 *    order: The order up to (excluding) which to compute the mosaics.
 */
// TODO: Refactor into smaller methods.
OrderKDelaunay_3::OrderKDelaunay_3(
    const std::vector<Point>& bpoints, int order) {
  for (auto const& p : bpoints) {
    squared_lengths.push_back(Vector(CGAL::ORIGIN, p).squared_length());
  }

  // Step 1: Compute order-1 Delaunay mosaics for k >= 2
  // Make combinatorial vertices of the first-order Delaunay mosaic.
  diagrams_vertices.push_back(compute_first_order_vertices(bpoints));

  // Turn each point into a weighted point with weight 0.
  std::vector< std::pair<Weighted_point,PIndex> > points;
  for (PIndex i = 0; i < bpoints.size(); ++i) {
    points.push_back(std::make_pair(Weighted_point(bpoints[i], 0), i));
  }

  // Compute first-order Delaunay triangulation as a regular triangulation
  // with weights 0.
  Reg_Tri T(points.begin(), points.end());

  // Turn the output into two vectors: of Cells and ICells.
  std::vector<Cell> cells;
  std::vector<ICell> icells;
  Reg_Tri::Finite_cells_iterator cit;
  for (cit = T.finite_cells_begin(); cit != T.finite_cells_end(); ++cit) {
    Cell cell;
    cell.k = 1;
    cell.anchor_vertex = std::vector<PIndex>();
    ICell icell;
    for (int i = 0; i <= dimension; ++i) {
      VIndex vindex = cit->vertex(i)->info();
      cell.subvertices.push_back(diagrams_vertices[0][vindex]);
      icell.push_back(vindex);
    }
    cells.push_back(cell);
    std::sort(icell.begin(), icell.end());
    icells.push_back(icell);
  }

  diagrams_cells.push_back(cells);
  diagrams_simplices.push_back(icells);

  // Queue of first-generation cells from the order-(k-1) and
  // order-(k-2) Delaunay mosaics from which we will obtain the
  // 2nd and 3rd generation cells of the order-k Delaunay mosaic.
  std::vector<Cell>& cell_queue = cells;

  // Step 2: Compute order-k Delaunay mosaics for k >= 2
  for (int k = 2; k <= order; ++k) {
    // In the process, we will make a new queue of first-generation cells,
    // for the next iteration of the loop.
    std::vector<Cell> cell_queue_new;

    // Vector to store the new cells. In step 2.1 these are the higher
    // generation cells, and in 2.2 we also add the first-generation cells.
    std::vector<Cell> new_cells;
    // Vector to store the vertex set of the order-k Delaunay mosaic.
    std::set<CVertex> new_vertices;

    /* Step 2.1: Compute the vertices and higher-generation cells of the
    order-k Delaunay mosaic. */
    for (auto const& cell : cell_queue) {
      // The cell queue contains first-generation vertices. Dependingon whether
      // the cell is from the order-(k-1) or order-(k-2) Delaunay mosaic, we
      // obtain a 2nd or 3rd generation cell for the order-k mosaic from it.
      int generation = k - cell.k + 1;
      // Create the higher-generation cell.
      // The new cell will have the same anchor vertex as the first-generation
      // cell it is derived from. Its subvertices will be pairs or triplets of
      // X_on of the first-generation cell it was derived from.
      Cell new_cell;
      new_cell.k = k;
      new_cell.anchor_vertex = cell.anchor_vertex;
      if (generation == 2) {
        // Take all 2-tuples of subvertices of the first-generation cell.
        // These will be the subvertices of the new cell.
        for (auto const& pair : combinatorial_pairs) {
          // cell is a first-gen cell, thus its subvertices are singletons.
          CSubVertex new_subvertex = {cell.subvertices[pair[0]][0], 
                                      cell.subvertices[pair[1]][0]};
          std::sort(new_subvertex.begin(), new_subvertex.end());
          new_cell.subvertices.push_back(new_subvertex);
          // Make the full new combinatorial vertex, and add it to the order-k
          // vertex set. Need to be sorted for uniqueness within the set.
          auto cp = std::vector<PIndex>(cell.anchor_vertex);
          cp.insert(cp.end(), new_subvertex.begin(), new_subvertex.end());
          std::sort(cp.begin(), cp.end());
          new_vertices.insert(cp);
        }

        // Put original first-generation cell back into the queue because we
        // will also derive a third-generation cell from it for the
        // order-(k+1) mosaic.
        cell_queue_new.push_back(cell);

      } else if (generation == 3) {
        // Take all 3-tuples of subvertices of the first-generation cell.
        // These will be the subvertices of the new cell.
        for (auto const& triplet : combinatorial_triplets) {
          // cell is a first-gen cell, thus its subvertices are singletons.
          CSubVertex new_subvertex = {cell.subvertices[triplet[0]][0], 
                                      cell.subvertices[triplet[1]][0], 
                                      cell.subvertices[triplet[2]][0]};
          std::sort(new_subvertex.begin(), new_subvertex.end());
          new_cell.subvertices.push_back(new_subvertex);
          // Make the full new combinatorial vertex, and add it to the order-k
          // vertex set. Need to be sorted for uniqueness within the set.
          auto cp = std::vector<PIndex>(cell.anchor_vertex);
          cp.insert(cp.end(), new_subvertex.begin(), new_subvertex.end());
          std::sort(cp.begin(), cp.end());
          new_vertices.insert(cp);
        }
      }
      new_cells.push_back(new_cell);
    }

    // Step 2.2: Compute the remaining first-generation cells of
    // the order-k Delaunay mosaic via a regular triangulation.

    // Store the new vertex set.
    std::vector<CVertex> new_vertices_vector(new_vertices.begin(), 
                                             new_vertices.end());
    diagrams_vertices.push_back(new_vertices_vector);

    // Step 2.2.1: Construct the geometric vertices and their weights.
    std::vector< std::pair<Weighted_point,PIndex> > new_points;
    int index = 0;
    for (auto const& cv: new_vertices) {
      // We compute the mean of the (lifted) points that are part of the vertex
      // mean: will be the mean of the 3D coordinates.
      Vector mean = CGAL::NULL_VECTOR;
      // mean_sq_length: will be the mean of the 4-th coordinate (i.e. "height")
      K::FT mean_sq_length = 0;
      for (auto const& pt: cv) {
        mean = mean + Vector(CGAL::ORIGIN, bpoints[pt]);
        mean_sq_length += squared_lengths[pt];
      }
      mean = mean / k;
      mean_sq_length = mean_sq_length / k;
      // mean.squared_length() is the "height" of the mean point in 3D
      // mean_sq_length is the "height" of the mean of the 4D (lifted) points
      // so the weight is how much the mean of the 4D points lies above the
      // paraboloid, because mean.squared_length() is the height of the
      // paraboloid at the mean point in 3D.
      auto weight = mean.squared_length() - mean_sq_length;

      new_points.push_back(std::make_pair(
          Weighted_point(Point(mean.x(), mean.y(), mean.z()), weight),
          index));
      ++index;
    }

    // Vector to store the triangulated cells as 4-tuples of vertex indices.
    std::vector<ICell> icells;

    // Step 2.2.2: Get weighted Delaunay triangulation and indentify its
    // first-generation cell.
    Reg_Tri T(new_points.begin(), new_points.end());
    Reg_Tri::Finite_cells_iterator cit;
    for (cit = T.finite_cells_begin(); cit != T.finite_cells_end(); ++cit) {
      // TODO: Use pointers here instead of copies?
      std::vector<CVertex> vertices;

      // Store all triangulated cells as k-tuples of vertex indices.
      ICell icell;
      for (int i = 0; i < 4; ++i) {
        VIndex vindex = cit->vertex(i)->info();
        vertices.push_back(new_vertices_vector[vindex]);
        icell.push_back(vindex);
      }
      std::sort(icell.begin(), icell.end());
      icells.push_back(icell);

      // Check whether the simplex is a first-generation cell.
      CSubVertex intersec = intersection(vertices);
      if (intersec.size() == k-1) {
        // If yes, turn it into a cell struct and store it.
        Cell cell;
        cell.k = k;
        cell.anchor_vertex = intersec;
        for (int i = 0; i < 4; ++i) {
          CSubVertex difference;
          std::set_difference(
              new_vertices_vector[cit->vertex(i)->info()].begin(), 
              new_vertices_vector[cit->vertex(i)->info()].end(),
              intersec.begin(), intersec.end(),
              std::inserter(difference, difference.end()));
          cell.subvertices.push_back(difference);
        }
        new_cells.push_back(cell);

        // Add it to cell queue to get 2nd and 3rd generation cell later.
        cell_queue_new.push_back(cell);
      }
    }

    // Add all the combinatorial cells and simplices to the output.
    diagrams_cells.push_back(new_cells);
    diagrams_simplices.push_back(icells);

    // The cell queue has been fully processed and a new one has been created.
    // For computing higher-generation cells for the next order,
    // we use the new cell queue.
    cell_queue = cell_queue_new;

  }

}


std::vector<CVertex> OrderKDelaunay_3::get_vertices(int order) {
  return diagrams_vertices[order-1];
}


// Each cell is a Cell struct.
std::vector<Cell> OrderKDelaunay_3::get_cells(int order) {
  return diagrams_cells[order-1];
}


// Note: the triangulation is not always unique.
// Each cell is a set of VIndices.
std::vector<ICell> OrderKDelaunay_3::get_triangulated_cells(int order) {
  return diagrams_simplices[order-1];
}


// The canonical representation is a list of cells, where
// each cell is a list of combinatorial vertices, with each
// combinatorial vertex being a list of indices into the original
// point set. Each vertex is sorted in ascending order, the vertices in
// each cell are sorted lexicographically, and the whole
// list of cells is sorted lexicographically. Thus the canonical
// representation is unique, and can be used to compare for equality
// when testing the output.
std::vector<std::vector<std::vector<unsigned>>>
    OrderKDelaunay_3::get_canonical_representation(int order) {
  std::vector<std::vector<std::vector<unsigned>>> cells;
  for (const auto& cell : get_cells(order)) {
    std::vector<std::vector<unsigned>> cell_vertices;
    // make a copy of the anchor vertex
    std::vector<unsigned> anchor_vertex(cell.anchor_vertex.begin(),
                                        cell.anchor_vertex.end());
    for (const auto& subvertex : cell.subvertices) {
      // each vertex shares the points of the anchor vertex, so we make a copy
      std::vector<unsigned> full_vertex(anchor_vertex);
      // and add the remaining point indices
      full_vertex.insert(full_vertex.end(),
                         subvertex.begin(), subvertex.end());
      std::sort(full_vertex.begin(), full_vertex.end());
      cell_vertices.push_back(full_vertex);
    }
    std::sort(cell_vertices.begin(), cell_vertices.end());
    cells.push_back(cell_vertices);
  }
  std::sort(cells.begin(), cells.end());
  return cells;
}
