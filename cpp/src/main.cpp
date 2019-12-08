/*
 * Copyright (c) 2019 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#include "orderk_delaunay.h"

#include <iostream>
#include <fstream>
#include <string>


/*
 * Read points from a file and return a vector of points.
 */
std::vector<Point> read_points(std::istream &stream) {
    //file should consists only of points in R^3
    std::vector<Point> points;
    for (std::string s; std::getline(stream, s);) {
        if (s.size() == 0 || s[0] == '#')
            continue;

        std::stringstream ss(s);

        std::vector<double> p;

        double v;
        for (int i = 0; i < dimension; i++) {
            ss >> v;
            p.push_back(v);
        }
        points.push_back(Point(p[0], p[1], p[2]));
    }
    return points;
}


/* 
 * Give nice string representation of vectors when streaming to output. 
 */
template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  if (!v.empty()) {
    out << '[';
    for (unsigned i = 0; i < v.size() - 1; ++i) {
      out << v[i] << ", ";
    }
    out << v[v.size() - 1] << "]";
  } else {
    out << "[]";
  }

  return out;
}


/*
 * Commandline tool for order-k Delaunay mosaics.
 */
int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cout << "Usage: order infile outfile [order]" << std::endl;
        std::cout << "infile: Text file with 3 space separated coordinates per line." << std::endl;
        std::cout << "outfile: Output filename." << std::endl;
        std::cout << "order: Order k up to which to compute order-k Delaunay mosaics." << std::endl;
        exit(0);
    }

    // TODO: Proper error handling.
    std::string infile = std::string(argv[1]);
    std::string outfile = std::string(argv[2]);
    int max_order = std::atoi(argv[3]);

    std::ifstream pfile(infile.c_str());
    auto points = read_points(pfile);
    pfile.close();
    std::cout << "Points loaded (" << infile << ").\n" << std::endl;

    std::ofstream ofile(outfile.c_str());

    auto orderkdelaunay = OrderKDelaunay_3(points, max_order+1);
    for (int order = 1; order <= max_order; ++order) {

        auto vertices = orderkdelaunay.get_vertices(order);
        auto cells = orderkdelaunay.get_triangulated_cells(order);
        auto canon = orderkdelaunay.get_canonical_representation(order);

        std::cout << "Order " << order << ": " 
                  << vertices.size() << " vertices, " 
                  << cells.size() << " 3-dimensional cells." << std::endl;

        ofile << order << " " << canon << std::endl;
    }

    ofile.close();
    std::cout << "DONE." << std::endl;
}
