/*
 * Copyright (c) 2019 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#include "orderk_delaunay.h"

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel                    K;

typedef OrderKDelaunay_3<K>::Point                                       Point;


TEST_CASE("Standard toy example", "[orderk_delaunay]") {

    Point p0(0, 0, 0);
    Point p1(0, 4, 4);
    Point p2(4, 4, 0);
    Point p3(4, 0, 4);
    Point p4(-10, 2, 2);

    std::vector<Point> points = {p0, p1, p2, p3, p4};

    std::vector<std::vector<std::vector<unsigned>>> o1del_expected = {
        {{0}, {1}, {2}, {3}},
        {{0}, {1}, {2}, {4}},
        {{0}, {1}, {3}, {4}}
    };

    std::vector<std::vector<std::vector<unsigned>>> o2del_expected = {
        {{0, 1}, {0, 2}, {0, 3}, {0, 4}},
        {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}},
        {{0, 1}, {0, 2}, {0, 4}, {1, 2}, {1, 4}, {2, 4}},
        {{0, 1}, {0, 3}, {0, 4}, {1, 3}, {1, 4}, {3, 4}},
        {{0, 1}, {1, 2}, {1, 3}, {1, 4}}
    };

    std::vector<std::vector<std::vector<unsigned>>> o3del_expected = {
        {{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {0, 2, 3}, {0, 2, 4}, {0, 3, 4}},
        {{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {1, 2, 3}, {1, 2, 4}, {1, 3, 4}},
        {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}},
        {{0, 1, 2}, {0, 1, 4}, {0, 2, 4}, {1, 2, 4}},
        {{0, 1, 3}, {0, 1, 4}, {0, 3, 4}, {1, 3, 4}}
    };

    std::vector<std::vector<std::vector<unsigned>>> o4del_expected = {
        {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}},
        {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {1, 2, 3, 4}}
    };

    auto orderkdelaunay = OrderKDelaunay_3<K>(points, 4);

    auto o1del = orderkdelaunay.get_canonical_representation(1);
    REQUIRE(o1del == o1del_expected);

    auto o2del = orderkdelaunay.get_canonical_representation(2);
    REQUIRE(o2del == o2del_expected);

    auto o3del = orderkdelaunay.get_canonical_representation(3);
    REQUIRE(o3del == o3del_expected);

    auto o4del = orderkdelaunay.get_canonical_representation(4);
    REQUIRE(o4del == o4del_expected);
}


TEST_CASE("Minimal example with a single cell", "[orderk_delaunay]") {

    Point p0(0, 0, 0);
    Point p1(0, 4, 4);
    Point p2(4, 4, 0);
    Point p3(4, 0, 4);

    std::vector<Point> points = {p0, p1, p2, p3};

    std::vector<std::vector<std::vector<unsigned>>> o1del_expected = {
        {{0,}, {1,}, {2,}, {3,}}
    };

    std::vector<std::vector<std::vector<unsigned>>> o2del_expected = {
        {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}
    };

    std::vector<std::vector<std::vector<unsigned>>> o3del_expected = {
        {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}
    };

    auto orderkdelaunay = OrderKDelaunay_3<K>(points, 3);

    auto o1del = orderkdelaunay.get_canonical_representation(1);
    REQUIRE(o1del == o1del_expected);

    auto o2del = orderkdelaunay.get_canonical_representation(2);
    REQUIRE(o2del == o2del_expected);

    auto o3del = orderkdelaunay.get_canonical_representation(3);
    REQUIRE(o3del == o3del_expected);
}


TEST_CASE("Non-convex cluster", "[orderk_delaunay]") {

    Point p0(0,0,3);
    Point p1(0,-1.5,-1.8);
    Point p2(-0.07,3.67,-2.03);
    Point p3(-2.37,3.08,2.49);
    Point p4(2.32,4.37,0.4);
    Point p5(0,-1.5,0);

    std::vector<Point> points = {p0, p1, p2, p3, p4, p5};

    // The cells
    // {{0, 5}, {1, 5}, {2, 5}, {3, 5}} and
    // {{0, 5}, {1, 5}, {2, 5}, {4, 5}}
    // form a non-convex cluster.
    std::vector<std::vector<std::vector<unsigned>>> o2del_expected = {
        {{0, 2}, {0, 3}, {0, 4}, {0, 5}},
        {{0, 2}, {0, 3}, {0, 4}, {2, 3}, {2, 4}, {3, 4}},
        {{0, 2}, {0, 3}, {0, 5}, {2, 3}, {2, 5}, {3, 5}},
        {{0, 2}, {0, 4}, {0, 5}, {2, 4}, {2, 5}, {4, 5}},
        {{0, 2}, {2, 3}, {2, 4}, {2, 5}},
        {{0, 4}, {1, 4}, {2, 4}, {4, 5}},
        {{0, 5}, {1, 5}, {2, 5}, {3, 5}},
        {{0, 5}, {1, 5}, {2, 5}, {4, 5}},
        {{1, 2}, {1, 3}, {1, 5}, {2, 3}, {2, 5}, {3, 5}},
        {{1, 2}, {1, 4}, {1, 5}, {2, 4}, {2, 5}, {4, 5}},
        {{1, 2}, {2, 3}, {2, 4}, {2, 5}}
    };

    auto orderkdelaunay = OrderKDelaunay_3<K>(points, 2);

    auto o2del = orderkdelaunay.get_canonical_representation(2);
    REQUIRE(o2del == o2del_expected);
}
