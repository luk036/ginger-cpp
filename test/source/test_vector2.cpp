// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <vector>

#include "fmt/format.h"        // for print
#include "ginger/vector2.hpp"  // for Vector2

TEST_CASE("test vector2") {
    auto h = std::vector<double>{1.0, 2.0, 3.0, 4.0};
    auto x = 1.0;
    auto y = 2.0;
    auto v = ginger::Vector2<double&, double&>(x, y);
    CHECK(v.x() == 1.0);
    CHECK(v.y() == 2.0);

    auto v2 = ginger::Vector2<double&, double&>(x, y);
    v2 *= 2.0;
    CHECK(v.y() == 4.0);

    // auto vA = horner(h, 3, v);
    // CHECK(vA.x() == 8.0);
    // CHECK(vA.y() == 10.0);
}


TEST_CASE("Vector2 Construction") {
    SUBCASE("Default construction") {
        ginger::Vector2<int> v;
        CHECK(v.x() == 0);
        CHECK(v.y() == 0);
    }

    SUBCASE("Value construction") {
        ginger::Vector2<int> v(5, 10);
        CHECK(v.x() == 5);
        CHECK(v.y() == 10);
    }

    SUBCASE("Copy construction") {
        ginger::Vector2<int> v1(5, 10);
        ginger::Vector2<double> v2(v1);
        CHECK(v2.x() == doctest::Approx(5.0));
        CHECK(v2.y() == doctest::Approx(10.0));
    }
}

TEST_CASE("Vector2 Operations") {
    ginger::Vector2<int> v1(3, 4);
    ginger::Vector2<int> v2(5, 6);
    ginger::Vector2<double> v3(1.5, 2.5);

    SUBCASE("Dot product") {
        CHECK(v1.dot(v2) == doctest::Approx(3*5 + 4*6));
        CHECK(v1.dot(v3) == doctest::Approx(3*1.5 + 4*2.5));
    }

    SUBCASE("Cross product") {
        CHECK(v1.cross(v2) == doctest::Approx(3*6 - 4*5));
        CHECK(v1.cross(v3) == doctest::Approx(3*2.5 - 1.5*4));
    }

    SUBCASE("Negation") {
        auto v = -v1;
        CHECK(v.x() == -3);
        CHECK(v.y() == -4);
    }
}

TEST_CASE("Vector2 Arithmetic Operators") {
    ginger::Vector2<int> v1(3, 4);
    ginger::Vector2<int> v2(5, 6);
    // ginger::Vector2<double> v3(1.5, 2.5);

    SUBCASE("Addition") {
        auto v = v1 + v2;
        CHECK(v.x() == 8);
        CHECK(v.y() == 10);

        v1 += v2;
        CHECK(v1.x() == 8);
        CHECK(v1.y() == 10);
    }

    SUBCASE("Subtraction") {
        auto v = v2 - v1;
        CHECK(v.x() == 2);
        CHECK(v.y() == 2);

        v2 -= v1;
        CHECK(v2.x() == 2);
        CHECK(v2.y() == 2);
    }

    SUBCASE("Scalar multiplication") {
        auto v = v1 * 2;
        CHECK(v.x() == 6);
        CHECK(v.y() == 8);

        v = 2 * v1;
        CHECK(v.x() == 6);
        CHECK(v.y() == 8);

        v1 *= 2;
        CHECK(v1.x() == 6);
        CHECK(v1.y() == 8);
    }

    SUBCASE("Scalar division") {
        ginger::Vector2<double> v(6.0, 9.0);
        auto v4 = v / 3.0;
        CHECK(v4.x() == doctest::Approx(2.0));
        CHECK(v4.y() == doctest::Approx(3.0));

        v /= 3.0;
        CHECK(v.x() == doctest::Approx(2.0));
        CHECK(v.y() == doctest::Approx(3.0));
    }
}

// TEST_CASE("Vector2 Mixed Type Operations") {
//     ginger::Vector2<int> vi(3, 4);
//     ginger::Vector2<double> vd(1.5, 2.5);

//     SUBCASE("Mixed type addition") {
//         auto v = vi + vd;
//         CHECK(v.x() == doctest::Approx(4.5));
//         CHECK(v.y() == doctest::Approx(6.5));
//     }

//     SUBCASE("Mixed type subtraction") {
//         auto v = vd - vi;
//         CHECK(v.x() == doctest::Approx(-1.5));
//         CHECK(v.y() == doctest::Approx(-1.5));
//     }

//     SUBCASE("Mixed type multiplication") {
//         auto v = vi * 1.5;
//         CHECK(v.x() == doctest::Approx(4.5));
//         CHECK(v.y() == doctest::Approx(6.0));
//     }
// }

// TEST_CASE("Vector2 Output Stream") {
//     ginger::Vector2<int> v(3, 4);
//     std::stringstream ss;
//     ss << v;
//     CHECK(ss.str() == "{3, 4}");
// }