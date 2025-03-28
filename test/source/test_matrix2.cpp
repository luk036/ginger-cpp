// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

// #include <vector>

// #include "fmt/format.h"        // for print
#include <ginger/vector2.hpp>  // for Vector2
#include <ginger/matrix2.hpp>  // for Matrix2

TEST_CASE("Matrix2 Construction") {
    SUBCASE("Value construction") {
        ginger::Vector2<int> row1(1, 2);
        ginger::Vector2<int> row2(3, 4);
        ginger::Matrix2<ginger::Vector2<int>> m(std::move(row1), std::move(row2));
        
        CHECK(m.x().x() == 1);
        CHECK(m.x().y() == 2);
        CHECK(m.y().x() == 3);
        CHECK(m.y().y() == 4);
    }
}

TEST_CASE("Matrix2 Operations") {
    ginger::Vector2<int> row1(1, 2);
    ginger::Vector2<int> row2(3, 4);
    ginger::Matrix2<ginger::Vector2<int>> m1(std::move(row1), std::move(row2));
    
    ginger::Vector2<int> row3(5, 6);
    ginger::Vector2<int> row4(7, 8);
    ginger::Matrix2<ginger::Vector2<int>> m2(std::move(row3), std::move(row4));

    SUBCASE("Negation") {
        auto m = -m1;
        CHECK(m.x().x() == -1);
        CHECK(m.x().y() == -2);
        CHECK(m.y().x() == -3);
        CHECK(m.y().y() == -4);
    }

    SUBCASE("Matrix addition") {
        auto m = m1 + m2;
        CHECK(m.x().x() == 6);
        CHECK(m.x().y() == 8);
        CHECK(m.y().x() == 10);
        CHECK(m.y().y() == 12);

        m1 += m2;
        CHECK(m1.x().x() == 6);
        CHECK(m1.x().y() == 8);
        CHECK(m1.y().x() == 10);
        CHECK(m1.y().y() == 12);
    }

    SUBCASE("Matrix subtraction") {
        auto m = m2 - m1;
        CHECK(m.x().x() == 4);
        CHECK(m.x().y() == 4);
        CHECK(m.y().x() == 4);
        CHECK(m.y().y() == 4);

        m2 -= m1;
        CHECK(m2.x().x() == 4);
        CHECK(m2.x().y() == 4);
        CHECK(m2.y().x() == 4);
        CHECK(m2.y().y() == 4);
    }

    SUBCASE("Scalar multiplication") {
        auto m = m1 * 2;
        CHECK(m.x().x() == 2);
        CHECK(m.x().y() == 4);
        CHECK(m.y().x() == 6);
        CHECK(m.y().y() == 8);

        m = 2 * m1;
        CHECK(m.x().x() == 2);
        CHECK(m.x().y() == 4);
        CHECK(m.y().x() == 6);
        CHECK(m.y().y() == 8);

        m1 *= 2;
        CHECK(m1.x().x() == 2);
        CHECK(m1.x().y() == 4);
        CHECK(m1.y().x() == 6);
        CHECK(m1.y().y() == 8);
    }

    SUBCASE("Scalar division") {
        ginger::Vector2<double> lrow1(2.0, 4.0);
        ginger::Vector2<double> lrow2(6.0, 8.0);
        ginger::Matrix2<ginger::Vector2<double>> m(std::move(lrow1), std::move(lrow2));

        auto lm2 = m / 2.0;
        CHECK(lm2.x().x() == doctest::Approx(1.0));
        CHECK(lm2.x().y() == doctest::Approx(2.0));
        
        CHECK(lm2.y().x() == doctest::Approx(3.0));
        CHECK(lm2.y().y() == doctest::Approx(4.0));

        m /= 2.0;
        CHECK(m.x().x() == doctest::Approx(1.0));
        CHECK(m.x().y() == doctest::Approx(2.0));
        CHECK(m.y().x() == doctest::Approx(3.0));
        CHECK(m.y().y() == doctest::Approx(4.0));
    }
}

TEST_CASE("Matrix2 Vector Operations") {
    ginger::Vector2<double> row1(1, 2);
    ginger::Vector2<double> row2(3, 4);
    ginger::Matrix2<ginger::Vector2<double>> m(std::move(row1), std::move(row2));

    ginger::Vector2<double> v(5, 6);

    SUBCASE("Matrix-vector multiplication") {
        auto result = m.mdot(v);
        CHECK(result.x() == 1*5 + 2*6);
        CHECK(result.y() == 3*5 + 4*6);
    }
}

TEST_CASE("Matrix2 Determinant") {
    SUBCASE("Integer matrix") {
        ginger::Vector2<int> row1(3, 4);
        ginger::Vector2<int> row2(5, 6);
        ginger::Matrix2<ginger::Vector2<int>> m(std::move(row1), std::move(row2));
        
        CHECK(m.det() == doctest::Approx(3*6 - 4*5));
    }

    SUBCASE("Double matrix") {
        ginger::Vector2<double> row1(1.5, 2.5);
        ginger::Vector2<double> row2(3.5, 4.5);
        ginger::Matrix2<ginger::Vector2<double>> m(std::move(row1), std::move(row2));
        
        CHECK(m.det() == doctest::Approx(1.5*4.5 - 2.5*3.5));
    }
}

// TEST_CASE("Matrix2 Mixed Type Operations") {
//     ginger::Vector2<int> row1(1, 2);
//     ginger::Vector2<int> row2(3, 4);
//     ginger::Matrix2<ginger::Vector2<int>> mi(std::move(row1), std::move(row2));

//     ginger::Vector2<double> row3(1.5, 2.5);
//     ginger::Vector2<double> row4(3.5, 4.5);
//     ginger::Matrix2<ginger::Vector2<double>> md(std::move(row3), std::move(row4));

//     SUBCASE("Mixed type addition") {
//         auto m = mi + md;
//         CHECK(m.x().x() == doctest::Approx(2.5));
//         CHECK(m.x().y() == doctest::Approx(4.5));
//         CHECK(m.y().x() == doctest::Approx(6.5));
//         CHECK(m.y().y() == doctest::Approx(8.5));
//     }

//     SUBCASE("Mixed type subtraction") {
//         auto m = md - mi;
//         CHECK(m.x().x() == doctest::Approx(0.5));
//         CHECK(m.x().y() == doctest::Approx(0.5));
//         CHECK(m.y().x() == doctest::Approx(0.5));
//         CHECK(m.y().y() == doctest::Approx(0.5));
//     }
// }