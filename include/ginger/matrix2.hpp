#pragma once

#include <utility>  // import std::move

#include "vector2.hpp"

#if __cpp_constexpr >= 201304
#    define CONSTEXPR14 constexpr
#else
#    define CONSTEXPR14 inline
#endif

namespace ginger {
    /**
     * @brief Matrix2
     *
     * A template class representing a 2x2 matrix with column vectors.
     * The matrix is stored as two column vectors (x and y).
     *
     * @tparam T1 The type of the first column vector
     * @tparam T2 The type of the second column vector (defaults to T1)
     */
    template <typename T1, typename T2 = T1> class Matrix2 {
      private:
        T1 _x;
        T2 _y;

      public:
        /**
         * @brief Construct a new Matrix2 object
         *
         * @param[in] x The first column vector
         * @param[in] y The second column vector
         */
        constexpr Matrix2(T1&& x, T2&& y) noexcept : _x{std::move(x)}, _y{std::move(y)} {}

        /**
         * @brief Get the first column vector
         *
         * @return constexpr const T1& Reference to the first column vector
         */
        constexpr auto x() const -> const T1& { return this->_x; }

        /**
         * @brief Get the second column vector
         *
         * @return constexpr const T2& Reference to the second column vector
         */
        constexpr auto y() const -> const T2& { return this->_y; }

        // /**
        //  * @brief Construct a new Matrix2 object
        //  *
        //  * @param[in] x
        //  * @param[in] y
        //  */
        // constexpr Matrix2(const T1& x, const T2& y) : Vector2<T1, T2>{x, y} {}

        /** @name Arithmetic operators
         *  definie +, -, *, /, +=, -=, *=, /=, etc.
         */
        ///@{

        /**
         * @brief Negate the matrix
         *
         * @return Matrix2<T1, T2> A new matrix with negated column vectors
         */
        constexpr auto operator-() const -> Matrix2<T1, T2> { return {-this->x(), -this->y()}; }

        /**
         * @brief Add another matrix to this one
         *
         * @tparam U1 Type of the first column of the other matrix
         * @tparam U2 Type of the second column of the other matrix
         * @param[in] other The matrix to add
         * @return Matrix2<T1, T2>& Reference to this matrix
         */
        template <typename U1, typename U2>
        CONSTEXPR14 auto operator+=(const Matrix2<U1, U2>& other) -> Matrix2<T1, T2>& {
            this->_x += other.x();
            this->_y += other.y();
            return *this;
        }

        /**
         * @brief Subtract another matrix from this one
         *
         * @tparam U1 Type of the first column of the other matrix
         * @tparam U2 Type of the second column of the other matrix
         * @param[in] other The matrix to subtract
         * @return Matrix2<T1, T2>& Reference to this matrix
         */
        template <typename U1, typename U2>  //
        CONSTEXPR14 auto operator-=(const Matrix2<U1, U2>& other) -> Matrix2<T1, T2>& {
            this->_x -= other.x();
            this->_y -= other.y();
            return *this;
        }

        /**
         * @brief Multiply
         *
         * @tparam R
         * @param[in] alpha
         * @return Matrix2<T1, T2>&
         */
        template <typename R> CONSTEXPR14 auto operator*=(const R& alpha) -> Matrix2<T1, T2>& {
            this->_x *= alpha;
            this->_y *= alpha;
            return *this;
        }

        /**
         * @brief Divide matrix by a scalar
         *
         * @tparam R The scalar type
         * @param[in] alpha The scalar value
         * @return Matrix2<T1, T2>& Reference to this matrix
         */
        template <typename R> CONSTEXPR14 auto operator/=(const R& alpha) -> Matrix2<T1, T2>& {
            this->_x /= alpha;
            this->_y /= alpha;
            return *this;
        }

        /**
         * @brief Add two matrices
         *
         * @tparam U1 Type of the first column of the first matrix
         * @tparam U2 Type of the second column of the first matrix
         * @param[in] x The first matrix
         * @param[in] y The second matrix
         * @return Matrix2<T1, T2> The result matrix
         */
        template <typename U1, typename U2>  //
        friend constexpr auto operator+(Matrix2<T1, T2> x, const Matrix2<U1, U2>& y)
            -> Matrix2<T1, T2> {
            return std::move(x) += y;
        }

        /**
         * @brief Subtract two matrices
         *
         * @tparam U1 Type of the first column of the first matrix
         * @tparam U2 Type of the second column of the first matrix
         * @param[in] x The first matrix
         * @param[in] y The second matrix
         * @return Matrix2<T1, T2> The result matrix
         */
        template <typename U1, typename U2>  //
        friend constexpr auto operator-(Matrix2<T1, T2> x, const Matrix2<U1, U2>& y)
            -> Matrix2<T1, T2> {
            return std::move(x) -= y;
        }

        /**
         * @brief Multiply matrix by a scalar
         *
         * @tparam R The scalar type
         * @param[in] x The matrix
         * @param[in] alpha The scalar value
         * @return Matrix2<T1, T2> The result matrix
         */
        template <typename R> friend constexpr auto operator*(Matrix2<T1, T2> x, const R& alpha)
            -> Matrix2<T1, T2> {
            return x *= alpha;
        }

        /**
         * @brief Multiply a scalar by a matrix
         *
         * @tparam R The scalar type
         * @param[in] alpha The scalar value
         * @param[in] x The matrix
         * @return Matrix2<T1, T2> The result matrix
         */
        template <typename R> friend constexpr auto operator*(const R& alpha, Matrix2<T1, T2> x)
            -> Matrix2<T1, T2> {
            return x *= alpha;
        }

        /**
         * @brief Divide matrix by a scalar
         *
         * @tparam R The scalar type
         * @param[in] x The matrix
         * @param[in] alpha The scalar value
         * @return Matrix2<T1, T2> The result matrix
         */
        template <typename R> friend constexpr auto operator/(Matrix2<T1, T2> x, const R& alpha)
            -> Matrix2<T1, T2> {
            return x /= alpha;
        }

        /**
         * @brief Multiply matrix with vector
         *
         * Performs matrix-vector multiplication: this * other.
         *
         * @tparam U1 Type of the first component of the vector
         * @tparam U2 Type of the second component of the vector
         * @param[in] other The vector to multiply
         * @return T1 The resulting vector (as column vector type)
         */
        template <typename U1, typename U2>  //
        constexpr auto mdot(const Vector2<U1, U2>& other) const -> T1 {
            return T1{this->_x.dot(other), this->_y.dot(other)};
        }

        /**
         * @brief Calculate the determinant
         *
         * Computes the determinant of the 2x2 matrix.
         *
         * @return double The determinant value
         */
        constexpr auto det() const -> double {
            return this->x().x() * this->y().y() - this->x().y() * this->y().x();
        }

        ///@}
    };
}  // namespace ginger
