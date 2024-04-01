#include "MatrixUtils.hpp"
#include "StaticMatrix.hpp"
#include <cassert>

int main()
{
    constexpr StaticMatrix<double> m1{
        {1, 2, 3,
         4, 5, 6,
         7, 8, 9}};
    constexpr StaticMatrix<double> m2{{1, 0, 0,
                                       0, 1, 0,
                                       0, 0, 1}};

    constexpr StaticMatrix<double> m3{
        {2, 2, 3,
         4, 6, 6,
         7, 8, 10}};

    constexpr StaticMatrix<double> m4{
        {2, 4, 7,
         2, 6, 8,
         3, 6, 10}};

    constexpr auto identity3x3 = StaticMatrix<double, 3, 3>::identity();
    static_assert(identity3x3 * m3 == m3);

    static_assert(m1 * m2 == m1);
    static_assert(m1 + m2 == m3);

    constexpr auto zeros = StaticMatrix<double, 3, 3>::zeros();
    static_assert(zeros * m3 == zeros);

    static_assert(m3.getTransposed() == m4);

    constexpr StaticMatrix<double> A{
        {2, 4, 7,
         2, 6, 8,
         3, 6, 10}};

    constexpr auto PLU = DecompPALU(A);
    constexpr auto P = PLU.P;
    constexpr auto L = PLU.L;
    constexpr auto U = PLU.U;

    static_assert(P*A == L*U);

    constexpr StaticMatrix<double, 3, 1> v{
        {2, 4, 7}};

    static_assert((StaticMatrix<double, 3, 3>::identity() * v) == v);

    static_assert((StaticMatrix<double, 3, 3>::identity() * v) == v);

    static_assert(EqualWithTolerance(CalcInvMatriz(m4)*m4, m4.identity()));

    static_assert(StaticMatrix<double, 2, 2>::identity()*5 == StaticMatrix<double, 2, 2>{{5.0, 0.0, 0.0, 5.0 }});

    showMatrix(m3.getTransposed());
    showMatrix(m4);

    return 0;
}