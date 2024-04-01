#pragma once
#include <iostream>
#include <limits>
#include "StaticMatrix.hpp"

template<typename T>
constexpr T absVal(const T v)
{
    if(v > T{0})
    {
        return v;
    }

    return -v;
}

template <typename T, size_t nlines, size_t ncols>
void showMatrix(const StaticMatrix<T, nlines, ncols> &M)
{
    std::cout << "= {\n";
    for (size_t i = 0; i < M.nlines_; i++)
    {
        std::cout << "[";

        for (size_t j = 0; j < M.ncols_; j++)
        {
            std::cout << M.getVal(i, j) << ",\t ";
        }

        std::cout << "]," << std::endl;
    }

    std::cout << "};" << std::endl;
}

template<typename T>
constexpr T TrocaLinha(const T& matriz, size_t l1, size_t l2)
{
    T res = matriz;
	for (size_t j = 0; j < matriz.ncols_; j++)
	{
		typename T::InternalType tmp = res.data[j + (l2*res.ncols_)];
		res.data[j + (l2*res.ncols_)] = res.data[j + (l1*res.ncols_)];
		res.data[j + (l1*res.ncols_)] = tmp;
	}
    return res;
}

template<typename T>
constexpr size_t AchaIndicePivo(const T& m, size_t nColuna, size_t nLinInicial)
{
	size_t res = 0;
	typename T::InternalType maior = {0.0};
	for (size_t i = nLinInicial; i < m.nlines_; ++i)
	{
		const auto val = absVal(m.data[nColuna + (i * m.ncols_)]);

		if (val > maior)
		{
			res = i;
			maior = val;
		}
	}

	return res;
}

template<typename L_T, typename T_b>
constexpr T_b SubstSucessivas(const L_T& L, const T_b& b)
{
	T_b res;
	
	res.data[0] = b.data[0] / L.data[0 + (0 * L.ncols_)];

	for (size_t index = 1; index < res.data.size(); index++)
	{
		typename T_b::InternalType soma = 0;
		for (size_t j = 0; j < index; j++)
		{
			soma += L.data[j + (index * L.ncols_)] * res.data[j];
		}
		res.data[index] = (b.data[index] - soma) / L.data[index + (index * L.ncols_)];
	}

	return res;
}

template<typename U_T, typename T_b>
constexpr T_b SubstRetroativas(const U_T& U, const T_b& b)
{
	T_b res;
	res.data[b.data.size() - 1] = b.data.back() / U.data[(U.ncols_ - 1) + ((U.nlines_ - 1) * U.ncols_)];

	for (size_t index = res.data.size() - 2; index > 0; index--)
	{
		double soma = 0.0;
		for (size_t j = U.ncols_ - 1; j > index; j--)
		{
			soma += U.data[j + (index * U.ncols_)] * res.data[j];
		}
		res.data[index] = (b.data[index] - soma) / U.data[index + (index * U.ncols_)];
	}

	double soma = 0.0;
	for (size_t j = U.ncols_ - 1; j > 0; j--)
	{
		soma += U.data[j + (0 * U.ncols_)] * res.data[j];
	}
	res.data[0] = (b.data[0] - soma) / U.data[0 + (0 * U.ncols_)];

	return res;
}

template<typename P_T, typename L_T, typename U_T>
struct DecompPALUResult
{
    P_T P;
    L_T L;
    U_T U;
};

template<typename T, size_t nlines, size_t ncols> 
constexpr DecompPALUResult<StaticMatrix<T, ncols, ncols>, StaticMatrix<T, ncols, ncols>, StaticMatrix<T, nlines, ncols>> DecompPALU(const StaticMatrix<T, nlines, ncols>& A) noexcept
{
	/*
	L = {[1 0 0 .. 0],
		 [M 1 0 .. 0],
		 [M M 1 .. 0],
		 [M M M .. 1]};

	U = {[ Linha Piv� 1],
		 [ Linha Piv� 2],
		 [ Linha Piv� n],};

	P = {[ 1 na coluna = linha piv� 1],
		 [ 1 na coluna = linha piv� 2],
		 [ 1 na coluna = linha piv� n]};
	*/

	auto ML = StaticMatrix<T, ncols, ncols>::zeros();
	auto MP = StaticMatrix<T, ncols, ncols>::identity();
	auto MU = A;

	auto& L = ML.data;
	auto& U = MU.data;

	for (size_t k = 0; k < A.ncols_; k++)
	{
		const size_t ipiv = AchaIndicePivo(MU, k, k);
		MU = TrocaLinha(MU, k, ipiv);
		MP = TrocaLinha(MP, k, ipiv);
		ML = TrocaLinha(ML, k, ipiv);

		for (size_t j = k + 1; j < A.ncols_; j++)
		{
			L[k + (j * ML.ncols_)] = U[k + (j * MU.ncols_)] / U[k + (k * MU.ncols_)];
			for (size_t index = k; index < A.ncols_; index++)
			{
				U[index + (j * MU.ncols_)] -= L[k + (j * MU.ncols_)] * U[index + (k * MU.ncols_)];
			}
		}
	}

	ML += StaticMatrix<T, ncols, ncols>::ThisMatrixType::identity();

    return DecompPALUResult{.P=MP, .L=ML, .U=MU};
}

template<typename T, size_t nlines, size_t ncols>
constexpr StaticMatrix<T, nlines, ncols> CalcInvMatriz(const StaticMatrix<T, nlines, ncols>& A)
{
	// 49 -> frederico
	const auto [P, L, U] = DecompPALU(A);

    StaticMatrix<T, nlines, ncols> res;

	for(size_t j = 0; j < A.ncols_; j++)
	{
		auto v = StaticMatrix<T, nlines, 1>::zeros();

		v.data[j] = 1.0;

		v = P*v;
		auto K = SubstSucessivas(L, v);

		// K = UX

		auto x = SubstRetroativas(U, K);

		for(size_t i = 0; i < A.nlines_; i++)
		{
			res.data[j + (i * A.ncols_)] = x.data[i];
		}
	}

    return res;
}

template<typename T>
constexpr bool EqualWithTolerance(const T& A, const T& B) noexcept
{
    using InternalType = typename T::InternalType;
    constexpr auto min_val = std::numeric_limits<InternalType>::epsilon();

    for(size_t i = 0; i < A.data.size(); i++)
    {
        if(absVal(A.data[i] - B.data[i]) > min_val* InternalType{100})
        {
            return false;
        }
    }

    return true;
}
