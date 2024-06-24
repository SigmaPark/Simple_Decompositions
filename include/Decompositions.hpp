#pragma once
#ifndef _PRAC_DECOMPOSITIONS_
#define _PRAC_DECOMPOSITIONS_


#include <iostream>


namespace prac
{

	template<class T, int R, int C>
	struct Matrix;

	template<class T, int R>
	static auto LDLT_Decomposition(Matrix<T, R, R> const& A) noexcept;

}


template<class T, int R, int C>
struct prac::Matrix
{
	T data[R*C];

	auto operator()(int const i, int const j) noexcept-> T&
	{
		return data[i + j*R];
	}

	auto operator()(int const i, int const j) const noexcept-> T const&
	{
		return data[i + j*R];
	}

	template<int C2>
	auto operator*(Matrix<T, C, C2> const& m) const noexcept-> Matrix<T, R, C2>
	{
		Matrix<T, R, C2> res{0, };

		for(int i = 0;  i < R;  ++i)
			for(int j = 0;  j < C2;  ++j)
			{
				auto& e = res(i, j);

				for(int k = 0;  k < C;  ++k)
					e += (*this)(i, k)*m(k, j);
			}

		return res;
	}

	auto transpose() const noexcept-> Matrix<T, C, R>
	{
		Matrix<T, C, R> res;

		for(int i = 0;  i < R;  ++i)
			for(int j = 0;  j < C;  ++j)
				res(j, i) = (*this)(i, j);

		return res;
	}

	auto diag() const noexcept
	{
		if constexpr(R == 1 || C == 1)
		{
			Matrix<T, R*C, R*C> res{0, };

			for(int i = 0;  i < R*C;  ++i)
				res(i, i) = this->data[i];

			return res;
		}
		else if constexpr(R == C && R > 1)
		{
			Matrix<T, R, 1> res{0, };

			for(int i = 0;  i < R;  ++i)
				res(i, 0) = (*this)(i, i);

			return res;
		}
	}
};


template<class T, int R, int C>
static decltype(auto) operator<<(std::ostream& os, prac::Matrix<T, R, C> const& m)
{
	for(int i = 0;  i < R;  ++i,  os << std::endl)
		for(int j = 0;  j < C;  ++j,  os << ", ")
			os << m(i, j);

	return os;
}
//========//========//========//========//=======#//========//========//========//========//=======#


template<class T, int R>
auto prac::LDLT_Decomposition(Matrix<T, R, R> const& A) noexcept
{
	struct res_t
	{
		Matrix<T, R, R> L{0, };
		Matrix<T, R, 1> d{0, };
	}	res;

	auto& L = res.L;
	auto& d = res.d;

	auto LLD_f
	=	[&L, &d](int const i, int const j) noexcept-> T
		{
			T res = 0;

			for(int k = 0;  k < j;  ++k)
				res += L(i, k)*L(j, k)*d(k, 0);

			return res;
		};

	for(int i = 0;  i < R;  ++i)
	{
		for(int j = 0;  j < i;  ++j)
			L(i, j) = ( A(i, j) - LLD_f(i, j) ) / d(j, 0);
		
		L(i, i) = 1;
		d(i, 0) = A(i, i) - LLD_f(i, i);
	}

	return res;
}


#endif // end of #ifndef _PRAC_DECOMPOSITIONS_