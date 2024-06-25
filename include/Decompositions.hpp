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

	template<class T, int R>
	static auto Solve_Linear_Equation_by_LDLT(Matrix<T, R, R> const& A, Matrix<T, R, 1> const& b)
	noexcept-> Matrix<T, R, 1>;

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
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


namespace prac::_LDLT_solver
{
	
	template<class T, int R>
	static auto Forward_solve(Matrix<T, R, R> const& L, Matrix<T, R, 1> const& b) noexcept
	->	Matrix<T, R, 1>;

	template<class T, int R>
	static auto Middle_solve(Matrix<T, R, 1> const& d, Matrix<T, R, 1> const& y) noexcept
	->	Matrix<T, R, 1>;

	template<class T, int R>
	static auto Backward_solve(Matrix<T, R, R> const& LT, Matrix<T, R, 1> const& z) noexcept
	->	Matrix<T, R, 1>;

}


template<class T, int R>
auto prac::_LDLT_solver::Forward_solve(Matrix<T, R, R> const& L, Matrix<T, R, 1> const& b) 
noexcept-> Matrix<T, R, 1>
{
	Matrix<T, R, 1> y;

	auto Ly_f
	=	[&L, &y](int const i) noexcept-> T
		{
			T res = 0;

			for(int k = 0;  k < i;  ++k)
				res += L(i, k)*y(k, 0);

			return res;
		};

	for(int i = 0;  i < R;  ++i)
		y(i, 0) = b(i, 0) - Ly_f(i);
	
	return y;
}


template<class T, int R>
auto prac::_LDLT_solver::Middle_solve(Matrix<T, R, 1> const& d, Matrix<T, R, 1> const& y) 
noexcept-> Matrix<T, R, 1>
{
	auto z = y;

	for(int i = 0;  i < R;  ++i)
		z(i, 0) /= d(i, 0);

	return z;
}


template<class T, int R>
auto prac::_LDLT_solver::Backward_solve(Matrix<T, R, R> const& LT, Matrix<T, R, 1> const& z) 
noexcept-> Matrix<T, R, 1>
{
	Matrix<T, R, 1> x;

	auto LTx_f
	=	[&LT, &x](int const r) noexcept-> T
		{
			T res = 0;

			for(int k = R - 1;  k > r;  --k)
				res += LT(r, k)*x(k, 0);

			return res;
		};

	for(int r = R - 1;  r >= 0;  --r)
		x(r, 0) = z(r, 0) - LTx_f(r);

	return x;
}
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


template<class T, int R>
auto prac::Solve_Linear_Equation_by_LDLT(Matrix<T, R, R> const& A, Matrix<T, R, 1> const& b)
noexcept-> Matrix<T, R, 1>
{
	using namespace _LDLT_solver;

	auto const [L, d] = LDLT_Decomposition(A);

	return
	Backward_solve
	(	L.transpose()
	,	Middle_solve
		(	d
		,	Forward_solve(L, b)
		)
	);
}


#endif // end of #ifndef _PRAC_DECOMPOSITIONS_