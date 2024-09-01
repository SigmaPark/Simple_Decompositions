#include "Decompositions.hpp"

int main()
{
	prac::Matrix<float, 3, 3> const A
	{	4, 12, -16
	,	12, 37, -43
	,	-16, -43, 98
	};

	{
		auto const [L, d] = prac::LDLT_Decomposition(A);

		std::cout << L << std::endl;
		std::cout << d.transpose() << std::endl;

		std::cout << "L*D*L^T = " << std::endl << L*d.diag()*L.transpose() << std::endl;
	}
	{
		prac::Matrix<float, 3, 1> const 
			x{-2, 4, -8}, 
			b = A*x,
			x_sol = prac::Solve_Linear_Equation_by_LDLT(A, b);

		std::cout << "x_sol = " << x_sol.transpose() << std::endl;
	}

	return 0;
}