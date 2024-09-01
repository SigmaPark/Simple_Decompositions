#include "Decompositions.hpp"


using std::cout;
using std::endl;


static void Test_LDLT_decomposition()
{
	cout << "Test_LDLT_decomposition" << endl << endl;

	prac::Matrix<float, 3, 3> const A
	{	4, 12, -16
	,	12, 37, -43
	,	-16, -43, 98
	};
	
	cout << "A = " << endl << A << endl << endl;
	
	auto const [L, d] = prac::LDLT_Decomposition(A);

	cout 
	<<	"L = " << endl << L << endl << endl
	<<	"d^T = " << d.transpose() << endl << endl
	<<	"L * diag(d) * L^T = " << endl << L*d.diag()*L.transpose() << endl;
}


static void Test_LDLT_linear_equation_solver()
{
	cout << "Test_LDLT_linear_equation_solver" << endl << endl;

	prac::Matrix<float, 3, 3> const A
	{	4, 12, -16
	,	12, 37, -43
	,	-16, -43, 98
	};

	prac::Matrix<float, 3, 1> const x{-2, 4, -8},  b = A*x;
	
	cout 
	<<	"A = " << endl << A << endl << endl
	<<	"x^T = " << x.transpose() << endl
	<<	"b^T = (A*x)^T = " << b.transpose() << endl << endl;
	
	auto const x_sol = prac::Solve_Linear_Equation_by_LDLT(A, b);
	
	cout << "x_sol^T = " << x_sol.transpose() << endl;
}


int main()
{
	::Test_LDLT_decomposition();
	cout << "--------  --------  --------  --------  --------  " << endl;

	::Test_LDLT_linear_equation_solver();
	cout << "--------  --------  --------  --------  --------  " << endl;

	return 0;
}