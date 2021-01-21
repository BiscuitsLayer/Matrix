#include <fstream>
#include <algorithm>

#include "MatrixLib.hpp"
#include "Solver/SolverLib.hpp"

int main () {
	Linear::Matrix <double> m1 {}, m2 {};
	std::cin >> m1 >> m2;
	
	//std::cout << m1 << std::endl;
	//std::cerr << m1.Rank () << std::endl;
	//std::cout << m2 << std::endl;
	//std::cerr << m2.Rank () << std::endl;
	
	Solver solver { m1, m2 };
	solver.Execute ();
	return 0;
}