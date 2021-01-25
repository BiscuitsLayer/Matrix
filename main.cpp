#include <fstream>
#include <algorithm>

#include "MatrixLib.hpp"
#include "Solver/SolverLib.hpp"
#include "Circuit/Circuit.hpp"

int main () {
	//Linear::Matrix <double> m1 {}, m2 {};
	//std::cin >> m1 >> m2;
	
	//std::cout << m1 << std::endl;
	//std::cerr << m1.Rank () << std::endl;
	//std::cout << m2 << std::endl;
	//std::cerr << m2.Rank () << std::endl;
	
	//Solver solver { m1, m2 };
	//solver.Execute ();

	/*
	std::vector <std::vector <RV> > adjTable = {
		{ { -1, -1 }, { 2, 0 }, { -1, -1 }, { 5, 0 } },
		{ { 2, 0 }, { -1, -1 }, { 2, 0 }, { -1, -1 } },
		{ { -1, -1 }, { 2, 0 }, { -1, -1 }, { 1, 10 } },
		{ { 5, 0 }, { -1, -1 }, { 1, 10 }, { -1, -1 } },
	};
	*/

	std::vector <std::vector <RV> > adjTable = {
		{ { -1, -1 }, { 2, 0 }, { -1, -1 }, { 5, 0 }, { -1, -1 }, { 2, 0 } },
		{ { 2, 0 }, { -1, -1 }, { 2, 0 }, { -1, -1 }, { -1, -1 }, { -1, -1 } },
		{ { -1, -1 }, { 2, 0 }, { -1, -1 }, { 1, 10 }, { -1, -1 }, { -1, -1 } },
		{ { 5, 0 }, { -1, -1 }, { 1, -10 }, { -1, -1 }, { 1, -10 }, { -1, -1 } },
		{ { -1, -1 }, { -1, -1 }, { -1, -1 }, { 1, 10 }, { -1, -1 }, { 2, 0 } },
		{ { 2, 0 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { 2, 0 }, { -1, -1 } },
	};
	
	Circuit circuit { adjTable };
	circuit.GetCycles ();
	std::cout << "Second law:" << std::endl;
	circuit.SecondKhLaw ();
	std::cout << "First law:" << std::endl;
	circuit.FirstKhLaw ();
	return 0;
}