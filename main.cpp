#include <fstream>
#include <algorithm>

#include "MatrixLib.hpp"
#include "Solver/SolverLib.hpp"
#include "Circuit/Circuit.hpp"

int main () {
    std::vector <std::vector <RV> > adjTable = {
		{ { -1, -1 }, { 2, 0 }, { -1, -1 }, { 5, 0 }, { -1, -1 }, { 2, 0 } },
		{ { 2, 0 }, { -1, -1 }, { 2, 0 }, { -1, -1 }, { -1, -1 }, { -1, -1 } },
		{ { -1, -1 }, { 2, 0 }, { -1, -1 }, { 1, 10 }, { -1, -1 }, { -1, -1 } },
		{ { 5, 0 }, { -1, -1 }, { 1, -10 }, { -1, -1 }, { 1, -10 }, { -1, -1 } },
		{ { -1, -1 }, { -1, -1 }, { -1, -1 }, { 1, 10 }, { -1, -1 }, { 2, 0 } },
		{ { 2, 0 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { 2, 0 }, { -1, -1 } },
	};
	
	Circuit circuit { adjTable };
	std::cout << "Second law:" << std::endl;
	circuit.SecondKhLaw ();
	std::cout << "First law:" << std::endl;
	circuit.FirstKhLaw ();
	return 0;
}