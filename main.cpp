#include <fstream>
#include <algorithm>

#include "MatrixLib.hpp"
#include "Solver/SolverLib.hpp"
#include "Circuit/Circuit.hpp"

int main () {
    Linear::Matrix <double> m1 {}, m2 {};
	std::cin >> m1 >> m2;
	Solver solver { m1, m2 };
	solver.Execute ();
}