#include <fstream>
#include <algorithm>

#include "Reader/Language/driver.hpp"
#include "Matrix/Matrix.hpp"
#include "Solver/Solver.hpp"
#include "Circuit/Circuit.hpp"

int main (int argc, char** argv) {
	std::ifstream infile { argv[1] };
	if (!infile) {
		ERRSTREAM << "Error opening file!" << std::endl;
		return 0;
	}
	
	yy::LangDriver driver { infile };
	if (driver.parse ()) {
		driver.execute ();
	}
	return 0;
}