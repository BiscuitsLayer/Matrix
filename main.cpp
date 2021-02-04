#include <fstream>
#include <algorithm>

//#include "Reader/Language/driver.hpp"
#include "Matrix/Matrix.hpp"
//#include "Solver/Solver.hpp"
//#include "Circuit/Circuit.hpp"

#include "Generator/Generator.hpp"

int main (int argc, char** argv) {
	/*
	std::ifstream infile { argv[1] };
	if (!infile) {
		ERRSTREAM << "Error opening file!" << std::endl;
		return 0;
	}
	*/
	Generator test {};
	test.Execute ();
	/*
	yy::LangDriver driver { infile };
	if (driver.parse ()) {
		driver.execute ();
	}
	
	Linear::Matrix <double> m1 {};
	infile >> m1;
	std::cout << m1 << std::endl << "Determinant = " << m1.Determinant (Linear::Determinant::Type::GAUSS) << std::endl;
	*/
	return 0;
}