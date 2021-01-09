#include <fstream>

#include "MatrixLib.hpp"

#define DEBUG

int main () {
#ifdef DEBUG
	Matrix <double> m1 {};
	std::cin >> m1;
	std::cout << m1.Determinant (Determinant::Type::GAUSS) <<std::endl;
#else
	Matrix <double> m2 {};
	std::ifstream infile { "tests/test_mul.txt" };
	if (!infile) {
		std::cerr << "Error opening file" << std::endl;
	}
	else {
		infile >> m2;
		std::cout << m2 << std::endl << m2.Transpose () 
		<< std::endl << m2.Determinant (Determinant::Type::GAUSS) << std::endl;
	}
#endif
	return 0;
}