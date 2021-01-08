#include <fstream>

#include "MatrixLib.hpp"

int main () {
	Matrix <double> m1 {}, m2 {};
	std::ifstream infile { "tests/test_mul.txt" };
	if (infile) {
		infile >> m1 >> m2;
		std::cout << std::boolalpha << (m1 != m2) << std::endl << (m1 *= m2).Clear () << std::endl;
	}
	else {
		std::cerr << "Error opening file" << std::endl;
	}
	return 0;
}