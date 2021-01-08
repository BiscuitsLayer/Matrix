#include <fstream>

#include "MatrixLib.hpp"

int main () {
	Matrix <double> m1 {};
	std::cin >> m1;
	std::cout << m1.Determinant (Determinant::Type::GAUSS) <<std::endl;
	return 0;
}