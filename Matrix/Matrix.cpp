#include "Matrix.hpp"

double Linear::Determinant::Gauss (const Linear::Matrix <double>& matrix) {
	auto shape = matrix.Shape ();
	int nRows = shape.first, nCols = shape.second;
	if (nRows != nCols) {
		throw (std::invalid_argument ("Trying to calcute non-square matrix determinant."));
	}
	if (nRows == 1) {
		return matrix.At (0, 0);
	}
	else {
		double ans = 1.0;
		int gaussFactor = 1;
		Linear::Matrix <double> temp = matrix;
		temp.DirectGauss (&gaussFactor);
		ans *= gaussFactor;

		for (int i = 0; i < nRows; ++i) {
			ans *= temp.At (i, i);
		}
		ans = (std::fabs (ans) < EPS ? 0 : ans);
		//return std::round (ans);
		return ans;
	}
}