#pragma once

#include "../MatrixLib.hpp"

class Solver {
    private:
        Linear::Matrix <double> main_ {};
        Linear::Matrix <double> additional_ {};
        Linear::Matrix <double> ansGeneral_ {};
        Linear::Matrix <double> ansParticular_ {};
    public:
        Solver (Linear::Matrix <double> main, Linear::Matrix <double> additional):
            main_ (main),
            additional_ (additional),
            ansGeneral_ ({}),
            ansParticular_ ({})
            {}
        void Execute () {
            auto shape = main_.Shape ();
            auto mainVec = static_cast <std::vector <double> > (main_.Transpose ());
	        auto additionalVec = static_cast <std::vector <double> > (additional_);
	        mainVec.insert (mainVec.end(), additionalVec.begin (), additionalVec.end ());
	        Linear::Matrix <double> result { shape.second + 1, shape.first, mainVec };
            /* можно после if */
            result.Transpose ();
            std::cerr << main_ << std::endl;
            if (main_.Rank () != result.Rank ()) {
                std::cout << "No solutions" << std::endl;
                return;
            }
            //result.Transpose ();
            shape = result.Shape ();
	        size_t nRows = shape.first, nCols = shape.second;
	        for (size_t i = 0; i < std::min (nRows - 1, nCols); ++i) {
	        	double maxElement = {};
	        	size_t maxIdx = i + 1;
	        	for (size_t j = i + 1; j < nRows; ++j) {
	        		if (result.At (j, i) > maxElement) {
	        			maxElement = result.At (j, i);
	        			maxIdx = j;
	        		}
	        	}
	        	if (maxElement < Linear::EPS) {
	        		//	Matrix has a zero-column
	        		continue;
	        	}
	        	else {
	        		result.SwapRows (i, maxIdx);
	        	}
	        	for (size_t j = i + 1; j < nRows; ++j) {
	        		result.AddRows (i, j, (- 1) * (result.At (j, i) / result.At (i, i)));
	        	}
	        }
            std::cout << result << std::endl;
        }
};