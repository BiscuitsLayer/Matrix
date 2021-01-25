#pragma once

#include "../MatrixLib.hpp"

class Solver {
    private:
        Linear::Matrix <double> main_ {};
        Linear::Matrix <double> additional_ {};
        Linear::Matrix <double> ansFundamental_ {};
        Linear::Matrix <double> ansParticular_ {};
    public:
        Solver (Linear::Matrix <double> main, Linear::Matrix <double> additional):
            main_ (main),
            additional_ (additional),
            ansFundamental_ ({}),
            ansParticular_ ({})
            {}
        void Execute () {
            auto shape = main_.Shape ();
            auto mainVec = static_cast <std::vector <double> > (main_.Transpose ());
	        auto additionalVec = static_cast <std::vector <double> > (additional_);
	        mainVec.insert (mainVec.end(), additionalVec.begin (), additionalVec.end ());
	        Linear::Matrix <double> result { shape.second + 1, shape.first, mainVec };
            
            int mainRank = -1, resultRank = -1;
            main_.Transpose ();
            result.Transpose ();
            main_.Diagonalize (mainRank, false);
            result.Diagonalize (resultRank, true);
            
            if (mainRank != resultRank) {
                std::cout << "No solutions: mainRank = " << mainRank << ", resultRank = " << resultRank << std::endl;
                return;
            }
            
            shape = result.Shape ();
            int nRows = shape.first, nCols = shape.second;
            std::vector <double> ansFundamentalVec {};
            std::vector <double> ansParticularVec {};

            //  FUNDAMENTAL
            for (int i = 0; i < std::min <int> (nRows, nCols - 1); ++i) {
                for (int j = mainRank; j < nCols - 1; ++j) {
                    ansFundamentalVec.push_back (result.At (i, j));
                }
            }
            ansFundamentalVec.resize ((nCols - 1) * ((nCols - 1) - mainRank));
            //  PARTICULAR
            for (int i = 0; i < std::min <int> (nRows, nCols - 1); ++i) {
                ansParticularVec.push_back (result.At (i, nCols - 1));
            }
            ansParticularVec.resize ((nCols - 1) * 1);

            ansFundamental_ = { nCols - 1, (nCols - 1) - mainRank, ansFundamentalVec};
            for (int i = 0; i < (nCols - 1) - mainRank; ++i) {
                ansFundamental_.At (mainRank + i, i) = -1;
            }
            ansFundamental_.Negate ();
            //DEBUG (ansFundamental_);
            ansParticular_ = { nCols - 1, 1, ansParticularVec };
            //DEBUG (ansParticular_);
        }
};