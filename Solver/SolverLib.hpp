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

            
            int mainRank = main_.Rank (), resultRank = result.Rank ();
            if (mainRank != resultRank) {
                std::cout << "No solutions: mainRank = " << mainRank << ", resultRank = " << resultRank << std::endl;
                return;
            }
        
            result.Transpose ();
            main_.Diagonalize ();
            result.Diagonalize (true);
            DEBUG (result);
            std::vector <double> ansFundamentalVec {};
            std::vector <double> ansParticularVec {};
            //  FUNDAMENTAL
            for (int i = 0; i < nRows; ++i) {
                for (int j = mainRank; j < nCols - 1; ++j) {
                    ansFundamentalVec.push_back (result.At (i, j));
                }
            }
            for (int i = nRows; i < nCols - 1; ++i) {
                for (int j = mainRank; j < nCols - 1; ++j) {
                    ansFundamentalVec.push_back (0);
                }
            }
            //  PARTICULAR
            for (int i = 0; i < nRows; ++i) {
                ansParticularVec.push_back (result.At (i, nCols - 1));
            }
            for (int i = nRows; i < nCols - 1; ++i) {
                ansParticularVec.push_back (0);
            }
            ansFundamental_ = { nCols - 1, nCols - mainRank - 1, ansFundamentalVec};
            ansParticular_ = { nCols - 1, 1, ansParticularVec };
            for (int i = 0; i < nCols - mainRank - 1; ++i) {
                ansFundamental_.At (nCols - mainRank - 1 + i, i) = -1;
            }
            DEBUG (ansFundamental_);
            DEBUG (ansParticular_);
        }
};