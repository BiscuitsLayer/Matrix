#include "Solver.hpp"

int Solver::CheckRank () {
    int mainRank = main_.Rank (), resultRank = result_.Rank ();
    if (mainRank != resultRank) {
        std::stringstream strstream {};
        strstream << "No solutions: mainRank = " << mainRank << ", resultRank = " << resultRank;
        throw std::invalid_argument (strstream.str ());
    }
    return mainRank;
}

void Solver::CreateFundamental (int generalRank) {
    int nRows = result_.Shape ().first, nCols = result_.Shape ().second;
    std::vector <double> ansFundamentalVec {};
    for (int i = 0; i < std::min <int> (nRows, nCols - 1); ++i) {
        for (int j = generalRank; j < nCols - 1; ++j) {
            ansFundamentalVec.push_back (result_.At (i, j));
        }
    }
    ansFundamentalVec.resize ((nCols - 1) * ((nCols - 1) - generalRank));
    ansFundamental_ = { nCols - 1, (nCols - 1) - generalRank, ansFundamentalVec};
    for (int i = 0; i < (nCols - 1) - generalRank; ++i) {
        ansFundamental_.At (generalRank + i, i) = -1;
    }
    ansFundamental_.Negate ();
}

void Solver::CreateParticular () {
    int nRows = result_.Shape ().first, nCols = result_.Shape ().second;
    std::vector <double> ansParticularVec {};
    for (int i = 0; i < std::min <int> (nRows, nCols - 1); ++i) {
        ansParticularVec.push_back (result_.At (i, nCols - 1));
    }
    ansParticularVec.resize ((nCols - 1) * 1);
    ansParticular_ = { nCols - 1, 1, ansParticularVec };
}

PairMatrix Solver::Execute () {
    result_ = main_;
    result_.AppendCols (additional_);
    main_.MakeEye (false);
    result_.MakeEye (true);
    int generalRank = CheckRank ();
    CreateFundamental (generalRank);
    CreateParticular ();
    return { ansFundamental_, ansParticular_ };
}