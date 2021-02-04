#pragma once

//  MATRIX
#include "../Matrix/Matrix.hpp"

//  TYPEDEFS
using PairMatrix = std::pair <Linear::Matrix <double>, Linear::Matrix <double>>;

class Solver final {
    private:
        //  GIVEN
        Linear::Matrix <double> main_ {};
        Linear::Matrix <double> additional_ {};
        Linear::Matrix <double> result_ {};

        //  COMPUTATIONS
        Linear::Matrix <double> ansFundamental_ {};
        Linear::Matrix <double> ansParticular_ {};
    public:
        //  CTOR
        Solver (Linear::Matrix <double> main, Linear::Matrix <double> additional):
            main_ (main),
            additional_ (additional),
            result_ ({}),
            ansFundamental_ ({}),
            ansParticular_ ({})
            {}

        //  CHECK RANK EQUALITY
        int CheckRank ();

        //  SOLVE
        void CreateFundamental (int generalRank);
        void CreateParticular ();

        //  EXECUTE
        PairMatrix Execute ();
};