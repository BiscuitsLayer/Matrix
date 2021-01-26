#pragma once

//TODO:
#define OUTSTREAM std::cout
#define ERRSTREAM std::cerr

//TODO: may be in paracl they also should be smwhr in settings
//  ERROR CODES
enum ErrorCodes {
    ERROR_INV_ARG = 1,
    ERROR_OVF = 2,
    ERROR_SYNTAX = 3
};

//  SYSTEM
#include <cstring>
#include <fstream>

//	BISON AND FLEX
#include "../Build/lang.tab.hh"
#include "../Language/SyntaxCheck.hpp"

//  MATRIX
#include "../../Matrix/Matrix.hpp"
#include "../../Solver/Solver.hpp"

namespace yy {

    class LangDriver {
        private:
            SyntaxCheck* lexer_ {};
            Linear::Matrix <RV> adjTable_ {};
        public:
            //  METHODS
            parser::token_type yylex (parser::semantic_type* yylval, parser::location_type* location);
            bool parse ();
            void execute ();

            RV& TableAt (int i, int j);
            void PrintErrorAndExit (yy::location location, const std::string& message) const;
            std::string GetCurrentString () const;

            //  CTOR AND DTOR
            LangDriver (std::ifstream& infile);
            ~LangDriver ();
    };

}