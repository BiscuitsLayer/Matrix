#pragma once

//  SYSTEM
#include <cstring>
#include <fstream>

//	BISON AND FLEX
#include "../Build/lang.tab.hh"
#include "../Language/SyntaxCheck.hpp"

//  SOLVER
#include "../../Solver/Solver.hpp"

//  SETTINGS
#include "../../Settings/Settings.hpp"

namespace yy {
    class LangDriver {
        private:
            //  LEXER
            SyntaxCheck* lexer_ {};

            //  CIRCUIT STUFF
            Linear::Matrix <RV> adjTable_ {};
            std::vector <Edge> givenEdges_ {};
        public:
            //  METHODS
            parser::token_type yylex (parser::semantic_type* yylval, parser::location_type* location);
            bool parse ();
            void execute ();

            //  CIRCUIT METHODS
            RV& TableAt (int i, int j);
            void PushGivenEdge (Edge edge);

            //  ERROR HANDLING METHODS
            void PrintErrorAndExit (yy::location location, const std::string& message) const;
            std::string GetCurrentString () const;

            //  CTOR
            LangDriver (std::ifstream& infile):
                lexer_ (new SyntaxCheck)
                {
                    lexer_->switch_streams (infile, OUTSTREAM);
                }

            //  DTOR
            ~LangDriver () { delete lexer_; }
    };
}