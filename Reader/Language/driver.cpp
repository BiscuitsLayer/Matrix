#include "driver.hpp"

yy::LangDriver::LangDriver (std::ifstream& infile):
    lexer_ (new SyntaxCheck)
    {
        lexer_->switch_streams (infile, OUTSTREAM);
    }

yy::LangDriver::~LangDriver () {
    delete lexer_;
}

yy::parser::token_type yy::LangDriver::yylex (yy::parser::semantic_type* yylval, parser::location_type* location) {
    yy::parser::token_type tokenType = static_cast <yy::parser::token_type> (lexer_->yylex ());
    switch (tokenType) {
        case yy::parser::token_type::UINT: {
            //  Getting the number itself
            yylval->as <unsigned int> () = std::stod (lexer_->YYText ());
            break;
        }
        case yy::parser::token_type::DOUBLE: {
            //  Getting the number itself
            yylval->as <double> () = std::stod (lexer_->YYText ());
            break;
        }
        default: {
            break;
        }
    }
    *location = lexer_->GetLocation ();
    return tokenType;
}

bool yy::LangDriver::parse () {
    yy::parser parser (this);
    bool failure = parser.parse ();   
    if (!failure) {
        std::cout << "Parsed successfully!" << std::endl;
    }
    return !failure;
}

void yy::LangDriver::execute () {
    Circuit circuit { adjTable_ };
    PairMatrix temp = circuit.Execute ();
    Solver solver { temp.first, temp.second };
    temp = solver.Execute ();
    DEBUG (temp.first);
    DEBUG (temp.second);
}

RV& yy::LangDriver::TableAt (int i, int j) {
    auto shape = adjTable_.Shape ();
    shape.first = std::max <int> (i + 1, shape.first);
    shape.second = std::max <int> (j + 1, shape.second);
    if (shape != adjTable_.Shape ()) {
        adjTable_.Resize (shape);
    }
    return adjTable_.At (i, j);
}

void yy::LangDriver::PrintErrorAndExit (yy::location location, const std::string& message) const {
    std::string wholeString = lexer_->GetCurrentString (), 
    trueString = wholeString.substr (0, location.begin.column - 1),
    falseString = wholeString.substr (location.begin.column - 1, location.end.column - 1);
    
    ERRSTREAM << message << std::endl;
    OUTSTREAM << "Line: " << location.begin.line << ", Columns: " << location.begin.column << " - " << location.end.column << ":" << std::endl;
    OUTSTREAM << trueString;
    ERRSTREAM << falseString;
    OUTSTREAM << " ..." << std::endl;

    exit (ErrorCodes::ERROR_SYNTAX);
}

std::string yy::LangDriver::GetCurrentString () const {
    return lexer_->GetCurrentString ();
}