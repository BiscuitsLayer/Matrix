%option c++
%option noyywrap nounput noinput batch debug
%option yylineno
%option yyclass="SyntaxCheck"

%{

//	BISON AND FLEX
#include "../Build/lang.tab.hh"

//  LOCATION UPDATER
#include "../Language/SyntaxCheck.hpp"
#define YY_USER_ACTION SetLocation ();

%}


WS              [ \t]+
CR              [\n\v]+
DOT             [.]
SIGN            [+-]
DIGIT           [0-9]
DIGIT_NON_ZERO  [1-9]
DOUBLE          {SIGN}?{DIGIT}+{DOT}?{DIGIT}*
UINT             {DIGIT}+
NUMBER          {DOUBLE}

%%
                            /* Сопоставляем каждой лексеме её константу для парсера */
{WS}                        /* skip blanks and tabs */
{CR}                        /* next line or vertical tab */
"R"                         return yy::parser::token_type::RESISTANCE;
"V"                         return yy::parser::token_type::VOLTAGE;
"--"                        return yy::parser::token_type::DOUBLEDASH;
","                         return yy::parser::token_type::COMMA;
";"                         return yy::parser::token_type::SEMICOLON;
{UINT}                      return yy::parser::token_type::UINT;
{DOUBLE}                    return yy::parser::token_type::DOUBLE;
.			                return yy::parser::token_type::ERROR;

%%