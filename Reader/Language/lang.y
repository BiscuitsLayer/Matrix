%language "c++"

%skeleton "lalr1.cc"
%defines
%define api.value.type variant
%define parse.error custom

%param {LangDriver* driver}
%locations

%code requires
{
	//	SYSTEM
	#include <algorithm>
	#include <string>
	#include <vector>

	//	CIRCUIT
	#include "../../Circuit/Circuit.hpp"

	// forward declaration of argument to parser
	namespace yy { class LangDriver; }

}

%code
{
	//	BISON AND FLEX
	#include "../Language/driver.hpp"

	//	CIRCUIT
	#include "../../Circuit/Circuit.hpp"

	namespace yy {

		parser::token_type yylex (parser::semantic_type* yylval, parser::location_type* location, LangDriver* driver);

	}

}

%token
/* Список токенов */
  RESISTANCE    "R"
  VOLTAGE       "V"
  DOUBLEDASH    "--"
  COMMA         ","
  SEMICOLON     ";"
  ERROR
;

%token <unsigned int> UINT
%token <double> DOUBLE

/* Объявление нетерминалов */
%nterm <Edge> edge
%nterm <RV> values

/* Левоассоциативные и правоассоциативные лексемы */
/* empty */

%start command

%%

command:
    command edge COMMA values					{ 
													driver->TableAt ($2.first, $2.second) = $4;
													driver->TableAt ($2.second, $2.first) = RV { $4.Resistance (), -1 * $4.Voltage () };
												}
|
;

edge:
    UINT DOUBLEDASH UINT						{ $$ = Edge { $1, $3 }; }
;

values:
    DOUBLE resistance SEMICOLON					{ $$ = RV { $1, 0 }; } 
|   DOUBLE resistance SEMICOLON DOUBLE voltage	{ $$ = RV { $1, $4 }; } 
;

resistance:
    RESISTANCE
|
;

voltage:
    VOLTAGE
|
;

%%

namespace yy {

	parser::token_type yylex (parser::semantic_type* yylval, parser::location_type* location, LangDriver* driver) {
		return driver->yylex (yylval, location);
	}

	void parser::error (const parser::location_type& location, const std::string& what) {
		/* empty */
	}

	void parser::report_syntax_error (parser::context const& context) const {
		driver->PrintErrorAndExit (context.location (), "Syntax error!");
	}

}