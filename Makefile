all:	fb b r
fb:
		$(MAKE) -C Reader/Build
b:
		g++ main.cpp Reader/Language/driver.cpp Reader/Language/SyntaxCheck.cpp \
		Matrix/Matrix.cpp Solver/Solver.cpp Circuit/Circuit.cpp \
		Reader/Build/lex.yy.cc Reader/Build/lang.tab.cc -ggdb3 -o main
r:
		./main Test/Input/Determinant/1