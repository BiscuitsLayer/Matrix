all:	fb b r
fb:
		$(MAKE) -C Reader/Build
b:
		g++ main.cpp Reader/Language/driver.cpp Matrix/Matrix.cpp Reader/Build/lex.yy.cc Reader/Build/lang.tab.cc -ggdb3 -o main
r:
		./main Test/Circuit/1
