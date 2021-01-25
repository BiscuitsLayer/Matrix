all:	b r
b:
		g++ main.cpp -ggdb3 -o main
r:
		./main < Test/test_det1.txt
