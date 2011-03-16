main:	main.o nahs_planar.o makefile
		g++ -o main main.o nahs_planar.o -lfftw3 -lm

main.o:	main.cpp
		g++ -c -Wall main.cpp

nahs_planar.o:	nahs_planar.cpp phi_code.c
		g++ -O3 -c -Wall nahs_planar.cpp
