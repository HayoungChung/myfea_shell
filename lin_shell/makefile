# do not use -O3 (inline)

all:
	g++ -O3 -g -fopenmp -W -std=c++14 -c FEA_hy.cpp
	g++ -O3 -g -fopenmp -W -std=c++14 -c f_DKT_OPT.cpp
	g++ -O3 -g -fopenmp -W -std=c++14 -c f_core_element.cpp
	g++ -O3 -g -fopenmp -W -std=c++14 -c f_lin_shell.cpp
	g++ -O3 -g -fopenmp -W -std=c++14  -o s_linshell s_linshell.cpp FEA_hy.o f_core_element.o  f_DKT_OPT.o  f_lin_shell.o -I../../../../eigen/Eigen -lblas -llapack 

build:
	g++ -O3 -g -fopenmp -W -std=c++14  -o s_linshell s_linshell.cpp FEA_hy.o f_core_element.o  f_DKT_OPT.o  f_lin_shell.o -I../../../../eigen/Eigen -lblas -llapack 