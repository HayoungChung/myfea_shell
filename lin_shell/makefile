#do not use -O3 (inline)

all:
	g++ -O3 -g -fopenmp -W -std=c++14 -c FEA_hy.cpp
	g++ -O3 -g -fopenmp -W -std=c++14 -c f_DKT_OPT.cpp
	#g++ -O3 -g -fopenmp -W -std=c++14 -c Core_element.cpp // NEED DEBUG
	g++ -O3 -g -fopenmp -W -std=c++14 -c f_core_element.cpp
	g++ -O3 -g -fopenmp -W -std=c++14 -c lin_shell.cpp
	g++ -O3 -g -fopenmp -W -std=c++14  -o s_linshell s_linshell.cpp FEA_hy.o f_core_element.o  f_DKT_OPT.o  lin_shell.o -I../../../../eigen/Eigen -lblas -llapack 

build:
	g++ -O3 -g -fopenmp -W -std=c++14  -o s_linshell s_linshell.cpp FEA_hy.o f_core_element.o  f_DKT_OPT.o  lin_shell.o -I../../../../eigen/Eigen -lblas -llapack 

main:
	g++ -O3 -g -fopenmp -W -std=c++14 -c main_lin.cpp
	g++ -fopenmp -o main_lin main_lin.o lin_shell.o FEA_hy.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt