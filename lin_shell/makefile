#do not use -O3 (inline)

all:
	g++ -O3 -fopenmp -W -std=c++14 -c FEA_hy.cpp
	g++ -O3 -fopenmp -W -std=c++14 -c f_DKT_OPT.cpp
	#g++ -O3 -fopenmp -W -std=c++14 -c Core_element.cpp // NEED DEBUG
	g++ -O3 -fopenmp -W -std=c++14 -c f_core_element.cpp
	g++ -O3 -fopenmp -W -std=c++14 -c lin_shell.cpp
	g++ -O3 -fopenmp -W -std=c++14 -c LinSensitivity.cpp

fea:
	++ -O3 -fopenmp -W -std=c++14 -c s_linshell.cpp
	g++ -fopenmp -o s_linshell s_linshell.o lin_shell.o FEA_hy.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt

main:
	g++ -O3 -fopenmp -W -std=c++14 -c main_lin.cpp
	g++ -fopenmp -o main_lin main_lin.o lin_shell.o FEA_hy.o LinSensitivity.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt

m2doSens:
	g++ -O3 -fopenmp -W -std=c++14 -c m2do_lin_sensitivity.cpp
	g++ -O3 -fopenmp -W -std=c++14 -c main_testSens.cpp
	g++ -fopenmp -o main_testSens main_testSens.o lin_shell.o m2do_lin_sensitivity.o FEA_hy.o LinSensitivity.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt