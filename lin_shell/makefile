#do not use -g -O3 (inline)

all:
	g++ -g -O3 -fopenmp -W -std=c++14 -c FEA_hy.cpp
	g++ -g -O3 -fopenmp -W -std=c++14 -c f_DKT_OPT.cpp
	#g++ -g -O3 -fopenmp -W -std=c++14 -c Core_element.cpp // NEED DEBUG
	g++ -g -O3 -fopenmp -W -std=c++14 -c f_core_element.cpp
	g++ -g -O3 -fopenmp -W -std=c++14 -c lin_shell.cpp
	#g++ -g -O3 -fopenmp -W -std=c++14 -c LinSensitivity.cpp
	g++ -g -O3 -fopenmp -W -std=c++14 -c m2do_lin_sensitivity.cpp

fea:
	g++ -g -O3 -fopenmp -W -std=c++14 -c s_linshell.cpp
	g++ -fopenmp -o s_linshell s_linshell.o lin_shell.o FEA_hy.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt

main:
	g++ -g -O3 -fopenmp -W -std=c++14 -c main_lin.cpp
	g++ -fopenmp -o main_lin main_lin.o lin_shell.o FEA_hy.o LinSensitivity.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt

m2doSens:
	g++ -g -O3 -fopenmp -W -std=c++14 -c main_testSens.cpp
	g++ -fopenmp -o main_testSens main_testSens.o lin_shell.o m2do_lin_sensitivity.o FEA_hy.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt

nonIsOverlaid:
	g++ -g -O3 -fopenmp -W -std=c++14 -c nonIsOverlaid.cpp
	g++ -fopenmp -o nonIsOverlaid nonIsOverlaid.o lin_shell.o m2do_lin_sensitivity.o FEA_hy.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt

checkRadMesh:
	g++ -g -O3 -fopenmp -W -std=c++14 -c mesh_test.cpp
	g++ -fopenmp -o mesh_test mesh_test.o FEA_hy.o -I/home/hac210/00.Working/0.Working/eigen3/Eigen 

nzCurvFea:
	g++ -g -O3 -fopenmp -W -std=c++14 -c main_nzCurv.cpp
	g++ -fopenmp -o main_nzCurv main_nzCurv.o lin_shell.o FEA_hy.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt -lmetis -lma57 -lgfortran -larpack -lblas -lmetis
	

ma57:
	g++ -g -O3 -fopenmp -W -std=c++14 -c s_linshell.cpp
	g++ -fopenmp -o s_linshell s_linshell.o lin_shell.o FEA_hy.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt -lmetis -lma57 -lgfortran -larpack -lblas -lmetis

nzCurvLSTO:
	g++ -g -O3 -fopenmp -W -std=c++14 -c LSTO_nzCurv.cpp
	g++ -fopenmp -o LSTO_nzCurv LSTO_nzCurv.o lin_shell.o m2do_lin_sensitivity.o FEA_hy.o f_core_element.o f_DKT_OPT.o -I/home/hac210/00.Working/0.Working/13.GeometricNonlinear/M2DO_LSM/include -I/home/hac210/00.Working/0.Working/eigen3/Eigen -lblas -llapack -L /usr/local/lib -lnlopt -lmetis -lma57 -lgfortran -larpack -lblas -lmetis
