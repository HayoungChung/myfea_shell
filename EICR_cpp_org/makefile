# do not use -O3 (inline)

all:
#	g++ -g -W -std=c++11 s_shell_static.cpp FEA_hy.cpp f_core_element.cpp f_DKT_OPT.cpp f_EICR_shell.cpp f_nlgeom.cpp ma57_solver.tpp -I./../../../../eigen3/Eigen -I./
#	g++ -g -W -std=c++14 -o qq s_shell_static.cpp FEA_hy.cpp -I../../../../eigen/Eigen -I. -l blas -l lapack -l gfortran -l metis -l ma57 -L /usr/local/lib -l metis -l nlopt -l arpack -L /usr/local/lib -L /usr/local/gfortran/lib 
# #	g++ -g -W -std=c++14 -o qq s_shell_static.cpp FEA_hy.cpp f_core_element.cpp f_DKT_OPT.cpp f_EICR_shell.cpp f_nlgeom.cpp ma57_solver.tpp -I../../../../eigen/Eigen -I. -l blas -l lapack -l gfortran -l metis -l ma57 -L /usr/local/lib -l metis -l nlopt -l arpack -L /usr/local/lib -L /usr/local/gfortran/lib 
	# # icpc -O3 -fopenmp -W -std=c++14 -c ma57_solver.cpp
	g++ -O3 -W -std=c++14 -c FEA_hy.cpp
	g++ -O3 -W -std=c++14 -c f_DKT_OPT.cpp
	g++ -O3 -W -std=c++14 -c f_core_element.cpp
	g++ -O3 -W -std=c++14 -c f_EICR_shell.cpp 
	g++ -O3 -W -std=c++14 -c f_nlgeom.cpp
	g++ -O3 -W -std=c++14  -o s_shell_static s_shell_static.cpp f_nlgeom.o FEA_hy.o f_core_element.o  f_EICR_shell.o f_DKT_OPT.o  -I../../../../eigen/Eigen -lblas -llapack #-L /opt/intel/compilers_and_libraries_2017.2.174/linux/compiler/lib/intel64/ -liomp5 -L/opt/intel/mkl/lib/intel64/ # -lmkl_core -lmkl_intel_thread -lpthread
	#-l blas -l lapack -l gfortran -l metis -l ma57 -L /usr/local/lib -l metis -l arpack -L /usr/local/lib -L /usr/local/gfortran/lib  ma57_solver.o -l blas -l lapack -l ma57 -lmetis -llapack -larpack -lblas

win10_bash:
	g++ -O3 -fopenmp -W -std=c++11 -c FEA_hy.cpp
	g++ -O3 -fopenmp -W -std=c++11 -c f_DKT_OPT.cpp
	g++ -O3 -fopenmp -W -std=c++11 -c f_core_element.cpp
	g++ -O3 -fopenmp -W -std=c++11 -c f_EICR_shell.cpp 
	g++ -O3 -fopenmp -W -std=c++11 -c f_nlgeom.cpp
	g++ -O3 -fopenmp -W -std=c++11 -o win10_bash s_shell_static.cpp f_nlgeom.o FEA_hy.o f_core_element.o  f_EICR_shell.o f_DKT_OPT.o  -I../../../../eigen/Eigen -lblas -llapack 

debug:
	g++ -g -fopenmp -W -std=c++14 -c FEA_hy.cpp
	g++ -g -fopenmp -W -std=c++14 -c f_DKT_OPT.cpp
	g++ -g -fopenmp -W -std=c++14 -c f_core_element.cpp
	g++ -g -fopenmp -W -std=c++14 -c f_EICR_shell.cpp 
	g++ -g -fopenmp -W -std=c++14 -c f_nlgeom.cpp
	g++ -g -fopenmp -W -std=c++14 -o dbs s_shell_static.cpp f_nlgeom.o FEA_hy.o f_core_element.o  f_EICR_shell.o f_DKT_OPT.o  -I../../../../eigen/Eigen -lblas -llapack -L /opt/intel/compilers_and_libraries_2017.2.174/linux/compiler/lib/intel64/ -liomp5 -L/opt/intel/mkl/lib/intel64/ # -lmkl_core -lmkl_intel_thread -lpthread


clear:
	rm *.o
	rm log*
	rm suite*

test_FEA_cpp:
	g++ -g -W -std=c++11 -c FEA_hy.cpp
	g++ -g -W -std=c++11 -o suite_FEA_cpp test_FEA_hy.cpp  FEA_hy.o
	./suite_FEA_cpp > log_FEA_hy

test_DKT_OPT_cpp:
	g++ -g -W -std=c++11 -c FEA_hy.cpp
	g++ -g -W -std=c++11 -c f_DKT_OPT.cpp
	g++ -g -W -std=c++11 -o suite_DKT_OPT test_DKT_OPT.cpp  FEA_hy.o f_DKT_OPT.o
	./suite_DKT_OPT > log_DKT_OPT
	
test_core_elem:
	g++ -g -W -std=c++11 -c FEA_hy.cpp
	g++ -g -W -std=c++11 -c f_DKT_OPT.cpp
	g++ -g -W -std=c++11 -c f_core_element.cpp
	g++ -g -W -std=c++11 -o suite_core_elem test_core_elem.cpp  FEA_hy.o f_DKT_OPT.o f_core_element.o
	./suite_core_elem > log_core_elem

test_EICR:
	g++ -g -W -std=c++11 -c FEA_hy.cpp
	g++ -g -W -std=c++11 -c f_DKT_OPT.cpp
	g++ -g -W -std=c++11 -c f_core_element.cpp
	g++ -g -W -std=c++11 -c f_EICR_shell.cpp
	g++ -g -W -std=c++11 -o suite_EICR test_EICR.cpp f_EICR_shell.o f_core_element.o f_DKT_OPT.o  FEA_hy.o
	./suite_EICR > log_EICR

test_nlgeom:
