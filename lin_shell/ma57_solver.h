#ifndef M2DO_FEA_MA57_SOLVER_H
#define M2DO_FEA_MA57_SOLVER_H

// #include <iostream>
// #include "./../../eigen3/Eigen/Dense"
// #include "./../../eigen3/Eigen/Sparse"

using namespace std ;
using namespace Eigen ;

/* 
	MA57 Functions
*/

extern "C" void ma57id_(double *cntl, int *icntl) ;

extern "C" void ma57ad_(int *n, int *nz, int *irn, int *jcn, int *lkeep, int *keep, int *iw, int *icntl, int *info, double *rinfo) ;
    
extern "C" void ma57bd_(int *n, int *nz, double *a, double *fact, int *lfact, int *ifact, int *lifact, int *lkeep, int *keep, int *iw, int *icntl, double *cntl, int *info, double *rinfo) ;
    
extern "C" void ma57cd_(int *job, int *n, double *fact, int *lfact, int *ifact, int *lifact, int *nrhs, double *rhs, int *lrhs, double *w, int *lw, int *iw, int *icntl, int *info) ;

class MA57Solver {
	
	private:
		//

	public:
		// Properties:
		VectorXd cntl, fact ;
		VectorXi icntl, ifact, iw, info ;
		int n, lfact, lifact ;

		// Methods:
		MA57Solver (bool use_metis, bool print_output) ;
		void print () ;
		void compute (SparseMatrix<double> & A) ;
		VectorXd solve (VectorXd & b) ;

} ;

#include "./ma57_solver.tpp"

#endif