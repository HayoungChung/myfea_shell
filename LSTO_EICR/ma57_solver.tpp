#include "./ma57_solver.h"
#include <iostream>
MA57Solver :: MA57Solver (bool use_metis, bool print_output) {
	
	/*

		Initialize MA57:

	*/

	cntl.resize(5) ;
	icntl.resize(20) ;

	ma57id_(cntl.data(), icntl.data()) ;

	// Ask for full printing:
	// diagonistic printing: if larger than 3, it prints all
	icntl(4) = -1 ; 
	if (print_output) {
		icntl(4) = 4 ; 
	}
	
	// Use Metis ordering:
	if (use_metis) {
		icntl(5) = 4 ;
	}

}

void MA57Solver :: print () {

	cout << "MA57 Solver." ;

}

void MA57Solver :: compute (SparseMatrix<double> & A) {

	/*
		Check dimensions:
	*/

	if ( A.rows() != A.cols() ) {
		
		throw invalid_argument( "\n\n*****\nInput matrix dimension problem in MA57Solver.compute()\n*****\n\n" );
	
	}

	/*
		Analyze the matrix A:
	*/

	n  = A.rows() ; // Matrix order.
	int nz = A.nonZeros() ;

    VectorXi irn (nz) ;
    VectorXi jcn (nz) ;
    VectorXd a (nz) ;

    int counter = 0 ;

    for (int k=0; k < A.outerSize(); ++k) {

        for (SparseMatrix<double>::InnerIterator it(A,k); it; ++it) {
            
            irn (counter) = 1+it.row() ; // row index
            jcn (counter) = 1+it.col() ; // col index (here it is equal to k)
            a (counter) = it.value() ;

            ++counter ;

        }
    }

	int lkeep = 5*n + nz + max(n, nz) + 42 + n ;
	VectorXi keep(lkeep) ;
    iw.resize(10*n) ;
    info.resize(40) ;       // integer array of length 40
    VectorXd rinfo(20) ;   // double array of length 20 

	ma57ad_(&n, &nz, irn.data(), jcn.data(), &lkeep, keep.data(), iw.data(), icntl.data(), info.data(), rinfo.data()) ;

	/*
		Factorize the matrix A:
	*/

    lfact = 2 * info(8) ;
    fact.resize(lfact) ;
    lifact = 2 * info(9) ;
    ifact.resize(lifact) ;

    ma57bd_(&n, &nz, a.data(), fact.data(), &lfact, ifact.data(), &lifact, &lkeep, keep.data(), iw.data(), icntl.data(), cntl.data(), info.data(), rinfo.data()) ;

}

VectorXd MA57Solver :: solve (VectorXd & b) {

	/*
		Check dimensions:
	*/

	if ( b.size() != n ) {
		
		throw invalid_argument( "\n\n*****\nInput vector dimension problem in MA57Solver.solve()\n*****\n\n" );
	
	}

	/*
		Solve that mofo:
	*/

	int job = 1 ; 	  // if job <= 1, solves Ax = b.
	int nrhs = 1 ;    // number of right hand side being solved
    int lw = n*nrhs ; // length of w; lw>=n*nrhs
    VectorXd w (lw) ; // double workspace
    int lrhs = n ;    // integer, length of rhs
    
    VectorXd x = b ;

    ma57cd_(&job, &n, fact.data(), &lfact, ifact.data(), &lifact, &nrhs, x.data(), &lrhs, w.data(), &lw, iw.data(), icntl.data(), info.data()) ;

	/*
		Check for errors:
	*/

	if (info(0) < 0) {

		cerr << "\n\n*****\nMA57 Solver ma57cd_ failed with info = " << info(0) << " \n*****\n\n" ;
		throw ;

	}
    
    return x ;

}

