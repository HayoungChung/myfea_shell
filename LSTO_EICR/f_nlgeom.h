#ifndef F_NLGEOM_H
#define F_NLGEOM_H

// #include <iostream>
// #include <cmath>
// #include <algorithm>
// #include "./EICR.h"
// #include "./f_EICR_shell.h"
#include "./EICR_shell.h"
#include <stdio.h>
#include <fstream>
// #include "./f_nlgeom.h"
// #include "./ma57_solver.h"
// #include "./../../../../eigen3/Eigen/OrderingMethods"
using namespace Eigen;

struct OPTION{
    int logflag = 1;
    double initper = 10;
    int saveflag = 1;
};

VectorXd f_nlgeom(std::vector<Material_ABD> &, Force &, FEAMesh &, OPTION &, 
    MatrixXd &, VectorXd &);

// MatrixXd fstore2mat(const MatrixXd xm_Rv);
void fRupdate(VectorXd & del_xm_w, int nNODE , MatrixXd & xm_R);

#endif
