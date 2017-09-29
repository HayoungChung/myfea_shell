#ifndef F_DKT_OPT_H
#define F_DKT_OPT_H

// #include <iostream>
// #include <cmath>
// #include <algorithm>
// #include "./EICR.h"
#include "./FEA_hy.h"

// #include "../../../../eigen3/Eigen/Dense"
// #include "../../../../eigen3/Eigen/Sparse"
using namespace Eigen;

struct DKT_OPT{
    MatrixXd M_thetax = MatrixXd::Zero(6,9), M_thetay = MatrixXd::Zero(6,9);
    MatrixXd M_etax = MatrixXd::Zero(3,9), M_etay = MatrixXd::Zero(3,9), 
                        M_gamxy = MatrixXd::Zero(3,9), M_psi = MatrixXd::Zero(3,9);
};

struct OPT_AUX{
    MatrixXd B = MatrixXd::Zero(3,9), Ttu = MatrixXd::Zero(3,9), 
                Te = MatrixXd::Zero(3,3), Q1 = MatrixXd::Zero(3,3),
                 Q2 = MatrixXd::Zero(3,3), Q3 = MatrixXd::Zero(3,3);
    double A;
}; 

DKT_OPT f_DKT_OPT(const MatrixXd & xycoord, double const nu);
void OPT_shape_function(const MatrixXd & xycoord, double const nu, MatrixXd& M_etax, MatrixXd& M_etay, MatrixXd& M_gamxy, MatrixXd& M_psi);
OPT_AUX f_OPT_aux(std::vector<double> fpars, const MatrixXd & xycoord);
void DKT_shape_functions(const VectorXd x, const VectorXd y, MatrixXd& M_thetax, MatrixXd& M_thetay);

#endif 