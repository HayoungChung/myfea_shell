#ifndef F_EICR_SHELL_H
#define F_EICR_SHELL_H

// #include <iostream>
// #include <cmath>
// #include <cassert>
// #include "./EICR.h"
// #include "./../../eigen3/Eigen/Dense"
// #include "./../../eigen3/Eigen/Sparse"
// #include "./../../../../eigen3/unsupported/Eigen/SparseExtra" 

#include "./FEA_hy.h"
#include "./f_core_element.h"

using namespace Eigen;

struct AuxMat{
    MatrixXd Te, P_, G_, H_, Fh_nm, F_n, F_nm_ext, M_ext, M_;
};

void f_EICR_shell(class FEAMesh & feaMesh, MatrixXd & GU_u0, VectorXd & GU_Rv0, std::vector<Material_ABD> & material, 
                    struct Force & force, SparseMatrix<double> & sGKT, VectorXd& Res);
AuxMat f_EICR(Matrix3d & x_R, Matrix<double,9,1> & th_d, VectorXd & f_el, VectorXd & f_el_f, Matrix3d & T);

MatrixXd getT(MatrixXd x);
double getArea(MatrixXd x);
Vector3d rot2vec(MatrixXd R);
MatrixXd fstore2mat(const VectorXd Rv);
VectorXd fmat2store(const MatrixXd xm_R);


#endif