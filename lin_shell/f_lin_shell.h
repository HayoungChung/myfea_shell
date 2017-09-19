#ifndef F_EICR_SHELL_H
#define F_EICR_SHELL_H

#include "./FEA_hy.h"
#include "./f_core_element.h"

using namespace Eigen;

void f_lin_shell(class FEAMesh & feaMesh, std::vector<Material_ABD> & material, 
                    struct Force & force, SparseMatrix<double> & sGKT, VectorXd& Res);

MatrixXd getT(MatrixXd x);
double getArea(MatrixXd x);


#endif