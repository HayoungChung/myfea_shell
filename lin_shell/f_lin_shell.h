#ifndef F_LIN_SHELL
#define F_LIN_SHELL

#include "./FEA_hy.h"
#include "./f_core_element.h"

using namespace Eigen;

void f_lin_shell(class FEAMesh & feaMesh, std::vector<Material_ABD> & material, 
                    struct Force & force, SparseMatrix<double> & sGKT, VectorXd& Res);

double getArea(MatrixXd x);


#endif