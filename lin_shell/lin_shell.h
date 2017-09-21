#ifndef LIN_SHELL
#define LIN_SHELL

#include "./FEA_hy.h"
#include "./f_core_element.h"

using namespace Eigen;

class LinShell
{
public:
  LinShell(class FEAMesh &feaMesh, std::vector<Material_ABD> &material,
           struct Force &force);
  void compute();
  SparseMatrix<double> sGKT;
  VectorXd Res;

private:
  FEAMesh &feaMesh;
  std::vector<Material_ABD> &material;
  Force &force;
  double getArea(MatrixXd x);
  int dpn, dpe, npe, nNODE, nELEM, nDOF;
  bool isOverlaid;

  MatrixXd FNM, Ffix;
  MatrixXd NODE;
  MatrixXi ELEM;

  // std::vector<Triplet<double>> GKT;
};

#endif