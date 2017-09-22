#ifndef LIN_SHELL
#define LIN_SHELL

#include "./FEA_hy.h"
#include "./f_core_element.h"
// #include "./Core_element.h"
using namespace Eigen;

struct GptsCompl
{
  double x, y, z;
  double sens;
};

class LinShell
{
public:
  LinShell(class FEAMesh &feaMesh, std::vector<Material_ABD> &material,
           struct Force &force);
  void compute();
  std::vector<GptsCompl> get_GaussCompl(MatrixXd &GU_u6);

  SparseMatrix<double> sGKT;
  VectorXd Res;

  std::vector<GptsCompl> gptsSens;

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
