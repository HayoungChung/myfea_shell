#ifndef EICR_SHELL_H
#define EICR_SHELL_H

#include "./FEA_hy.h"
#include "./CoreElement.h"

using namespace Eigen;

struct AuxMat
{
    MatrixXd Te, P_, G_, H_, Fh_nm, F_n, F_nm_ext, M_ext, M_;
};

struct GptsCompl
{
  // it comes from lin_shell..
  double x, y, z;
  double sens;
  Matrix3d N_tilde, M_tilde; // resultants // WIP
};

class EICR_SHELL
{
  public:
    EICR_SHELL(class FEAMesh &feaMesh, int ngp = 3);
    void assembly(MatrixXd &GU_u0, VectorXd &GU_Rv0, std::vector<Material_ABD> &material, struct Force &force); //, bool issens = false);
    // void sensitivity(MatrixXd & GU_u0, VectorXd & GU_Rv0, std::vector<Material_ABD> & material, struct Force & force);

    // void push_forward(Matrix3d X, Matrix3d u, Matrix3d Ra, Material_ABD & material, Matrix3d & X_R, Matrix3d & x_R, Matrix3d & u_d, Matrix<double,9,1> & th_d, Material_ABD & Mater_e);
    void Filter_Def(Matrix3d X, Matrix3d u, Matrix<double, 9, 3> Ra, Matrix3d &X_R, Matrix3d &x_R, Matrix3d &u_d, Matrix<double, 9, 1> &th_d, Matrix3d &T0, Matrix3d &T);

    // void push_forward(Matrix3d X, Matrix3d u, Matrix3d Ra, Material_ABD & material, Matrix3d & X_R, Matrix3d & x_R, Matrix3d & u_d, Matrix<double,9,1> & th_d, Material_ABD & Mater_e);
    SparseMatrix<double> sGKT;
    VectorXd Res;

    static MatrixXd fstore2mat(const VectorXd Rv);
    static VectorXd fmat2store(const MatrixXd xm_R);

    std::vector<GptsCompl> get_GaussCompl(const MatrixXd &GU_u3,const MatrixXd& GU_R, std::vector<Material_ABD> &material, struct Force &force);
    double compliance = 0;

  protected:
    std::vector<GptsCompl> gptsSens;
    std::vector<double> ri, si, wi;
    Matrix<int, 9, 1> mdof, bdof;
    MatrixXi elem_id0, elem_order;

    void get_AuxMat(Matrix3d &x_R, Matrix<double, 9, 1> &th_d, VectorXd &f_el, VectorXd &f_el_f, Matrix3d &T);
    AuxMat aux;
    MatrixXd getT(MatrixXd x);
    double getArea(MatrixXd x);
    Vector3d rot2vec(MatrixXd R);

  private:
    FEAMesh &feaMesh;
    std::vector<Triplet<double>> GKT;
    // MatrixXd & GU_u0;
    // VectorXd & GU_Rv0;
    // std::vector<Material_ABD> & material;
    // Force & force;
};

#endif