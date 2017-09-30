#ifndef EICR_CORE_ELEMENT_H
#define EICR_CORE_ELEMENT_H
// this is an oop version of f_core_element.cpp

// #include "./EICR.h"
#include "./f_DKT_OPT.h"
#include "./FEA_hy.h"

using namespace Eigen;

// struct FilteredMat{
//     // MatrixXd Km(9,9), Kb(9,9), Kmb(9,9);
//     // MatrixXd Fm(9,1), Fb(9,1), Fmf(9,1), Fbf(9,1);
//     MatrixXd Km = MatrixXd::Zero(9,9), Kb = MatrixXd::Zero(9,9), Kmb = MatrixXd::Zero(9,9);
//     VectorXd Fm = VectorXd::Zero(9), Fb = VectorXd::Zero(9), Fmf = VectorXd::Zero(9), Fbf = VectorXd::Zero(9);
// };

struct FilteredP
{
    Matrix<double, 9, 1> membrane = Matrix<double, 9, 1>::Zero(9), bending = Matrix<double, 9, 1>::Zero(9);
};

class CoreElement
{
  public:
    CoreElement(Matrix<double, 2, 3> &xycoord_, Material_ABD &Mater_e_); //Matrix<double, 6, 1> & FNM_r_, FilteredP & p_d_);

    Vector2d get_StressStrain(FilteredP &, Matrix<double, 6, 1> &, double, double, Matrix<double, 1, 6> &, Matrix<double, 1, 6> &, Matrix3d &); //MatrixXd & );
    // Vector2d get_StressStrain_noVoigt(FilteredP &, Matrix<double, 6, 1> &, double, double, Matrix<double, 1, 6> &, Matrix<double, 2, 4> &, Matrix3d &); //MatrixXd & );
    void get_gradStress(Matrix<double, 1, 6> &, double, double, Matrix<double, 1, 6> &);
    DKT_OPT dkt_opt;

    Matrix<double, 9, 9> Km, Kb, Kmb;
    Matrix<double, 9, 1> Fm, Fb, Fmf, Fbf;

    void get_filteredmat(FilteredP &p_d, Matrix<double, 6, 1> &FNM_r, std::vector<double> ri, std::vector<double> si, std::vector<double> wi);
    double set_Bmat(Vector3d x_, Vector3d y_, double r, double s);
    double set_Bmat(Vector3d x_, Vector3d y_, double r, double s, Matrix3d &globalXY0, Vector2d &globalgptsxy);
    // void get_stress_strain(VectorXd & r_pos, VectorXd & s_pos, VectorXd & stress);

  private:
    Matrix<double, 2, 3> &xycoord;
    Material_ABD &Mater_e;
    // Matrix<double, 6, 1> & FNM_r;
    // FilteredP & p_d;

    // Bmatrices
    MatrixXd N_kappa;      //= MatrixXd::Zero(3,9);
    MatrixXd N_psi, N_eta; //= MatrixXd::Zero(3,9);
};

#endif
