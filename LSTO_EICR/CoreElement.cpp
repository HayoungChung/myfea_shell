#include "./CoreElement.h"
// this is an example
using namespace Eigen;

CoreElement::CoreElement(Matrix<double, 2, 3> &xycoord_, Material_ABD &Mater_e_) : xycoord(xycoord_), Mater_e(Mater_e_)
{

    // DKT_OPT matrices
    double nu = std::min(Mater_e.Amat(1, 0) / Mater_e.Amat(0, 1), Mater_e.Amat(1, 0) / Mater_e.Amat(0, 0));
    dkt_opt = f_DKT_OPT(xycoord, nu);

    Km.fill(0), Kb.fill(0), Kmb.fill(0);
    Fm.fill(0), Fb.fill(0), Fmf.fill(0), Fbf.fill(0);

    N_kappa.resize(3, 9); // = MatrixXd::Zero(3,9);
    N_psi.resize(3, 9);   // = MatrixXd::Zero(3,9);
    N_eta.resize(3, 9);   // = MatrixXd::Zero(3,9);
                          //     MatrixXd N_kappa; //= MatrixXd::Zero(3,9);
                          // MatrixXd N_psi, N_eta; //= MatrixXd::Zero(3,9);
}

Vector2d CoreElement::get_StressStrain(FilteredP &p_d, Matrix<double, 6, 1> &FNM_r, double r, double s, Matrix<double, 1, 6> &stress, Matrix<double, 1, 6> &strain)
{
    // computes stress and strain
    Matrix<double, 9, 1> u_m = p_d.membrane, u_b = p_d.bending;
    Vector3d Nph = FNM_r.topRows(3), Mph = FNM_r.bottomRows(3);
    Vector3d x_ = xycoord.row(0), y_ = xycoord.row(1);
    Matrix3d Amat = Mater_e.Amat, Dmat = Mater_e.Dmat, Bmat = Mater_e.Bmat;

    Vector2d globalGpts;
    double Area = set_Bmat(x_, y_, r, s);

    // strain << N_eta*u_m, N_kappa*u_b;
    // stress << Amat*N_eta*u_m + Bmat*N_kappa*u_b, Bmat*N_eta*u_m + Dmat*N_kappa*u_b;
    strain.leftCols(3) = (N_eta * u_m).transpose();
    strain.rightCols(3) = (N_kappa * u_b).transpose();
    stress.leftCols(3) = (Amat * N_eta * u_m + Bmat * N_kappa * u_b).transpose();
    stress.rightCols(3) = (Bmat * N_eta * u_m + Dmat * N_kappa * u_b).transpose();

    return globalGpts;
}

Vector2d CoreElement::get_StressStrain(FilteredP &p_d, Matrix<double, 6, 1> &FNM_r, double r, double s, Matrix<double, 1, 6> &stress, Matrix<double, 1, 6> &strain, Matrix3d &globalXY0)
{ // computes stress and strain
    Matrix<double, 9, 1> u_m = p_d.membrane, u_b = p_d.bending;
    Vector3d Nph = FNM_r.topRows(3), Mph = FNM_r.bottomRows(3);
    Vector3d x_ = xycoord.row(0), y_ = xycoord.row(1);
    Matrix3d Amat = Mater_e.Amat, Dmat = Mater_e.Dmat, Bmat = Mater_e.Bmat;

    Vector2d globalGpts;
    double Area = set_Bmat(x_, y_, r, s, globalXY0, globalGpts);

    // strain << N_eta*u_m, N_kappa*u_b;
    // stress << Amat*N_eta*u_m + Bmat*N_kappa*u_b, Bmat*N_eta*u_m + Dmat*N_kappa*u_b;
    strain.leftCols(3) = (N_eta * u_m).transpose();
    strain.rightCols(3) = (N_kappa * u_b).transpose();
    stress.leftCols(3) = (Amat * N_eta * u_m + Bmat * N_kappa * u_b).transpose();
    stress.rightCols(3) = (Bmat * N_eta * u_m + Dmat * N_kappa * u_b).transpose();

    return globalGpts;
}

void CoreElement::get_gradStress(Matrix<double, 1, 6> &stress, double r, double s, Matrix<double, 1, 6> &GradStress)
{ // computes stress and strain
    Vector3d x_ = xycoord.row(0), y_ = xycoord.row(1);
    Matrix3d Amat = Mater_e.Amat, Dmat = Mater_e.Dmat, Bmat = Mater_e.Bmat;

    double Area = set_Bmat(x_, y_, r, s);

    Vector3d Nij_j, Mij_j;
    Nij_j.fill(0);
    Mij_j.fill(0);

    // TOFIX: 귀찮아서 걍 nodal force를 average 하기로 함. 나중에 수정할 것
    MatrixXd N_trans_, M_trans_; // 3 x 9
    N_trans_ = N_eta.transpose();
    M_trans_ = N_kappa.transpose();

    MatrixXd n_ = N_trans_ * stress.leftCols(3).transpose();
    MatrixXd m_ = M_trans_ * stress.rightCols(3).transpose(); // 9 x 1

    for (int mmm = 0; mmm < 3; ++mmm)
    {
        GradStress(mmm) = (n_(mmm) + n_(3 + mmm) + n_(6 + mmm)) / 3.0;
        GradStress(3 + mmm) = (m_(mmm) + m_(3 + mmm) + m_(6 + mmm)) / 3.0;
    }
}

void CoreElement::get_filteredmat(const FilteredP &p_d, const Matrix<double, 6, 1> &FNM_r, const std::vector<double> ri, const std::vector<double> si, const std::vector<double> wi)
{
    Matrix<double, 9, 1> u_m = p_d.membrane, u_b = p_d.bending;
    Vector3d Nph = FNM_r.topRows(3), Mph = FNM_r.bottomRows(3);
    Vector3d x_ = xycoord.row(0), y_ = xycoord.row(1);
    Matrix3d Amat = Mater_e.Amat, Dmat = Mater_e.Dmat, Bmat = Mater_e.Bmat;

    for (unsigned int gg = 0; gg < ri.size(); ++gg)
    {
        double r = ri[gg], s = si[gg], w = wi[gg];

        double Area = set_Bmat(x_, y_, r, s);

        // Outputs
        Km.noalias() += (N_eta.transpose() * Amat * N_eta) * w * Area;
        Kb.noalias() += (N_kappa.transpose() * Dmat * N_kappa) * w * Area;
        Kmb.noalias() += (N_eta.transpose() * Bmat * N_kappa) * w * Area;
        // Km = Km + (N_psi.transpose()*Gmat(9)*N_psi)*w*Area; // if phiflag

        Fm.noalias() += (N_eta.transpose() * Amat * N_eta * u_m) * w * Area;
        Fb.noalias() += (N_kappa.transpose() * Dmat * N_kappa * u_b) * w * Area;

        Fmf.noalias() += (N_eta.transpose() * Nph) * w * Area;
        Fbf.noalias() += (N_kappa.transpose() * Mph) * w * Area;
    }
}

std::vector<GaussProp> CoreElement::get_GaussProp(FilteredP &p_d, Matrix<double, 6, 1> &FNM_r, int ngp)
{
    std::vector<GaussProp> gaussProps;
    Matrix3d Amat = Mater_e.Amat, Dmat = Mater_e.Dmat, Bmat = Mater_e.Bmat;
    Matrix<double, 9, 1> u_m = p_d.membrane, u_b = p_d.bending;
    Vector3d Nph = FNM_r.topRows(3), Mph = FNM_r.bottomRows(3);
    Vector3d x_ = xycoord.row(0), y_ = xycoord.row(1);

    gaussProps.resize(ngp);
    std::vector<double> ri(ngp), si(ngp), wi(ngp);
    if (ngp == 3)
    {
        ri = {0, 0.5, 0.5}, si = {0.5, 0, 0.5}, wi = {0.333333333, 0.333333333, 0.333333333};
    }
    else if (ngp == 1)
    {
        ri[0] = 1.0 / 3, si[0] = 1.0 / 3, wi[0] = 1;
    }
    
    for (unsigned int gg = 0; gg < ngp; ++gg)
    {
        double r = ri[gg], s = si[gg], w = wi[gg];
        double Area = set_Bmat(x_, y_, r, s);
        // compute gaussProp
        gaussProps[gg].eps0 = N_eta * u_m;
        gaussProps[gg].kappa = N_kappa * u_b;
        gaussProps[gg].resN = Amat * gaussProps[gg].eps0 + Bmat * gaussProps[gg].kappa;
        gaussProps[gg].resM = Dmat * gaussProps[gg].kappa + Bmat * gaussProps[gg].eps0;
        gaussProps[gg].sensLin = gaussProps[gg].resN.dot(gaussProps[gg].eps0) +
                              gaussProps[gg].resM.dot(gaussProps[gg].kappa);
    }
    return gaussProps;
}

double CoreElement::set_Bmat(Vector3d x_, Vector3d y_, double r, double s, Matrix3d &globalXY, Vector2d &GptsXY)
{
    // // Jacobian (dx/dXi)
    // Matrix<double, 2, 2> Jmat, invJ;
    // Jmat << x_(1) - x_(0), x_(2) - x_(0), y_(1) - y_(0), y_(2) - y_(0);
    // double Area = Jmat.determinant() / 2;
    // invJ = Jmat.inverse();

    // //  areal coordinates
    // // MatrixXd L = MatrixXd(1,3); //(1,3);
    // RowVector3d L;
    // Matrix<double, 1, 6> Lsq;
    // Matrix<double, 2, 6> d_Lsq_dxi;
    // // MatrixXd Lsq(1,6), d_Lsq_dxi(2,6);
    // L << r, s, 1 - r - s;
    // // shape function: quadratic monomials of areal coordinates
    // Lsq << r * r, s * s, (1 - r - s) * (1 - r - s), s * (1 - r - s), (1 - r - s) * r, r * s;
    // // derivatives (wrt xi,eta) of quadratic monomials of areal coord.
    // d_Lsq_dxi << -2 * r, 2 * s, 0, (1 - r - s), -(1 - r - s), r - s,
    //     -2 * r, 0, 2 * (1 - r - s), s, r - (1 - r - s), -s;

    // // Bmatrices ================================================
    // MatrixXd d_N_thetax_dxi, d_N_thetax_dx, d_N_thetay_dxi, d_N_thetay_dx;
    // // rotations
    // // N_thetax = Lsq*dkt_opt.M_thetax;
    // // N_thetay = Lsq*dkt_opt.M_thetay;
    // // shape functions of derivative of thetax wrt xi,eta
    // d_N_thetax_dxi = d_Lsq_dxi * dkt_opt.M_thetax;
    // // shape functions of derivative of thetax wrt x,y
    // d_N_thetax_dx = invJ.transpose() * d_N_thetax_dxi;

    // // shape functions of derivative of thetay wrt xi,eta
    // d_N_thetay_dxi = d_Lsq_dxi * dkt_opt.M_thetay;
    // // shape functions of derivative of thetay wrt x,y
    // d_N_thetay_dx = invJ.transpose() * d_N_thetay_dxi;

    // // shape functions of curvatures kappax, kappay, kappaxy
    // N_kappa.row(0) = d_N_thetay_dx.row(0);
    // N_kappa.row(1) = -d_N_thetax_dx.row(1);
    // N_kappa.row(2) = d_N_thetay_dx.row(1) - d_N_thetax_dx.row(0);

    // // OPT shape function --------------
    // /// shape functions of in-plane displacements and drilling rotation
    // //MatrixXd N_psi, N_eta = MatrixXd::Zero(3,9);
    // N_psi = L * dkt_opt.M_psi;

    // // shape functions of in-plane deformations etax, etay, gamxy
    // N_eta << L * dkt_opt.M_etax, L * dkt_opt.M_etay, L * dkt_opt.M_gamxy;
    // // ========================================================
    Matrix<double, 2, 2> Jmat;
    Jmat << x_(1) - x_(0), x_(2) - x_(0), y_(1) - y_(0), y_(2) - y_(0);
    double Area = Jmat.determinant() / 2;
    
    set_Bmat(x_, y_, r, s);
    // gpt coordinates ========================================
    RowVector3d L;
    L << r, s, 1 - r - s;
    
    GptsXY << L.dot(globalXY.row(0)), L.dot(globalXY.row(1));
    return Area;
    // ========================================================
}

double CoreElement::set_Bmat(Vector3d x_, Vector3d y_, double r, double s)
{
    // Jacobian (dx/dXi)
    Matrix<double, 2, 2> Jmat, invJ;
    Jmat << x_(1) - x_(0), x_(2) - x_(0), y_(1) - y_(0), y_(2) - y_(0);
    double Area = Jmat.determinant() / 2;
    invJ = Jmat.inverse();

    //  areal coordinates
    // MatrixXd L = MatrixXd(1,3); //(1,3);
    RowVector3d L;
    Matrix<double, 1, 6> Lsq;
    Matrix<double, 2, 6> d_Lsq_dxi;
    // MatrixXd Lsq(1,6), d_Lsq_dxi(2,6);
    L << r, s, 1 - r - s;
    // shape function: quadratic monomials of areal coordinates
    Lsq << r * r, s * s, (1 - r - s) * (1 - r - s), s * (1 - r - s), (1 - r - s) * r, r * s;
    // derivatives (wrt xi,eta) of quadratic monomials of areal coord.
    d_Lsq_dxi << -2 * r, 2 * s, 0, (1 - r - s), -(1 - r - s), r - s,
        -2 * r, 0, 2 * (1 - r - s), s, r - (1 - r - s), -s;

    // Bmatrices ================================================
    MatrixXd d_N_thetax_dxi, d_N_thetax_dx, d_N_thetay_dxi, d_N_thetay_dx;
    // rotations
    // N_thetax = Lsq*dkt_opt.M_thetax;
    // N_thetay = Lsq*dkt_opt.M_thetay;
    // shape functions of derivative of thetax wrt xi,eta
    d_N_thetax_dxi = d_Lsq_dxi * dkt_opt.M_thetax;
    // shape functions of derivative of thetax wrt x,y
    d_N_thetax_dx = invJ.transpose() * d_N_thetax_dxi;

    // shape functions of derivative of thetay wrt xi,eta
    d_N_thetay_dxi = d_Lsq_dxi * dkt_opt.M_thetay;
    // shape functions of derivative of thetay wrt x,y
    d_N_thetay_dx = invJ.transpose() * d_N_thetay_dxi;

    // shape functions of curvatures kappax, kappay, kappaxy
    N_kappa.row(0) = d_N_thetay_dx.row(0);
    N_kappa.row(1) = -d_N_thetax_dx.row(1);
    N_kappa.row(2) = d_N_thetay_dx.row(1) - d_N_thetax_dx.row(0);

    // OPT shape function --------------
    /// shape functions of in-plane displacements and drilling rotation
    //MatrixXd N_psi, N_eta = MatrixXd::Zero(3,9);
    N_psi = L * dkt_opt.M_psi;

    // shape functions of in-plane deformations etax, etay, gamxy
    N_eta << L * dkt_opt.M_etax, L * dkt_opt.M_etay, L * dkt_opt.M_gamxy;
    // ========================================================

    // gpt coordinates ========================================
    return Area;
    // ========================================================
}
