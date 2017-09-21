#include "./f_lin_shell.h"
#define PI 3.14159265359
void f_lin_shell(class FEAMesh &feaMesh, std::vector<Material_ABD> &material,
                 struct Force &force, SparseMatrix<double> &sGKT, VectorXd &Res)
{

    Matrix3d zero3 = Matrix3d::Zero();

    // int dpn = 6, dpe = 18, npe = 3;
    const int dpn = feaMesh.dpn, dpe = feaMesh.dpe, npe = feaMesh.npe;
    const bool isOverlaid = feaMesh.isOverlaid;
    VectorXi BCid = feaMesh.BCid;

    // parsing inputs
    MatrixXd FNM = force.NM, Ffix = force.fix;
    MatrixXd NODE = feaMesh.NODE;
    MatrixXi ELEM = feaMesh.ELEM;
    const int nNODE = NODE.rows(), nELEM = ELEM.rows();
    const int nDOF = nNODE * dpn;

    // assert(FNM.rows() == nELEM, "size of F_nm does not match to nELEM")

    std::vector<Triplet<double>> GKT;

    GKT.reserve(nELEM * dpe * dpe);

    // MatrixXd GKT_mat = MatrixXd::Zero(nDOF,nDOF);
    VectorXd GFi = VectorXd::Zero(nDOF);
    MatrixXd GF_dof = Ffix.transpose();
    MatrixXd GF_dof_map = Map<VectorXd>(GF_dof.data(), nDOF);

    // ids for membrane and bending
    Matrix<int, 9, 1> mdof, bdof;
    mdof << 0, 1, 5, 6, 7, 11, 12, 13, 17;
    bdof << 2, 3, 4, 8, 9, 10, 14, 15, 16;

    // Overlaid case

    MatrixXi elem_id0(3, 1), elem_order(1, 3);
    elem_order << 0, 1, 2;
    // Vector3i elem_id0;
    // Matrix<int, 1, 3> elem_order; elem_order << 0, 1, 2;
    if (isOverlaid == true)
    {
        elem_id0 = MatrixXi::Zero(4, 1);
        elem_order = MatrixXi::Zero(4, 3);
        elem_order << 0, 1, 2, 1, 2, 3, 2, 3, 0, 3, 0, 1;
        // Matrix<int, 4, 1> elem_id0;
        // Matrix<int, 4, 3> elem_order; elem_order << 0,1,2, 1,2,3, 2,3,0, 3,0,1;
    }

    for (int ee = 0; ee < nELEM; ++ee)
    {
        elem_id0 = ELEM.row(ee);

        Matrix<double, 3, 3> Amat, Dmat, Bmat;
        Amat = material[ee].Amat;
        Dmat = material[ee].Dmat;
        Bmat = material[ee].Bmat;

        for (unsigned int kk = 0; kk < elem_order.rows(); ++kk)
        {
            Vector3i elem_id;
            for (unsigned int ppp = 0; ppp < 3; ++ppp)
            {
                elem_id(ppp) = elem_id0(elem_order(kk, ppp));
            }

            Matrix<double, 3, 3> X;
            for (int mmm = 0; mmm < 3; ++mmm){
                X.col(mmm) = NODE.row(elem_id(mmm)).transpose();
            }

            // Transpose matrix [e1; e2; e3]
            // x axis of T0 is x1-x2 line
            // TOFIX: x axis of Tmat (material) should be computed analytically
            Matrix<double, 3, 3> T0, Tmat;
            T0.row(0) = X.col(1) - X.col(0);
            T0.row(0) /= T0.row(0).norm();
            T0.row(2) = T0.row(0).cross(X.col(2) - X.col(0));
            T0.row(2) /= T0.row(2).norm();
            T0.row(1) = T0.row(2).cross(T0.row(0));

            // TOFIX: X must be flat currently
            Matrix3d X_R;
            X_R = T0 * X;
            Matrix<double, 2, 3> xycoord = X_R.topRows(2);

            // Material of local element
            // coordinate trasfer from global to initial (T0)
            Matrix<double, 3, 3> Q;
            Q << T0(0, 0) * T0(0, 0), T0(0, 1) * T0(0, 1), T0(0, 0) * T0(0, 1),
                T0(1, 0) * T0(1, 0), T0(1, 1) * T0(1, 1), T0(1, 0) * T0(1, 1),
                2 * T0(0, 0) * T0(1, 0), 2 * T0(0, 1) * T0(1, 1), T0(0, 0) * T0(1, 1) + T0(0, 1) * T0(1, 0);

            Matrix<double, 6, 1> FNM_r;
            FNM_r.topRows(3) = (FNM.block(ee, 0, 1, 3) * Q.transpose()).transpose();
            FNM_r.bottomRows(3) = (FNM.block(ee, 3, 1, 3) * Q.transpose()).transpose();
            Material_ABD Mater_e;
            Mater_e.Amat = Q * Amat * Q.transpose();
            Mater_e.Dmat = Q * Dmat * Q.transpose();
            Mater_e.Bmat = Q * Bmat * Q.transpose();

            CoreElement coreelem = f_core_element(xycoord, Mater_e, FNM_r);

            MatrixXd K_el = MatrixXd::Zero(18, 18);
            VectorXd f_el = VectorXd::Zero(18), f_el_f = VectorXd::Zero(18);

            if (isOverlaid == false)
            {
                for (int mmm = 0; mmm < 9; ++mmm)
                {
                    for (int nnn = 0; nnn < 9; ++nnn)
                    {
                        K_el(mdof(mmm), mdof(nnn)) += coreelem.Km(mmm, nnn);
                        K_el(bdof(mmm), bdof(nnn)) += coreelem.Kb(mmm, nnn);
                        K_el(mdof(mmm), bdof(nnn)) += coreelem.Kmb(mmm, nnn);
                        K_el(bdof(mmm), mdof(nnn)) += coreelem.Kmb(nnn, mmm);
                    }
                    f_el(mdof(mmm)) += coreelem.Fm(mmm);
                    f_el(bdof(mmm)) += coreelem.Fb(mmm);

                    f_el_f(mdof(mmm)) += coreelem.Fmf(mmm);
                    f_el_f(bdof(mmm)) += coreelem.Fbf(mmm);
                }
            }
            else
            {
                for (int mmm = 0; mmm < 9; ++mmm)
                {
                    for (int nnn = 0; nnn < 9; ++nnn)
                    {
                        K_el(mdof(mmm), mdof(nnn)) += coreelem.Km(mmm, nnn) / 2;
                        K_el(bdof(mmm), bdof(nnn)) += coreelem.Kb(mmm, nnn) / 2;
                        K_el(mdof(mmm), bdof(nnn)) += coreelem.Kmb(mmm, nnn) / 2;
                        K_el(bdof(mmm), mdof(nnn)) += coreelem.Kmb(nnn, mmm) / 2;
                    }
                    f_el(mdof(mmm)) += coreelem.Fm(mmm) / 2;
                    f_el(bdof(mmm)) += coreelem.Fb(mmm) / 2;

                    f_el_f(mdof(mmm)) += coreelem.Fmf(mmm) / 2;
                    f_el_f(bdof(mmm)) += coreelem.Fbf(mmm) / 2;
                }
            }

            MatrixXd LKT;
            VectorXd Fi;

            MatrixXd Te = MatrixXd::Zero(18, 18);
            Te << T0, zero3, zero3, zero3, zero3, zero3,
                zero3, T0, zero3, zero3, zero3, zero3,
                zero3, zero3, T0, zero3, zero3, zero3,
                zero3, zero3, zero3, T0, zero3, zero3,
                zero3, zero3, zero3, zero3, T0, zero3,
                zero3, zero3, zero3, zero3, zero3, T0;

            LKT = Te.transpose() * K_el * Te;
            Fi = Te.transpose() * (f_el - f_el_f);

            // assembly
            Matrix<int, 18, 1> dof_sum(18);
            dof_sum << VectorXi::LinSpaced(6, elem_id(0) * 6, elem_id(0) * 6 + 5),
                VectorXi::LinSpaced(6, elem_id(1) * 6, elem_id(1) * 6 + 5),
                VectorXi::LinSpaced(6, elem_id(2) * 6, elem_id(2) * 6 + 5);

            for (int mmm = 0; mmm < 18; ++mmm)
            {
                if (BCid.cwiseEqual(dof_sum(mmm)).any() > 0)
                {
                    GF_dof_map(dof_sum(mmm)) = 0;
                    continue;
                }
                GFi(dof_sum(mmm)) += Fi(mmm);
                for (int nnn = 0; nnn < 18; ++nnn)
                {
                    if (BCid.cwiseEqual(dof_sum(nnn)).any() > 0)
                        continue;
                    GKT.push_back(Triplet<double>(dof_sum(mmm), dof_sum(nnn), LKT(mmm, nnn)));
                }
            }
        }
    }
    Res = GF_dof_map - GFi;
    sGKT.setFromTriplets(GKT.begin(), GKT.end());
}

double getArea(MatrixXd x)
{
    // calculate Area of triangular element
    MatrixXd x21 = (x.col(1) - x.col(0)), x31 = (x.col(2) - x.col(0));
    Vector3d crossx2x3;
    crossx2x3 << x21(1) * x31(2) - x21(2) * x31(1), x21(2) * x31(0) - x21(0) * x31(2), x21(0) * x31(1) - x21(1) * x31(0);
    double Area = 0.5 * crossx2x3.norm();
    return Area;
}
