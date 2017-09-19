#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "./linshell.h"

using namespace Eigen;

int main()
{
    Matrix3d eye3 = Matrix3d::Identity();

    int npe = 3, dpn = 6, dpe = 18;
    bool isOverlaid = true;

    // Mesh generation
    const double Lxy[2] = {12, 2};
    const int exy[2] = {12, 2};
    const double h = 0.5; // the h/L < 0.1

    FEAMesh feaMesh(Lxy, exy, isOverlaid);
    const unsigned int nELEM = feaMesh.ELEM.rows();
    const unsigned int nNODE = feaMesh.NODE.rows();
    const unsigned nDOF = nNODE * dpn;

    // Material

    const double E = 1.2e6, v = 0.3;
    MatrixXd Cijkl(3, 3), Amat(3, 3), Dmat(3, 3), Bmat(3, 3);
    Cijkl << 1., v, 0., v, 1, 0., 0., 0., 0.5 * (1 - v);
    Cijkl *= E / (1 - v * v);
    Amat = Cijkl * h;
    Dmat = Cijkl * std::pow(h, 3) / 12.0;
    Bmat.setZero(3, 3);

    Material_ABD material0;
    material0.Amat = Amat;
    material0.Dmat = Dmat;
    material0.Bmat = Bmat;

    std::vector<Material_ABD> material;

    for (int ee = 0; ee < nELEM; ++ee)
    {
        material.push_back(material0);
    }

    // BC
    double posX = 0, posY = 0, Xtol = 1e-3, Ytol = 1e3; // Dirichlet.
    std::vector<int> Xlo = feaMesh.get_nodeID(posX, posY, Xtol, Ytol);
    std::vector<int> BCtmp;
    for (unsigned int dd = 0; dd < feaMesh.dpn; dd++)
    {
        std::vector<int> tmp = feaMesh.get_dof(dd, Xlo);
        for (unsigned int mm = 0; mm < tmp.size(); mm++)
        {
            BCtmp.push_back(tmp[mm]);
        }
    }

    feaMesh.get_BCid(BCtmp);
    /*
    std::cout << "NODE: \n"
              << feaMesh.NODE << "\n============" << std::endl;
    std::cout << "ELEM: \n"
              << feaMesh.ELEM << "\n============" << std::endl;
    std::cout << "BCid: \n"
              << feaMesh.BCid << "\n============" << std::endl;
    */
    MatrixXd Force_NM(nELEM, feaMesh.dpn);
    MatrixXd Force_Fix = MatrixXd::Zero(nNODE, feaMesh.dpn);
    Force_NM.fill(0.0);

    std::vector<int> Xtip = feaMesh.get_nodeID(Lxy[0], Lxy[1], 1e-3, 1e-3);

    Force_Fix(Xtip[0], 2) = -10;

    Force force;
    force.NM = Force_NM;
    force.fix = Force_Fix;

    SparseMatrix<double> GKT(nDOF, nDOF);
    VectorXd Res(nDOF);

    Res = VectorXd::Zero(nDOF);
    f_lin_shell(feaMesh, material, force, GKT, Res);

    SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> spQR;
    GKT.makeCompressed();
    spQR.compute(GKT);

    MatrixXd u;
    MatrixXd u6;

    u = spQR.solve(Res);
    u6 = Map<MatrixXd>(u.data(), 6, nNODE).transpose();

    // Post-process
    feaMesh.to_vtk(u6);
    std::cout << "(" << u.rows() << ", " << u.cols() << ")" << std::endl;
    // std::cout << u << std::endl;

    std::ofstream Dfile("disp_test.txt");
    std::ofstream Ffile("force_test.txt");
    if (Dfile.is_open())
    {
        Dfile << u6 << std::endl;
    }
    if (Ffile.is_open())
    {
        Ffile << force.fix << std::endl;
        Ffile << "Res: \n";
        Ffile << Res << std::endl; // TOFIX: too many nans
    }
}
