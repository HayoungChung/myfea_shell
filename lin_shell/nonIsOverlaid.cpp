#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <numeric>

#include "./lin_shell.h"

using namespace Eigen;

int main()
{
    Matrix3d eye3 = Matrix3d::Identity();

    int npe = 3, dpn = 6, dpe = 18;
    bool isOverlaid = false;

    // Mesh generation
    const double Lxy[2] = {80., 40.};
    const int exy[2] = {80, 40};
    const double h = 2; // the h/L < 0.1

    FEAMesh feaMesh(Lxy, exy, isOverlaid);
    const unsigned int nELEM = feaMesh.ELEM.rows();
    const unsigned int nNODE = feaMesh.NODE.rows();
    const unsigned nDOF = nNODE * dpn;
    feaMesh.areafraction.assign(nELEM,1);
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

    std::vector<int> Xtip = feaMesh.get_nodeID(Lxy[0], Lxy[1]/2, 1e-3, 1e-3);

    for (int tt = 0; tt < Xtip.size(); tt++)
    {
        Force_Fix(Xtip[tt], 1) = -1;
    }

    Force force;
    force.NM = Force_NM;
    force.fix = Force_Fix;

    // Res = VectorXd::Zero(nDOF);

    // f_lin_shell(feaMesh, material, force, GKT, Res);
    LinShell lin_shell(feaMesh, material, force);
    lin_shell.compute();

    //SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> spsolver;
    ConjugateGradient<SparseMatrix<double>> spsolver;
    lin_shell.sGKT.makeCompressed(); // this is time-consuming!!!! (DONNO WHY)
    spsolver.compute(lin_shell.sGKT);

    MatrixXd u;
    MatrixXd u6;
    std::cout << "starting lin solve\n";
    
    u = spsolver.solve(lin_shell.Res);
    u6 = Map<MatrixXd>(u.data(), 6, nNODE).transpose();
    std::vector<GptsCompl> gptsCompl = lin_shell.get_GaussCompl(u6);
    std::cout << "ending lin solve\n";
    
    std::ofstream sfile("sens_test.txt");
    if (sfile.is_open())
    {
        for (int ii = 0; ii < gptsCompl.size() ; ++ii)
        {
            // sfile << "(" << gptsCompl[ii].x << ", " << gptsCompl[ii].y <<
            // ", " << gptsCompl[ii].z << ") " << gptsCompl[ii].sens << std::endl;
            sfile << gptsCompl[ii].x << " " << gptsCompl[ii].y <<
            " " << gptsCompl[ii].z << " " << gptsCompl[ii].sens << std::endl;
        }
    }

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
        Ffile << lin_shell.Res << std::endl; 
    }
    std::cout << "end" << std::endl;
    return 0;
}
