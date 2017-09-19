
#ifndef FEA_HY_H
#define FEA_HY_H

#include <iostream>
#include <vector>
// #include "./EICR.h"
#include <cmath>
#include <algorithm>
#include "./../../eigen3/Eigen/Dense"
#include "./../../eigen3/Eigen/Sparse"

using namespace Eigen;

struct Material_ABD{
    Matrix<double,3,3> Amat, Dmat, Bmat;
};

struct Force{
    MatrixXd NM;
    MatrixXd fix;
};

class FEAMesh{
    public:
    MatrixXd NODE;
    MatrixXi ELEM; 
    VectorXi BCid;
    FEAMesh(const double* _Lxy,const int* _exy, bool isOverlaid = false);
    
    void get_Mesh(int nVertices = 3);
    // static void get_N(const double r, const double s, const int order = 2, MatrixXd& Ni);
    // static void get_N_rs(const double r, const double s, const int order = 2, MatrixXd& Ni_rs);
    // static void get_gpts(int IntegrationPoints =2, MatrixXd& ri, MatrixXd& si, MatrixXd& wi);
    std::vector<int> get_dof(int direction, std::vector<int> & ID_nodes);
    std::vector<int> get_nodeID(double posX, double posY, double Xtol, double Ytol);
    void get_BCid(std::vector<int> Dofs);

    void set_Force(int direction, std::vector<int> & ID_nodes, double val, MatrixXd & force_fix); //fix
    void set_Force(std::vector<Material_ABD> & material, MatrixXd eps0, MatrixXd kappa0, MatrixXd & force_NM); // NM

    void to_vtk(MatrixXd & u);
    // void Compute_Area();
    unsigned int dpn = 6, npe = 3, dpe = 18;
    bool isOverlaid;

    private:
    std::vector<double> AreaFraction;
    
    int nNODE, nELEM;
    const double* Lxy;
    const int* exy;

};

#endif 


