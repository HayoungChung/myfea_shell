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
    
    std::vector<int> get_dof(int direction, std::vector<int> & ID_nodes);
    std::vector<int> get_nodeID(double posX, double posY, double Xtol, double Ytol);
    void get_BCid(std::vector<int> Dofs);

    void set_Force(int direction, std::vector<int> & ID_nodes, double val, MatrixXd & force_fix); //fix
    void set_Force(std::vector<Material_ABD> & material, MatrixXd eps0, MatrixXd kappa0, MatrixXd & force_NM); // NM

    // comptue centeroids (for least-square analysis)
    MatrixXd Centeroids;

    void to_vtk(MatrixXd & u, const char* str);
    void to_vtk(MatrixXd & u);
    void to_vtk();
    // void Compute_Area();
    unsigned int dpn = 6, npe = 3, dpe = 18;
    bool isOverlaid;
    std::vector<double> areafraction;
    std::vector<double> ElemArea;

    private:
    
    void ComputeCentroids();
    void set_ElemArea();
    unsigned int nNODE, nELEM;
    const double* Lxy;
    const int* exy;
};

#endif


