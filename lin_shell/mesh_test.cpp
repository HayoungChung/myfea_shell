#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <numeric>

#include "./FEA_hy.h"
#define PI 3.141592

using namespace std;
using namespace Eigen;

int main()
{
    const double max_theta = PI / 2; 
    const double Length = 80.;
    const double R0 = Length/max_theta;

    const double Lxy[2] = {max_theta, 40.0 }; // angle[radian], y dimension
    const int exy[2] = { 80, 40 };
    const bool isOverlaid = true;
    FEAMesh feaMesh(Lxy, exy, isOverlaid);

    // generate to radian
    MatrixXd Rtheta = feaMesh.NODE;
    for (unsigned int nn = 0; nn < feaMesh.NODE.rows(); nn++)
    {
        feaMesh.NODE(nn,0) = R0 * sin(Rtheta(nn,0));
        feaMesh.NODE(nn,1) = Rtheta(nn,1);
        feaMesh.NODE(nn,2) = R0 * cos(Rtheta(nn,0));
    }
    feaMesh.to_vtk();

    std::ofstream eFile("elem.txt");
    eFile << feaMesh.ELEM << std::endl;
    eFile.close();

    std::ofstream nFile("node.txt");
    nFile << feaMesh.NODE << std::endl;
    nFile.close();

    return 0;
}