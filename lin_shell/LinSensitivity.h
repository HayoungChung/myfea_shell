#ifndef LIN_SENS_H
#define LIN_SENS_H

#define PI 3.141592

#include "./FEA_hy.h"
#include "./lin_shell.h"

using namespace Eigen;

class LinSensitivity
{
  public:
    LinSensitivity(class FEAMesh&, std::vector<GptsCompl>&);
    ~LinSensitivity(){};

    std::vector<double> get_bptsCompl(MatrixXd &BptsCoord, double radius = 2);
    double ComputeBoundaryPointSensitivity2D(std::vector<double> &Pointxy, double Radius, unsigned int WeightFlag, double Tolerance = 0.001);
    
private:
    MatrixXd Gpts, GptsSensitivities; 
    std::vector<double> BptsSensitivities;

    std::vector<double> areafraction;
    std::vector<GptsCompl> gaussComplSens;
    FEAMesh& feaMesh;
    void leastsquares();
    unsigned int nGpts, nBpts;    

};

#endif
