#ifndef EICR_SENS_H
#define EICR_SENS_H

#define PI 3.141592

#include "./FEA_hy.h"
#include "./CoreElement.h"
#include "./EICR_shell.h"
#include "cblas.h"
#include <fstream>

using namespace Eigen;

struct GptsCompl
{
    double x, y, z;
    double sens;
};

class Sensitivity : public EICR_SHELL
{
  public:
    Sensitivity(class FEAMesh &feaMEsh, std::vector<Material_ABD> &material, MatrixXd &GU_u, VectorXd &GU_Rv, VectorXd &p_Adjoint);
    ~Sensitivity()
    {
    }

    // VectorXd BptsSensitivities;
    MatrixXd GptsSensitivities;

    //   void ComputeComplianceSensitivities(double areamin);
    void ComputeComplianceSensitivities(double Tolerance = 0.001); // Gausspoint sensitivities
    // void ComputeBoundarySensitivities(double radius, double * bPoints);
    // double ComputeBoundaryPointSensitivity(std::vector<double> & Pointxy, std::vector<double> & Sensitivities, double Radius, unsigned int WeightFlag, double Tolerance = 0.001);
    // double ComputeBoundaryPointSensitivity(std::vector<double> & Pointxy, double Radius, unsigned int WeightFlag, bool isLinear = false, double Tolerance = 0.001);
    // double compliance;
    MatrixXd Gpts; // global position of gausspoints
    void to_gptSens(bool isLinear = false);
    std::vector<GptsCompl> gptsSens;

  private:
    std::vector<double> ri, si, wi;
    FEAMesh &feaMesh;
    std::vector<Material_ABD> &material;
    MatrixXd GU_u, GU_R;
    VectorXd GU_Rv;

    VectorXd p_Adjoint;
    MatrixXd nabla_0; // Bmatrix of undeformed condition
    MatrixXd result_cr;
    MatrixXd strain_cr;

    MatrixXd pull_back(MatrixXd param_CR);

    // void LeastSquare(unsigned int N, unsigned int M, unsigned int nRhs, std::vector <double>& A, std::vector <double>& B);
};

#endif
/*
#ifdef __cplusplus
extern "C" {
#endif
	extern  int dgels_(const char*, unsigned int*, unsigned int*, unsigned int*, double*, unsigned int*, double*, unsigned int*, double*, unsigned int*, int*); 

#ifdef __cplusplus
}
#endif*/