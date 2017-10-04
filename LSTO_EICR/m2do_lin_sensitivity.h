#ifndef M2DO_LIN_SENSITIVITY_H
#define M2DO_LIN_SENSITIVITY_H

#include "./FEA_hy.h"
#include "./EICR_shell.h"
// #include "./Sensitivity.h"

// using namespace std ;
using namespace Eigen;

// -------------------------------------------------------------------------------- //
// LEAST SQUARES EXTERN FUNCTION
extern "C" int dgels_(const char *trans, int *m, int *n, int *nrhs,
        double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);

// -------------------------------------------------------------------------------- //
// SENSITIVITIES CLASS
class M2DOSensitivity {

    private:
        //

    public:
        // Sensitivities at prescribed point (is a vector so we can compute at Gauss points when handy)
        std::vector<double> dx;
        std::vector<double> dx1;
        std::vector<double> dx2;
        double dx_average; // Average
        // Global coordinates of the sensitivities points
        std::vector< std::vector<double> > dxcoord;

        // Is element excluded from the optimisation?
        bool isExcluded;

};

// -------------------------------------------------------------------------------- //
// STRESS, STRAIN AND ADJOINT STRAIN CLASS
class StressStrain {

    private:
        //

    public:
        // von Mises stress.
        std::vector<double> Tvm; // For all Gauss points.
        double Tvm_average; // Average for the element.

        // Strain.
        std::vector<double> strain; // For all Gauss points.

};

// -------------------------------------------------------------------------------- //
// LEAST SQUARES CLASS
class LeastSquares {
    public:
        // Distances from gauss point to boundary point
        double distGauss;
        // Area fraction at the Gauss point
        double areaFractionGauss;
        // Element number
        int elementNumber;
        // Node number (if applicable).
        int nodeNumber;
        // Gauss point local number
        int gaussPointNumber;
        // Coordinates of the Gauss point
        std::vector<double> coord;

};
/*
// -------------------------------------------------------------------------------- //
// SENSITIVITY Analysis CLASS
template <int spacedim, int dim, int order, class Mesh, class Physics, class Study, class LSM, int sens_order>
class SensitivityAnalysis {

    private:
        //

    public:
        // Properties:
        int material_id; // Material id number.
        // Vector of strains.
        std::vector<StressStrain> strains;
        // Vector of sensitivities.
        std::vector<Sensitivity> sensitivities;
        // Least squares information class.
        std::vector<LeastSquares> leastsquares;
        // Boundary Sensitivities.
        std::vector<double> boundarySens;

        // Attaching classes.
        Study & study;
        LSM & lsm;

        // For stress analysis.
        double objective; // Objective function value.
        // Maximum von Mises stress.
        double Tvm_max;

        // General methods:
        SensitivityAnalysis(Mesh &, Study &, LSM &); // Constructor
        void ExcludeElements (std::vector<int>);
        void ComputeSensitivitiesCoordinates(Mesh &);
        std::vector<double> ComputeSignedDistanceAtGaussPoints (Mesh &, int);

        // Compliance problem
        void ComputeComplianceSensitivities(Mesh &, Physics &, double);
        void ComputeThermomechanicalComplianceSensitivities(Mesh &, Physics &, MatrixXd, double, double, double);

        // Stress problem
        void ComputeStressShapeSensitivities(Mesh &, Physics &, HomogeneousDirichletBC, double, double);
        void ComputeStressShapeSensitivities_multiple_loads(Mesh &, Physics &, HomogeneousDirichletBC, double, double, double, double, double, int);
        void ComputeStressShapeSensitivitiesAllaire(Mesh &, Physics &, HomogeneousDirichletBC, double, double);

        // Eigenfrequency problems.
        void ComputePlateEigenfrequencySensitivities(Mesh &, Physics &, double, int);

        // Least squares functionalities.
        void ComputeBoundarySensitivities(Mesh &, Physics &, double, double, std::vector<double>);
        void ComputeBoundarySensitivities_stress(Mesh &, Physics &, double, double, std::vector<double>, double);
        double SolveLeastSquares(std::vector<LeastSquares>, std::vector<double>);
        double SolveLeastSquares_stress(std::vector<LeastSquares>, std::vector<double>, int);
        
        // Printing
        void WriteSensitivitiesTXT(int, int, int);
        void WriteStressTXT(int, int, int);
        void WriteSolutionTXT(Mesh &, int, int, int);

};
*/

class SensitivityAnalysis {    
public:
    SensitivityAnalysis(FEAMesh&, std::vector<GptsCompl>&);
    ~SensitivityAnalysis(){};
    void ComputeBoundarySensitivities(double radius, std::vector<double> bPoint, double areamin = 1e-2);
    
    std::vector<LeastSquares> leastsquares;
    std::vector<double> boundarySens;
    
private:
    FEAMesh& feaMesh;
    std::vector<GptsCompl> gptsCompls;
    void GptsCompl2sensitivities();
    double SolveLeastSquares(std::vector<LeastSquares>&, std::vector<double>&);
    
    int spacedim = 2;
    std::vector<M2DOSensitivity> sensitivities;
    

};

#endif
