#ifndef EICR_SENSITIVITY_H
#define EICR_SENSITIVITY_H

#include "./EICR.h"

class SensitivityAnalysis{
    public: 

    SensitivityAnalysis(FEAMesh & feaMesh_, MatrixXd & GU_u_, MatrixXd & GU_Rv_,
                std::vector<Material_ABD> & materials_, const unsigned int SpaceDimension_);

    protected:

    FEAMesh & feaMesh;
    MatrixXd & GU_u, GU_Rv;
    std::vector<Material_ABD> & material;
    const unsigned int SpaceDimension;

    MatrixXd IntegrationPoints;
    void CopmuteBoundaryPointSEnsitivities(std::vector<double> pos, std::vector <double> Sensitivities, std::vector <double> PointSensitivities, 
        unsigned int nDual, double Radius, unsigned int WeightFlag, double Tolerance);
};

class EICR_Sensitivities : public SensitivityAnalysis{
    public:

    EICR_Sensitivities(FEAMesh &, MatrixXd &, MatrixXd &, std::vector<Material_ABD>&, double AllowedAreaFraction = 0.001);
    void Compliance(MatrixXd & BoundaryPoints, MatrixXd& BoundarySensitivities, MatrixXd& Weights, double Radius, unsigned int WeightFlag, double Tolerance);

    private:

    AuxMat & auxmat;
    
};

#endif