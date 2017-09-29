#include "./EICR.h"
#include "./../M2DO_LSM/include/M2DO_LSM.h"
#include <fstream>

namespace LSM = M2DO_LSM ;
using namespace Eigen;

int main(int argc, char *argv[]){
    Matrix3d eye3 = Matrix3d::Identity(3,3);

    // mesh generation
    const int exy[2] = {80, 40};
    const double Lxy[2] = {80., 40.};
    const double h = 1;

    int npe = 3, dpn = 6, dpe = 18;
    FEAMesh feaMesh(Lxy, exy, true);
    uint nNODE = feaMesh.NODE.rows();
    uint nELEM = feaMesh.ELEM.rows();
    for (uint ee = 0; ee < nELEM; ++ee){
        feaMesh.areafraction[ee] = 1;
    }
    
    // Material
    double E = 1000, v = 0.3, rho = 1.0;
    MatrixXd Cijkl(3,3);//, Amat(3,3), Dmat(3,3), Bmat(3,3);
	Material_ABD material0;
    Cijkl << 1, v, 0, v, 1, 0, 0, 0, (1-v)*0.5;
    Cijkl *= E/(1-std::pow(v,2));
    material0.Amat = Cijkl * h;
    material0.Dmat = Cijkl * std::pow(h,3)/12.0;
    material0.Bmat.setZero(3,3);

	std::vector<Material_ABD> material;
    
    for (unsigned int ee = 0; ee < nELEM; ++ee){
        material.push_back(material0);        
    }


    // load file
    // const char* f1 = argv[1]; // configuration file
    const char* f2 = argv[2]; // result file
    const char* f1 = "./FEAConfiguration.txt"; // configuration file
    // const char* f2 = "./FEA_solution.txt"; // result file
    myreadTXT read(f1, f2);

    read.read_config();
    feaMesh.BCid = read.feaConfig.BCid;
    Force force;
    force.NM = read.feaConfig.Force_NM;
    force.fix = read.feaConfig.Force_fix;

    read.read_results();
    MatrixXd GU_u = read.feaField.GU_u;
    VectorXd GU_Rv = read.feaField.GU_Rv;
    // VectorXd p_Adjoint = read.feaField.p_Adjoint;
    VectorXd p_Adjoint;
    double lambdaR = read.feaField.lambdaR;
    force.NM *= lambdaR ;
    force.fix *= lambdaR ;

    OPTION option;
    option.saveflag = 0;
    option.initper = 1;
    p_Adjoint = f_nlgeom(material, force, feaMesh, option, GU_u, GU_Rv);

    // declare sensitivity class
    Sensitivity sens(feaMesh, material, GU_u, GU_Rv, p_Adjoint);
    sens.ComputeComplianceSensitivities(); 
    
    // Level-set mesh
    LSM::Mesh lsmMesh(Lxy[0], Lxy[1], false);
    std::vector<LSM::Hole> holes;
    LSM::LevelSet levelSet(lsmMesh, holes, 0.5, 6, false);

    levelSet.reinitialise();
    LSM::Boundary boundary(levelSet);

    boundary.discretise();

    MatrixXd bptsSens(boundary.points.size(),3);
    for (uint ii = 0; ii < boundary.points.size(); ii++){
        std::vector<double> bpts(2,0);
        bpts[0] = boundary.points[ii].coord.x;
        bpts[1] = boundary.points[ii].coord.y;
        double tmp_sensP = sens.ComputeBoundaryPointSensitivity(bpts, 5, 5);
        bptsSens(ii, 0) = bpts[0];
        bptsSens(ii, 1) = bpts[1];
        bptsSens(ii, 2) = tmp_sensP;
    }

    std::ofstream file;
    file.open("sens.data");
    file << bptsSens;
    file.close();

    file.open("Gpts.data");
    file << sens.Gpts;
    file.close();

    file.open("GptsSens.data");
    file << sens.GptsSensitivities;
    file.close();
            
}

