// this is a replica of prev. cantilever example, 
// infor 2D, triangular shell element

// #include "../../M2DO_FEA/include/M2DO_FEA.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>

#include "./../M2DO_LSM/include/M2DO_LSM.h"
#include "./EICR.h"
#include "./read_FEA_txt.h"
// #define PI 3.14159265359

using namespace std ;
using namespace Eigen;

// namespace FEA = M2DO_FEA ;
namespace LSM = M2DO_LSM ;

int main(int argc, char *argv[]){
	
    Matrix3d eye3 = Matrix3d::Identity(3,3);
    
	// 0.1. Mesh generation 
    const int exy[2] = {80, 40};
    const double Lxy[2] = {160.,80};
    const double h = 0.5;  
    int npe = 3, dpn = 6, dpe = 18;
    class FEAMesh feaMesh(Lxy, exy, true);

    // load file
    const char* f1 = argv[1];
    const char* f2 = argv[2];
    myreadTXT read(f1, f2);
    // myreadTXT read("./fea_80_40_0.1/FEAConfiguration.txt", "./fea_80_40_0.1/FEA_solution.txt");
    read.read_config(); 
    // MatrixXd NODE = read.feaConfig.NODE;
    // MatrixXi ELEM = read.feaConfig.ELEM;
    feaMesh.BCid = read.feaConfig.BCid;
    Force force; 
    force.NM = read.feaConfig.Force_NM;
    force.fix = read.feaConfig.Force_fix;
    	
    const uint nELEM = feaMesh.ELEM.rows();
    const uint nNODE = feaMesh.NODE.rows();
    const uint nDOF = nNODE*feaMesh.dpn;
	
    
    // Materials =================================
    double E = 1.0e6, v = 0.3, rho = 1.0;
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

    // OPTION =========================================
    struct OPTION option;
    option.initper = 1;

    // Read output fields =================================
    read.read_results();
    MatrixXd GU_u = read.feaField.GU_u;
    VectorXd GU_Rv = read.feaField.GU_Rv;
    VectorXd p_Adjoint = read.feaField.p_Adjoint;
    
    LSM::InputOutput io;
    double radius = 5;
    
    // setup levelset layer ================================
    LSM::Mesh lsmMesh(exy[0], exy[1], false);
    vector<LSM::Hole> holes;
    LSM::LevelSet levelSet(lsmMesh, holes, 0.5, 6, false);
    levelSet.reinitialise();

    // boundary discretization & bptss
    LSM::Boundary boundary(levelSet); 
    boundary.discretise();
    boundary.computeAreaFractions();
 
    // purturb areafraction =================================
    std::vector<double> areafraction0(nELEM);
    double h_eps = 0.01; // amount of perturbation
    for (uint ee = 0; ee < nELEM; ee++){
        if (lsmMesh.elements[ee].area < 1e-3) feaMesh.areafraction[ee] = 1e-3 ;
        else feaMesh.areafraction[ee] = lsmMesh.elements[ee].area ;
    }
    areafraction0 = feaMesh.areafraction;

    // find elements close to the boundary
    uint nbpts = boundary.points.size();
    vector<int> id_bndELEM;
    double p = radius*2; 
    double x_, y_, xc_, yc_;
    for (uint ee = 0; ee < nELEM; ++ee){
        xc_ = feaMesh.Centeroids(ee,0); 
        yc_ = feaMesh.Centeroids(ee,1); 
        for (uint bb = 0; bb < nbpts; ++bb){
            x_ = boundary.points[bb].coord.x;
            y_ = boundary.points[bb].coord.y;
            if (std::pow(xc_-x_,2) + std::pow(yc_-y_,2) < p*p){
                id_bndELEM.push_back(ee);
                continue;
            }
        }
    }

    // let's say in this case, gausspoints lies at the centeroid
    // and used for least square interpolation
    double compliance0, compliance;
    compliance0 = GU_u.cwiseProduct(force.fix.leftCols(3)).norm();
    MatrixXd GU_u_tmp; VectorXd GU_Rv_tmp;

    MatrixXd Sensitivity_center(id_bndELEM.size(),3);
    for (uint pp = 0; pp < id_bndELEM.size(); ++pp){
        int id_ = id_bndELEM[pp];
        // MatrixXd GU_u_tmp = MatrixXd::Zero(nNODE,3);
        // MatrixXd GU_R_tmp = eye3.replicate(nNODE,1);
        // VectorXd GU_Rv_tmp = EICR_SHELL::fmat2store(GU_R_tmp);
        GU_u_tmp = GU_u; GU_Rv_tmp = GU_Rv;
        feaMesh.areafraction = areafraction0;
        feaMesh.areafraction[id_] -= h_eps;
        x_ = feaMesh.Centeroids(id_,0);
        y_ = feaMesh.Centeroids(id_,1);
        // start from previous result
        f_nlgeom(material,force, feaMesh, option, GU_u_tmp, GU_Rv_tmp);
        compliance = GU_u_tmp.cwiseProduct(force.fix.leftCols(3)).norm();
        Sensitivity_center(pp,0) = x_;
        Sensitivity_center(pp,1) = y_;
        Sensitivity_center(pp,2) = -(compliance - compliance0)/h_eps;
    }

    // print sensitivity .. #TOFIX
    std::ofstream outfile;
	outfile.open("Num_sens.txt");
	
	outfile << "Num_sens: \n";
	outfile << Sensitivity_center;
	
	outfile.close();	

}