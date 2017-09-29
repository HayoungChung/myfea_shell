#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>

#include "EICR.h"

using namespace Eigen;

int main(){ // this function test restart capability
    // mesh generation
    const int exy[2] = {8, 4};
    const double Lxy[2] = {160., 80.};
    const double h = 0.5;
    int npe = 3, dpn = 6, dpe = 18;
    FEAMesh feaMesh(Lxy, exy, true);

    // read previous output and force
    const char* f1 = "./FEAConfiguration.txt";
    const char* f2 = "./FEA_solution.txt";

    myreadTXT read(f1, f2);
    read.read_config();
    feaMesh.BCid = read.feaConfig.BCid;
    Force force;

    force.NM = read.feaConfig.Force_NM;
    force.fix = read.feaConfig.Force_fix;

    const uint nELEM = feaMesh.ELEM.rows();
    const uint nNODE = feaMesh.NODE.rows();
    const uint nDOF = nNODE*feaMesh.dpn;
    
    read.read_results();
    MatrixXd GU_u = read.feaField.GU_u;
    VectorXd GU_Rv = read.feaField.GU_Rv;
    VectorXd p_Adjoint = read.feaField.p_Adjoint;

    // material
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

    for (uint ee = 0; ee < nELEM; ee++){
        feaMesh.areafraction[ee] = 1;
    }

    // OPTION =========================================
    struct OPTION option;
    option.initper = 1;

    f_nlgeom(material,force, feaMesh, option, GU_u, GU_Rv);


}