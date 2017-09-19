// #include <FEA_hy.h>
#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "../../eigen3/Eigen/Dense"
#include "../../eigen3/Eigen/Sparse"
#include "linshell.h"
// #include "FEA_hy.h"
// #include "f_nlgeom.h"

#define EIGEN_USE_MKL_ALL
#define PI 3.14159265359

using namespace Eigen;

int main(){
    setNbThreads(4);
    initParallel();
    // 
    Matrix3d eye3 = Matrix3d::Identity();
        
    int npe = 3, dpn = 6, dpe = 18;
    bool isOverlaid = true;
    // Mesh generation
    const double Lxy[2] = {12, 2};
    const int exy[2] = {12,2}; //{80,10}; //{2, 2};
    const double h = 0.5; 
    
    class FEAMesh feaMesh(Lxy, exy, isOverlaid); 
    int nELEM = feaMesh.ELEM.rows();
    int nNODE = feaMesh.NODE.rows();
    int nDOF = nNODE*dpn;

    // std::cout << "NODE: \n" << feaMesh.NODE << "\n============" << std::endl;
    // std::cout << "ELEM: \n" << feaMesh.ELEM << "\n============" << std::endl;
    // std::cout << "BCid: \n" << feaMesh.BCid << "\n============" << std::endl;

    struct OPTION option;
    option.initper = 10;

    // Material property
    double E = 1.2e6, v = 0.3; 
    MatrixXd Cijkl(3,3), Amat(3,3), Dmat(3,3), Bmat(3,3);
    Cijkl << 1, v, 0, v, 1, 0, 0, 0, (1-v)*0.5;
    Cijkl *= E/(1-std::pow(v,2));
    Amat = Cijkl * h;
    Dmat = Cijkl * std::pow(h,3)/12.0;
    Bmat.setZero(3,3);

    Material_ABD material0;
    material0.Amat = Amat;
    material0.Dmat = Dmat;
    material0.Bmat = Bmat;

    std::vector<Material_ABD> material;
    
    for (int ee = 0; ee < nELEM; ++ee){
        material.push_back(material0);        
    }
    // BC
    double posX = 0., posY = 0., Xtol = 1e-3, Ytol = 1e3;

	std::vector<int> Xlo = feaMesh.get_nodeID(posX, posY, Xtol, Ytol);
	std::vector<int> BCtmp;
	for (unsigned int dd = 0; dd < feaMesh.dpn; ++dd){ // clamp
		std::vector<int> tmp = feaMesh.get_dof(dd,Xlo);
		// std::cout << Map<VectorXi>(tmp.data(), tmp.size()).transpose() << " \n";
		for (unsigned int mm = 0; mm < tmp.size(); ++mm){
			BCtmp.push_back(tmp[mm]);
		}
	}
	for (unsigned int dd = 0; dd < nNODE; ++ dd){ // to 2D
		// BCtmp.push_back(dd*6+2); 
		// BCtmp.push_back(dd*6+3);
		// BCtmp.push_back(dd*6+4);
		BCtmp.push_back(dd*6+5); // no rotation
	}
	feaMesh.get_BCid(BCtmp);
    
    std::cout << "NODE: \n" << feaMesh.NODE << "\n============" << std::endl;
    std::cout << "ELEM: \n" << feaMesh.ELEM << "\n============" << std::endl;
    std::cout << "BCid: \n" << feaMesh.BCid << "\n============" << std::endl;

    // Resultant-driven rotation (2pi)
    double curv_ = -2*PI / Lxy[0]/4;
    Vector3d eps0, kappa0;
    eps0 << 0, 0, 0;
    kappa0 << curv_, 0, 0;

    MatrixXd NM = MatrixXd::Zero(feaMesh.dpn,1); 
    NM << Amat*eps0, Dmat*kappa0;

    MatrixXd Force_NM(nELEM, feaMesh.dpn); 
    MatrixXd Force_Fix = MatrixXd::Zero(nNODE, feaMesh.dpn);
    Force_NM = NM.transpose().replicate(nELEM,1);

    std::cout << "Force_NM[1]: \n" << Force_NM.row(1) << "\n============\n" << std::endl;
    std::cout << "Amat[2]: \n" << material[2].Amat << "\n============\n" << std::endl;
    std::cout << "Dmat[10]: \n" << material[10].Dmat << "\n============\n" << std::endl;
    
    Force force; 
    force.NM = Force_NM;
    force.fix = Force_Fix;

    // Solve! 
    MatrixXd GU_u = MatrixXd::Zero(nNODE,3);
    MatrixXd GU_R = eye3.replicate(nNODE,1);
    
    // VectorXd GU_Rv = Map<VectorXd>(GU_R.transpose().data(),nNODE*9,1);
    VectorXd GU_Rv = fmat2store(GU_R);
    // std::cout << "GU_Rv: \n" << GU_Rv << "\n============\n" << std::endl;
    
    f_nlgeom(material, force, feaMesh, option, GU_u, GU_Rv);
    // Post-process
    double sol_anal = 1.0/curv_*(1-std::cos(curv_*Lxy[0]));

    std::cout << GU_u.col(2).maxCoeff() << "\t" << sol_anal << std::endl;
    return 0;

}

// void get_Mesh
