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
// #define PI 3.14159265359

using namespace std ;
using namespace Eigen;

// namespace FEA = M2DO_FEA ;
namespace LSM = M2DO_LSM ;

int main(int argc, char *argv[]){ 
    int ex = std::stoi(argv[1]), ey = std::stoi(argv[2]);
    const int exy[2] = {ex,ey};//{80, 40};
    const double Lxy[2] = {160.,80.}; //{160.,80.};
    const double h = 0.5;  
	const double force_in = std::stod(argv[3]); 

	Matrix3d eye3 = Matrix3d::Identity(3,3);

	// 0.1. Mesh generation 
    int npe = 3, dpn = 6, dpe = 18;
    bool isOverlaid = true;
    class FEAMesh feaMesh(Lxy, exy, isOverlaid);
    const int nELEM = feaMesh.ELEM.rows();
    const int nNODE = feaMesh.NODE.rows();
    const int nDOF = nNODE*feaMesh.dpn;
	
	// 0.1.1. BC (Dirichlet)
	// X clamping & remove transverse directions (reduced to 2D)
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
		BCtmp.push_back(dd*6+2); 
		BCtmp.push_back(dd*6+3);
		BCtmp.push_back(dd*6+4);
		// BCtmp.push_back(dd*6+5); // no rotation
	}

	// 1. material
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

	// 2. force 
	std::vector<int> tipnode1, tipnode2, tipnode3;
	tipnode1 = feaMesh.get_nodeID(160,40, 1e-3, 1e-3);//(160.0, 40.0, 1e-3, 1e-3);
	tipnode2 = feaMesh.get_nodeID(160,39, 1e-3, 1e-3);//(160.0, 39.0, 1e-3, 1e-3);
	tipnode3 = feaMesh.get_nodeID(160,41, 1e-3, 1e-3);//(160.0, 41.0, 1e-3, 1e-3);

	Force force;
	force.NM  = MatrixXd::Zero(nELEM,feaMesh.dpn);
	force.fix = MatrixXd::Zero(nNODE,feaMesh.dpn);
	
	feaMesh.set_Force(1,tipnode1,-force_in, force.fix);
	feaMesh.set_Force(1,tipnode1,-force_in, force.fix);
	feaMesh.set_Force(1,tipnode1,-force_in, force.fix);
	//feaMesh.set_Force(1,tipnode1,-0.5, force.fix);
	//feaMesh.set_Force(1,tipnode2,-0.5, force.fix);
	//feaMesh.set_Force(1,tipnostod(argv[2])

	double curv_ = 0.0;//-2*PI / Lxy[0]/4;
    Vector3d eps0, kappa0;
    eps0 << 0, 0, 0;
    kappa0 << curv_, 0, 0;
	MatrixXd eps = eps0.transpose().replicate(nELEM,1);
	MatrixXd kappa = kappa0.transpose().replicate(nELEM,1);
	feaMesh.set_Force(material, eps, kappa, force.NM);
	feaMesh.get_BCid(BCtmp);

    // ======================================================================================================//
    // ===============================   FEA    =============================================================//
    // ======================================================================================================//

	// 3. compute
	struct OPTION option;
    option.initper = 1;

	// ===================================== //
	// Level set
	double moveLimit = 0.5;
	double maxTime = 150;
	int maxIter = 100;
	double sampleInterval = 50;
	double nextSample = 50;
	double maxArea = 0.5;

	LSM::Mesh lsmMesh(exy[0], exy[1], false);
	double meshArea = lsmMesh.width * lsmMesh.height;

	vector<LSM::Hole> holes;

	// holes.push_back(LSM::Hole(16, 14, 5)) ;
    // holes.push_back(LSM::Hole(32, 27, 5)) ;
    // holes.push_back(LSM::Hole(48, 14, 5)) ;
    // holes.push_back(LSM::Hole(64, 27, 5)) ;
    // holes.push_back(LSM::Hole(80, 14, 5)) ;
    // holes.push_back(LSM::Hole(96, 27, 5)) ;
    // holes.push_back(LSM::Hole(112, 14, 5)) ;
    // holes.push_back(LSM::Hole(128, 27, 5)) ;
    // holes.push_back(LSM::Hole(144, 14, 5)) ;
    // holes.push_back(LSM::Hole(16, 40, 5)) ;
    // holes.push_back(LSM::Hole(32, 53, 5)) ;
    // holes.push_back(LSM::Hole(48, 40, 5)) ;
    // holes.push_back(LSM::Hole(64, 53, 5)) ;
    // holes.push_back(LSM::Hole(80, 40, 5)) ;
    // holes.push_back(LSM::Hole(96, 53, 5)) ;
    // holes.push_back(LSM::Hole(112, 40, 5)) ;
    // holes.push_back(LSM::Hole(128, 53, 5)) ;
    // holes.push_back(LSM::Hole(144, 40, 5)) ;
    // holes.push_back(LSM::Hole(16, 66, 5)) ;
    // holes.push_back(LSM::Hole(48, 66, 5)) ;
    // holes.push_back(LSM::Hole(80, 66, 5)) ;
    // holes.push_back(LSM::Hole(112, 66, 5)) ;
    // holes.push_back(LSM::Hole(144, 66, 5)) ;

	LSM::LevelSet levelSet(lsmMesh, holes, moveLimit, 6, false); 

	LSM::InputOutput io;
	levelSet.reinitialise();

	LSM::Boundary boundary(levelSet); 
	// LSM::MersenneTwister rng ;
	unsigned int nReinit = 0 ;
	double time = 0 ;

    // Loging outputs
	vector<double> times ;
	vector<double> compliances ; 
	vector<double> areas ;

    vector<double> lambdas(2) ;
    // const int sSensitivityens_order = 2 ;

	double areamin = 0.5; // minimum element area to compute sensitivity
	double radius = 0.02; // radius for least square calculation

    // comptue element centeroids (for least square analysis)
    // feaMesh.ComputeCentroids();
    // Sensitivity coordinates
    // sens.ComputeSensitivitiesCoordinates(feaMesh); 


    for (unsigned int ee=0 ; ee< nELEM ; ee++)
    {
        // if (lsmMesh.elements[ee].area < 1e-3) feaMesh.areafraction[ee] = 1e-3 ;
        // else feaMesh.areafraction[ee] = lsmMesh.elements[ee].area ;
		 feaMesh.areafraction[ee] = 1.0; // if triangulated
    }


	// print -------------------------------------
	std::ofstream file;
	file.open("FEAConfiguration.txt");
	file << "NODE: \n";
	file << feaMesh.NODE;//Map<VectorXi>(BCtmp.data(), BCtmp.size()).transpose();
	file << "\n==================\n";

	file << "ELEM: \n";
	file << feaMesh.ELEM;
	file << "\n==================\n";

	file << "BCid: \n";
	file << feaMesh.BCid;
	file << "\n==================\n";

	file << "Force.NM: \n";
	file << force.NM;
	file << "\n==================\n";

	file << "Force.fix: \n";
	file << force.fix;
	file << "\n==================\n";

    file << "AreaFraction: \n";
	for (int ee = 0; ee < nELEM; ++ee){
		// std::cout << feaMesh.areafraction[ee] << "\n";
		file << feaMesh.areafraction[ee] << "\n";
	}
    file << "\n==================\n";

	file.close(); // VERIFIED

    feaMesh.to_vtk();

	// initial condition
    MatrixXd GU_u = MatrixXd::Zero(nNODE,3);
    MatrixXd GU_R = eye3.replicate(nNODE,1);
    VectorXd GU_Rv = EICR_SHELL::fmat2store(GU_R);
    // output: GU_u, GU_Rv
    VectorXd p_Adjoint;
    p_Adjoint = f_nlgeom(material, force, feaMesh, option, GU_u, GU_Rv);
    feaMesh.to_vtk(GU_u);

	std::ofstream outfile;
	outfile.open("FEA_solution.txt");
	outfile << "GU_u: \n";
	outfile << GU_u;
	outfile << "\n==================\n";
	
	outfile << "GU_Rv: \n";
	outfile << GU_Rv;
	outfile << "\n==================\n";

	outfile << "p_Adjoint: \n";
	outfile << p_Adjoint;
	outfile << "\n==================\n";

	outfile.close();	

    double compliance = (GU_u.cwiseProduct(force.fix.leftCols(3))).sum();

	return 0;
}
