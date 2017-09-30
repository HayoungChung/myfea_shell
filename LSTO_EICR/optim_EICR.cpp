// this is a replica of prev. cantilever example,  
// infor 2D, triangular shell element 
// 2004, Allaire, JCP benchmark
 
// #include "../../M2DO_FEA/include/M2DO_FEA.h" 
#include <cassert> 
#include <iostream> 
#include <cmath> 
#include <algorithm> 
#include <fstream> 
 
// #include "./../M2DO_LSM/include/M2DO_LSM.h" 
#include "M2DO_LSM.h" 
#include "EICR.h" 
#include "./m2do_lin_sensitivity.h"
// #define PI 3.14159265359 
 
using namespace std ; 
using namespace Eigen; 
 
// namespace FEA = M2DO_FEA ; 
namespace LSM = M2DO_LSM ; 
// using namespace EICR_FEA;
 
int main(int argc, char *argv[]){  
	// int ex = std::stoi(argv[1]), ey = std::stoi(argv[2]); 
	std::ofstream log;
	log.open("ComplianceLog.txt");
	std::ofstream file; 
	char fname[32]; 
 
	// const int exy[2] = {std::stoi(argv[1]), std::stoi(argv[2])};// ex,ey};//{80, 40}; 
	const int exy[2] = {80, 40}; 
	const double Lxy[2] = {80.0, 40.0};//{std::stod(argv[1]), std::stod(argv[2])}; // NOTE: 만약 lxy exy 다르면 현재는 island 발생...
	const double h = 40;   
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
	} 
 
	// 1. material 
	double E = 1.0/1.6, v = 0.3, rho = 1.0; 
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
	//tipnode1 = feaMesh.get_nodeID(std::stod(argv[1]),std::stod(argv[2])/2, 1e-3, 1e-3);//(160.0, 40.0, 1e-3, 1e-3); 
	tipnode1 = feaMesh.get_nodeID(80.0, 20.0, 1e-3, 1e-3); 
 
	Force force; 
	force.NM  = MatrixXd::Zero(nELEM,feaMesh.dpn); 
	force.fix = MatrixXd::Zero(nNODE,feaMesh.dpn); 
	 
	feaMesh.set_Force(1,tipnode1,-std::stod(argv[1]), force.fix); 
	// feaMesh.set_Force(1,tipnode1,-1, force.fix); 
 
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
	option.saveflag = true;
 
	// ===================================== // 
	// Level set 
	double moveLimit = 0.5; 
	int maxIter = 200; 
	double sampleInterval = 50; 
	double maxArea = 0.5; 
 
	LSM::Mesh lsmMesh(Lxy[0], Lxy[1], false); 
	double meshArea = lsmMesh.width * lsmMesh.height; 
 
	vector<LSM::Hole> holes; 
 
	/*
	holes.push_back(LSM::Hole(10,  10, 3)) ;
	holes.push_back(LSM::Hole(20,  10, 3)) ;
	holes.push_back(LSM::Hole(30,  10, 3)) ;
	holes.push_back(LSM::Hole(40,  10, 3)) ;
	holes.push_back(LSM::Hole(50,  10, 3)) ;
	holes.push_back(LSM::Hole(60,  10, 3)) ;
	holes.push_back(LSM::Hole(70,  10, 3)) ;

	holes.push_back(LSM::Hole(10,  30, 3)) ;
	holes.push_back(LSM::Hole(20,  30, 3)) ;
	holes.push_back(LSM::Hole(30,  30, 3)) ;
	holes.push_back(LSM::Hole(40,  30, 3)) ;
	holes.push_back(LSM::Hole(50,  30, 3)) ;
	holes.push_back(LSM::Hole(60,  30, 3)) ;
	holes.push_back(LSM::Hole(70,  30, 3)) ;

	holes.push_back(LSM::Hole(15,  20, 3)) ;
	holes.push_back(LSM::Hole(25,  20, 3)) ;
	holes.push_back(LSM::Hole(35,  20, 3)) ;
	holes.push_back(LSM::Hole(45,  20, 3)) ;
	holes.push_back(LSM::Hole(55,  20, 3)) ;
	holes.push_back(LSM::Hole(65,  20, 3)) ;
	*/

	    holes.push_back(LSM::Hole(8, 7, 2.5)) ; 
	    holes.push_back(LSM::Hole(16, 14, 2.5)) ; 
	    holes.push_back(LSM::Hole(24, 7, 2.5)) ; 
	    holes.push_back(LSM::Hole(32, 14, 2.5)) ; 
	    holes.push_back(LSM::Hole(40, 7, 2.5)) ; 
	    holes.push_back(LSM::Hole(48, 14, 2.5)) ; 
	    holes.push_back(LSM::Hole(56, 7, 2.5)) ; 
	    holes.push_back(LSM::Hole(64, 14, 2.5)) ; 
	    holes.push_back(LSM::Hole(72, 7, 2.5)) ; 
	    holes.push_back(LSM::Hole(8, 20, 2.5)) ; 
	    holes.push_back(LSM::Hole(16, 27, 2.5)) ; 
	    holes.push_back(LSM::Hole(24, 20, 2.5)) ; 
	    holes.push_back(LSM::Hole(32, 27, 2.5)) ; 
	    holes.push_back(LSM::Hole(40, 20, 2.5)) ; 
	    holes.push_back(LSM::Hole(48, 27, 2.5)) ; 
	    holes.push_back(LSM::Hole(56, 20, 2.5)) ; 
	    holes.push_back(LSM::Hole(64, 27, 2.5)) ; 
	    holes.push_back(LSM::Hole(72, 20, 2.5)) ; 
	    holes.push_back(LSM::Hole(8, 33, 2.5)) ; 
	    holes.push_back(LSM::Hole(24, 33, 2.5)) ; 
	    holes.push_back(LSM::Hole(40, 33, 2.5)) ; 
	    holes.push_back(LSM::Hole(56, 33, 2.5)) ; 
	    holes.push_back(LSM::Hole(72, 33, 2.5)) ; 

	LSM::LevelSet levelSet(lsmMesh, holes, moveLimit, 6, false);  
	LSM::InputOutput io; 
	levelSet.reinitialise(); 
 
	LSM::Boundary boundary(levelSet);  
	LSM::MersenneTwister rng ; 
 
	boundary.discretise() ; 
	boundary.computeAreaFractions() ; 
	
	#if __DEBUGFLAG__
		io.saveLevelSetVTK(9999, levelSet) ; 
		io.saveAreaFractionsVTK(9999, lsmMesh) ; 
		io.saveBoundarySegmentsTXT(9999, boundary) ; 
	#endif
 
	unsigned int nReinit = 0 ; 
	double time = 0 ; 
 
	// Loging outputs 
	vector<double> times ; 
	vector<double> compliances ;  
	vector<double> areas ; 
 
	vector<double> lambdas(2) ; 
	// const int sSensitivityens_order = 2 ; 
 
	double areamin = 0.01; // minimum element area to compute sensitivity 
	double radius = 4; // radius for least square calculation 
	double Tolerance = areamin; //1e-6;
 
	printf("--------------------------------\n"); 
	printf("%8s %12s %10s\n", "Iteration", "Compliance", "Area"); 
	printf("--------------------------------\n"); 

	int n_iterations = 0 ; 
	double compliance; 
	// string savefile0 = "FEAconfig"; 
	while (n_iterations < maxIter) { 
		++n_iterations ; 

			// Discretize and compute AreaFraction 
			boundary.discretise() ; 
			boundary.computeAreaFractions() ; 

			// Assign areafractions 
			for (unsigned int ee=0 ; ee< nELEM ; ee++) 
			{ 
					if (lsmMesh.elements[ee].area <= 1e-2) feaMesh.areafraction[ee] = 0 ; 
//					if (lsmMesh.elements[ee].area < 1e-6) feaMesh.areafraction[ee] = 0 ; // 이 경우에 주변에 areafraction=0 인 mesh로 둘러싸인 경우에만 essential BC 를 걸어줘여
					else feaMesh.areafraction[ee] = lsmMesh.elements[ee].area ; 
			} 

			// Write level set and boundary segments to file. 
			io.saveLevelSetVTK(n_iterations, levelSet) ; 
			io.saveAreaFractionsVTK(n_iterations, lsmMesh) ; 
		
			// initial condition 
			MatrixXd GU_u = MatrixXd::Zero(nNODE,3); 
			MatrixXd GU_R = eye3.replicate(nNODE,1); 
			VectorXd GU_Rv = EICR_SHELL::fmat2store(GU_R); 
			// output: GU_u, GU_Rv 
			VectorXd p_Adjoint; 
			p_Adjoint = f_nlgeom(material, force, feaMesh, option, 
				GU_u, GU_Rv); 
				
			compliance = (GU_u.cwiseProduct(force.fix.leftCols(3))).sum();

			Sensitivity sens(feaMesh, material, GU_u, GU_Rv, p_Adjoint); 
			// Compute compliance sensitivities (stress*strain) at the Gauss points. 
			sens.ComputeComplianceSensitivities(Tolerance) ; 
			
			bool isLinear; 
			if (n_iterations < 100) isLinear = true; 
			else isLinear = false; 
			
			sens.to_gptSens(true);
			
			SensitivityAnalysis m2doLeastsquare(feaMesh, sens.gptsSens);
			
			double max_sens = 0, min_sens = 0; 
			int max_i, min_i; 
			
			#if __DEBUGFLAG__
				snprintf(fname, sizeof(char)*32, "step_%i.bndsens.txt", n_iterations); 
				file.open(fname);
			#endif

			for (int i=0 ; i<boundary.points.size() ; i++) 
			{ 
					vector<double> bPoint (2, 0) ; 
					bPoint[0] = boundary.points[i].coord.x; 
					bPoint[1] = boundary.points[i].coord.y; 

					// double boundarySens = sens.ComputeBoundaryPointSensitivity(bPoint, radius, 5, isLinear, Tolerance); 
					m2doLeastsquare.ComputeBoundarySensitivities(2, bPoint, 0.01); // pushing back boundarySens
					
					// sens.ComputeBoundarySensitivities(radius, bPoint) ; 

					boundary.points[i].sensitivities[0] = m2doLeastsquare.boundarySens[i]; // -boundarySens ; 
					boundary.points[i].sensitivities[1] = -1 ; 
					
					#if __DEBUGFLAG__
						file << bPoint[0] << "\t" << bPoint[1] << "\t" << -boundarySens << "\n";
					#endif
			} 
			
			#if __DEBUGFLAG__
				file << "\n==================\n"; 
				file.close(); // VERIFIED  
			#endif

			double timeStep ; 

			vector<double> constraintDistances ; 
			constraintDistances.push_back(meshArea*maxArea - boundary.area) ; 

			/* Initialise the optimisation object. 

					The Optimise class is a lightweight object so there is no cost for 
					reinitialising at every iteration. A smart compiler will optimise 
					this anyway, i.e. the same memory space will be reused. It is better 
					to place objects in the correct scope in order to aid readability 
					and to avoid unintended name clashes, etc. 
				*/ 
			LSM::Optimise optimise(boundary.points, constraintDistances, 
					lambdas, timeStep, levelSet.moveLimit, false) ; 

			// Perform the optimisation. 
			optimise.solve() ; 

			// Extend boundary point velocities to all narrow band nodes. 
			levelSet.computeVelocities(boundary.points, timeStep, 0, rng) ; 

			// Compute gradient of the signed distance function within the narrow band. 
			levelSet.computeGradients() ; 

			// Update the level set function. 
			bool isReinitialised = levelSet.update(timeStep) ; 

			// Reinitialise the signed distance function, if necessary. 
				// if (!isReinitialised) 
				// { 
				// 		// Reinitialise at least every 20 iterations. 
				// 		if (nReinit == 20) 
				// 		{ 
							levelSet.reinitialise() ; 
				// 			nReinit = 0 ; 
				// 		} 
				// } 
				// else nReinit = 0 ; 

			// Increment the number of steps since reinitialisation. 
			nReinit++ ; 

			// Increment the time. 
			time += timeStep; 

			// Calculate current area fraction. 
			double area = boundary.area / meshArea ; 

			// Record the time, compliance, and area. 
			times.push_back(time) ; 
			areas.push_back(area) ; 

			// Print statistics. 
			printf("%8.1f %12.4f %10.4f\n", double (n_iterations), compliance, area) ; 
			log << n_iterations << "\t" << compliance << "\t" << area <<"\n"; 


			snprintf(fname, sizeof(char)*32, "step_%i.vtk", n_iterations); 
			feaMesh.to_vtk(GU_u, fname); //verified
			// print ------------------------------------- 
			// #if __DEBUGFLAG__	
				snprintf(fname, sizeof(char)*32, "step_%i.txt", n_iterations); 
				file.open(fname); 
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

				file << "GU_u: \n"; 
				file << GU_u; 
				file << "\n==================\n"; 

				file << "AreaFraction: \n"; 
				for (int ee = 0; ee < nELEM; ++ee){ 
						file << feaMesh.areafraction[ee] << "\n"; 
				} 
				file << "\n==================\n"; 

				file.close(); // VERIFIED 
			// #endif
		} 
	log.close();
	return 0; 
} 
