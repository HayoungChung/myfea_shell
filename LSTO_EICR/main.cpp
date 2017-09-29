#include "M2DO_FEA.h"
#include "M2DO_LSM.h"

using namespace std ;
using namespace Eigen ;

namespace FEA = M2DO_FEA ;
namespace LSM = M2DO_LSM ;

int main () {

	/*
		Dimensionality of problem:
	*/

	const int spacedim = 2, dim = spacedim, order = 2 ;

	/*
		Level-set mesh
	*/

    const unsigned int nelx = 160, nely = 80;

	/*
		Typedefs. These clean up the code a lot.
	*/

	typedef FEA::LinearElasticity<spacedim> linear_elasticity_physics ;
	typedef FEA::LinearElasticMaterial linear_elastic_material ;
	typedef FEA::Mesh<spacedim, order, dim, linear_elasticity_physics > linear_elasticity_mesh ;
	typedef FEA::StationaryStudy<linear_elasticity_mesh, linear_elasticity_physics> stationary_study ;

	/*
		Set-up the physics:
	*/

	linear_elasticity_physics physics ;

	/*
		Add material properties:
	*/

	double E = 1.0, nu = 0.3, rho = 1.0 ;
	physics.material.push_back(linear_elastic_material (E, nu, rho)) ;

	/*
		Create a mesh object.
	*/

	linear_elasticity_mesh mesh ;

	/*
		Mesh a hyper rectangle.
	*/

	MatrixXd x_rectangle(int(pow(2, spacedim)), spacedim) ;

	x_rectangle <<    0,  0,
			  		 160,  0,
			  		 160, 80,
			  		  0, 80 ;


	vector<int> nel = {nelx, nely} ;

	double duration ;
	clock_t start ;

	mesh.hyper_rectangle(nel, x_rectangle, false) ;
    mesh.is_structured = true ;

	/*
		Add a homogeneous Dirichlet boundary condition (fix some nodes).
	*/

	vector<double> coord, tol ;

	coord = {0.0, 0.0} ;
	tol   = {1e-12, 1e9} ;

	vector<int> selected_nodes = mesh.get_nodes_by_coordinate(coord, tol) ;
	vector<int> selected_dof   = mesh.dof(selected_nodes) ;

	FEA::HomogeneousDirichletBC homogeneous_dirichlet_bc (selected_dof, mesh.n_dof) ;

	/*
		Apply a point load.
	*/

	coord = {160.0, 40.0} ;
    tol   = {1e-12, 1e-12} ;

    selected_nodes = mesh.get_nodes_by_coordinate(coord, tol) ;

    vector<double> values (selected_nodes.size() * 2, 0) ;

    for (int i = 0 ; i < selected_nodes.size() ; ++i) {
        values[2*i]   = 0;
        values[2*i+1] = -0.5;
    }

    selected_dof   = mesh.dof(selected_nodes) ;

	FEA::PointValues point_values (selected_dof, values) ;

	/*
		Next we specify that we will undertake a stationary study, which takes the form Ku = F.
	*/

	stationary_study study (mesh, physics) ;

    study.add_boundary_condition (homogeneous_dirichlet_bc) ;

	/*
		Tread carefully; these be Renato's additions.
	*/

    // Maximum displacement per iteration, in units of the mesh spacing.
    // This is the CFL limit.
    double moveLimit = 0.5 ;

    // Set maximum running time.
    double maxTime = 6000 ;
    int maxit = 100;

    // Set sampling interval.
    double sampleInterval = 50 ;

    // Set time of the next sample.
    double nextSample = 50 ;

    // Maximum material area.
    double maxArea = 0.5 ;

	// Default temperature of the thermal bath.
    double temperature = 0 ;

    // Initialise the level set mesh (same resolution as the FE mesh).
    LSM::Mesh lsmMesh(nelx, nely, false) ;

    double meshArea = lsmMesh.width * lsmMesh.height ;

    // Create two horizontal rows with four equally space holes.
    vector<LSM::Hole> holes ;

    holes.push_back(LSM::Hole(16, 14, 5)) ;
    holes.push_back(LSM::Hole(32, 27, 5)) ;
    holes.push_back(LSM::Hole(48, 14, 5)) ;
    holes.push_back(LSM::Hole(64, 27, 5)) ;
    holes.push_back(LSM::Hole(80, 14, 5)) ;
    holes.push_back(LSM::Hole(96, 27, 5)) ;
    holes.push_back(LSM::Hole(112, 14, 5)) ;
    holes.push_back(LSM::Hole(128, 27, 5)) ;
    holes.push_back(LSM::Hole(144, 14, 5)) ;
    holes.push_back(LSM::Hole(16, 40, 5)) ;
    holes.push_back(LSM::Hole(32, 53, 5)) ;
    holes.push_back(LSM::Hole(48, 40, 5)) ;
    holes.push_back(LSM::Hole(64, 53, 5)) ;
    holes.push_back(LSM::Hole(80, 40, 5)) ;
    holes.push_back(LSM::Hole(96, 53, 5)) ;
    holes.push_back(LSM::Hole(112, 40, 5)) ;
    holes.push_back(LSM::Hole(128, 53, 5)) ;
    holes.push_back(LSM::Hole(144, 40, 5)) ;
    holes.push_back(LSM::Hole(16, 66, 5)) ;
    holes.push_back(LSM::Hole(48, 66, 5)) ;
    holes.push_back(LSM::Hole(80, 66, 5)) ;
    holes.push_back(LSM::Hole(112, 66, 5)) ;
    holes.push_back(LSM::Hole(144, 66, 5)) ;

   	// Initialise the level set object (from the hole vector).
    LSM::LevelSet levelSet(lsmMesh, holes, moveLimit, 6, false) ;

    // Initialise io object.
    LSM::InputOutput io ;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise() ;

    // Initialise the boundary object.
    LSM::Boundary boundary(levelSet) ;

    // Initialise random number generator.
    LSM::MersenneTwister rng ;

    // Number of cycles since signed distance reinitialisation.
    unsigned int nReinit = 0 ;

    // Running time.
    double time = 0 ;

    // Time measurements.
    vector<double> times ;

    // Compliance measurements.
    vector<double> compliances ;

    // Boundary curvature measurements.
    vector<double> areas ;

    /* Lambda values for the optimiser.
       These are reused, i.e. the solution from the current iteration is
       used as an estimate for the next, hence we declare the vector
       outside of the main loop.
     */
    vector<double> lambdas(2) ;

    // Order for sensitivity points calculation at Gauss points (Similar to integration order: sens_order = 2 -> 2x2 sensitivity points per element)
    const int sens_order = 2 ;

    // Create sensitivity analysis instance.
    FEA::SensitivityAnalysis<spacedim, dim, order, linear_elasticity_mesh, linear_elasticity_physics, stationary_study, LSM::LevelSet, sens_order> sens(mesh, study, levelSet) ;

    // Sensitivities
    double areamin = 0.1 ; // Minimum element area to compute sensitivity.
    double radius = 2 ; // Radius for least squares calculation.

    // Compute element centroids (needed for least squares analysis).
    mesh.ComputeCentroids() ;
    // Compute Sensitivity coordinates.
    sens.ComputeSensitivitiesCoordinates(mesh) ;

    cout << "\nStarting compliance minimisation demo...\n\n" ;

    // Print output header.
    printf("--------------------------------\n");
    printf("%8s %12s %10s\n", "Iteration", "Compliance", "Area");
    printf("--------------------------------\n");

    // Integrate until we exceed the maximum time.
    int n_iterations = 0 ;

    while (n_iterations < maxit) {

    	++n_iterations ;

        // Perform boundary discretisation.
        boundary.discretise() ;

        // Compute element area fractions.
        boundary.computeAreaFractions() ;

        // Assign area fractions.
        for (unsigned int i=0 ; i<mesh.elements.size() ; i++)
        {
            if (lsmMesh.elements[i].area < 1e-3) mesh.elements[i].areafraction = 1e-3 ;
            else mesh.elements[i].areafraction = lsmMesh.elements[i].area ;
        }

        // Assemblying K and F.
        study.assemble_K_with_area_fractions (false) ;
        study.assemble_f (point_values, false) ;

        // Solve equation
        study.solve_with_HSL_MA57(false, false, false) ;

        // Compute compliance sensitivities (stress*strain) at the Gauss points.
        sens.ComputeComplianceSensitivities(mesh, physics, areamin) ;

        double max_sens = 0, min_sens = 0;
        int max_i, min_i;
        for (int i=0 ; i<boundary.points.size() ; i++)
        {
            vector<double> bPoint (2, 0) ;
            bPoint[0] = boundary.points[i].coord.x;
            bPoint[1] = boundary.points[i].coord.y;

            // Interpolate Guass point sensitivities by least squares.
            sens.ComputeBoundarySensitivities(mesh, physics, areamin, radius, bPoint) ;

            // std::cout << "(" << bPoint[0] << "; " << bPoint[1] << ") - sens:  " << sens.boundarySens[i] << std::endl;

            // Assign sensitivities.
            boundary.points[i].sensitivities[0] = -sens.boundarySens[i] ;
            boundary.points[i].sensitivities[1] = -1 ;

            // if (sens.boundarySens[i] > max_sens)
            // {
            //     max_sens = sens.boundarySens[i];
            //     max_i = i;
            // }
            // if (sens.boundarySens[i] < min_sens)
            // {
            //     min_sens = sens.boundarySens[i];
            //     min_i = i;
            // }
        }
        // clearing sens.boundarysens vector.
        sens.boundarySens.clear() ;

        // Time step associated with the iteration.
        double timeStep ;

        // Constraint distance vector.
        vector<double> constraintDistances ;

        // Push current distance from constraint violation into vector.
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
        levelSet.computeVelocities(boundary.points, timeStep, temperature, rng) ;

        // Compute gradient of the signed distance function within the narrow band.
        levelSet.computeGradients() ;

        // Update the level set function.
        bool isReinitialised = levelSet.update(timeStep) ;

        // Reinitialise the signed distance function, if necessary.
        if (!isReinitialised)
        {
            // Reinitialise at least every 20 iterations.
            if (nReinit == 20)
            {
                levelSet.reinitialise() ;
                nReinit = 0 ;
            }
        }
        else nReinit = 0 ;

        // Increment the number of steps since reinitialisation.
        nReinit++ ;

        // Increment the time.
        time += timeStep;

        // Check if the next sample time has been reached.
        // while (time >= nextSample)
        // {
            // Calculate current area fraction.
            double area = boundary.area / meshArea ;

            // Record the time, compliance, and area.
            times.push_back(time) ;
            //compliances.push_back(study.compliance);
            areas.push_back(area) ;

            // Update the time of the next sample.
            nextSample += sampleInterval ;

            // Print statistics.
            printf("%8.1f %12.4f %10.4f\n", double (n_iterations), sens.objective, area) ;

            // Write level set and boundary segments to file.
            io.saveLevelSetVTK(n_iterations, levelSet) ;
            io.saveAreaFractionsVTK(n_iterations, lsmMesh) ;
            io.saveBoundarySegmentsTXT(n_iterations, boundary) ;
        // }
    }

	/*
		Aaaaaand that's all, folks!
	*/

	cout << "\nProgram complete.\n\n" ;

	return 0 ;
}
