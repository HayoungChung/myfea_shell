#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "./lin_shell.h"
// #include "./LinSensitivity.h"
#include "./m2do_lin_sensitivity.h"
#include "./../../01.M2DO_code/M2DO_LSM/include/M2DO_LSM.h"

using namespace Eigen;
namespace LSM = M2DO_LSM;

int main()
{
    Matrix3d eye3 = Matrix3d::Identity();
    bool isOverlaid = true;

    int npe = 3, dpn = 6, dpe = 18;

    // Mesh generation
    const double Lxy[2] = {80., 40.};
    const int exy[2] = {80, 40};
    const double h = 2; // the h/L < 0.1

    LSM::Mesh lsmMesh(exy[0], exy[1], false);
    double meshArea = lsmMesh.width * lsmMesh.height;

    std::vector<LSM::Hole> holes;
    // holes.push_back(LSM::Hole(10, 10, 3));
    // holes.push_back(LSM::Hole(20, 10, 3));
    // holes.push_back(LSM::Hole(30, 10, 3));
    // holes.push_back(LSM::Hole(15, 20, 3));
    // holes.push_back(LSM::Hole(25, 20, 3));
    // holes.push_back(LSM::Hole(15, 0, 3));
    // holes.push_back(LSM::Hole(25, 0, 3));
    holes.push_back(LSM::Hole(8, 7, 2.5));
    holes.push_back(LSM::Hole(16, 13.5, 2.5));
    holes.push_back(LSM::Hole(24, 7, 2.5));
    holes.push_back(LSM::Hole(32, 13.5, 2.5));
    holes.push_back(LSM::Hole(40, 7, 2.5));
    holes.push_back(LSM::Hole(48, 13.5, 2.5));
    holes.push_back(LSM::Hole(56, 7, 2.5));
    holes.push_back(LSM::Hole(64, 13.5, 2.5));
    holes.push_back(LSM::Hole(72, 7, 2.5));
    holes.push_back(LSM::Hole(8, 20, 2.5));
    holes.push_back(LSM::Hole(16, 26.5, 2.5));
    holes.push_back(LSM::Hole(24, 20, 2.5));
    holes.push_back(LSM::Hole(32, 26.5, 2.5));
    holes.push_back(LSM::Hole(40, 20, 2.5));
    holes.push_back(LSM::Hole(48, 26.5, 2.5));
    holes.push_back(LSM::Hole(56, 20, 2.5));
    holes.push_back(LSM::Hole(64, 26.5, 2.5));
    holes.push_back(LSM::Hole(72, 20, 2.5));
    holes.push_back(LSM::Hole(8, 33, 2.5));
    holes.push_back(LSM::Hole(24, 33, 2.5));
    holes.push_back(LSM::Hole(40, 33, 2.5));
    holes.push_back(LSM::Hole(56, 33, 2.5));
    holes.push_back(LSM::Hole(72, 33, 2.5));

    double maxArea = 0.5;
    double temperature = 0;
    LSM::MersenneTwister rng;
    std::vector<double> lambdas(2);
    std::vector<double> times;
    std::vector<double> compliances;
    std::vector<double> areas;

    double time = 0;

    LSM::LevelSet levelSet(lsmMesh, holes);
    levelSet.reinitialise();
    LSM::InputOutput io;
    LSM::Boundary boundary(levelSet);

    std::cout << "meshing fea \n";

    FEAMesh feaMesh(Lxy, exy, isOverlaid);
    const unsigned int nELEM = feaMesh.ELEM.rows();
    const unsigned int nNODE = feaMesh.NODE.rows();
    const unsigned nDOF = nNODE * dpn;

    std::ofstream eFile("elem.txt");
    eFile << feaMesh.ELEM << std::endl;
    eFile.close();

    std::ofstream nFile("node.txt");
    nFile << feaMesh.NODE << std::endl;
    nFile.close();

    // Material

    const double E = 1.2e6, v = 0.3;
    MatrixXd Cijkl(3, 3), Amat(3, 3), Dmat(3, 3), Bmat(3, 3);
    Cijkl << 1., v, 0., v, 1, 0., 0., 0., 0.5 * (1 - v);
    Cijkl *= E / (1 - v * v);
    Amat = Cijkl * h;
    Dmat = Cijkl * std::pow(h, 3) / 12.0;
    Bmat.setZero(3, 3);

    Material_ABD material0;
    material0.Amat = Amat;
    material0.Dmat = Dmat;
    material0.Bmat = Bmat;

    std::vector<Material_ABD> material;

    for (int ee = 0; ee < nELEM; ++ee)
    {
        material.push_back(material0);
    }

    // BC
    std::cout << "BCs in fea \n";
    double posX = 0, posY = 0, Xtol = 1e-3, Ytol = 1e3; // Dirichlet.
    std::vector<int> Xlo = feaMesh.get_nodeID(posX, posY, Xtol, Ytol);
    std::vector<int> BCtmp;
    for (unsigned int dd = 0; dd < feaMesh.dpn; dd++)
    {
        std::vector<int> tmp = feaMesh.get_dof(dd, Xlo);
        for (unsigned int mm = 0; mm < tmp.size(); mm++)
        {
            BCtmp.push_back(tmp[mm]);
        }
    }

    feaMesh.get_BCid(BCtmp);
    /*
    std::cout << "NODE: \n"
              << feaMesh.NODE << "\n============" << std::endl;
    std::cout << "ELEM: \n"
              << feaMesh.ELEM << "\n============" << std::endl;
    std::cout << "BCid: \n"
              << feaMesh.BCid << "\n============" << std::endl;
    */
    MatrixXd Force_NM(nELEM, feaMesh.dpn);
    MatrixXd Force_Fix = MatrixXd::Zero(nNODE, feaMesh.dpn);
    Force_NM.fill(0.0);

    std::vector<int> Xtip = feaMesh.get_nodeID(Lxy[0], Lxy[1]/2, 1e-3, 1e-3);

    for (int tt = 0; tt < Xtip.size(); tt++)
    {
        Force_Fix(Xtip[tt], 1) = -1;
        Force_Fix(Xtip[tt], 4) = -1;
    }

    // std::vector<int> Xtip2 = feaMesh.get_nodeID(Lxy[0], 0, 1e-3, 1e-3);

    // for (int tt = 0; tt < Xtip2.size(); tt++)
    // {
    //     Force_Fix(Xtip2[tt], 2) = -1;
    // }

    Force force;
    force.NM = Force_NM;
    force.fix = Force_Fix;
    io.saveLevelSetVTK(0, levelSet);
    unsigned int MAXITER = 200;
    for (unsigned int n_iterations = 0; n_iterations < MAXITER; n_iterations++)
    {
        // =========== SOLVE ======================== //
        std::cout << "discretizing in LS \n";
        boundary.discretise();
        boundary.computeAreaFractions();

        for (unsigned int ee = 0; ee < nELEM; ee++)
        {
            if (lsmMesh.elements[ee].area < 1e-3)
                feaMesh.areafraction[ee] = 1e-3;
            else
                feaMesh.areafraction[ee] = lsmMesh.elements[ee].area;
            // feaMesh.areafraction[ee] = 1.0; // if triangulated
        }

        std::cout << "fea assemblying in LS \n";
        // f_lin_shell(feaMesh, material, force, GKT, Res);
        LinShell lin_shell(feaMesh, material, force);
        lin_shell.compute();

        std::cout << "making compressed (sparse) \n";
        //SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> spsolver;
        ConjugateGradient<SparseMatrix<double>> spsolver;
        lin_shell.sGKT.makeCompressed(); // this is time-consuming!!!! (DONNO WHY)
        spsolver.compute(lin_shell.sGKT);

        MatrixXd u;
        MatrixXd u6;

        std::cout << "starting lin solve\n";
        u = spsolver.solve(lin_shell.Res);
        std::cout << "end of lin solve\n";

        u6 = Map<MatrixXd>(u.data(), 6, nNODE).transpose();
        std::vector<GptsCompl> gptsCompl = lin_shell.get_GaussCompl(u6);

        std::ofstream dispFile("displacement.txt");
        dispFile << u6 << std::endl;
        dispFile.close();
        feaMesh.to_vtk(u6);

        /*
        // WIP: checking if the locations of the Gpts are the source of the unsymm.
        // this leads to rank dificiency... 
        std::vector<GptsCompl> elemCompl(feaMesh.ELEM.rows()); // elem-wise
        double x_temp, y_temp, s_temp;
        for (unsigned int ee = 0; ee < feaMesh.ELEM.rows(); ee++)
        {
            x_temp = 0; y_temp = 0; s_temp = 0;
            
            x_temp = gptsCompl[ee*4].x + gptsCompl[ee*4+1].x + gptsCompl[ee*4+2].x + gptsCompl[ee*4+3].x;
            x_temp /= 4;
            y_temp = gptsCompl[ee*4].y + gptsCompl[ee*4+1].y + gptsCompl[ee*4+2].y + gptsCompl[ee*4+3].y;
            y_temp /= 4;
            s_temp = gptsCompl[ee*4].sens + gptsCompl[ee*4+1].sens + gptsCompl[ee*4+2].sens + gptsCompl[ee*4+3].sens;
            s_temp /= 4;
            
            elemCompl[ee].x = x_temp;
            elemCompl[ee].y = y_temp;
            elemCompl[ee].sens = s_temp;
        }
        feaMesh.isOverlaid = false;
        LinSensitivity linSens(feaMesh, elemCompl); 
        feaMesh.isOverlaid = true;
        // end of WIP
        */
        // LinSensitivity linSens(feaMesh, gptsCompl);
        SensitivityAnalysis m2doSens(feaMesh, gptsCompl);

        std::ofstream gsfile("gpts_sens_test.txt");
        for (unsigned int ii = 0; ii < gptsCompl.size(); ii++)
        {
            gsfile << gptsCompl[ii].x << " " << gptsCompl[ii].y << " " << gptsCompl[ii].sens << std::endl;
        }
        gsfile.close();

        std::ofstream bsfile("bpts_sens_test.txt");

        std::cout << "starting least square interpolation \n";
        for (int i = 0; i < boundary.points.size(); i++)
        {
            std::vector<double> bPoint(2, 0);
            bPoint[0] = boundary.points[i].coord.x;
            bPoint[1] = boundary.points[i].coord.y;

            // Interpolate Guass point sensitivities by least squares.
            // TOFIX: currently 2D
            // boundary.points[i].sensitivities[0] = -linSens.ComputeBoundaryPointSensitivity2D(bPoint, 4, 5, 0.01);
            m2doSens.ComputeBoundarySensitivities(2, bPoint, 0.01);
        }
        for (int i = 0; i < boundary.points.size(); i++)
        {
            std::vector<double> bPoint(2, 0);
            bPoint[0] = boundary.points[i].coord.x;
            bPoint[1] = boundary.points[i].coord.y;

            boundary.points[i].sensitivities[0] = m2doSens.boundarySens[i];
            if (bsfile.is_open())
            {
                bsfile << bPoint[0] << " " << bPoint[1] << " " << boundary.points[i].sensitivities[0] << std::endl;
            }

            // Assign sensitivities.
            boundary.points[i].sensitivities[1] = -1;

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
        bsfile.close();
        std::cout << "end of the least square interpolation \n";

        double timeStep;

        // Constraint distance vector.
        std::vector<double> constraintDistances;

        // Push current distance from constraint violation into vector.
        constraintDistances.push_back(meshArea * maxArea - boundary.area);

        /* Initialise the optimisation object.
    
               The Optimise class is a lightweight object so there is no cost for
               reinitialising at every iteration. A smart compiler will optimise
               this anyway, i.e. the same memory space will be reused. It is better
               to place objects in the correct scope in order to aid readability
               and to avoid unintended name clashes, etc.
             */
        LSM::Optimise optimise(boundary.points, constraintDistances,
                               lambdas, timeStep, levelSet.moveLimit, false);

        // Perform the optimisation.
        optimise.solve();

        // Extend boundary point velocities to all narrow band nodes.
        levelSet.computeVelocities(boundary.points, timeStep, temperature, rng);

        // Compute gradient of the signed distance function within the narrow band.
        levelSet.computeGradients();

        // Update the level set function.
        bool isReinitialised = levelSet.update(timeStep);

        // // Reinitialise the signed distance function, if necessary.
        // if (!isReinitialised)
        // {
        //     // Reinitialise at least every 20 iterations.
        //     if (nReinit == 20)
        //     {
        levelSet.reinitialise();
        //         nReinit = 0;
        //     }
        // }
        // else
        //     nReinit = 0;

        // // Increment the number of steps since reinitialisation.
        // nReinit++;

        // Increment the time.
        time += timeStep;

        double area = boundary.area / meshArea;

        times.push_back(time);
        areas.push_back(area);

        printf("%8.1f %12.4f %10.4f\n", double(n_iterations), lin_shell.compliance, area);

        lin_shell.compliance = 0;

        io.saveLevelSetVTK(n_iterations, levelSet);
        io.saveAreaFractionsVTK(n_iterations, lsmMesh);
        io.saveBoundarySegmentsTXT(n_iterations, boundary);
    }
    return 0;
}
