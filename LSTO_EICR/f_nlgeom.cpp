#include "./f_nlgeom.h"
#include "./ma57_solver.h" // NOTE: EICR seems unsymmtric (a speed of ma57 ~ iterative solver.)

// VectorXd f_nlgeom(std::vector<Material_ABD> & material, struct Force & force,class FEAMesh & feaMesh, 
//                     struct OPTION & option, MatrixXd & GU_u, VectorXd & GU_Rv){
    std::vector<GptsCompl> f_nlgeom(std::vector<Material_ABD> & material, struct Force & force,class FEAMesh & feaMesh, 
                    struct OPTION & option, MatrixXd & GU_u, VectorXd & GU_Rv, MatrixXd& p_Adjoint6){
    int nNODE = feaMesh.NODE.rows();
    int nELEM = feaMesh.ELEM.rows();
    int nDOF  = nNODE*feaMesh.dpn;

    MatrixXd GF_dof = force.fix.transpose();
    VectorXd GF_dof_map = Map<VectorXd>(GF_dof.data(),nDOF);

    // initialize deformations
    MatrixXd xm_u = MatrixXd::Zero(nNODE*3,1);
    MatrixXd xm_R = MatrixXd::Identity(3,3).replicate(nNODE,1);
    MatrixXd xm_Rv = EICR_SHELL::fmat2store(xm_R);

    // options for N-R solver
    int MAX_ITER = 32;
    int MAX_STEP = 1e3;
    double RTOL = 1e-3;
    double XTOL = 1e-2;
    double CTOL = 1e-4;
    int initper = option.initper;

    // dof for u
    // VectorXi id_DOF = VectorXi::LinSpaced(nDOF,0,nDOF-1);    
    MatrixXd id_all = MatrixXd::Zero(nNODE,feaMesh.dpn);//Map<MatrixXi>(id_DOF.data(),feaMesh.dpn,nNODE);
    // VectorXi id_u   = Map<VectorXi>(id_all.topRows(3).data(),3*nNODE);
    // VectorXi id_w   = Map<VectorXi>(id_all.bottomRows(3).data(),3*nNODE);

    // std::cout << id_all << "\n=====\n" << std::endl;
    // std::cout << id_u << "\n=====\n" << std::endl;
    // std::cout << id_w << "\n=====\n" << std::endl;

    // initialize a loop
    int istep = 0;
    double del_lambdaR = 1.0/initper;
    double lambdaR = 0;
    double lambdaR0 = 0;
    VectorXd del_xm = VectorXd::Zero(nDOF);
    static double R0 = 0;
    MatrixXd GU_R = xm_R;

    if (option.logflag){
        std::cout << "====================================================" << std::endl;
        std::cout << "istep     iITER   R_residual  X_del       load_ratio" << std::endl;
        std::cout << "====================================================" << std::endl;
    }

    char fname[32];
    std::ofstream file;

    Force F_target;
    
    // SparseMatrix<double> GKT(nDOF,nDOF); 
    // VectorXd Res(nDOF); 

    VectorXd del_xm_u = VectorXd::Zero(3*nNODE), del_xm_w = VectorXd::Zero(3*nNODE);
    EICR_SHELL eicr_shell(feaMesh);
    // SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > spsolver;
    ConjugateGradient<SparseMatrix<double> > spsolver; 
    /*
    // WIP: ma57
    MA57Solver ma57(true, false);
    */
    while (istep < MAX_STEP){

        MatrixXd xm_u3 = GU_u; //Map<MatrixXd>(GU_u.transpose(),3,nNODE).transpose();
        // MatrixXd xm_u  = Map<MatrixXd>(xm_u3.transpose().data(),nNODE*3,1); //SEVERE ERROR
        MatrixXd xm_u3_tmp = GU_u.transpose();
        VectorXd xm_u = Map<VectorXd>(xm_u3_tmp.data(),nNODE*3,1);
        
        VectorXd xm_Rv = GU_Rv;
        
        xm_R  = EICR_SHELL::fstore2mat(xm_Rv); //error 

        if (std::abs(lambdaR-1) < 1e-3){ //RTOL){
            std::cout << "calculating Adjoint..." << std::endl;
            eicr_shell.assembly(GU_u, GU_Rv, material, force);
            eicr_shell.sGKT.makeCompressed();
            spsolver.compute(eicr_shell.sGKT);
            /*
            WIP: ma57
            // SparseMatrix<double> K_tri = eicr_shell.sGKT.triangularView<Lower>();
            // ma57.compute(K_tri);
            // VectorXd p_Adjoint = ma57.solve(GF_dof_map);
            */
            VectorXd p_Adjoint = spsolver.solve(GF_dof_map); 

            // std::cout << "calculating sensitivity..." << std::endl;
            std::cout << "job finished" << std::endl;
            // std::vector<GptsCompl> gaussSens = eicr_shell.get_GaussCompl(GU_u, GU_R, material, force);

            p_Adjoint6 = Map<MatrixXd>(p_Adjoint.data(),6,nNODE).transpose();  
            MatrixXd GU_u_adj = p_Adjoint6.leftCols(3);            
            MatrixXd GU_w_adj = p_Adjoint6.rightCols(3);
            
            VectorXd xm_R_adj = VectorXd::Zero(3*nNODE);
            // Map<MatrixXd> xm_R_adj(GU_w_adj.data(), nNODE * 3, 1);
            
            MatrixXd R_identity = MatrixXd::Identity(3,3).replicate(nNODE,1);

            for (int uu = 0; uu < nNODE; ++uu){
                xm_R_adj[3*uu]      = GU_w_adj(uu,0);
                xm_R_adj[3*uu+1]    = GU_w_adj(uu,1);
                xm_R_adj[3*uu+2]    = GU_w_adj(uu,2);
            }

            fRupdate(xm_R_adj, nNODE, R_identity);
            std::vector<GptsCompl> gaussSens = eicr_shell.get_GaussCompl(GU_u_adj, R_identity, material, force);
            
            return gaussSens;
            // return p_Adjoint;
        }

        if ((istep > 0) && (del_lambdaR < CTOL)){
            std::cout << "step increment is smaller than CTOL" << std::endl;
        }

        if ((istep > 0) && (lambdaR0 + del_lambdaR >= 1)){ // last step
            lambdaR = 1;
            del_lambdaR = 1 - lambdaR0;
        }
        else{
            lambdaR = lambdaR0 + del_lambdaR;
        }

        for (int iITER = 0; iITER < MAX_ITER; iITER++){
            
            F_target.NM = force.NM*lambdaR;
            F_target.fix = force.fix*lambdaR;
            // f_EICR_shell(feaMesh, xm_u3, xm_Rv, material, F_target, GKT, Res); 
            eicr_shell.assembly(xm_u3, xm_Rv, material, F_target);

            if (iITER == 0) R0 = eicr_shell.Res.norm();

            if ((iITER > 0) && ((eicr_shell.Res.norm()/R0 < RTOL) || (del_xm.norm() < XTOL))){
                if (option.logflag) std::cout << "=====================================================" << std::endl;

                lambdaR0 = lambdaR;
                GU_u = xm_u3; GU_Rv = xm_Rv; GU_R = xm_R;
                
                if (option.saveflag){
                    snprintf(fname, sizeof(char)*32, "step_%i.txt", istep);
                    file.open(fname);
                    file << "lambdaR = " << lambdaR0 << "\n";
                    file << "GU_u: \n";
                    file << GU_u;
                    file << "\n=================\n";

                    file << "GU_Rv: \n";
                    file << GU_Rv;
                    file << "\n=================\n";

                    file << "Force.NM: \n";
                    file << F_target.NM;
                    file << "\n=================\n";

                    file << "Force.fix: \n";
                    file << F_target.fix;
                    file << "\n=================\n";

                    file.close();
                }

                if (iITER <= 4) del_lambdaR = del_lambdaR*2;

                istep += 1;
                break;
            }

            if (iITER == MAX_ITER-1){ // too many iterations
                lambdaR = lambdaR0;
                del_lambdaR *= 0.5;
                if (option.logflag){
                    std::cout << " ----------------- MAX. ITER. NUM ------------------" << std::endl;
                }
            break;
            }
            eicr_shell.sGKT.makeCompressed();
            spsolver.compute(eicr_shell.sGKT);
            /*
            WIP: ma57
            // SparseMatrix<double> K_tri = eicr_shell.sGKT.triangularView<Lower>();
            // ma57.compute(K_tri);
            */
            try{
                // del_xm = ma57.solve(eicr_shell.Res);
                del_xm = spsolver.solve(eicr_shell.Res); 
            }
            catch (...) {
                if (option.logflag){
                    std::cout << "--------------------- BAD MAT. CONDITION ------------------------" << std::endl;
                }
                lambdaR = lambdaR0;
                del_lambdaR = del_lambdaR/2;
                break;                
            }

            for (int uu = 0; uu < nNODE; ++uu){
                del_xm_u[3*uu]      = del_xm(6*uu);
                del_xm_u[3*uu+1]    = del_xm(6*uu+1);
                del_xm_u[3*uu+2]    = del_xm(6*uu+2);

                del_xm_w[3*uu]      = del_xm(6*uu+3);
                del_xm_w[3*uu+1]    = del_xm(6*uu+4);
                del_xm_w[3*uu+2]    = del_xm(6*uu+5);
            }
            xm_u += del_xm_u;
            fRupdate(del_xm_w, nNODE, xm_R);

            xm_Rv = EICR_SHELL::fmat2store(xm_R);
            xm_u3 = Map<MatrixXd>(xm_u.data(),3,nNODE).transpose();
            if (option.logflag){
                std::cout << istep  << "\t" << iITER << "\t" << eicr_shell.Res.norm() << "\t" << lambdaR << "\n";
            }
            feaMesh.to_vtk(xm_u3);
        }   
    }
}

void fRupdate(VectorXd & xm_w, int nNODE, MatrixXd& Rnew){
        
    MatrixXd spin; 

    for (int ii = 0; ii < nNODE ; ++ii){
        // Vector3i id_dof; id_dof << ii*3, ii*3+1, ii*3+2;
        VectorXd we = xm_w.middleRows(ii*3,3);
        double n_we = we.norm();
        Matrix3d delR;
        Matrix3d spin;
        spin << 0, -we(2), we(1), we(2), 0, -we(0), -we(1), we(0), 0;
        if (n_we < 0.05){
            delR.noalias() = MatrixXd::Identity(3,3) + (1 - n_we*n_we/6 + pow(n_we,4)/120 - pow(n_we,6)/5040.)*spin
                + 1./2*(1-n_we*n_we/12 + pow(n_we,4)/360)*spin*spin;
        }
        else{
            delR.noalias() = MatrixXd::Identity(3,3) + (std::sin(n_we)/n_we)*spin
            + 2*std::pow((std::sin(n_we/2)/n_we),2)*spin*spin;
        }
        Rnew.middleRows(ii*3,3) = delR*Rnew.middleRows(ii*3,3);
    }
}