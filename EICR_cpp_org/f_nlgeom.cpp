#include "./f_nlgeom.h"

void f_nlgeom(std::vector<Material_ABD> & material, struct Force & force,class FEAMesh & feaMesh, 
                    struct OPTION & option, MatrixXd & GU_u, VectorXd & GU_Rv){
    int nNODE = feaMesh.NODE.rows();
    int nELEM = feaMesh.ELEM.rows();
    int nDOF  = nNODE*feaMesh.dpn;

    // initialize deformations
    MatrixXd xm_u = MatrixXd::Zero(nNODE*3,1);
    MatrixXd xm_R = MatrixXd::Identity(3,3).replicate(nNODE,1);
    MatrixXd xm_Rv = fmat2store(xm_R);

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

    Force F_target;
    
    SparseMatrix<double> GKT(nDOF,nDOF); 
    VectorXd Res(nDOF); 

    VectorXd del_xm_u = VectorXd::Zero(3*nNODE), del_xm_w = VectorXd::Zero(3*nNODE);
    
    while (istep < MAX_STEP){

        MatrixXd xm_u3 = GU_u; //Map<MatrixXd>(GU_u.transpose(),3,nNODE).transpose();
        // MatrixXd xm_u  = Map<MatrixXd>(xm_u3.transpose().data(),nNODE*3,1); //SEVERE ERROR
        MatrixXd xm_u3_tmp = GU_u.transpose();
        VectorXd xm_u = Map<VectorXd>(xm_u3_tmp.data(),nNODE*3,1);
        
        VectorXd xm_Rv = GU_Rv;
        
        xm_R  = fstore2mat(xm_Rv); //error 

        if (std::abs(lambdaR-1) < 1e-3){ //RTOL){
            // MatrixXd GKT_mat = f_EICR_shell(feaMesh, xm_u3, xm_Rv, material, F_target, GKT, Res); 
            std::cout << "job finished" << std::endl;
            return;
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
            Res = VectorXd::Zero(nDOF);
            f_EICR_shell(feaMesh, xm_u3, xm_Rv, material, F_target, GKT, Res); 

            if (iITER == 0) R0 = Res.norm();

            if ((iITER > 0) && (Res.norm()/R0 < RTOL) && (del_xm.norm() < XTOL)){
                if (option.logflag) std::cout << "=====================================================" << std::endl;

                //update
                lambdaR0 = lambdaR;
                GU_u = xm_u3; GU_Rv = xm_Rv; GU_R = xm_R;
                if (option.saveflag){
                    //TODO
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

            // continued manuscript
            // std::cout << "====GKT====" << MatrixXd(GKT) << "\n ============\n" << std::endl;
            // std::cout << "====Res====" << Res << "\n ============\n" << std::endl;

            // MatrixXd K_test = MatrixXd(GKT);
            // std::cout << "K_test : \n " << K_test << "\n==========\n"<< std::endl;
            // std::cout << "K_dense : \n " << GKT_mat << "\n==========\n"<< std::endl;
            // MA57Solver ma57 (false,false);
            // SparseMatrix<double> K_tri = GKT.triangularView<Lower>();
            // ma57.compute(K_tri);
            // LLT<MatrixXd> llt;
            // llt.compute(GKT_mat);
            // GKT_mat = (GKT_mat + GKT_mat.transpose())/2.;
            // FullPivHouseholderQR<MatrixXd> QR;
            // QR.compute(GKT_mat);
            // // TEMP! this is very slow ============
            SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > spQR;
            GKT.makeCompressed();
            spQR.compute(GKT);
            // // ==================
            // VectorXd Res_tmp = VectorXd::Zero(nDOF); Res_tmp.bottomRows(6) = VectorXd::Ones(6);
            // VectorXd tu = ma57.solve (Res_tmp); // Error : 
            // LeastSquaresConjugateGradient<SparseMatrix<double> > lscg(GKT);
            try{
                // del_xm = ma57.solve (Res); // ERROR
                // del_xm = GKT_mat.colPivHouseholderQr().solve(Res);
                // del_xm = QR.solve(Res);
                // del_xm = lscg.solve(Res);
                // // SLOW =========
                del_xm = spQR.solve(Res); 
                // // =============+
                // del_xm = llt.solve(Res); // SOLUTION DIVERGES
                // for (int kk = 0; kk < 13 ; ++kk){
                //     std::cout << "del_xm(" << kk << ") : \n " << del_xm.middleRows(kk*6,6).transpose() << "\n" << std::endl;    
                // }
                // std::cout << "del_xm(5) : \n " << del_xm.middleRows(24,6).transpose() << std::endl;
                // std::cout << "del_xm(8) : \n " << del_xm.middleRows(42,6).transpose() << std::endl;
                // std::cout << "del_xm(9) : \n " << del_xm.middleRows(48,6).transpose() << "\n==========\n"<< std::endl;
            }
            catch (...) {
                if (option.logflag){
                    std::cout << "--------------------- BAD MAT. CONDITION ------------------------" << std::endl;
                }
                lambdaR = lambdaR0;
                del_lambdaR = del_lambdaR/2;
                break;                
            }

            // del_xm_u.setZero(3*nNODE); del_xm_w.setZero(3*nNODE);
            // id_all = Map<MatrixXd>(del_xm.data(),6,nNODE);
            // del_xm_u = Map<VectorXd>(id_all.topRows(3).data(),nNODE*3); // TOFIX: SEVERE ERROR
            // del_xm_w = Map<VectorXd>(id_all.bottomRows(3).data(),nNODE*3); // ERROR SEVERE
            for (int uu = 0; uu < nNODE; ++uu){
                del_xm_u[3*uu]      = del_xm(6*uu);
                del_xm_u[3*uu+1]    = del_xm(6*uu+1);
                del_xm_u[3*uu+2]    = del_xm(6*uu+2);

                del_xm_w[3*uu]      = del_xm(6*uu+3);
                del_xm_w[3*uu+1]    = del_xm(6*uu+4);
                del_xm_w[3*uu+2]    = del_xm(6*uu+5);
            }
            xm_u += del_xm_u;
            // for (int kk = 0; kk < 13 ; ++kk){
            //         std::cout << "xm(" << kk << ") : \n " << del_xm.middleRows(kk*3,3).transpose() << "\n" << std::endl;    
            // }
            fRupdate(del_xm_w, nNODE, xm_R);

            xm_Rv = fmat2store(xm_R);
            xm_u3 = Map<MatrixXd>(xm_u.data(),3,nNODE).transpose();
            // for (unsigned int ii =0; ii < 13; ++ii){
            //     for (unsigned int jj =0; jj < 13; ++jj){
            //         if (GKT_mat(6*ii,6*jj) > 1e-3){
            //             std::cout << "GKT( " << ii << ", " << jj <<
            //                     "):\n" << GKT_mat.block(6*ii,6*jj,6,6) << "\n\n" << std::endl;
            //         }
            //     }
            // }
            // for (int kk = 0; kk < 13 ; ++kk){
            //         std::cout << "Res(" << kk << ") : \n " << Res.middleRows(kk*6,6).transpose() << "\n" << std::endl;    
            // }
            // for (int kk = 0; kk < 13 ; ++kk){
            //         std::cout << "xm_u3(" << kk << ") : \n " << xm_u3.row(kk) << "\n" << std::endl;    
            // }
            // for (int kk = 0; kk < 13 ; ++kk){
            //         std::cout << "xm_Rv(" << kk << ") : \n " << xm_Rv.middleRows(kk*9,9).transpose() << "\n" << std::endl;    
            // }
            // for (int kk = 0; kk < 13 ; ++kk){
            //         std::cout << "del_xm_w(" << kk << ") : \n " << del_xm_w.middleRows(kk*3,3).transpose() << "\n" << std::endl;    
            // }
            // if (option.logflag){
                std::cout << istep  << "\t" << iITER << "\t" << Res.norm() << "\t" << lambdaR << "\n";
            // }

        }   
    }
}
// MatrixXd fstore2mat(MatrixXd Rv){
//     int num_el = Rv.rows()/9;
//     MatrixXd store2mat; store2mat.setZero(num_el*3,3);
//     for (int ii = 0; ii < num_el; ++ ii){
//         store2mat.middleRows(ii*3,3) = Map<MatrixXd>(Rv.middleRows(ii*9,9).transpose().data(),3,3);
//     }
// }

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