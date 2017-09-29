// this script uses del_p and 1st PK, which are nonsymm.
#include "./Sensitivity.h"
#include "./f_DKT_OPT.h"
Sensitivity::Sensitivity(FEAMesh & feaMesh_, std::vector<Material_ABD> & material_, MatrixXd & GU_u_, VectorXd & GU_Rv_, VectorXd & p_Adjoint_) 
: feaMesh(feaMesh_), material(material_), p_Adjoint(p_Adjoint_), GU_u(GU_u_), GU_Rv(GU_Rv_), EICR_SHELL(feaMesh_, 0) {
    // this->GU_u = GU_u_;
    this->GU_R = fstore2mat(GU_Rv);

    unsigned int nGPTS;
    nGPTS = feaMesh.isOverlaid == true ? feaMesh.ELEM.rows()*4 : feaMesh.ELEM.rows();
    
    Gpts.resize(nGPTS,2); // TODO: assuming global domain is flat; 
    GptsSensitivities.resize(nGPTS,4); // assume (nth term, sum)
    result_cr.resize(nGPTS,6); // thickness integration (N, M)
    strain_cr.resize(nGPTS,6); // eps, kappa

}

void Sensitivity::ComputeComplianceSensitivities(double Tolerance){ 
    // this->GU_u = GU_u;
    // this->GU_R = fstore2mat(GU_Rv);

    int nGPTS = Gpts.rows();
    DKT_OPT dkt_opt;

    int counter = 0; // counter for each element 
    for (unsigned int ee = 0; ee < feaMesh.ELEM.rows(); ++ee){
        elem_id0 = feaMesh.ELEM.row(ee); 

        if (feaMesh.areafraction[ee] < Tolerance) continue;       
        
        Vector3i elem_id; 
        for (unsigned int kk = 0; kk < elem_order.rows(); ++kk){
            for (unsigned int ppp = 0; ppp < 3; ++ppp ){
                elem_id(ppp) = elem_id0(elem_order(kk,ppp));
            }

            Matrix<double,3,3> X, u;//, x; 
            for (int mmm = 0; mmm < 3; ++mmm){
                X.col(mmm) = feaMesh.NODE.row(elem_id(mmm)).transpose();
                u.col(mmm) = GU_u.row(elem_id(mmm)).transpose();
            }

            Matrix<double,9,3> Ra;
            for (int nnn = 0; nnn < 3; ++nnn){
                for (int mmm = 0; mmm < 3; ++mmm){
                    Ra.row(3*nnn+mmm) = GU_R.row(elem_id(nnn)*3+mmm);
                }    
            }            

            Matrix3d X_R, x_R;
            Matrix3d u_d;
            Matrix<double,9,1> th_d;
            Matrix<double,3,3> Q;
            Matrix3d T0, T;
            
            Filter_Def(X, u, Ra, X_R, x_R, u_d, th_d, T0, T);

            // Material of local element (stress and strain being voigt notation: note that shear = 2*eps12)
            Q << T0(0,0)*T0(0,0), T0(0,1)*T0(0,1), T0(0,0)*T0(0,1),
            T0(1,0)*T0(1,0), T0(1,1)*T0(1,1), T0(1,0)*T0(1,1),
            2*T0(0,0)*T0(1,0), 2*T0(0,1)*T0(1,1), T0(0,0)*T0(1,1)+T0(0,1)*T0(1,0);

            // from corotational to initial (not global)
            Matrix3d Rdot = T.transpose()*T0; // JUN10th. Final (R 계산)
            // Matrix3d Rdot = T*T0.transpose(); // JUN9th 
            Matrix3d Rotator;
            // Rotator <<  Rdot(0,0)*Rdot(0,0), Rdot(1,0)*Rdot(1,0), 2*Rdot(0,0)*Rdot(1,0),
            //             Rdot(0,1)*Rdot(0,1), Rdot(1,1)*Rdot(1,1), 2*Rdot(0,1)*Rdot(1,1), 
            //             Rdot(0,0)*Rdot(0,1), Rdot(1,0)*Rdot(1,1), Rdot(0,0)*Rdot(1,1)+Rdot(0,1)*Rdot(1,0);
            Rotator = Rdot;

            MatrixXd Q2(6,6);
            Q2 << Rotator, MatrixXd::Identity(3,3), MatrixXd::Identity(3,3), Rotator; 

            // rotator R을 쓰는게 맞는데 양 변에 곱해줘야 한다는 사실을 잊고 있었다. .// JUN8
            // 1st PK로 formulation되어있기 때문에 .. JUN9 rotated 랑 reference 간 mapping 이므로 일단 Jacobian은 1로, F는 R로 근사.

            Matrix<double,6,1> Fnm_zeros;
            Fnm_zeros.fill(0);
            
            Material_ABD Mater_e;
            Mater_e.Amat = Q*material[ee].Amat*Q.transpose();
            Mater_e.Dmat = Q*material[ee].Dmat*Q.transpose();
            Mater_e.Bmat = Q*material[ee].Bmat*Q.transpose();

            FilteredP p_d; // corotational coordinate
            p_d.membrane << u_d(0), u_d(1), th_d(2), u_d(3), u_d(4), th_d(5), u_d(6), u_d(7), th_d(8);
            p_d.bending  << u_d(2), th_d(0), th_d(1), u_d(5), th_d(3), th_d(4), u_d(8), th_d(6), th_d(7);

            FilteredP adj_e; // Global coordinate -> TODO: to initial coord. T0
            Matrix<double,9,1> adj_undf1_, adj_undf2_;
            // u, v, w
            adj_undf1_.middleRows(0,3) = T0*p_Adjoint.middleRows(elem_id(0)*6,3); 
            adj_undf1_.middleRows(3,3) = T0*p_Adjoint.middleRows(elem_id(1)*6,3); 
            adj_undf1_.middleRows(6,3) = T0*p_Adjoint.middleRows(elem_id(2)*6,3);
            // thetax thetay thetaz
            adj_undf2_.middleRows(0,3) = T0*p_Adjoint.middleRows(elem_id(0)*6+3,3);
            adj_undf2_.middleRows(3,3) = T0*p_Adjoint.middleRows(elem_id(1)*6+3,3);
            adj_undf2_.middleRows(6,3) = T0*p_Adjoint.middleRows(elem_id(2)*6+3,3);

            // adj_e.membrane << adj_undf1_(0), adj_undf1_(1), adj_undf2_(2), adj_undf1_(3), adj_undf1_(4), adj_undf2_(5), adj_undf1_(6), adj_undf1_(7), adj_undf2_(8);
            // adj_e.bending  << adj_undf1_(2), adj_undf2_(0), adj_undf2_(1), adj_undf1_(5), adj_undf2_(3), adj_undf2_(4), adj_undf1_(8), adj_undf2_(6), adj_undf2_(7);

            Matrix<double,2,3> xycoord0 = X_R.topRows(2);
            Matrix<double,2,3> xycoord  = x_R.topRows(2);

            CoreElement coreelem0(xycoord0, Mater_e); // NOTE: this is not global coord, and must be consider material anisotropy direction 
            CoreElement coreelem (xycoord,  Mater_e);

            Matrix<double,1,6>  stress_tmp, strain_tmp;

            // std::cout << coreelem.get_StressStrain(p_d, Fnm_zeros, 1./3.0, 1./3.0, stress_tmp, strain_tmp, X) << std::endl;
            Gpts.row(counter) = coreelem.get_StressStrain(p_d, Fnm_zeros, 1./3.0, 1./3.0, stress_tmp, strain_tmp, X);
            //result_cr.row(counter) = stress_tmp * feaMesh.areafraction[ee];
            //strain_cr.row(counter) = strain_tmp;

            // TODO :: formulation 이 부분은 1st PK 구하듯이 해야함
            // Matrix<double,1,6>  stress_initial;
            // stress_initial.noalias() = (stress_tmp*Q2).transpose(); // JUN9
            // stress_initial.noalias() = (stress_tmp*Q2.transpose()).transpose();
            
            Matrix<double, 2,2 > N_1PK, M_1PK, del_p_mem, del_p_ben;
            N_1PK << stress_tmp[0], stress_tmp[2], stress_tmp[2], stress_tmp[1];
            M_1PK << stress_tmp[3], stress_tmp[5], stress_tmp[5], stress_tmp[4];
            N_1PK *= Rdot.topLeftCorner(2,2); M_1PK *= Rdot.topLeftCorner(2,2);
            
            del_p_mem.fill(0); del_p_ben.fill(0);
            // F_ij = x,i/X,j 여야 rotation 정의랑 맞음. (1st PK 정의는? 그러므로 del_p 도 pi,j인것으로..)
            Matrix<double,1,6>  stress_adjoint, del_adjoint;//, gradstress; // gradstress는 isoparametric 가정으로 사라짐
            
            // nabla_p 계산을 위해서 따로따로 구분해서 Bmat 얻는 과정. theta_z = u,y - v,x이고, w_x&w_y는 theta랑 연관이 있으니까 무시하자. 여기서는 strain이랑 curvature만 감안 (CST, KT로 reduce)
            // only u, theta_x
            adj_e.membrane << adj_undf1_(0), 0, 0, adj_undf1_(3), 0, 0, adj_undf1_(6), 0, 0;
            adj_e.bending  << adj_undf1_(2), 0, 0, adj_undf1_(5), 0, 0, adj_undf1_(8), 0, 0;
            coreelem0.get_StressStrain(adj_e, Fnm_zeros, 1./3.0, 1./3.0, stress_adjoint, del_adjoint, X);
            del_p_mem(0,0) = del_adjoint(0); del_p_mem(0,1) = del_adjoint(2); 
            del_p_ben(0,0) = del_adjoint(3); del_p_ben(0,1) = del_adjoint(5); 

            // only v, theta_y
            adj_e.membrane << 0, adj_undf1_(1), 0, 0, adj_undf1_(4), 0, 0, adj_undf1_(7), 0;
            adj_e.bending  << 0, adj_undf2_(0), 0, 0, adj_undf2_(3), 0, 0, adj_undf2_(6), 0;
            coreelem0.get_StressStrain(adj_e, Fnm_zeros, 1./3.0, 1./3.0, stress_adjoint, del_adjoint, X);

            del_p_mem(1,0) = del_adjoint(0); del_p_mem(1,1) = del_adjoint(2); 
            del_p_ben(1,0) = del_adjoint(3); del_p_ben(1,1) = del_adjoint(5); 
            
            // coreelem0.get_StressStrain(adj_e, Fnm_zeros, 1./3.0, 1./3.0, stress_adjoint, del_adjoint, X);
            N_1PK *= del_p_mem.transpose(); M_1PK *= del_p_ben.transpose();
            //coreelem0.get_gradStress(stress_tmp, 1./3.0, 1./3.0, gradstress);
            double sum_cnt2 = N_1PK.trace() + M_1PK.trace();

            GptsSensitivities(counter,0) = stress_tmp.dot(strain_tmp) * feaMesh.areafraction[ee];
            GptsSensitivities(counter,1) = -sum_cnt2 * feaMesh.areafraction[ee]; //-stress_initial.dot(strain_adjoint) * feaMesh.areafraction[ee];
            GptsSensitivities(counter,2) = GptsSensitivities(counter,1);//-gradstress.dot(strain_tmp) * feaMesh.areafraction[ee];  // this should be similar to [1]
            GptsSensitivities(counter,3) = GptsSensitivities(counter,0) + GptsSensitivities(counter,1) + GptsSensitivities(counter,2); // ,2);
            
            // increase counter
            counter += 1;  
        } 
    }
}

double Sensitivity::ComputeBoundaryPointSensitivity(std::vector<double> & Pointxy, double Radius, unsigned int WeightFlag, bool isLinear, double Tolerance){
    unsigned int Counter = 0;
    unsigned int p;
    double el_dist;
    int nGpts;
    double temp;

    if (feaMesh.isOverlaid == true){
        p = ((PI * Radius * Radius) / feaMesh.ElemArea[0]) * 16; 
        nGpts = 4;
    }        
    else{
        p = ((PI * Radius * Radius) / feaMesh.ElemArea[0]); 
        nGpts = 1;
    } 
    p *= 1.25; // for a conservative estimate
    std::vector<double> Distances(p);
    std::vector<unsigned int> ElementIndices(p);
    std::vector<unsigned int> Indices(p);
    
    double PointSensitivity = 0.0;
    
    int CntPoints = 0;
    
    Vector2d el_cood, gg_cood;
    double gg_dist;
    for (unsigned int ee = 0; ee < feaMesh.ELEM.rows(); ++ee){
        // initial domain
        el_cood << feaMesh.Centeroids(ee,0), feaMesh.Centeroids(ee,1);

        if (feaMesh.areafraction[ee] > Tolerance){
            el_dist = std::sqrt(std::pow(Pointxy[0] - el_cood[0],2)+std::pow(Pointxy[1] - el_cood[1],2));
            if (el_dist < 1.5*Radius){
                for (int gg = 0; gg < nGpts; ++gg){
                    gg_cood = Gpts.row(nGpts*ee + gg);
                    gg_dist = std::sqrt(std::pow(Pointxy[0] - gg_cood[0],2)+std::pow(Pointxy[1] - gg_cood[1],2));
                    if (gg_dist < Radius){
                        Distances[CntPoints] = gg_dist;
                        ElementIndices[CntPoints] = ee;
                        Indices[CntPoints] = gg;
                        CntPoints += 1;
                    }
                }
            }            
        }
    }

    if (CntPoints < 10){
        std::cout << "a very small island is found" << std::endl;
        // std::fill(PointSensitivities.begin(), PointSensitivities.end(), 0.0);
        PointSensitivity = 0.0;
        return PointSensitivity;
    }

	// std::vector <double> A (6*CntPoints);
	// std::vector <double> B (1 * CntPoints);
	// std::vector <double> BMax (1);
	// std::vector <double> BMin (1);
    MatrixXd A(CntPoints, 6);
    VectorXd b_sens(CntPoints);
    double BMax = 1e20, BMin = -1e20;
	std::vector <double> RelativeCoordinate (2);

    // std::fill (BMax.begin (), BMax.end (), 1e20);
	// std::fill (BMin.begin (), BMin.end (), -1e20);
    // TOFIX: spacedim = 2, nDual = 1
    for (int ii = 0; ii < CntPoints; ++ii){
        p = ElementIndices[ii]*nGpts + Indices[ii]; //Indices[ii];
        switch (WeightFlag){
            case 1:
                temp = 1.0; // least squares
                break;
            case 2:
                temp = 1.0 / Distances[ii];
                break;
            case 3:
                temp = feaMesh.areafraction[ii];
                break;
            case 4:
                temp = feaMesh.areafraction[ii]/Distances[ii];
                break;
            case 5:
                temp = std::sqrt(feaMesh.areafraction[ii]/Distances[ii]);
                break;
            default:
                temp = 1.0;
				std::cout << "Weight Flag should lie in [1, 5]. Using Least Squares.\n";
        }
        for (unsigned int jj = 0 ; jj < 2; ++ jj){
            RelativeCoordinate[jj] = Gpts(p,jj) - Pointxy[jj];
        }
        A (ii,0) = temp;
        A (ii,1) = RelativeCoordinate [0] * temp;
        A (ii,2) = RelativeCoordinate [1] * temp;
        A (ii,3) = RelativeCoordinate [0] * RelativeCoordinate [1] * temp;
        A (ii,4) = RelativeCoordinate [0] * RelativeCoordinate [0] * temp;
        A (ii,5) = RelativeCoordinate [1] * RelativeCoordinate [1] * temp;    
        
        if (isLinear == true) b_sens(ii) = -GptsSensitivities (p,0)*temp; // linear case test
        else b_sens(ii) =  GptsSensitivities (p,3)*temp;

    }
    MatrixXd B_tmp = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b_sens);
    double B = B_tmp(0);
    if ((B > BMax*10) || (B < BMin*10)){
        B = 0.0; 
        temp = 0.0;

        for (int nn = 0; nn < CntPoints; ++nn){
            B += GptsSensitivities(ElementIndices[nn]*nGpts + Indices[nn],3)*feaMesh.areafraction[ElementIndices[nn]];
            temp += feaMesh.areafraction[ElementIndices[nn]];
        }
        PointSensitivity = B/temp;
    }
    else if (B > BMax) PointSensitivity = BMax;
    else if (B < BMin) PointSensitivity = BMin;
    else PointSensitivity = B;
    
    return PointSensitivity;
}
/*
void Sensitivity::LeastSquare(unsigned int N, unsigned int M, unsigned int nRhs, std::vector <double>& A, std::vector <double>& B){
    int Info;
	char trans = 'N';
	unsigned int nBlock = 48;
	unsigned int lwork;
	if (N < M)
	{
		lwork = N;
	}
	else
	{
		lwork = M;
	}
	if (lwork > nRhs)
	{
		lwork = lwork + lwork + lwork*nBlock;
	}
	else
	{
		lwork = lwork + nRhs*nBlock;
	}
	std::vector <double> work (lwork, 0.0);
	dgels_ (&trans, &M, &N, &nRhs, &A [0], &M, &B [0], &M, &work [0], &lwork, &Info);
	if (Info != 0)
	{
		if (Info < 0)
		{
			std::cout << "Error in DGELS: " << -Info << " argument had an illegal value.\n";
		}
		else
		{
			std::cout << "Error in DGELS: " << Info << "the diagonal element of the triangular factor of A is zero. So, A doesn't have full rank. The least squares solution is terminated.\n";
		}
	}
}*/
