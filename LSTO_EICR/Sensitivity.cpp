// this uses cauchy stress and p-induced strain
// Jun21: add jacobian computation
// Jun21_both_Jacob: add jacobian to linear... & x_R2 remove
// Jun22 : only 2nd
#include "./Sensitivity.h"
#include "./f_DKT_OPT.h"
Sensitivity::Sensitivity(FEAMesh &feaMesh_, std::vector<Material_ABD> &material_, MatrixXd &GU_u_, VectorXd &GU_Rv_, VectorXd &p_Adjoint_)
    : feaMesh(feaMesh_), material(material_), p_Adjoint(p_Adjoint_), GU_u(GU_u_), GU_Rv(GU_Rv_), EICR_SHELL(feaMesh_, 0)
{
    // this->GU_u = GU_u_;
    this->GU_R = fstore2mat(GU_Rv);

    unsigned int nGPTS;
    nGPTS = feaMesh.isOverlaid == true ? feaMesh.ELEM.rows() * 4 : feaMesh.ELEM.rows();

    Gpts.resize(nGPTS, 2);              // TODO: assuming global domain is flat;
    GptsSensitivities.resize(nGPTS, 4); // assume (nth term, sum)
    result_cr.resize(nGPTS, 6);         // thickness integration (N, M)
    strain_cr.resize(nGPTS, 6);         // eps, kappa
}

void Sensitivity::ComputeComplianceSensitivities(double Tolerance)
{
    // this->GU_u = GU_u;
    // this->GU_R = fstore2mat(GU_Rv);

    int nGPTS = Gpts.rows();
    DKT_OPT dkt_opt;

#if __DEBUGFLAG__
    std::ofstream sigma_log;
    sigma_log.open("stresses.txt");
    sigma_log << "position x \t position y \t sigma_ \t sigma_at_T0 \n";

    std::ofstream delP_log;
    delP_log.open("strain.txt");
    delP_log << "position x \t position y \t strain_ \t del_P_T0 \n";
#endif

    int counter = 0; // counter for each element
    for (unsigned int ee = 0; ee < feaMesh.ELEM.rows(); ++ee)
    {
        elem_id0 = feaMesh.ELEM.row(ee);

        // if (feaMesh.areafraction[ee] < Tolerance) continue;

        Vector3i elem_id;
        for (unsigned int kk = 0; kk < elem_order.rows(); ++kk)
        {
            for (unsigned int ppp = 0; ppp < 3; ++ppp)
            {
                elem_id(ppp) = elem_id0(elem_order(kk, ppp));
            }

            Matrix<double, 3, 3> X, u, x;
            for (int mmm = 0; mmm < 3; ++mmm)
            {
                X.col(mmm) = feaMesh.NODE.row(elem_id(mmm)).transpose();
                u.col(mmm) = GU_u.row(elem_id(mmm)).transpose();
            }
            x = X + u;
            double Jacobian;
            Jacobian = getArea(x) / getArea(X);
            Matrix<double, 9, 3> Ra;
            for (int nnn = 0; nnn < 3; ++nnn)
            {
                for (int mmm = 0; mmm < 3; ++mmm)
                {
                    Ra.row(3 * nnn + mmm) = GU_R.row(elem_id(nnn) * 3 + mmm);
                }
            }

            Matrix3d X_R, x_R;
            Matrix3d u_d;
            Matrix<double, 9, 1> th_d;
            Matrix<double, 3, 3> Q;
            Matrix3d T0, T;

            Filter_Def(X, u, Ra, X_R, x_R, u_d, th_d, T0, T);

            // Material of local element (stress and strain being voigt notation: note that shear = 2*eps12)
            Q << T0(0, 0) * T0(0, 0), T0(0, 1) * T0(0, 1), T0(0, 0) * T0(0, 1),
                T0(1, 0) * T0(1, 0), T0(1, 1) * T0(1, 1), T0(1, 0) * T0(1, 1),
                2 * T0(0, 0) * T0(1, 0), 2 * T0(0, 1) * T0(1, 1), T0(0, 0) * T0(1, 1) + T0(0, 1) * T0(1, 0);

            // from corotational to initial (not global)
            Matrix3d Rdot = T.transpose() * T0; // JUN10th. Final (R 계산)
            Matrix3d Rotator;
            Rotator << Rdot(0, 0) * Rdot(0, 0), Rdot(1, 0) * Rdot(1, 0), 2 * Rdot(0, 0) * Rdot(1, 0),
                Rdot(0, 1) * Rdot(0, 1), Rdot(1, 1) * Rdot(1, 1), 2 * Rdot(0, 1) * Rdot(1, 1),
                Rdot(0, 0) * Rdot(0, 1), Rdot(1, 0) * Rdot(1, 1), Rdot(0, 0) * Rdot(1, 1) + Rdot(0, 1) * Rdot(1, 0);
            // Rotator = Rdot; // Jun9 + noRR

            MatrixXd Q2(6, 6);
            Q2 << Rotator, MatrixXd::Zero(3, 3), MatrixXd::Zero(3, 3), Rotator; // cauchy stress = R'*sigma_CR*R // Zero: jun13_Zero (before=Iden)

            Matrix<double, 6, 1> Fnm_zeros;
            Fnm_zeros.fill(0);

            Material_ABD Mater_e;
            Mater_e.Amat = Q * material[ee].Amat * Q.transpose();
            Mater_e.Dmat = Q * material[ee].Dmat * Q.transpose();
            Mater_e.Bmat = Q * material[ee].Bmat * Q.transpose();

            FilteredP p_d; // corotational coordinate
            p_d.membrane << u_d(0), u_d(1), th_d(2), u_d(3), u_d(4), th_d(5), u_d(6), u_d(7), th_d(8);
            p_d.bending << u_d(2), th_d(0), th_d(1), u_d(5), th_d(3), th_d(4), u_d(8), th_d(6), th_d(7);

            FilteredP adj_e; // Global coordinate -> TODO: to initial coord. T0
            Matrix<double, 9, 1> adj_undf1_, adj_undf2_;
            // // u, v, w
            // adj_undf1_.middleRows(0,3) = T0*p_Adjoint.middleRows(elem_id(0)*6,3);
            // adj_undf1_.middleRows(3,3) = T0*p_Adjoint.middleRows(elem_id(1)*6,3);
            // adj_undf1_.middleRows(6,3) = T0*p_Adjoint.middleRows(elem_id(2)*6,3);
            // // thetax thetay thetaz
            // adj_undf2_.middleRows(0,3) = T0*p_Adjoint.middleRows(elem_id(0)*6+3,3);
            // adj_undf2_.middleRows(3,3) = T0*p_Adjoint.middleRows(elem_id(1)*6+3,3);
            // adj_undf2_.middleRows(6,3) = T0*p_Adjoint.middleRows(elem_id(2)*6+3,3);

            // u, v, w
            adj_undf1_.middleRows(0, 3) = T * p_Adjoint.middleRows(elem_id(0) * 6, 3); // TEST for JUN15_n
            adj_undf1_.middleRows(3, 3) = T * p_Adjoint.middleRows(elem_id(1) * 6, 3);
            adj_undf1_.middleRows(6, 3) = T * p_Adjoint.middleRows(elem_id(2) * 6, 3);
            // thetax thetay thetaz
            adj_undf2_.middleRows(0, 3) = T * p_Adjoint.middleRows(elem_id(0) * 6 + 3, 3);
            adj_undf2_.middleRows(3, 3) = T * p_Adjoint.middleRows(elem_id(1) * 6 + 3, 3);
            adj_undf2_.middleRows(6, 3) = T * p_Adjoint.middleRows(elem_id(2) * 6 + 3, 3);

            adj_e.membrane << adj_undf1_(0), adj_undf1_(1), adj_undf2_(2), adj_undf1_(3), adj_undf1_(4), adj_undf2_(5), adj_undf1_(6), adj_undf1_(7), adj_undf2_(8);
            adj_e.bending << adj_undf1_(2), adj_undf2_(0), adj_undf2_(1), adj_undf1_(5), adj_undf2_(3), adj_undf2_(4), adj_undf1_(8), adj_undf2_(6), adj_undf2_(7);

            Matrix3d X_R2, X_2;
            X_R2.noalias() = T * T0.transpose() * X_R; // Jun22_filterX in {T} coordinate (no deformation)
            X_2.noalias() = T * X;                     // non-centered X
            Matrix<double, 2, 3> xycoord = X_R2.topRows(2);
            Matrix<double, 2, 3> xycoord0 = X_2.topRows(2);

            // NOTE: this is not global coord, and must consider material anisotropy direction
            CoreElement coreelem0(xycoord0, Mater_e);
            CoreElement coreelem(xycoord, Mater_e);

            Matrix<double, 1, 6> stress_tmp, strain_tmp;

            Gpts.row(counter) = coreelem.get_StressStrain(p_d, Fnm_zeros, 1. / 3.0, 1. / 3.0, stress_tmp, strain_tmp, X); // Deformed coordinate, filtered p_d
            //result_cr.row(counter) = stress_tmp * feaMesh.areafraction[ee];
            //strain_cr.row(counter) = strain_tmp;

            Matrix<double, 1, 6> stress_adjoint, del_adjoint; //, gradstress; // gradstress는 isoparametric 가정으로 사라짐
            // coreelem0.get_StressStrain(adj_e, Fnm_zeros, 1./3.0, 1./3.0, stress_adjoint, del_adjoint, X); // Undeformed coordinate
            coreelem0.get_StressStrain(adj_e, Fnm_zeros, 1. / 3.0, 1. / 3.0, stress_adjoint, del_adjoint, X); // double_jacob

            Matrix<double, 1, 6> stress_initial;
            // stress_initial.noalias() = (stress_tmp*Q2.transpose()).transpose(); // Q2*stress = reference state
            // stress_initial = stress_tmp; // Sens1
            stress_initial = Jacobian * stress_tmp; // Jun21_Jac

#if __DEBUGFLAG__
            sigma_log << Gpts.row(counter) << ": " << stress_tmp << "\t ||" << stress_initial << "\n";
            delP_log << Gpts.row(counter) << ": " << strain_tmp << "\t ||" << del_adjoint << "\n";
#endif

            GptsSensitivities(counter, 0) = stress_initial.dot(strain_tmp) * feaMesh.areafraction[ee]; // Jacobian - double_jacob //stress_tmp.dot(strain_tmp) * feaMesh.areafraction[ee];
            GptsSensitivities(counter, 1) = stress_initial.dot(del_adjoint) * feaMesh.areafraction[ee];
            GptsSensitivities(counter, 2) = GptsSensitivities(counter, 1);
            GptsSensitivities(counter, 3) = -GptsSensitivities(counter, 1);

            // increase counter
            counter += 1;
        }
    }
#if __DEBUGFLAG__
    sigma_log.close();
    delP_log.close();
#endif
}

void Sensitivity::to_gptSens(bool isLinear)
{
    MatrixXi elem_id0, elem_order;

    if (feaMesh.isOverlaid == true)
    {
        elem_id0 = MatrixXi::Zero(4, 1);
        elem_order = MatrixXi::Zero(4, 3);
        elem_order << 0, 1, 2, 1, 2, 3, 2, 3, 0, 3, 0, 1;

        gptsSens.resize(feaMesh.ELEM.rows() * 4);
    }
    else
    {
        elem_order << 0, 1, 2;
        gptsSens.resize(feaMesh.ELEM.rows());
    }

    for (unsigned int ee = 0; ee < feaMesh.ELEM.rows(); ++ee)
    {
        elem_id0 = feaMesh.ELEM.row(ee);
        for (unsigned int kk = 0; kk < elem_order.rows(); ++kk)
        {
            Vector3i elem_id;
            for (unsigned int ppp = 0; ppp < 3; ++ppp)
            {
                elem_id(ppp) = elem_id0(elem_order(kk, ppp));
            }

            Matrix<double, 3, 3> X;
            for (int mmm = 0; mmm < 3; ++mmm)
            {
                X.col(mmm) = feaMesh.NODE.row(elem_id(mmm)).transpose();
            }

            gptsSens[ee * elem_order.rows() + kk].x = (X(0, 0) + X(0, 1) + X(0, 2)) / 3.0;
            gptsSens[ee * elem_order.rows() + kk].y = (X(1, 0) + X(1, 1) + X(1, 2)) / 3.0;
            gptsSens[ee * elem_order.rows() + kk].z = (X(2, 0) + X(2, 1) + X(2, 2)) / 3.0;
            if (isLinear == false)
            {
                gptsSens[ee * elem_order.rows() + kk].sens = GptsSensitivities(ee * elem_order.rows() + kk, 2);
            }
            else
            {
                gptsSens[ee * elem_order.rows() + kk].sens = GptsSensitivities(ee * elem_order.rows() + kk, 0);
            }
        }
    }
}

/*
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
        // initial domaincounter
        el_cood << feaMesh.Centeroids(ee,0), feaMesh.Centeroids(ee,1);

        if (feaMesh.areafraction[ee] > Tolerance){
            el_dist = std::sqrt(std::pow(Pointxy[0] - el_cood[0],2)+std::pow(Pointxy[1] - el_cood[1],2));
            if (el_dist < 1.5*Radius){
                for (int gg = 0; gg < nGpts; ++gg){
                    unsigned int gg_indx = nGpts*ee +  gg; // ee *= nGpts 로 인식하는 것 같다.. eigen lib. 문제 Jun20
                    gg_cood = Gpts.row(gg_indx);
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
        std::cout << "a very small island is found at (" << Pointxy[0] << ", " << Pointxy[1] << "): npts = " << CntPoints << std::endl;
        PointSensitivity = 0.0;
        return PointSensitivity;
    }

    MatrixXd A(CntPoints, 6);
    VectorXd b_sens(CntPoints);
    double BMax = 1e20, BMin = -1e20;
	std::vector <double> RelativeCoordinate (2);

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
*/