#include "./f_core_element.h"

using namespace Eigen;

CoreElement f_core_element(const MatrixXd & xycoord,const Material_ABD & Mater_e, 
    const VectorXd & FNM_r){
        FilteredP p_d;
        CoreElement Core_E;
        Core_E = f_core_element(xycoord, Mater_e, FNM_r, p_d);
        return Core_E;
    }


CoreElement f_core_element(const MatrixXd & xycoord,const Material_ABD & Mater_e, 
                        const VectorXd & FNM_r,const FilteredP & p_d){
    MatrixXd Amat = Mater_e.Amat, Dmat = Mater_e.Dmat, Bmat = Mater_e.Bmat;
    VectorXd u_m = p_d.membrane, u_b = p_d.bending;
    VectorXd Nph = FNM_r.topRows(3), Mph = FNM_r.bottomRows(3);
    MatrixXd x_ = xycoord.row(0), y_ = xycoord.row(1);
    double nu = std::min(Mater_e.Amat(1,0)/Mater_e.Amat(0,1),Mater_e.Amat(1,0)/Mater_e.Amat(0,0));
    DKT_OPT dkt_opt = f_DKT_OPT(xycoord, nu); // NOTE:opt_only
    CoreElement Core_E;

    int ngp = 3;    // num of gpts
    std::vector<double> ri(ngp), si(ngp), wi(ngp);
    if (ngp == 3){
        ri = {0, 0.5, 0.5}, si = {0.5, 0, 0.5}, wi = {0.333333333, 0.333333333, 0.333333333};
    }
    else if(ngp == 1){
        ri[0] = 1.0/3, si[0] = 1.0/3, wi[0] = 1;
    }
    
    for (int gg = 0; gg < ngp; ++gg){
        double r = ri[gg], s = si[gg], w = wi[gg];
        // Jacobian (dx/dXi)
        Matrix<double,2,2> Jmat, invJ;
        Jmat << x_(1)-x_(0), x_(2)-x_(0), y_(1)-y_(0), y_(2)-y_(0);
        double Area = Jmat.determinant()/2;
        invJ = Jmat.inverse();
        
        //  areal coordinates
        // MatrixXd L = MatrixXd(1,3); //(1,3);
        RowVector3d L;
        Matrix<double,1,6> Lsq;
        Matrix<double,2,6>d_Lsq_dxi;
        // MatrixXd Lsq(1,6), d_Lsq_dxi(2,6);
        L << r,s,1-r-s;
        // shape function: quadratic monomials of areal coordinates
        Lsq << r*r, s*s, (1-r-s)*(1-r-s), s*(1-r-s), (1-r-s)*r, r*s;
        // derivatives (wrt xi,eta) of quadratic monomials of areal coord.
        d_Lsq_dxi << -2*r, 2*s, 0,      (1-r-s), -(1-r-s),     r-s, 
                    -2*r, 0,      2*(1-r-s), s, r-(1-r-s), -s;
        
        //  DKT shape function -----------
        MatrixXd N_thetax, N_thetay, d_N_thetax_dxi, d_N_thetax_dx, d_N_thetay_dxi, d_N_thetay_dx;
        N_thetax = Lsq*dkt_opt.M_thetax;
        N_thetay = Lsq*dkt_opt.M_thetay;
        // shape functions of derivative of thetax wrt xi,eta
        d_N_thetax_dxi = d_Lsq_dxi*dkt_opt.M_thetax;
        // shape functions of derivative of thetax wrt x,y
        d_N_thetax_dx = invJ.transpose()*d_N_thetax_dxi;
        
        // shape functions of derivative of thetay wrt xi,eta
        d_N_thetay_dxi = d_Lsq_dxi*dkt_opt.M_thetay;
        // shape functions of derivative of thetay wrt x,y
        d_N_thetay_dx = invJ.transpose()*d_N_thetay_dxi;
        
        // shape functions of curvatures kappax, kappay, kappaxy
        MatrixXd N_kappa(3,9);
        N_kappa.row(0) = d_N_thetay_dx.row(0); 
        N_kappa.row(1) =-d_N_thetax_dx.row(1);
        N_kappa.row(2) = d_N_thetay_dx.row(1) - d_N_thetax_dx.row(0);
        
        // OPT shape function --------------
        /// shape functions of in-plane displacements and drilling rotation
        MatrixXd N_psi, N_eta = MatrixXd::Zero(3,9);
        N_psi = L*dkt_opt.M_psi;
        
        // shape functions of in-plane deformations etax, etay, gamxy
        N_eta << L*dkt_opt.M_etax,  L*dkt_opt.M_etay,  L*dkt_opt.M_gamxy;
        
        // Outputs
        Core_E.Km.noalias() +=  (N_eta.transpose()*Amat*N_eta)*w*Area;
        Core_E.Kb.noalias() +=  (N_kappa.transpose()*Dmat*N_kappa)*w*Area;
        Core_E.Kmb.noalias() += (N_eta.transpose()*Bmat*N_kappa)*w*Area;
        // Core_E.Km = Core_E.Km + (N_psi.transpose()*Gmat(9)*N_psi)*w*Area; // if phiflag
        
        Core_E.Fm.noalias() += (N_eta.transpose()*Amat*N_eta*u_m)*w*Area;
        Core_E.Fb.noalias() += (N_kappa.transpose()*Dmat*N_kappa*u_b)*w*Area;
        
        Core_E.Fmf.noalias() += (N_eta.transpose()*Nph)*w*Area;
        Core_E.Fbf.noalias() += (N_kappa.transpose()*Mph)*w*Area;
    }
    return Core_E;
}
