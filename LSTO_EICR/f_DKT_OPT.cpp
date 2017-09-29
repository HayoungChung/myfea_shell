#include "./f_DKT_OPT.h"

DKT_OPT f_DKT_OPT(const MatrixXd & xycoord, double const nu){
    // Vector3d x, y;
    // x = xycoord.row(0);
    // y = xycoord.row(1);
    
    // MatrixXd M_thetax, M_thetay;// = MatrixXd::Zero(6,9), M_thetay = MatrixXd::Zero(6,9);
    // MatrixXd M_etax, M_etay, M_gamxy, M_psi; // = MatrixXd::Zero(3,9), M_etay = MatrixXd::Zero(3,9),
         // M_gamxy = MatrixXd::Zero(3,9), M_psi = MatrixXd::Zero(3,9);

    DKT_OPT dkt_opt;


    DKT_shape_functions(xycoord.row(0), xycoord.row(1), dkt_opt.M_thetax, dkt_opt.M_thetay);
    OPT_shape_function(xycoord, nu, dkt_opt.M_etax, dkt_opt.M_etay, dkt_opt.M_gamxy, dkt_opt.M_psi);

    // return value
    // dkt_opt.M_thetax = M_thetax;
    // dkt_opt.M_thetay = M_thetay;
    
    // dkt_opt.M_etax   = M_etax;
    // dkt_opt.M_etay   = M_etay;
    // dkt_opt.M_gamxy  = M_gamxy;
    // dkt_opt.M_psi    = M_psi;

    return dkt_opt;
}


void OPT_shape_function(const MatrixXd & xycoord, double const nu,
                MatrixXd& M_etax, MatrixXd& M_etay, MatrixXd& M_gamxy, MatrixXd& M_psi){
    //   Auxiliary matrices are required for OPT:
    //   beta parameters for OPT
    // double nu = std::min(Amat(1,0)/Amat(0,1),Amat(1,0)/Amat(0,0)), alpha_b = 3.0/2, beta0 = std::max((1.-4.*(std::pow(nu,2)))/2.,0.01);
    double alpha_b = 1.5, beta0 = std::max((1.-4.*(std::pow(nu,2)))/2.,0.01);

    std::vector<double> beta_pars;
    beta_pars = {alpha_b, beta0, 1.,2.,1.,0.,1.,-1.,-1.,-1.,-2.};
    
    // compute auxiliary matrices
    OPT_AUX opt_aux = f_OPT_aux(beta_pars,xycoord);

    // [etax;etay;gamxy] = M1*r+M2*s+M(3)*(1-r-s), eq. (49)
    // Mi(3x_9), Metai(3,9)
    // M_etai and M_psi are related to N(shape tensor) counterpart
    // by linear shape function
    MatrixXd M1(3,9), M2(3,9), M3(3,9);
    M1.noalias() = opt_aux.B + (1.5*std::sqrt(beta0)*opt_aux.Te*opt_aux.Q1*opt_aux.Ttu);
    M2.noalias() = opt_aux.B + (1.5*std::sqrt(beta0)*opt_aux.Te*opt_aux.Q2*opt_aux.Ttu);
    M3.noalias() = opt_aux.B + (1.5*std::sqrt(beta0)*opt_aux.Te*opt_aux.Q3*opt_aux.Ttu);

    // MatrixXd M_etax, M_etay, M_gamxy, M_psi;
    M_etax.row(0) = M1.row(0);
    M_etax.row(1) = M2.row(0);
    M_etax.row(2) = M3.row(0);

    M_etay.row(0) = M1.row(1);
    M_etay.row(1) = M2.row(1);
    M_etay.row(2) = M3.row(1);

    M_gamxy.row(0) = M1.row(2);
    M_gamxy.row(1) = M2.row(2);
    M_gamxy.row(2) = M3.row(2);

    // shape function of drilling rotation
    M_psi.setZero(3,9);
    M_psi(0,2) = -1; M_psi(1,5) = -1; M_psi(2,8) = -1;    
}

OPT_AUX f_OPT_aux(std::vector<double> fpars, const MatrixXd & xycoord){

    // parse the inputs
    MatrixXd x123 = xycoord.row(0), y123 = xycoord.row(1);
    double alp_b = fpars[0], beta0 = fpars[1];
    std::vector<double> beta;
    for (int ii = 0; ii < 9; ++ii){
        beta.push_back(fpars[ii+2]);
    }

    // compute inner variables (x, Area)
    std::vector<int> id1, id2;
    id1 = {0, 1, 2, 1, 2, 0};
    id2 = {1, 2, 0, 0, 1, 2};

    Matrix<double,3,3> x, y;
    for (int ii = 0; ii < 6; ++ii){ 
        x(id1[ii],id2[ii]) = x123(id1[ii])-x123(id2[ii]);
        y(id1[ii],id2[ii]) = y123(id1[ii])-y123(id2[ii]);
    }

    double A = (y(1,0)*x(0,2)-x(1,0)*y(0,2))/2., A2 = 2*A, A4 = 4*A;

    // 1. CST part B(3x9) matrix

    // Load-lumping matrix L (= B'*(h*Area))
    // L below is actually L*2/h
    MatrixXd L(9,3), B(3,9); 
    L << y(1,2),0,x(2,1), 
        0,x(2,1),y(1,2),
        y(1,2)*(y(0,2)-y(1,0)) , x(2,1)*(x(2,0)-x(0,1)), (x(2,0)*y(0,2)-x(0,1)*y(1,0))*2,
        y(2,0),0,x(0,2),
        0,x(0,2),y(2,0),
        y(2,0)*(y(1,0)-y(2,1)) , x(0,2)*(x(0,1)-x(1,2)), (x(0,1)*y(1,0)-x(1,2)*y(2,1))*2,
        y(0,1),0,x(1,0),
        0,x(1,0),y(0,1),
        y(0,1)*(y(2,1)-y(0,2)),x(1,0)*(x(1,2)-x(2,0)),(x(1,2)*y(2,1)-x(2,0)*y(0,2))*2;

    L.row(2) *= alp_b/6.; 
    L.row(5) *= alp_b/6.; 
    L.row(8) *= alp_b/6.; 
    B = L.transpose()/A2;

    // 2. Ttu(3x9) that extract the hierachical corner rotation from total one
    // done by substracting mean CST rotation
    MatrixXd Ttu(3,9);
    Ttu << x(2,1),y(2,1),A4,x(0,2),y(0,2),0 ,x(1,0),y(1,0),0,
           x(2,1),y(2,1),0 ,x(0,2),y(0,2),A4,x(1,0),y(1,0),0,
           x(2,1),y(2,1),0 ,x(0,2),y(0,2),0 ,x(1,0),y(1,0),A4;
    Ttu /= A4;

    // 3. Te(3,3) that transforms natural strain into global strain

    // squared side lengths
    double LL21, LL32, LL13;
    LL21 = std::pow(x(1,0),2)+std::pow(y(1,0),2);
    LL32 = std::pow(x(2,1),2)+std::pow(y(2,1),2); 
    LL13 = std::pow(x(0,2),2)+std::pow(y(0,2),2);

    MatrixXd Te(3,3);
    Te <<   y(1,2)*y(0,2)*LL21, y(2,0)*y(1,0)*LL32, y(0,1)*y(2,1)*LL13,
            x(1,2)*x(0,2)*LL21, x(2,0)*x(1,0)*LL32, x(0,1)*x(2,1)*LL13,
            (y(1,2)*x(2,0)+x(2,1)*y(0,2))*LL21, (y(2,0)*x(0,1)+x(0,2)*y(1,0))*LL32, (y(0,1)*x(1,2)+x(1,0)*y(2,1))*LL13;
    Te /= (A*A4);

    // 4. Qi(3x3), i=1,2,3: matrix relating the natural strains eps_i at corner i
    // Q1 + Q2 + Q3 = 0 (vanishment of high order K at the centeroid
    MatrixXd Q1(3,3), Q2(3,3), Q3(3,3);
    Q1 <<   beta[0]/LL21, beta[1]/LL21, beta[2]/LL21, 
            beta[3]/LL32, beta[4]/LL32, beta[5]/LL32,
            beta[6]/LL13, beta[7]/LL13, beta[8]/LL13;
    Q1 *= A2/3.;

    Q2 <<   beta[8]/LL21, beta[6]/LL21, beta[7]/LL21, 
            beta[2]/LL32, beta[0]/LL32, beta[1]/LL32,
            beta[5]/LL13, beta[3]/LL13, beta[4]/LL13;
    Q2 *= A2/3.;

    Q3 <<   beta[4]/LL21, beta[5]/LL21, beta[3]/LL21, 
            beta[7]/LL32, beta[8]/LL32, beta[6]/LL32,
            beta[1]/LL13, beta[2]/LL13, beta[0]/LL13;
    Q3 *= A2/3.;

    OPT_AUX opt_aux;
    opt_aux.B = B; opt_aux.Ttu = Ttu; opt_aux.Te = Te;
    opt_aux.Q1 = Q1; opt_aux.Q2 = Q2; opt_aux.Q3 = Q3;
    opt_aux.A = A;

    return opt_aux;
}


void DKT_shape_functions(const VectorXd x, const VectorXd y, MatrixXd& M_thetax, MatrixXd& M_thetay){
    // This function computes the shape functions of DKT plate
    //
    // input:
    //
    // x(2), y(2): nodal coordinates (the triangle is located in the x,y plane)
    //
    // output:
    //
    // M_thetax(6,9); M_thetay(6,9); M_w(15,9)
    //
    // the shape functions of bending rotations thetax and thetay
    // are quadratic homogeneous polynomials of areal coordinates L=[1-xi-eta, xi, eta];
    // N_thetax (shape function of thetax) is given by Lsq*M_thetax, where:
    // Lsq=[L(1)^2, L(2)^2, L(3)^2, L(2)*L(3), L(3)*L(1), L(1)*L(2)]
    // analogously, N_thetay=Lsq*M_thetay
    //
    // the shape function of transversal displacement w is a quartic homogeneous polynomial of areal coordinates
    // N_w (shape function of w) is given by Lsqsq*M_w, where:
    // Lsqsq = [ ...
    //     L(1) ^ 4, ...
    //     L(2) ^ 4, ...
    //     L(3) ^ 4, ...
    //     L(2) ^ 3 * L(3), ...
    //     L(2) * L(3) ^ 3, ...
    //     L(3) ^ 3 * L(1), ...
    //     L(3) * L(1) ^ 3, ...
    //     L(1) ^ 3 * L(2), ...
    //     L(1) * L(2) ^ 3, ...
    //     L(2) ^ 2 * L(3) ^ 2, ...
    //     L(3) ^ 2 * L(1) ^ 2, ...
    //     L(1) ^ 2 * L(2) ^ 2, ...
    //     L(1) ^ 2 * L(2) * L(3), ...
    //     L(2) ^ 2 * L(3) * L(1), ...
    //     L(3) ^ 2 * L(1) * L(2)];
    //
    // the 9 dofs are ordered as follows:
    //   w, thetax, thetay at node 1, then at node 2, then at node 3
    //
    // This code is part of a Matlab toolkit distributed as supplementary material of the paper:
    // Caselli F, Bisegna P. Polar decomposition based corotational framework
    // for triangular shell elements with distributed loads.
    // International Journal for Numerical Methods in Engineering, 2013
    // DOI: 10.1002/nme.4528
    //
    // Authors' e-mail addresses:
    // caselli@ing.uniroma2.it (Federica Caselli)
    // bisegna@uniroma2.it (Paolo Bisegna)
    //
    // (C) 2010-2013 Paolo Bisegna and Federica Caselli. License: GNU General Public License (GPLv3)

    // M_thetax, M_thetay: used to compute the shape functions of thetax and thetay

    // M_thetax = [0 1 0 0 0 0 0 0 0; ...
    //     0 0 0 0 1 0 0 0 0; ...
    //     0 0 0 0 0 0 0 1 0; ...
    //     0 0 0 6 * (-y(2) + y(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) (x(1) ^ 2 - 2 * x(2) * x(1) + x(2) ^ 2 - 2 * y(2) ^ 2 - 2 * y(1) ^ 2 + 4 * y(2) * y(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) 3 * (-y(2) + y(1)) * (-x(2) + x(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) -6 * (-y(2) + y(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) (x(1) ^ 2 - 2 * x(2) * x(1) + x(2) ^ 2 - 2 * y(2) ^ 2 - 2 * y(1) ^ 2 + 4 * y(2) * y(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) 3 * (-y(2) + y(1)) * (-x(2) + x(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2); ...
    //     6 * (y(0) - y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) (x(0) ^ 2 - 2 * x(0) * x(2) - 2 * y(0) ^ 2 + x(2) ^ 2 - 2 * y(2) ^ 2 + 4 * y(0) * y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) 3 * (x(0) - x(2)) * (y(0) - y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) 0 0 0 -6 * (y(0) - y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) (x(0) ^ 2 - 2 * x(0) * x(2) - 2 * y(0) ^ 2 + x(2) ^ 2 - 2 * y(2) ^ 2 + 4 * y(0) * y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) 3 * (x(0) - x(2)) * (y(0) - y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2); ...
    //     6 * (-y(1) + y(0)) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) (x(0) ^ 2 - 2 * x(1) * x(0) + 4 * y(1) * y(0) + x(1) ^ 2 - 2 * y(1) ^ 2 - 2 * y(0) ^ 2) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) 3 * (-y(1) + y(0)) * (-x(1) + x(0)) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) -6 * (-y(1) + y(0)) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) (x(0) ^ 2 - 2 * x(1) * x(0) + 4 * y(1) * y(0) + x(1) ^ 2 - 2 * y(1) ^ 2 - 2 * y(0) ^ 2) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) 3 * (-y(1) + y(0)) * (-x(1) + x(0)) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) 0 0 0;];

    // M_thetay = [0 0 1 0 0 0 0 0 0; ...
    //     0 0 0 0 0 1 0 0 0; ...
    //     0 0 0 0 0 0 0 0 1; ...
    //     0 0 0 -6 * (-x(2) + x(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) 3 * (-y(2) + y(1)) * (-x(2) + x(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) -(2 * x(1) ^ 2 - 4 * x(2) * x(1) + 2 * x(2) ^ 2 - y(2) ^ 2 - y(1) ^ 2 + 2 * y(2) * y(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) 6 * (-x(2) + x(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) 3 * (-y(2) + y(1)) * (-x(2) + x(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2) -(2 * x(1) ^ 2 - 4 * x(2) * x(1) + 2 * x(2) ^ 2 - y(2) ^ 2 - y(1) ^ 2 + 2 * y(2) * y(1)) / (x(2) ^ 2 - 2 * x(2) * x(1) + x(1) ^ 2 + y(2) ^ 2 - 2 * y(2) * y(1) + y(1) ^ 2); ...
    //     -6 * (x(0) - x(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) 3 * (x(0) - x(2)) * (y(0) - y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) -(2 * x(0) ^ 2 - 4 * x(0) * x(2) - y(0) ^ 2 + 2 * x(2) ^ 2 - y(2) ^ 2 + 2 * y(0) * y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) 0 0 0 6 * (x(0) - x(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) 3 * (x(0) - x(2)) * (y(0) - y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2) -(2 * x(0) ^ 2 - 4 * x(0) * x(2) - y(0) ^ 2 + 2 * x(2) ^ 2 - y(2) ^ 2 + 2 * y(0) * y(2)) / (x(0) ^ 2 - 2 * x(0) * x(2) + x(2) ^ 2 + y(0) ^ 2 - 2 * y(0) * y(2) + y(2) ^ 2); ...
    //     -6 * (-x(1) + x(0)) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) 3 * (-y(1) + y(0)) * (-x(1) + x(0)) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) -(2 * x(0) ^ 2 - 4 * x(1) * x(0) + 2 * y(1) * y(0) + 2 * x(1) ^ 2 - y(1) ^ 2 - y(0) ^ 2) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) 6 * (-x(1) + x(0)) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) 3 * (-y(1) + y(0)) * (-x(1) + x(0)) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) -(2 * x(0) ^ 2 - 4 * x(1) * x(0) + 2 * y(1) * y(0) + 2 * x(1) ^ 2 - y(1) ^ 2 - y(0) ^ 2) / (x(1) ^ 2 - 2 * x(1) * x(0) + x(0) ^ 2 + y(1) ^ 2 - 2 * y(1) * y(0) + y(0) ^ 2) 0 0 0;];

    // optimized expressions

    VectorXd t(76);
            
    t[1 ]= (-y(2) + y(1));
    t[2 ]= (x(2)*x(2));
    t[3 ]= (x(2) * x(1));
    t[4 ]= 2 * t[3 ];
    t[5 ]= (x(1)*x(1));
    t[6 ]= (y(2)*y(2));
    t[7 ]= (y(2) * y(1));
    t[8 ]= 2 * t[7 ];
    t[9 ]= (y(1)*y(1));
    t[11] = 1 / (t[2 ]- t[4 ]+ t[5 ]+ t[6 ]- t[8 ]+ t[9]);
    t[13] = 6 * t[11] * t[1 ];
    t[15] = 2 * t[9];
    t[16] = 2 * t[6];
    t[18] = t[11] * (t[5 ]- t[4 ]+ 4 * t[7 ]+ t[2 ]- t[15] - t[16]);
    t[19] = (-x(2) + x(1));
    t[22] = 3 * t[11] * t[19] * t[1 ]; 
    t[23] = (y(0) - y(2));
    t[24] = (x(0)*x(0));
    t[25] = (x(0) * x(2));
    t[26] = 2 * t[25];
    t[27] = (y(0) * y(0));
    t[28] = (y(0) * y(2));
    t[29] = 2 * t[28];
    t[31] = 1 / (t[24] - t[26] + t[2 ]+ t[27] - t[29] + t[6]);
    t[33] = 6 * t[31] * t[23];
    t[35] = 2 * t[27];
    t[37] = t[31] * (t[24] - t[26] - t[16] + 4 * t[28] - t[35] + t[2]);
    t[38] = (x(0) - x(2));
    t[41] = 3 * t[31] * t[23] * t[38];
    t[42] = (-y(1) + y(0));
    t[43] = (x(1) * x(0));
    t[44] = 2 * t[43];
    t[45] = (y(1) * y(0));
    t[46] = 2 * t[45];
    t[48] = 1 / (t[5 ]- t[44] + t[24] + t[9 ]- t[46] + t[27]);
    t[50] = 6 * t[48] * t[42];
    t[53] = t[48] * (t[24] - t[44] + t[5 ]- t[35] + 4 * t[45] - t[15]);
    t[54] = (-x(1) + x(0));
    t[57] = 3 * t[48] * t[54] * t[42];
    t[59] = 6 * t[11] * t[19];
    t[60] = 2 * t[5];
    t[62] = 2 * t[2];
    t[64] = t[11] * (t[60] - 4 * t[3 ]+ t[8 ]+ t[62] - t[9 ]- t[6]);
    t[66] = 6 * t[31] * t[38];
    t[67] = 2 * t[24];
    t[70] = t[31] * (t[67] - 4 * t[25] - t[6 ]+ t[29] - t[27] + t[62]);
    t[72] = 6 * t[48] * t[54];
    t[75] = t[48] * (t[67] - 4 * t[43] + t[60] - t[27] + t[46] - t[9]);

    M_thetax(0,1) = 1;
    M_thetax(1,4) = 1;
    M_thetax(2,7) = 1;
    M_thetax(3,3) = t[13];
    M_thetax(3,4) = t[18];
    M_thetax(3,5) = t[22];
    M_thetax(3,6) = -t[13];
    M_thetax(3,7) = t[18];
    M_thetax(3,8) = t[22];
    M_thetax(4,0) = t[33];
    M_thetax(4,1) = t[37];
    M_thetax(4,2) = t[41];
    M_thetax(4,6) = -t[33];
    M_thetax(4,7) = t[37];
    M_thetax(4,8) = t[41];
    M_thetax(5,0) = t[50];
    M_thetax(5,1) = t[53];
    M_thetax(5,2) = t[57];
    M_thetax(5,3) = -t[50];
    M_thetax(5,4) = t[53];
    M_thetax(5,5) = t[57];

    M_thetay(0,2) = 1;
    M_thetay(1,5) = 1;
    M_thetay(2,8) = 1;
    M_thetay(3,3) = -t[59];
    M_thetay(3,4) = t[22];
    M_thetay(3,5) = -t[64];
    M_thetay(3,6) = t[59];
    M_thetay(3,7) = t[22];
    M_thetay(3,8) = -t[64];
    M_thetay(4,0) = -t[66];
    M_thetay(4,1) = t[41];
    M_thetay(4,2) = -t[70];
    M_thetay(4,6) = t[66];
    M_thetay(4,7) = t[41];
    M_thetay(4,8) = -t[70];
    M_thetay(5,0) = -t[72];
    M_thetay(5,1) = t[57];
    M_thetay(5,2) = -t[75];
    M_thetay(5,3) = t[72];
    M_thetay(5,4) = t[57];
    M_thetay(5,5) = -t[75];
}