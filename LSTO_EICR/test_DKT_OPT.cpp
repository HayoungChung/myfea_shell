// #include <FEA_hy.h>
// rephrase of matlab toolkit (see DKT_OPT_SHELL)
#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "../../../../eigen3/Eigen/Dense"
#include "../../../../eigen3/Eigen/Sparse"
#include "EICR.h"
#include "FEA_hy.h"
#include "f_nlgeom.h"

#define PI 3.141592

int main(){

    // Nodal coordinate: core_el.col(i) contains the coord. of node V
    MatrixXd coor(2,3);
    coor <<  -5.7871e+00,   7.3542e+00,  -1.5671e+00,
                -3.1558e+00,  -3.1558e+00,   6.3116e+00;

    VectorXd bar_a(18);
    bar_a <<  
    // displacement of node 1 
    -2.9059e-01, -1.6536e+00, 0,
    // rotation vector of node 1
    -6.6781e-02, -4.4860e-03, 6.2155e-02,
    // displacement of node 2
    -1.4838e+00, 1.7447e+00, -8.8818e-16,
    // rotation vector of node 2
    -4.1467e-02, 4.8584e-02, 5.7000e-02,
    // displacement of node 3
    1.7744e+00, -9.1088e-02, -4.4409e-16,
    // rotation vector of node 3
    -1.0751e-01, 8.4418e-03, 1.1402e-01;
    
    VectorXd FNM(6);
    FNM << 0, 0, 0, -0.0017, -0.0058, 0;

    /*
    M_thetax =

            0   1.0000e+00            0            0            0            0            0            0            0
            0            0            0            0   1.0000e+00            0            0            0            0
            0            0            0            0            0            0            0   1.0000e+00            0
            0            0            0  -3.3568e-01  -5.8901e-01  -1.4974e+00   3.3568e-01  -5.8901e-01  -1.4974e+00
  -5.2871e-01  -1.5028e+00   1.1156e+00            0            0            0   5.2871e-01  -1.5028e+00   1.1156e+00
            0   1.0000e+00            0            0   1.0000e+00            0            0            0            0

    M_thetay =

            0            0   1.0000e+00            0            0            0            0            0            0
            0            0            0            0            0   1.0000e+00            0            0            0
            0            0            0            0            0            0            0            0   1.0000e+00
            0            0            0  -3.1632e-01  -1.4974e+00  -4.1099e-01   3.1632e-01  -1.4974e+00  -4.1099e-01
   2.3566e-01   1.1156e+00   5.0275e-01            0            0            0  -2.3566e-01   1.1156e+00   5.0275e-01
   4.5657e-01            0  -2.0000e+00  -4.5657e-01            0  -2.0000e+00            0            0            0

    M_etax =

  -1.0532e-01   3.1012e-02   3.8388e-01   6.2273e-02  -3.1012e-02   2.2743e-01   4.3046e-02            0   2.0377e-01
  -4.6873e-02  -3.1012e-02  -2.2743e-01   8.9919e-02   3.1012e-02  -3.8388e-01  -4.3046e-02            0  -2.0377e-01
  -7.6096e-02  -1.5702e-19   3.8388e-01   7.6096e-02   1.5702e-19  -3.8388e-01  -7.1560e-19            0            0


    M_etay =

   5.0498e-02  -1.2530e-01  -6.6335e-01   2.3887e-02   1.9671e-02  -1.2455e-01  -7.4385e-02   1.0563e-01  -6.2057e-01
  -3.0355e-02  -3.9494e-02  -2.8051e-01  -1.4359e-02  -6.6132e-02   3.9875e-01   4.4714e-02   1.0563e-01   7.2841e-01
  -2.0143e-02  -5.0331e-02   1.0166e-02  -9.5281e-03  -5.5295e-02   2.8706e-01   2.9671e-02   1.0563e-01   2.6460e-01

    M_gamxy =

  -4.5655e-02  -1.0374e-01  -3.4222e-01  -2.1596e-02   1.0374e-01  -7.0275e-01   6.7251e-02            0   3.1835e-01
  -1.6631e-02  -1.3454e-01  -9.2862e-01  -7.8670e-03   1.3454e-01  -7.2347e-01   2.4499e-02            0   1.1597e-01
  -1.5283e-01   9.9968e-03   7.8916e-01  -7.2293e-02  -9.9968e-03   4.0790e-01   2.2513e-01            0   1.0657e+00

    */
    DKT_OPT dkt_opt = f_DKT_OPT(coor,0.3);

    std::cout << "M_thetax: " << "\n\n" << dkt_opt.M_thetax << "\n\n" << std::endl;
    std::cout << "M_thetay: " << "\n\n" << dkt_opt.M_thetay << "\n\n" << std::endl;
    std::cout << "M_etax: " << "\n\n" << dkt_opt.M_etax << "\n\n" << std::endl; 
    std::cout << "M_etay: " << "\n\n" << dkt_opt.M_etay << "\n\n" << std::endl;
    std::cout << "M_gamxy: " << "\n\n" << dkt_opt.M_gamxy << "\n\n" << std::endl;

    return 0;
}