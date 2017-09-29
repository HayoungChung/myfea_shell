#include "./f_core_element.h"

int main(){
    double E = 1.2e6, v = 0.3; 
    Matrix<double,3,3> Cijkl, Amat, Dmat, Bmat;
    Cijkl << 1, v, 0, v, 1, 0, 0, 0, (1-v)/2;
    Cijkl *= E/(1-std::pow(v,2));
    Amat = Cijkl * 0.01;
    Dmat = Cijkl * std::pow(0.01,3)/12;
    Bmat.fill(0);

    Material_ABD Mater_e;
    Mater_e.Amat = Amat;
    Mater_e.Dmat = Dmat;
    Mater_e.Bmat = Bmat;

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
    FNM << 0.01, 0, 0, -0.0017, -0.0058, 0; 

    FilteredP p_d;
    p_d.membrane << bar_a.middleRows(0,3), bar_a.middleRows(6,3), bar_a.middleRows(12,3);
    p_d.bending  << bar_a.middleRows(3,3), bar_a.middleRows(9,3), bar_a.middleRows(15,3);
    
    CoreElement coreelem = f_core_element(coor,Mater_e,FNM,p_d);
    
    std::cout << "Km: \n\n" << coreelem.Km << "\n\n" << std::endl;
    std::cout << "Kb: \n\n" << coreelem.Kb << "\n\n" << std::endl;
    std::cout << "Kmb: \n\n" << coreelem.Kmb << "\n\n" << std::endl;

    std::cout << "Fm: \n\n" << coreelem.Fm.transpose() << "\n\n" << std::endl;
    std::cout << "Fmf: \n\n" << coreelem.Fmf.transpose() << "\n\n" << std::endl;
    std::cout << "Fb: \n\n" << coreelem.Fb.transpose() << "\n\n" << std::endl;
    std::cout << "Fbf: \n\n" << coreelem.Fbf.transpose() << "\n\n" << std::endl;

    return 0;    

}
