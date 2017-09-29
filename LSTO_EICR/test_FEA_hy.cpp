#include "./FEA_hy.h"

int main(){
    MatrixXd eye3 = MatrixXd::Identity(3,3);
        
    int npe = 3, dpn = 6, dpe = 18;

    // Mesh generation
    const double Lxy[2] = {12, 1};
    const int exy[2] = {2,2}; //{40, 2};
    const double h = 0.01; 
    
    class FEAMesh feaMesh(Lxy, exy);

    std::cout << "NODE: (" << feaMesh.NODE.rows() << ")" << std::endl;
    std::cout << feaMesh.NODE << std::endl;

    std::cout << "ELEM: (" << feaMesh.ELEM.rows() << ")" << std::endl;
    std::cout << feaMesh.ELEM << std::endl;

    return 0;

}