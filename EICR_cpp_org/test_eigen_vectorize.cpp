#include <iostream>
#include "./../../../../eigen3/Eigen/Dense"

#include "./../../../../eigen3/Eigen/Core"

using namespace Eigen;
int main(){
    Matrix3d A;
    A << 1, 4, 3, 8, 2, 7, 1, 5, 2;
    double k = 3;

    Matrix3d isEqMat;
    isEqMat = (A.array() - k).square();
    std::cout << isEqMat.cast<int>().cwiseEqual(0) << std::endl;

    Matrix3d B;
    B.fill(2);

    // Matrix<int,9,1> idi; // doesn't work
    Matrix3i idi;
    // MatrixXi idi; //OK
    idi = (A.array() == B.array()).cast<int>() + 1;
    std::cout << idi.norm() << std::endl;

    Vector3i q;
    q << 2, 3, 7;
    std::cout << q.cwiseEqual(7) << std::endl;
    return 0;

    
}