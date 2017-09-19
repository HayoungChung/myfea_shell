#include "./f_EICR_shell.h"
#define PI 3.14159265359

int main(){
    double Lxy[2] = {12,1};  
    int exy[2] = {2,2};   

    FEAMesh feaMesh(Lxy, exy);
    double h = 0.01;
    int nELEM = feaMesh.ELEM.rows(), nNODE = feaMesh.NODE.rows();

    double E = 1.2e6, v = 0.3;
    MatrixXd Cijkl(3,3), Amat(3,3), Dmat(3,3), Bmat(3,3);
    Cijkl << 1, v, 0, v, 1, 0, 0, 0, (1-v)*0.5;
    Cijkl *= E/(1-std::pow(v,2));
    Amat = Cijkl*h;
    Dmat = Cijkl*std::pow(h,3)/12.0;
    Bmat.setZero(3,3);

    std::vector<Material_ABD> material;
    Material_ABD material_0;
    material_0.Amat = Amat;
    material_0.Dmat = Dmat;
    material_0.Bmat = Bmat;
    for (unsigned int ee = 0; ee < nELEM; ++ee){
        material.push_back(material_0);
    }
    MatrixXd GU_u0  = MatrixXd::Zero(nNODE,3);
    MatrixXd eye3 = MatrixXd::Identity(3,3);
    MatrixXd GU_R0 = eye3.replicate(nNODE,1);
    VectorXd GU_Rv0 = fmat2store(GU_R0);

    std::cout << "GU_Rv0: \n\n" << GU_Rv0.transpose() << std::endl;

    double curv_ = -2*PI / Lxy[0];
    Vector3d eps0, kappa0;
    eps0 << 0, 0, 0;
    kappa0 << curv_, 0, 0;

    MatrixXd NM = MatrixXd::Zero(feaMesh.dpn,1); 
    NM << Amat*eps0, Dmat*kappa0;

    MatrixXd Force_NM(nELEM, feaMesh.dpn); 
    MatrixXd Force_Fix = MatrixXd::Zero(nNODE, feaMesh.dpn);
    Force_NM = NM.transpose().replicate(nELEM,1);
    Force_NM *= 0.1;

    std::cout << "Force_NM[1]: \n" << Force_NM.row(1) << "\n============\n" << std::endl;
    std::cout << "Amat[2]: \n" << material[2].Amat << "\n============\n" << std::endl;
    std::cout << "Dmat[10]: \n" << material[10].Dmat << "\n============\n" << std::endl;

    SparseMatrix<double> sGKT(nNODE*6,nNODE*6);

    VectorXd Res = VectorXd::Zero(nNODE*feaMesh.dpn);

    Force force;
    force.NM = Force_NM;
    force.fix = Force_Fix;
    MatrixXd GKT_m = f_EICR_shell(feaMesh, GU_u0, GU_Rv0, material, force, sGKT, Res);

    std::cout << "Res: \n" << Map<MatrixXd>(Res.data(),6,nNODE).transpose() << "\n\n" << std::endl;
    // std::cout << "GKT_m: \n" << GKT_m.transpose() << "\n\n" << std::endl;
    // std::cout << "sGKT: \n" << MatrixXd(sGKT) << "\n\n" << std::endl;

    // std::cout << "R_norm: \n" << Res.norm() << "\n\n" << std::endl;
    // std::cout << "GKT_m_norm: \n" << GKT_m.norm() << "\n\n" << std::endl;
    // std::cout << "sGKT_norm: \n" << sGKT.norm() << "\n\n" << std::endl;

    for (unsigned int ii =0; ii < 13; ++ii){
        for (unsigned int jj =0; jj < 13; ++jj){
            std::cout << "GKT( " << ii << ", " << jj <<
                    "):\n" << GKT_m.block(6*ii,6*jj,6,6) << "\n\n" << std::endl;
        }
    }

    // in matlab: getGKT_x = @(x,y) disp(full(GKT(x*6-5:x*6,y*6-5:y*6)))
    



    return 0;
}