#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <tuple>
#include "../../eigen3/Eigen/Dense"  
using namespace Eigen;

int main(int argc, char *argv[]){
    std::string line_;
    double e[6];
    int i[4];

    std::ifstream file;
	
	MatrixXd GU_u;
    VectorXd GU_Rv, p_Adjoint;
    // std::string tmp_ = "FEA_solution.txt";

    // const char* fname  = tmp_.c_str();
    const char* fname = argv[1];
    
    // first open to get NODE, ELEM numbs
    file.open(fname);
    uint cnt = 0;
    uint nNODE, nELEM, nBCid;

    if (file.is_open()){
       while (!file.eof()){
           file >> line_;
        //    std::cout << line_ << std::endl;
           if (line_.compare(0,5,"GU_u:") == 0 ){
               while (line_.compare(0,2,"==") != 0){
                    file >> line_;
                    cnt += 1;
               }             
               nNODE = cnt/3;
               cnt = 0;  
           }
        }
        file.close();
    }

    GU_u.resize(nNODE,3);
    GU_Rv.resize(nNODE*9);
    p_Adjoint.resize(nNODE*6);

    // reopen to get values   

    file.open("FEA_solution.txt");
    cnt = 0;
    if (file.is_open()){
        while (!file.eof()){
            file >> line_;
            // std::cout << line_ << std::endl;
            if (line_.compare(0,5,"GU_u:") == 0 ){
                while (cnt < nNODE){
                    file >> e[0] >> e[1] >> e[2];
                    GU_u(cnt,0) = e[0];
                    GU_u(cnt,1) = e[1]; 
                    GU_u(cnt,2) = e[2];
                    cnt += 1;
                }    
            cnt = 0;           
            std::cout << GU_u << std::endl;
            }
            if (line_.compare(0,6,"GU_Rv:") == 0 ){
                // std::cout << cnt << ", " << nELEM << std::endl;
                while (cnt < nNODE*9){
                    file >> e[0];
                    GU_Rv(cnt,0) = e[0];
                    cnt += 1;
                }   
            cnt = 0;            
            std::cout << GU_Rv << std::endl;  
            }
            if (line_.compare(0,10,"p_Adjoint:") == 0 ){
                // std::cout << cnt << ", " << nELEM << std::endl;
                while (cnt < nNODE*6){
                    file >> e[0];
                    p_Adjoint(cnt) = e[0];
                    cnt += 1;
                }   
            cnt = 0;            
            std::cout << p_Adjoint << std::endl;  
            }
        }
        file.close();
    }   
    return 0;
}
