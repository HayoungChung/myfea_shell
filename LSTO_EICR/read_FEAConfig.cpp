#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <tuple>
#include "../../eigen3/Eigen/Dense"  
using namespace Eigen;

int main(int argc, char *argv[]){
    std::string line_;
    double a, b, c, d; // temp
    double e[6];
    int i, j, k, l;

    std::ifstream file;
	
	MatrixXd NODE;
    MatrixXi ELEM;
    VectorXd BCid;
    MatrixXd Force_NM, Force_fix;

    // first open to get NODE, ELEM numbs
    file.open("FEAConfiguration.txt");
    uint cnt = 0;
    uint nNODE, nELEM, nBCid;

    if (file.is_open()){
       while (!file.eof()){
           file >> line_;
        //    std::cout << line_ << std::endl;
           if (line_.compare(0,5,"NODE:") == 0 ){
               while (line_.compare(0,2,"==") != 0){
                    file >> line_;
                    cnt += 1;
               }             
               nNODE = cnt/3;
               cnt = 0;  
           }
           if (line_.compare(0,5,"ELEM:") == 0 ){
               while (line_.compare(0,2,"==") != 0){
                    file >> line_;
                    cnt += 1;
               }             
               nELEM = cnt/4;
               cnt = 0;  
           }
            if (line_.compare(0,5,"BCid:") == 0 ){
               while (line_.compare(0,2,"==") != 0){
                    file >> line_;
                    cnt += 1;
               }             
               nBCid = cnt-1;
               cnt = 0;  
           }
        }
        file.close();
    }
    std::cout << "nNODE: " << nNODE << ", nELEM: " << nELEM << ", nBCid: "  << nBCid << std::endl;

    NODE.resize(nNODE,3);
    ELEM.resize(nELEM,4);
    BCid.resize(nBCid);
    Force_NM.resize(nELEM,6);
    Force_fix.resize(nNODE,6);
    // reopen to get values   

    file.open("FEAConfiguration.txt");
    cnt = 0;
    if (file.is_open()){
        while (!file.eof()){
            file >> line_;
            // std::cout << line_ << std::endl;
            if (line_.compare(0,5,"NODE:") == 0 ){
                while (cnt < nNODE){
                    file >> a >> b >> c;
                    NODE(cnt,0) = a;
                    NODE(cnt,1) = b; 
                    NODE(cnt,2) = c;
                    cnt += 1;
                }    
            cnt = 0;           
            // std::cout << NODE << std::endl;
            }
            if (line_.compare(0,5,"ELEM:") == 0 ){
                // std::cout << cnt << ", " << nELEM << std::endl;
                while (cnt < nELEM){
                    file >> i >> j >> k >> l;
                    ELEM(cnt,0) = i;
                    ELEM(cnt,1) = j; 
                    ELEM(cnt,2) = k;
                    ELEM(cnt,3) = l;
                    cnt += 1;
                }   
            cnt = 0;            
            // std::cout << ELEM << std::endl;  
            }
            if (line_.compare(0,5,"BCid:") == 0 ){
                // std::cout << cnt << ", " << nELEM << std::endl;
                while (cnt < nBCid){
                    file >> i;
                    BCid(cnt) = i;
                    cnt += 1;
                }   
            cnt = 0;            
            // std::cout << BCid << std::endl;  
            }
            if (line_.compare(0,9,"Force.NM:") == 0 ){
                // std::cout << cnt << ", " << nELEM << std::endl;
                while (cnt < nELEM){
                    file >> e[0] >> e[1] >> e[2] >> e[3] >> e[4] >> e[5];
                    for (uint mm = 0; mm < 6; ++mm){
                        Force_NM(cnt,mm) = e[mm];
                    }
                    cnt += 1;
                }   
            cnt = 0;            
            // std::cout << Force_NM << std::endl;  
            }
            if (line_.compare(0,10,"Force.fix:") == 0 ){
                // std::cout << cnt << ", " << nELEM << std::endl;
                while (cnt < nNODE){
                    file >> e[0] >> e[1] >> e[2] >> e[3] >> e[4] >> e[5];
                    for (uint mm = 0; mm < 6; ++mm){
                        Force_fix(cnt,mm) = e[mm];
                    }
                    cnt += 1;
                }   
            cnt = 0;            
            // std::cout << Force_fix << std::endl;  
            }
        }
        file.close();
    }   
    return 0;
}
