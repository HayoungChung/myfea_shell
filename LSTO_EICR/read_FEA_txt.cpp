#include "read_FEA_txt.h"

using namespace Eigen;

myreadTXT::myreadTXT(const char* fconfig_, const char* fresults_){
    fconfig = fconfig_;
    fresults = fresults_;

    std::ifstream file;
    file.open(fconfig);

    std::string line_;
    int cnt = 0;

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
    
    std::cout << "nNODE: " << nNODE << ", nELEM: " << nELEM << ", nBCid: "  << nBCid << std::endl;

    feaConfig.NODE.resize(nNODE,3);
    feaConfig.ELEM.resize(nELEM,4);
    feaConfig.BCid.resize(nBCid);
    feaConfig.Force_NM.resize(nELEM,6);
    feaConfig.Force_fix.resize(nNODE,6);

    feaField.GU_u.resize(nNODE,3);
    feaField.GU_Rv.resize(nNODE*9);
    feaField.p_Adjoint.resize(nNODE*6);
}

void myreadTXT::read_config(){
    std::ifstream file;
    file.open(fconfig);

    std::string line_;
    int cnt = 0;
     if (file.is_open()){
        while (!file.eof()){
            file >> line_;
            // std::cout << line_ << std::endl;
            if (line_.compare(0,5,"NODE:") == 0 ){
                while (cnt < nNODE){
                    for (uint mm = 0; mm < 3; ++mm){
                        file >> e[mm];
                        feaConfig.NODE(cnt,mm) = e[mm];
                    }
                    cnt += 1;
                }    
            cnt = 0;           
            // std::cout << NODE << std::endl;
            }
            if (line_.compare(0,5,"ELEM:") == 0 ){
                while (cnt < nELEM){
                    for (uint mm = 0; mm < 4; ++mm){
                        file >> i[mm];
                        feaConfig.ELEM(cnt,mm) = i[mm];
                    }
                    cnt += 1;
                }    
            cnt = 0;            
            // std::cout << ELEM << std::endl;  
            }
            if (line_.compare(0,5,"BCid:") == 0 ){
                while (cnt < nBCid){
                    for (uint mm = 0; mm < 1; ++mm){
                        file >> i[0];
                        feaConfig.BCid(cnt) = i[mm];
                    }
                    cnt += 1;
                }    
            cnt = 0;            
            // std::cout << BCid << std::endl;  
            }
            if (line_.compare(0,9,"Force.NM:") == 0 ){
                while (cnt < nELEM){
                    // file >> e[0] >> e[1] >> e[2] >> e[3] >> e[4] >> e[5];
                    for (uint mm = 0; mm < 6; ++mm){
                        file >> e[mm];
                        feaConfig.Force_NM(cnt,mm) = e[mm];
                    }
                    cnt += 1;
                }   
            cnt = 0;            
            // std::cout << Force_NM << std::endl;  
            }
            if (line_.compare(0,10,"Force.fix:") == 0 ){
                while (cnt < nNODE){
                    // file >> e[0] >> e[1] >> e[2] >> e[3] >> e[4] >> e[5];
                    for (uint mm = 0; mm < 6; ++mm){
                        file >> e[mm];
                        feaConfig.Force_fix(cnt,mm) = e[mm];
                    }
                    cnt += 1;
                }   
            cnt = 0;            
            // std::cout << feaConfig.Force_fix << std::endl;  
            }
        }
        file.close();
    }   
}

void myreadTXT::read_results(){
    std::ifstream file;
    file.open(fresults);

    std::string line_;
    int cnt = 0;
    if (file.is_open()){
        while (!file.eof()){
            file >> line_;
            // std::cout << line_ << std::endl;
            if (line_.compare(0,5,"GU_u:") == 0 ){
                while (cnt < nNODE){
                    // file >> e[0] >> e[1] >> e[2];
                    for (uint mm = 0; mm < 3; ++mm){
                        file >> e[mm];
                        feaField.GU_u(cnt,mm) = e[mm];
                    }
                    cnt += 1;
                }    
            cnt = 0;           
            // std::cout << GU_u << std::endl;
            }
            if (line_.compare(0,6,"GU_Rv:") == 0 ){
                // std::cout << cnt << ", " << nELEM << std::endl;
                while (cnt < nNODE*9){
                    file >> e[0];
                    feaField.GU_Rv(cnt,0) = e[0];
                    cnt += 1;
                }   
            cnt = 0;            
            // std::cout << GU_Rv << std::endl;  
            }
            if (line_.compare(0,10,"p_Adjoint:") == 0 ){
                // std::cout << cnt << ", " << nELEM << std::endl;
                while (cnt < nNODE*6){
                    file >> e[0];
                    feaField.p_Adjoint(cnt) = e[0];
                    cnt += 1;
                }   
            cnt = 0;            
            // std::cout << p_Adjoint << std::endl;  
            }
            if (line_.compare(0,7,"lambdaR") == 0 ){
                file >> line_;
                file >> e[0]; 
                
                feaField.lambdaR = e[0];
            }

        }
        file.close();
    }   
}
   
