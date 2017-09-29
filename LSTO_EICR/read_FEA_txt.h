#ifndef READ_FEA_TXT_H
#define READ_FEA_TXT_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <tuple>
#include "../../eigen3/Eigen/Dense" 

using namespace Eigen;

struct FEACONFIG{
    MatrixXd NODE, Force_NM, Force_fix;
    VectorXi BCid;
    MatrixXi ELEM;
};

struct FEAFIELD{
    MatrixXd GU_u;
    VectorXd p_Adjoint, GU_Rv;
    double lambdaR;
};


class myreadTXT{
    public:
    myreadTXT(const char* fconfig, const char* fresults); // config file: nNODE, nELEM, nBCid
    void read_config();
    void read_results();

    FEACONFIG feaConfig;
    FEAFIELD  feaField;

    private:
        
    double e[6];
    int i[4];

    const char* fconfig;
    const char* fresults;

    uint nNODE, nELEM, nBCid;
};

#endif