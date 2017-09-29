#include "./read_FEA_txt.h"

using namespace Eigen;

int main(){
    myreadTXT read("./fea_80_40_0.1/FEAConfiguration.txt", "./fea_80_40_0.1/FEA_solution.txt");
    read.read_config();
    std::cout << read.feaConfig.NODE << std::endl;
    read.read_results();
    std::cout << read.feaField.GU_u << std::endl;

    return 0;
}