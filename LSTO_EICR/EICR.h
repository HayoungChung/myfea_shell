#ifndef EICR_H
#define EICR_H

#include <iostream>
#include <iomanip>
#include <fstream>
// #include <stdexcept>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <cassert>
#include <set>

#include "../../eigen3/Eigen/Dense"  
#include "../../eigen3/Eigen/Sparse"
#include "../../eigen3/unsupported/Eigen/SparseExtra"

// namespace EICR_FEA {
	#include "FEA_hy.h"
	#include "read_FEA_txt.h"
	#include "f_DKT_OPT.h"
	#include "CoreElement.h"
	#include "EICR_shell.h"
	#include "f_nlgeom.h"
	#include "Sensitivity.h"
// }



#endif 
