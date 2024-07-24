#include "BC_Model.hpp"
#include "CV_Model.hpp"
#include "RS_Model.hpp"

#include <fstream>

int main() 
{
    double bc_et = BC_Model();
    double cv_et = CV_Model();
    double rs_et = RS_Model();

    // save the execution times to time_performance.csv
    std::ofstream file("dde_elapsed_times_cpp.csv");
    file << "model,execution_time\n";
    file << "BC," << bc_et << "\n";
    file << "CV," << cv_et << "\n";
    file << "RS," << rs_et << "\n";
    file.close();


    return 0;
}