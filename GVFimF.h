#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <cmath> 

#include "imDev.h"

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  GVFimF(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
double mu,
int Giter,
double sigma3,
int style
);
