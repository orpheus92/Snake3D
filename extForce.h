#define _USE_MATH_DEFINES

#include <Eigen/Core>
#include <string>
#include <cmath> 
#include <armadillo>
#include "imDev.h"

arma::cube extForce(
arma::cube Img,
double wl,
double we,
int sigma
);
