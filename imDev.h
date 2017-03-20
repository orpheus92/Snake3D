#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <cmath> 
#include <armadillo>

arma::cube imDev(
arma::cube Img,
int sigma,
int type
// type 1 = x; 2 = y; 3 = z; 4 = xx; 5 = yy; 6 = zz; 7 = xy; 8 = xz; 9 = yz
);
