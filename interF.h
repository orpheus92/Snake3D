#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <cmath> 
#include <armadillo>
#include <Eigen/SparseLU>
#include <Eigen/SparseCore>
#include <igl/viewer/Viewer.h>
//#include <igl/unique.h>
//#include <igl/triangle/unique.h>
#include "tutorial_shared_path.h"

Eigen::SparseMatrix<double>  interF(
Eigen::MatrixXd V, //vertices
Eigen::MatrixXi F, //Faces
double alpha,
double beta,
double gamma);

