#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <cmath> 
#include <armadillo>
	
#include <Eigen/SparseCore>


Eigen::SparseMatrix<double>  interF(
Eigen::MatrixXd V, //vertices
Eigen::MatrixXi F, //Faces
double alpha,
double beta,
double gamma);

