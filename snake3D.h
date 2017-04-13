#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <armadillo>	
#include <Eigen/SparseCore>

#include <unistd.h>
#include <string>
#include <cmath>
#include <iostream>
#include "interpcont.h"
#include "extForce.h"
#include "imDev.h"
#include "interF.h"
#include "GVFimF.h"

//#include <GLFW/glfw3.h>
//#include <igl/get_seconds.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/triangle/triangulate.h>
#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>
#include <igl/unique.h>

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snake3D(
Eigen::MatrixXd V, //vertices
Eigen::MatrixXi F, //initial mesh
arma::cube input, //input 3D image
double gamma, //time step, default 1
int iter, //# of iteration, default 100
//int npts, //# of pts to interpolate contours 
int sigma, // sigma to calculate img derivative default 10
double wl, //attraction to lines, < 0 to black line; > 0 to white line; default 0.04
double we, //attraction to edge, default 2.0
//double wt, //attraction to end points, default 0.01
double sigma2, //sigma to calculate gradient of edge energy image (give image force), default 20
double alpha, //membrane energy, default 0.2
double beta, //thin plate energy, default 0.2
double delta, //baloon force, default 0.1
double kappa, //weight of external img force, default 2
// the following is used for GVF snake
double mu, //tradeoff between real edge vectors and noise vectors, default 0.2
int Giter, //GVF iteration, default 0
double sigma3, //sigma used to calculate laplacian in GVF, default 1
double lamb, //Weight which changes the direction of the image potential force to the direction of the surface normal,
// (Keeps the surface from self intersecting)
igl::viewer::Viewer& viewer
);

