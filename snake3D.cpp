#include "snake3D.h"
#include <stb_image.h>
#include <iostream>
#include "boost/multi_array.hpp"
#include <cassert>



//snake movement
/*Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snakeMove(
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> intForce,
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P,
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
	double gamma,
	double kappa,
	double delta);
//baloon Force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baloonF(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P);
//bilinear interp
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interp2(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mypts,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P);

*/

//main function for snake
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snake3D(
Eigen::MatrixXd V, //vertices
Eigen::MatrixXi F, //Faces
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
)
{
// make clockwise contour  (always clockwise due to balloon force)
double volume = 0;
Eigen::Vector3d a;
Eigen::Vector3d b;
Eigen::Vector3d c;
Eigen::Vector3d k;
double v;
Eigen::MatrixXi Fin;

for(int j = 0; j<F.rows();j++){
	a = V.row(F(j,0)-1);
	b = V.row(F(j,1)-1);
	c = V.row(F(j,2)-1);
	k = b.cross(c);
	v = (a(0)*k(0)+a(1)*k(1)+a(2)*k(2))/6; 
	volume = volume + v; 
}
if (volume >0)
{
	Fin = F.rowwise().reverse();
}

//Calculate external force
//size = image size
arma::cube Eext;

Eext = extForce(input,wl,we,sigma);

//Make the external flow field
//size = image size
arma::cube Fx;
arma::cube Fy;
arma::cube Fz;

Fx = -imDev(Eext,sigma2,1)*2*sigma2*sigma2;
std::cout<<"after imDev 1"<<std::endl;
Fy = -imDev(Eext,sigma2,2)*2*sigma2*sigma2;
std::cout<<"after imDev 2"<<std::endl;
Fz = -imDev(Eext,sigma2,3)*2*sigma2*sigma2;
std::cout<<"after imDev 3"<<std::endl;
//Calcuate GVF Image Force  Might be needed later 
//Not used for now
/*
arma::cube Fx2;
arma::cube Fy2;
arma::cube Fz2;

//will be modified
Fx2=GVFimF(Fx,Fy,Fz, mu, Giter, sigma3, 1);
Fy2=GVFimF(Fx,Fy,Fz, mu, Giter, sigma3, 2);
Fz2=GVFimF(Fx,Fy,Fz, mu, Giter, sigma3, 3);
*/
//std::cout<<"after imDev"<<std::endl;
Eigen::SparseMatrix<double> intForce;


//internal Force for snake
//This is correct
intForce = interF(V,F,alpha,beta,gamma);


// Triangulated interior

//for (int it = 0;it<iter;it++){

//Vout = snakeMove(intForce,F,V,Eext,gamma,kappa,delta,lamb);

  //  }
  std::cout<<"end of snake"<<std::endl;
return V; 
}

//function to calculate baloon force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baloonF(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P){
//input n by 2
//output n by 2
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out; 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dx; 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dy;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> l;
dx.resize(P.rows(),1);
dy.resize(P.rows(),1);
l.resize(P.rows(),1);

Eigen::VectorXd f;
Eigen::VectorXd b;

f.setLinSpaced(P.rows(),1,P.rows()); 
b.setLinSpaced(P.rows(),1,P.rows()); 
f.array() += 4;                
b.array() -= 4; 
for(int i = 0; i<P.rows();i++){

	if (f(i)>P.rows())
	f(i) = f(i)-P.rows();
	if (b(i)<1)
	b(i) = b(i)+P.rows();

	dx(i,0) = P(f(i)-1, 0)-P(b(i)-1, 0);
	dy(i,0) = P(f(i)-1, 1)-P(b(i)-1, 1);
	
}
l = (dx.array().square()+dy.array().square()).cwiseSqrt() ;            

out.resize(P.rows(),2);
out.col(0) = -dy.array() / l.array();  
out.col(1) = dx.array() / l.array();  
return out;
}


// function for linear interpolation
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interp2(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mypic,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P)
{
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out;
out.resize(P.rows(),1);
double x;
double y;
int indx;
int indy;
double q11;
double q12;
double q21;
double q22;

double x1;
double x2;
double y1;
double y2;

//bilinear interpolation
for (int i = 0; i<P.rows();i++){

y = P(i,0);
x = P(i,1);

indx = (int)ceil(x);
indy = (int)ceil(y);

x1 = (double)(indx-1);
x2 = (double)(indx);
y1 = (double)(indy-1);
y2 = (double)(indy);

	if ((indx==0) & (indy==0))	
		out(i,0) = mypic(0,0);
	else if (indx == 0)
		out(i,0) = mypic(indy,0)*(y-y1) + mypic(indy-1,0)*(y2-y);
	else if (indy == 0)
		out(i,0) = mypic(0,indx)*(x-x1) + mypic(0,indx-1)*(x2-x);
	else{
		q22 = mypic(indy,indx);
		q12 = mypic(indy-1,indx);
		q21 = mypic(indy,indx-1);
		q11 = mypic(indy-1,indx-1);
		out(i,0) = q11*(x2-x)*(y2-y) + q21*(x-x1)*(y2-y) + q12*(x2-x)*(y-y1) + q22*(x-x1)*(y-y1);
	}

}

return out;
}

/*
//This function will calculate one iteration of contour Snake movement
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snakeMove(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> intForce, //internal force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P,//contour pts
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,//external vector field
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,//external vec field
double gamma, //time step
double kappa,//external field weight
double delta) //Balloon Force Weight)
{
//Clamp contour to boundary
//rows and cols might be checked later 
//image force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> imforce;
imforce.resize(P.rows(),2);

//baloon force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baforce;

//clamp boundary
P.col(0)=P.col(0).cwiseMax(0).cwiseMin(Fx.cols()-1);
P.col(1)=P.col(1).cwiseMax(0).cwiseMin(Fx.rows()-1);

//Get image force on the contour points
imforce.col(0)=kappa*interp2(Fx,P);
imforce.col(1)=kappa*interp2(Fy,P);

//Get baloon force on the contour points

baforce = delta * baloonF(P);

//Update contour positions

P.col(0) = intForce * (gamma*P.col(0) + imforce.col(0) + baforce.col(0));
P.col(1) = intForce * (gamma*P.col(1) + imforce.col(1) + baforce.col(1));

//Clamp contour to boundary
P.col(0)=P.col(0).cwiseMax(0).cwiseMin(Fx.cols()-1);
P.col(1)=P.col(1).cwiseMax(0).cwiseMin(Fx.rows()-1);


return P;
}
*/

