#include "snake3D.h"
#include <stb_image.h>
#include <iostream>
#include "boost/multi_array.hpp"
#include <cassert>


Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interp3(arma::cube Fx, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> V);
//baloon Force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baloonF(Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> F,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> V);

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snakeMove3D(
Eigen::SparseMatrix<double> intForce, //internal force
Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> F,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> V,//contour pts
arma::cube Fx,
arma::cube Fy,
arma::cube Fz,//external vec field
double gamma, //time step
double kappa,//external field weight
double delta,
double lamb);

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
){
//std::cout<<input.slice(9)<<std::endl;
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
Eext.save("myEext", arma::arma_ascii);
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
//std::cout<<intForce<<std::endl;

// Triangulated interior

for (int it = 0;it<0;it++){


V = snakeMove3D(intForce, F, V, Fx, Fy, Fz, gamma, kappa, delta, lamb);
  std::cout<<"V = "<<V<<std::endl;

  //  
  }
  std::cout<<"end of snake"<<std::endl;
return V; 
}

//function to calculate baloon force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baloonF(
Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> F,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> V)
{
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> e1;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> e2;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> e3;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> e1n;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> e2n;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> e3n;
e1.resize(F.rows(),F.cols());
e2.resize(F.rows(),F.cols());
e3.resize(F.rows(),F.cols());
e1n.resize(F.rows(),F.cols());
e2n.resize(F.rows(),F.cols());
e3n.resize(F.rows(),F.cols());
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Normal;
Normal.resize(F.rows(),F.cols());
for (int cure =0; cure<F.rows();cure++)
	{
	e1.row(cure) = V.row(F(cure,0))-V.row(F(cure,1));
	e2.row(cure) = V.row(F(cure,1))-V.row(F(cure,2));
	e3.row(cure) = V.row(F(cure,2))-V.row(F(cure,0));
	Normal.row(cure)<<e1(cure,1)*e3(cure,2)-e1(cure,2)*e3(cure,1),e1(cure,2)*e3(cure,0)-e1(cure,0)*e3(cure,2),e1(cure,0)*e3(cure,1)-e1(cure,1)*e3(cure,0);
	}
e1n = e1.array()/(((e1.col(0).array().square()+e1.col(1).array().square()+e1.col(2).array().square()).sqrt()).replicate(1,3));
e2n = e2.array()/(((e2.col(0).array().square()+e2.col(1).array().square()+e2.col(2).array().square()).sqrt()).replicate(1,3));
e3n = e3.array()/(((e3.col(0).array().square()+e3.col(1).array().square()+e3.col(2).array().square()).sqrt()).replicate(1,3));

//Eigen::MatrixXd e1nt = e1n.adjoint();
//Eigen::MatrixXd e2nt = e2n.adjoint();
//Eigen::MatrixXd e3nt = e3n.adjoint();
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> ang1 = (-e1n.array()*e3n.array()).rowwise().sum().acos();//.array().acos();
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> ang2 = (-e2n.array()*e1n.array()).rowwise().sum().acos();//.array().acos();
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> ang3 = (-e3n.array()*e2n.array()).rowwise().sum().acos();//.array().acos();
std::cout<<ang3.rows()<<ang3.cols()<<" here "<<std::endl;


Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> vNormal = Eigen::MatrixXd::Zero(V.rows(), V.cols());

for (int myi; myi<F.rows();myi++)
	{
	vNormal.row(F(myi,0))=vNormal.row(F(myi,0))+Normal.row(myi)*ang1(myi);
    vNormal.row(F(myi,1))=vNormal.row(F(myi,1))+Normal.row(myi)*ang2(myi);
    vNormal.row(F(myi,2))=vNormal.row(F(myi,2))+Normal.row(myi)*ang3(myi);
	}
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> vNorm;
vNorm = (vNormal.col(0).array().square()+vNormal.col(1).array().square()+vNormal.col(2).array().square()).sqrt()+0.000000000000001;
//Eigen::NumTraits::epsilon()
vNormal = vNormal.array()/vNorm.replicate(1,3).array();


return vNormal;
}

/*
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

*/
// Problem for boundary condition, Fix later
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interp3(arma::cube Fx, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> V)
{
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out;
//int r = Fx.n_rows;
//int c = Fx.n_cols;
//int l = Fx.n_slices;

out.resize(V.rows(),1);
double x;
double y;
double z;

double xd;
double yd;
double zd;

int x0;
int x1;
int y0;
int y1;
int z0;
int z1;

double c00;
double c01;
double c10;
double c11;
double c0;
double c1;
//double q110;
//double q111;


//double x1;
//double x2;
//double y1;
//double y2;

//bilinear interpolation
for (int i = 0; i<V.rows();i++){

y = V(i,0);
x = V(i,1);
z = V(i,2);

x0 = (int)ceil(x)-1;
x1 = (int)ceil(x);
y0 = (int)ceil(y)-1;
y1 = (int)ceil(y);
z0 = (int)ceil(z)-1;
z1 = (int)ceil(z);


xd = (x-(double)x0);
yd = (y-(double)y0);
zd = (z-(double)z0);
	if(x0<0)
		x0 = 0;
	if(y0<0)
		y0 = 0;
	if(z0<0)
		z0 = 0;
	c00 = Fx(y0,x0,z0)*(1-xd)+Fx(y0,x1,z0)*xd;
	c01 = Fx(y0,x0,z1)*(1-xd)+Fx(y0,x1,z1)*xd;
	c10 = Fx(y1,x0,z0)*(1-xd)+Fx(y1,x1,z0)*xd;
	c11 = Fx(y1,x0,z1)*(1-xd)+Fx(y1,x1,z1)*xd; 	
	c0 = c00*(1-yd)+c10*yd;
	c1 = c01*(1-yd)+c11*yd;
	out(i,0) = c0*(1-zd)+c1*zd;
	}



return out;




}


//This function will calculate one iteration of contour Snake movement
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> snakeMove3D(
Eigen::SparseMatrix<double> intForce, //internal force
Eigen::Matrix<int, Eigen::Dynamic,Eigen::Dynamic> F,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> V,//contour pts
arma::cube Fx,
arma::cube Fy,
arma::cube Fz,//external vec field
double gamma, //time step
double kappa,//external field weight
double delta,
double lamb) //Balloon Force Weight)
{
//Clamp contour to boundary
//rows and cols might be checked later 
//image force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> imforce;
imforce.resize(V.rows(),V.cols());

//baloon force
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> baforce;
baforce.resize(V.rows(),V.cols());
//clamp boundary
V.col(0)=V.col(0).cwiseMax(0).cwiseMin(Fx.n_cols-1);
V.col(1)=V.col(1).cwiseMax(0).cwiseMin(Fx.n_rows-1);
V.col(2)=V.col(2).cwiseMax(0).cwiseMin(Fx.n_slices-1);
//std::cout<<V<<std::endl;
//std::cout<<"======= first V"<<std::endl;

//Get image force on the contour points
imforce.col(0)=kappa*interp3(Fx,V);
imforce.col(1)=kappa*interp3(Fy,V);
imforce.col(2)=kappa*interp3(Fz,V);
//std::cout<<"imforce = "<<imforce<<std::endl;
//Get baloon force on the contour points
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> N;
N.resize(V.rows(),V.cols());
N = baloonF(F,V);
baforce = delta * N; 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fext3;
Fext3.resize(V.rows(),V.cols());
Fext3 = (imforce.col(0).array()*N.col(0).array()+imforce.col(1).array()*N.col(1).array()+imforce.col(2).array()*N.col(2).array()).replicate(1,3)*N.array();

//Update contour positions
V.col(0) = intForce * (gamma*V.col(0)) + imforce.col(0)*(1-lamb) + Fext3.col(0)*lamb + baforce.col(0);
V.col(1) = intForce * (gamma*V.col(1)) + imforce.col(1)*(1-lamb) + Fext3.col(1)*lamb + baforce.col(1);
V.col(2) = intForce * (gamma*V.col(2)) + imforce.col(2)*(1-lamb) + Fext3.col(2)*lamb + baforce.col(2);
 

//Clamp contour to boundary
V.col(0)=V.col(0).cwiseMax(0).cwiseMin(Fx.n_cols-1);
V.col(1)=V.col(1).cwiseMax(0).cwiseMin(Fx.n_rows-1);
V.col(2)=V.col(2).cwiseMax(0).cwiseMin(Fx.n_slices-1);

std::cout<<"end of snake move"<<std::endl;
return V;
}


