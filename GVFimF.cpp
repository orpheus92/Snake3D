#include "GVFimF.h"
#include <iostream>

//Function to calculate Laplacian and update vector field
/*
void Gupdate(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &u, 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &v,
double sigma,
double mu,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> sMag);
*/
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  GVFimF(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
double mu,
int Giter,
double sigma3,
int style){

//calculate magnitude,

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> sMag;
sMag.resize(Fx.rows(),Fy.cols());
sMag = Fx.array().square()+Fy.array().square();

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> uu;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> vv;
uu.resize(Fx.rows(),Fy.cols());
vv.resize(Fx.rows(),Fy.cols());
uu = Fx;
vv = Fy;

//Function to calculate Laplacian and update vector field
for (int i = 0; i < Giter; i++){
uu = uu + mu*(imDev(uu, sigma3, 3) + imDev(uu, sigma3, 4))- sMag.cwiseProduct(uu-Fx);
vv = vv + mu*(imDev(vv, sigma3, 3) + imDev(vv, sigma3, 4))- sMag.cwiseProduct(vv-Fy);
}

if (style == 1)
return uu;
else 
return vv;
}

/*
void Gupdate(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fx,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fy,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &u, 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &v,
double sigma,
double mu,
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> sMag)
{
u = u + mu*(imDev(u, sigma, 3) + imDev(u, sigma, 4))- sMag.cwiseProduct(u-Fx);

v = v + mu*(imDev(v, sigma, 3) + imDev(v, sigma, 4))- sMag.cwiseProduct(v-Fy);

return;
}
*/
