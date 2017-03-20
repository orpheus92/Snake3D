#include "interpcont.h"
#include <iostream>

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> interpcont(
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P,
int npts
){
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> out;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dis;

out.resize(npts,2);

// May be modified later if npts is used as params to set up # of pts on the contour
int row;
row = P.rows();
dis.resize(row,1);
//
for (int ii =0; ii<row-1; ii++){
dis(ii,0) = sqrt(pow(P(ii+1,0)-P(ii,0),2) + pow(P(ii+1,1)-P(ii,1),2));

}

for (int x =1; x<row-1; x++){
dis(x,0) = dis(x,0)+dis(x-1,0);

}

dis(row-1,0)=dis(row-2,0)+sqrt(pow(P(0,0)-P(row-1,0),2) + pow(P(0,1)-P(row-1,1),2));

double temp;
double tempx;
double tempy;
double totald = dis(row-1,0);
out(0,0)=P(0,0);
out(0,1)=P(0,1);
for (int i = 1; i<npts; i++){

	temp = (double)i/npts*totald;

	for (int j = 0; j<row-1; j++){
		if(temp<dis(0,0)){
			tempx = temp/dis(0,0)*(P(1,0)-P(0,0))+P(0,0);
			tempy = temp/dis(0,0)*(P(1,1)-P(0,1))+P(0,1);

			out(i,0)=tempx;
			out(i,1)=tempy;

		}
		else if(temp>dis(j,0)&(temp<=dis(j+1,0))){

			if(j+2>=row)
			{
			tempx = (temp-dis(j,0))/(dis(j+1,0)-dis(j,0))*(P(0,0)-P(j+1,0))+P(j+1,0);
			tempy = (temp-dis(j,0))/(dis(j+1,0)-dis(j,0))*(P(0,1)-P(j+1,1))+P(j+1,1);

			out(i,0)=tempx;
			out(i,1)=tempy;
			}

			else{
			tempx = (temp-dis(j,0))/(dis(j+1,0)-dis(j,0))*(P(j+2,0)-P(j+1,0))+P(j+1,0);
			tempy = (temp-dis(j,0))/(dis(j+1,0)-dis(j,0))*(P(j+2,1)-P(j+1,1))+P(j+1,1);

			out(i,0)=tempx;
			out(i,1)=tempy;
			}

		}
	}

}

return out;
}
