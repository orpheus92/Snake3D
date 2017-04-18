#include "imDev.h"
#include <iostream>

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorXd;

arma::cube compute(arma::cube X, arma::cube W, bool flip, int stride_x, int stride_y, int stride_z);

void imfilter2a(Eigen::MatrixXd src_image,Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image);

arma::cube imDev(
arma::cube Img,
int sigma,
int type
// type 1 = x; 2 = y; 3 = z; 4 = xx; 5 = yy; 6 = zz; 7 = xy; 8 = xz; 9 = yz
){
int total = 2*3*sigma+1;

arma::cube x(total,1,1);
arma::cube y(1,total,1);
arma::cube z(1,1,total);
for (int i = 0; i< total;i++){

		x(i,0,0) = (double)i-sigma*3;
		y(0,i,0) = (double)i-sigma*3;
		z(0,0,i) = (double)i-sigma*3;
}

int coli = Img.n_cols;
int rowi = Img.n_rows;
int slii = Img.n_slices;

arma::cube out(rowi,coli,slii);

arma::cube Gx(total,1,1);
arma::cube Gy(1,total,1);
arma::cube Gz(1,1,total);

//Gaussian Filters
//arma?
std::cout<<"problem starts here"<<std::endl;
if (type == 1){ //x
Gx = -(exp((-square(x)/2/sigma/sigma)) % x)/(pow(2*M_PI,1.5)*pow(sigma,5));
//std::cout<<"0 = "<<-exp((-square(x)/2/sigma/sigma));
//std::cout<<"1st = "<<x;
//std::cout<<"0 = "<<pow(2*M_PI,1.5);
//std::cout<<"1 = "<<pow(sigma,5);
//std::cout<<"2nd = "<<1/(pow(2*M_PI,3/2)*pow(sigma,5));
//std::cout<<"Gx = "<<Gx<<std::endl;
Gy = exp(-square(y)/2/sigma/sigma);
Gz = exp(-square(z)/2/sigma/sigma);
}
else if (type == 2){ //y
Gy = -(exp((-square(y)/2/sigma/sigma)) % y)/(pow(2*M_PI,1.5)*pow(sigma,5));
Gz = exp(-square(z)/2/sigma/sigma);
Gx = exp(-square(x)/2/sigma/sigma);
}
else if (type == 3){//z
Gz = -(exp((-square(z)/2/sigma/sigma)) % z)/(pow(2*M_PI,1.5)*pow(sigma,5));
Gy = exp(-square(y)/2/sigma/sigma);
Gx = exp(-square(x)/2/sigma/sigma);
}
else if (type == 4){//xx
Gx = (exp((-square(x)/2/sigma/sigma)) % (square(x)/sigma/sigma-1))/(pow(2*M_PI,1.5)*pow(sigma,5));
Gy = exp(-square(y)/2/sigma/sigma);
Gz = exp(-square(z)/2/sigma/sigma);
}
else if (type == 5){//yy
Gy = (exp((-square(y)/2/sigma/sigma)) % (square(y)/sigma/sigma-1))/(pow(2*M_PI,1.5)*pow(sigma,5));
Gx = exp(-square(x)/2/sigma/sigma);
Gz = exp(-square(z)/2/sigma/sigma);
}
else if (type == 6){//zz
Gz = (exp((-square(z)/2/sigma/sigma)) % (square(z)/sigma/sigma-1))/(pow(2*M_PI,1.5)*pow(sigma,5));
Gy = exp(-square(y)/2/sigma/sigma);
Gx = exp(-square(x)/2/sigma/sigma);
}
else if (type == 7){//xy
Gx = (exp(-square(x)/2/sigma/sigma) % x)/(pow(2*M_PI,1.5)*pow(sigma,7));
Gy = (exp(-square(y)/2/sigma/sigma) % y);
Gz = exp(-square(z)/2/sigma/sigma);
}
else if (type == 8){//xz
Gx = (exp(-square(x)/2/sigma/sigma) % x)/(pow(2*M_PI,1.5)*pow(sigma,7));
Gz = (exp(-square(z)/2/sigma/sigma) % z );
Gy = exp(-square(y)/2/sigma/sigma);
}
else {//yz
Gy = (exp(-square(y)/2/sigma/sigma) % y)/(pow(2*M_PI,1.5)*pow(sigma,7));
Gz = (exp(-square(z)/2/sigma/sigma) % z);
Gx = exp(-square(x)/2/sigma/sigma);
}
//G(i,j)=1/(2*M_PI*pow(sigma,6))*(x(i,j)*y(i,j)) *exp(-(pow(x(i,j),2)+pow(y(i,j),2))/(2*pow(sigma,2)));

//}//end for j
//end for i
//std::cout<<"Gx = "<<Gx<<std::endl;
//std::cout<<"Gy = "<<Gy<<std::endl;
//std::cout<<"Gz = "<<Gz<<std::endl;
//use 2D convolution as the filter
//If true, use correlation 
std::cout<<"reach here"<<std::endl;
//std::cout<<"Img20 = "<< Img.slice(20) <<std::endl;
out = compute(Img, Gx, false, 1, 1, 1);
out = compute(out, Gy, false, 1, 1, 1);
out = compute(out, Gz, false, 1, 1, 1);
//std::cout<<"out20 = "<< out.slice(20) <<std::endl;
return out;

}


void imfilter2a(Eigen::MatrixXd src_image, Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image)
{	Eigen::MatrixXd filterconv;
	filterconv = filter.rowwise().reverse().colwise().reverse();
	dst_image = Eigen::MatrixXd::Zero(src_image.rows(),src_image.cols()) ;

	int start_row = int(filter.rows()/2) ;
	int start_col = int(filter.cols()/2) ;

	Eigen::MatrixXd Mid_Matrix = Eigen::MatrixXd::Zero(src_image.rows()+2*start_row,src_image.cols()+2*start_col) ;

	for (int i = 0; i < src_image.rows(); i++)
	{
		for (int j = 0; j < src_image.cols(); j++)
		{
			Mid_Matrix(i+start_row,j+start_col) = src_image(i,j) ;
		}
	}

	int end_row = Mid_Matrix.rows() -1 - start_row ;
	int end_col = Mid_Matrix.cols() -1 - start_col ;

	int filter_row = filter.rows();
	int filter_col = filter.cols() ;
	
	for (int i = start_row; i <= end_row; i++)
	{
		for (int j = start_col; j <= end_col; j++)
		{			
			int tmp_row = i - start_row  ;
			int tmp_col = j - start_col  ;
			for (int m = 0; m < filter_row; m++)
			{				
				for (int n = 0; n < filter_col; n++)
				{
					dst_image(tmp_row,tmp_col) += Mid_Matrix(tmp_row + m,tmp_col + n)*filterconv(m,n) ; 
				}
			}
		}
	}

	return ;
}



arma::cube compute(arma::cube X, arma::cube W, bool flip, int stride_x, int stride_y, int stride_z)
	{
	int width = X.n_cols;
	int height = X.n_rows;
	int slice = X.n_slices;
	
	int kx = W.n_cols;
	int ky = W.n_rows;
	int kz = W.n_slices;
	
	int rx = (kx-1)/2;
	int ry = (ky-1)/2;
	int rz = (kz-1)/2;
	arma::cube out(height,width,slice);
	//out.resize(height, width);
           for (int x=0; x<width; x+=stride_x)
           {
               int xout = x/stride_x;
   
               for (int y=0; y<height; y+=stride_y)
               {
                   int yout = y/stride_y;
                   for (int z = 0; z < slice; z += stride_z)
                   {
                   	
               		   int zout = z/stride_z;
   
		               double sum = 0;

				       for (int x1=x-rx; x1<=x+rx; x1++)
				       {
				           //int wx = flip ? x1-x+rx : rx-x1+x;
				           for (int y1=y-ry; y1<=y+ry; y1++)
				           {
				           		for (int z1 = z-rz; z1<=z+rz;z1++)
				           		{
				           		
						           if (x1>=0 && y1>=0 && z1>=0 && x1<width && y1<height && z1<slice)
						           {
						               if (flip)
						                   sum += W(y1-y+ry,x1-x+rx,z1-z+rz)*X(y1,x1,z1);
						               else
						                   sum += W(ry-y1+y,rx-x1+x,rz-z1+z)*X(y1,x1,z1);
						           }
						           
				           		}
				           }
				       }
		               out(yout,xout,zout) = sum;
		               
		               
		       		}        
               }
           }
	return out;       
	}



