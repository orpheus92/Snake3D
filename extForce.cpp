#include "extForce.h"
#include <iostream>
//function to calculate the external force
//modified version for 3D dataset
//function to filter an image with a designed filter 
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorXd;

//Imfilter function that uses correlation
MatrixXd compute2(MatrixXd X, MatrixXd W, bool flip, int stride_x, int stride_y);
void imfilter2(Eigen::MatrixXd src_image,Eigen::MatrixXd filter,Eigen::MatrixXd &dst_image);

//IMGAUSSIAN filters a 3D greyscale image with an Gaussian filter.
arma::cube imgaussian(
arma::cube Img, 
int Sigma);

arma::cube extForce(
arma::cube Img, //use armadillo for 3D data
double wl,
double we,
int sigma
//int npts
){ int total = 2*3*sigma+1;
int row = Img.n_rows;
int col = Img.n_cols;
int sli = Img.n_slices;

arma::cube Ix(row,col,sli);
arma::cube Iy(row,col,sli);
arma::cube Iz(row,col,sli);
arma::cube out(row,col,sli);


Ix = imDev(Img, sigma, 1);
Iy = imDev(Img, sigma, 2);
Iz = imDev(Img, sigma, 3);



//Line Energy
arma::cube Eline(row,col,sli);
Eline = imgaussian(Img,sigma);

//Edge Energy
arma::cube Eedge(row,col,sli);
Eedge = sqrt(square(Ix) + square(Iy) + square(Iz)); 

//Externa =l Energy
arma::cube Eext(row,col,sli);

Eext = wl*Eline - we*Eedge; 

return Eext;
}

//Imfilter function that uses correlation 
void imfilter2(arma::cube src_image, arma::cube filter,arma::cube &dst_image)
{
	dst_image = arma::cube(src_image.n_rows,src_image.n_cols, src_image.n_slices) ;	
	
	int start_row = int(filter.n_rows/2) ;
	int start_col = int(filter.n_cols/2) ;
	int start_sli = int(filter.n_slices/2) ;
	//std::cout<<"start_row = "<<start_row<<std::endl;
	arma::cube Mid_Matrix = arma::cube(src_image.n_rows+2*start_row, src_image.n_cols+2*start_col, src_image.n_slices+2*start_sli) ;

	for (int i = 0; i < src_image.n_rows; i++)
	{
		for (int j = 0; j < src_image.n_cols; j++)
		{
			for (int k = 0; k<src_image.n_slices; k++)
			{
				Mid_Matrix(i+start_row,j+start_col, k+start_sli) = src_image(i,j,k) ;
			}
		}
	}

	int end_row = Mid_Matrix.n_rows -1 - start_row ;
	int end_col = Mid_Matrix.n_cols -1 - start_col ;
	int end_sli = Mid_Matrix.n_slices -1 - start_sli ;

	int filter_row = filter.n_rows;
	int filter_col = filter.n_cols;
	int filter_sli = filter.n_slices;
	
	for (int i = start_row; i <= end_row; i++)
	{
		for (int j = start_col; j <= end_col; j++)
		{			
			for(int k = start_sli; k<=end_sli; k++)
			{
				int tmp_row = i - start_row;
				int tmp_col = j - start_col;
				int tmp_sli = k - start_sli;
				for (int m = 0; m < filter.n_rows; m++)
				{				
					for (int n = 0; n < filter.n_cols; n++)
					{
						for (int p = 0;p < filter.n_slices; p++)
						{
							dst_image(tmp_row,tmp_col,tmp_sli) -= Mid_Matrix(tmp_row + m,tmp_col + n, tmp_sli+p)*filter(m,n,p) ; 
						}					
					}
				}
			}
		}
	}
	return ;
}

arma::cube imgaussian(
arma::cube Img, 
int Sigma){
int siz = Sigma*6;

if (Sigma>0){
	// Make 1D Gaussian kernel
	    
	arma::cube x(siz+1,1,1);
	arma::cube Img1; 
	arma::cube Img2;
	arma::cube out;  
	for (int ii = 0; ii<siz+1; ii++){
	
	x(ii,1,1) = (double) ii- siz/2;
	
	}
	
	arma::cube Hx(siz+1,1,1);
	Hx = exp(-square(x)/2/Sigma/Sigma)/sum(exp(-square(x)/2/Sigma/Sigma));

	arma::cube Hy = resize(Hx, 1, siz+1, 1);
	arma::cube Hz = resize(Hx, 1, 1, siz+1);
	
	imfilter2(Img, Hx, Img1);
	imfilter2(Img1, Hy, Img2);
	imfilter2(Img2, Hz, out);

	return out;
	}
else
return Img;

}


MatrixXd compute2(MatrixXd X, MatrixXd W, bool flip, int stride_x, int stride_y)
	{
	int width = X.cols();
	int height = X.rows();
	int kx = W.cols();
	int ky = W.rows();
	int rx = (kx-1)/2;
	int ry = (ky-1)/2;
	MatrixXd out;
	out.resize(height, width);
           for (int x=0; x<width; x+=stride_x)
           {
               int xout = x/stride_x;
   
               for (int y=0; y<height; y+=stride_y)
               {
                   int yout = y/stride_y;
   
                   double sum = 0;

		           for (int x1=x-rx; x1<=x+rx; x1++)
		           {
		               int wx = flip ? x1-x+rx : rx-x1+x;
		               for (int y1=y-ry; y1<=y+ry; y1++)
		               {
		                   if (x1>=0 && y1>=0 && x1<width && y1<height)
		                   {
		                       if (flip)
		                           sum += W(y1-y+ry,wx)*X(y1,x1);
		                       else
		                           sum += W(ry-y1+y,wx)*X(y1,x1);
		                   }
		               }
		           }
                   out(yout,xout) = sum;
               }
           }
	return out;       
	}



