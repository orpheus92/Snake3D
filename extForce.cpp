#include "extForce.h"
#include <iostream>
//function to calculate the external force
//modified version for 3D dataset
//function to filter an image with a designed filter 
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorXd;

//Imfilter function that uses correlation
arma::cube compute2(arma::cube X, arma::cube W, bool flip, int stride_x, int stride_y, int stride_z);
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
//std::cout<<"Iy20 = "<<Iy.slice(20)<<std::endl;
//std::cout<<"Iy35 = "<<Iy.slice(35)<<std::endl;
Iz = imDev(Img, sigma, 3);
//Ix.save("myIx", arma::arma_ascii);

//std::cout<<" find the problem "<<"Ix= "<< Ix << std::endl;
//Line Energy
arma::cube Eline(row,col,sli);
Eline = imgaussian(Img,sigma);
//Eline.save("myEline", arma::arma_ascii);
//std::cout<<" Eline20 =  "<<Eline.slice(20) << std::endl;
//std::cout<<" Eline38 =  "<<Eline.slice(38) << std::endl;
//Edge Energy
arma::cube Eedge(row,col,sli);
Eedge = sqrt(square(Ix) + square(Iy) + square(Iz)); 
//std::cout<<" Edge20 =  "<<Eedge.slice(20) << std::endl;
//std::cout<<" Edge38 =  "<<Eedge.slice(38) << std::endl;
//Eedge.save("myEedge", arma::arma_ascii);
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
	arma::cube xxx(siz+1,1,1);
	arma::cube y(1,siz+1,1);
	arma::cube yyy(1,siz+1,1);
	arma::cube z(1,1,siz+1);
	arma::cube zzz(1,1,siz+1);
	
	//arma::cube Img1; 
	//arma::cube Img2;
	arma::cube out;  
	//std::cout<<"enter gaussian"<<std::endl;
	for (int ii = 0; ii<siz+1; ii++){
	
	x(ii,0,0) = (double) ii- siz/2;
	y(0,ii,0) = (double) ii- siz/2;
	z(0,0,ii) = (double) ii- siz/2;
	}
	
	arma::cube Hx(siz+1,1,1);
	arma::cube Hy(1,siz+1,1);
	arma::cube Hz(1,1,siz+1);
	
	xxx.fill(as_scalar(sum(exp(-square(x)/2/Sigma/Sigma))));
	yyy.fill(as_scalar(sum(exp(-square(x)/2/Sigma/Sigma))));
	zzz.fill(as_scalar(sum(exp(-square(x)/2/Sigma/Sigma))));
	Hx = exp(-square(x)/2/Sigma/Sigma)/xxx;
	Hy = exp(-square(y)/2/Sigma/Sigma)/yyy;
	Hz = exp(-square(z)/2/Sigma/Sigma)/zzz;
	
	//std::cout<<"Hx = "<<Hx<<std::endl;
	//std::cout<<"Hy = "<<Hy<<std::endl;
	//std::cout<<"Hz = "<<Hz<<std::endl;
	
	out = compute2(Img, Hx, false, 1, 1, 1);
	out = compute2(out, Hy, false, 1, 1, 1);
	out = compute2(out, Hz, false, 1, 1, 1);
	//imfilter2(Img, Hx, Img1);
	//imfilter2(Img1, Hy, Img2);
	//imfilter2(Img2, Hz, out);	
	//std::cout<<out.n_rows<<" "<<out.n_cols<<" "<<std::endl;
	return out;
	}
else
return Img;

}


arma::cube compute2(arma::cube X, arma::cube W, bool flip, int stride_x, int stride_y, int stride_z)
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

