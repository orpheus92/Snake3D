#include "interF.h"
#include <iostream>


arma::field<arma::vec> vNbnr(
Eigen::MatrixXd V, //vertices
Eigen::MatrixXi F  //Faces);
);

Eigen::SparseMatrix<double>  interF(
Eigen::MatrixXd V, //vertices
Eigen::MatrixXi F, //Faces
double alpha,
double beta,
double gamma)
{

int nV = V.rows();
arma::field<arma::vec> Nebr(nV);
Nebr = vNbnr(V, F);

Eigen::SparseMatrix<double> regmat(nV,nV); 

for (int vv = 0; vv<nV; vv++)
{    
// Add the neighbours

//iterator used for uncertain size of Nebr pts

arma::vec::iterator a = Nebr(vv).begin();
arma::vec::iterator b = Nebr(vv).end();

for(arma::mat::iterator i=a; i!=b; ++i)
  {
 // cout << *i << endl;
 regmat(vv,(int)*i)=1/size(Nebr(vv),1);
  }



// Add the vertex it self 
regmat(vv,vv)=-1;
}

Eigen::SparseMatrix<double> idmat(nV,nV); 

//return type is a sparse regularization matrix
return (gamma*idmat.setIdentity()-alpha*regmat+beta*regmat*regmat).inverse();


}

arma::field<arma::vec> vNbnr(
Eigen::MatrixXd V, //vertices
Eigen::MatrixXi F  //Faces);
)
{
int r = V.rows();
arma::field<mat> Ne(r); // Neighbor vertices, field is used for different # of vals
arma::field<vec> Neunique(r); //return type
for (int fi = 0; fi<F.rows(); fi++)
{
	Ne(F(fi,0)-1) = join_cols( Ne(F(fi,0)-1), F(fi,1) ); 
	Ne(F(fi,0)-1) = join_cols( Ne(F(fi,0)-1), F(fi,2) );
	
	Ne(F(fi,1)-1) = join_cols( Ne(F(fi,1)-1), F(fi,2) ); 
	Ne(F(fi,1)-1) = join_cols( Ne(F(fi,1)-1), F(fi,0) ); 

	Ne(F(fi,2)-1) = join_cols( Ne(F(fi,2)-1), F(fi,0) ); 
	Ne(F(fi,2)-1) = join_cols( Ne(F(fi,2)-1), F(fi,1) ); 	
}

for (int vi = 0; vi<r; vi++)
{	
	// unique function in armadillo, gives vec 
	Neunique(vi) = unique(Ne(vi));
	
}

return Neunique;
} 


