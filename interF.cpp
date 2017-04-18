#include "interF.h"
#include <iostream>
#include <stb_image.h>

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
//std::cout<<"in interF"<<std::endl;
int nV = V.rows();

arma::field<arma::vec> Nebr(nV);


Nebr = vNbnr(V, F);

Eigen::SparseMatrix<double> regmat(nV,nV); 
//std::cout<<"spmat check"<<std::endl;
  // Act like matlab's [C,IA,IC] = unique(X)
  //
  // Templates:
  //   T  comparable type T
  // Inputs:
  //   A  #A vector of type T
  // Outputs:
  //   C  #C vector of unique entries in A
  //   IA  #C index vector so that C = A(IA);
  //   IC  #A index vector so that A = C(IC);

for (int vv = 0; vv<nV; vv++)
{    
// Add the neighbours
	  //Eigen::VectorXi pre = Nebr.row(vv);
	  //Eigen::VectorXi curR;
	 // igl::viewer::Viewer viewer;
	 // igl::unique(pre,curR);
//iterator used for uncertain size of Nebr pts
	 arma::vec curR = arma::unique(Nebr(vv));
//std::cout<<"UniqueNe = "<<curR<<std::endl;
//arma::vec::iterator a = Nebr(vv).begin();
//arma::vec::iterator b = Nebr(vv).end();

for(int iterC = 0; iterC<curR.n_elem; iterC++)
  {
 // cout << *i << endl;
 regmat.insert(vv,curR(iterC))=(double)1/curR.n_elem;
  }



// Add the vertex it self 
regmat.insert(vv,vv)=-1;
}
//std::cout<<"sparmat = "<<regmat<<std::endl;
//std::cout<<"after for loop in interF"<<std::endl;
Eigen::SparseMatrix<double> idmat(nV,nV); 
idmat.setIdentity();
//std::cout<<"idmat = "<<idmat<<std::endl;
//      
//return type is a sparse regularization matrix
Eigen::SparseMatrix<double> solvmat(nV,nV); 

solvmat = gamma*idmat-alpha*regmat+beta*regmat*regmat;
//std::cout<<"part1 = "<<(gamma*idmat).rowwise().sum()<<std::endl;
//std::cout<<"part2 = "<<-(alpha*regmat).rowwise().sum()<<std::endl;
//std::cout<<"part3 = "<<(beta*regmat*regmat).rowwise().sum()<<std::endl;
//std::cout<<"part4 = "<<solvmat<<std::endl;

Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
solver.compute(solvmat);
Eigen::SparseMatrix<double> I(nV,nV);
I.setIdentity();
return solver.solve(I);
//return ;


}

arma::field<arma::vec> vNbnr(
Eigen::MatrixXd V, //vertices
Eigen::MatrixXi F  //Faces);
)
{
int r = V.rows();
//Eigen::MatrixXi Nb;
//arma::field<arma::mat> Ne(r); // Neighbor vertices, field is used for different # of vals
//arma::field<arma::vec> Neunique(r); //return type
arma::field<arma::vec> Nb(F.rows()); //field of column vectors
//Nb.resize(r, 6);
std::cout<<"prb here"<<std::endl;
arma::vec mytemp = arma::vec(1);
for (int fi = 0; fi<F.rows(); fi++)
{

//Nb.row(fi)<<F(fi,1), F(fi,2), F(fi,2), F(fi,0), F(fi,0), F(fi,1);
	Nb(F(fi,0)) = arma::join_cols( Nb(F(fi,0)), mytemp.fill(F(fi,1))  );
	Nb(F(fi,0)) = arma::join_cols( Nb(F(fi,0)), mytemp.fill(F(fi,2))  );
	Nb(F(fi,1)) = arma::join_cols( Nb(F(fi,1)), mytemp.fill(F(fi,2))  );
	Nb(F(fi,1)) = arma::join_cols( Nb(F(fi,1)), mytemp.fill(F(fi,0))  );
	Nb(F(fi,2)) = arma::join_cols( Nb(F(fi,2)), mytemp.fill(F(fi,0))  );
	Nb(F(fi,2)) = arma::join_cols( Nb(F(fi,2)), mytemp.fill(F(fi,1))  );
	//Ne(F(fi,0)-1) = arma::join_cols( Ne(F(fi,0)-1), F(fi,1) ); 
	//Ne(F(fi,0)-1) = arma::join_cols( Ne(F(fi,0)-1), F(fi,2) );
	
	//Ne(F(fi,1)-1) = arma::join_cols( Ne(F(fi,1)-1), F(fi,2) ); 
	//Ne(F(fi,1)-1) = arma::join_cols( Ne(F(fi,1)-1), F(fi,0) ); 

	//Ne(F(fi,2)-1) = arma::join_cols( Ne(F(fi,2)-1), F(fi,0) ); 
	//Ne(F(fi,2)-1) = arma::join_cols( Ne(F(fi,2)-1), F(fi,1) ); 	
}
//std::cout<<"Nb = "<<Nb<<std::endl;
/*for (int vi = 0; vi<r; vi++)
{	
	// unique function in armadillo, gives vec 
	Neunique(vi) = arma::unique(Ne(vi));
	
}
*/
return Nb;
} 


