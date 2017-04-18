#ifdef IGL_VIEWERWITH_NANOGUI
int main(){
  std::cerr<<
  "Error: recompile with LIBIGL_VIEWER_WITH_NANOGUI defined"<<std::endl;
  return EXIT_FAILURE;
}
#else
#include <iostream>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/triangle/triangulate.h>
#include <igl/png/writePNG.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <igl/png/readPNG.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "igl/unique.h"
#include "igl/readCSV.h"
#include <armadillo>


#include "interpcont.h"
#include "extForce.h"
#include "imDev.h"
#include "snake3D.h"


#include <string>

using namespace std;
// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi E;
Eigen::MatrixXd H;


// Input contour Vertices
Eigen::MatrixXd xy;
Eigen::MatrixXi xyedge;
//Eigen::MatrixXd y;

// Triangulated interior
Eigen::MatrixXd V2;
Eigen::MatrixXi F2;

Eigen::MatrixXd V3;
Eigen::MatrixXi F3;
igl::viewer::Viewer viewer2;

// GLOBAL VARIABLES
float alpha = 0.2;
float beta = 0.2;
float iterations = 20;
// Function to press key; called every time a keyboard button is pressed 
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    // Allocate temporary buffers
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

    // Draw the scene in the buffers
    viewer.core.draw_buffer(viewer.data,viewer.opengl,false,R,G,B,A);

    // Save it to a PNG
    igl::png::writePNG(R,G,B,A,"out.png");
  }

  if (key == '2')
  {
  /*
    // Allocate temporary buffers
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A,I;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> temp;
    // Read the PNG
    temp = igl::png::readPNG("myim2.png",R,G,B,A,temp);

    // Replace the mesh with a triangulated square
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P;
	P.resize(7,2);
	P <<
	38,25,
	20,59,
	39,97,
	81,105,
	109,84,
	112,39,
	73,24;

//output contour
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> outcont;	
*/


Eigen::MatrixXd V; //vertices
Eigen::MatrixXi F;
arma::cube input(64,64,64);
//arma::cube input;
igl::readCSV("newV",V);
igl::readCSV("newF",F);
//Eigen::MatrixXi F2 = F-Eigen::MatrixXd::Ones(F.rows(),F.cols());   
input.load("target");
input.reshape(64,64,64);
//std::cout<<F.rows()<<"  "<< F.cols()<<"  "<<V.rows()<<" "<<std::endl;
std::cout<<input.n_rows<<" "<<input.n_cols<<" "<<input.n_slices<<std::endl;

viewer2.data.set_mesh(V,F);
//viewer2.core.align_camera_center(V,F);
cout << "alpha = " << alpha<<endl;
cout << "beta = " << beta<<endl;
cout << "iteratopms = " << iterations <<endl;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> outV;
outV =  snake3D(
V, //vertices
F, //Faces
input, //input 3D image
0.1, //time step, default 1
iterations, //# of iteration, default 100
//int npts, //# of pts to interpolate contours 
2, // sigma to calculate img derivative default 10
-1, //attraction to lines, < 0 to black line; > 0 to white line; default 0.04
0, //attraction to edge, default 2.0
//double wt, //attraction to end points, default 0.01
2, //sigma to calculate gradient of edge energy image (give image force), default 20
alpha,//double alpha, //membrane energy, default 0.2
beta,//double beta, //thin plate energy, default 0.2
0.1,//double delta, //baloon force, default 0.1
0.5,//double kappa, //weight of external img force, default 2
// the following is used for GVF snake
0.2,//double mu, //tradeoff between real edge vectors and noise vectors, default 0.2
0,//int Giter, //GVF iteration, default 0
1,//double sigma3, //sigma used to calculate laplacian in GVF, default 1
0.8,//double lamb, //Weight which changes the direction of the image potential force to the direction of the surface normal,
// (Keeps the surface from self intersecting)
viewer
);

/*
outcont = snake3D(V,
F, //initial contour
input, //input 3D image
1, //time step, default 1
20, //# of iteration, default 100
//int npts, //# of pts to interpolate contours 
2, // sigma to calculate img derivative default 1
-1, //attraction to lines, < 0 to black line; > 0 to white line; default 0.04
0, //attraction to edge, default 2.0
//double wt, //attraction to end points, default 0.01
2, //sigma to calculate gradient of edge energy image (give image force), default 20
0.2, //membrane energy, default 0.2
0.2, //thin plate energy, default 0.2
0.1, //baloon force, default 0.1
0.5, //weight of external img force, default 2
// the following is used for GVF snake
0.2, //tradeoff between real edge vectors and noise vectors, default 0.2
0, //GVF iteration, default 0
1, //sigma used to calculate laplacian in GVF, default 1
0.8, //Weight which changes the direction of the image potential force to the direction of the surface normal,
// (Keeps the surface from self intersecting)
viewer
);

std::cout<<"out contour = "<<outcont<<std::endl;

double rsize =(double)temp.rows();
double csize =(double)temp.cols();
Eigen::MatrixXd Vout;
Eigen::MatrixXi Eout;
Eigen::MatrixXd Hout;
Vout.resize(outcont.rows()+P.rows()+4,2);
Eout.resize(outcont.rows()+P.rows()+4,2);
for (int i = 0; i<outcont.rows();i++){
Vout(i,0) = outcont(i,0)-csize/2;
Vout(i,1) = outcont(i,1)-rsize/2;
Eout(i,0) = i;
Eout(i,1) = i+1;
}
Eout(outcont.rows()-1,1) = 0;
for (int j = 0; j<P.rows();j++){
Vout(j+outcont.rows(),0) = (double)(P(j,0)-csize/2);
Vout(j+outcont.rows(),1) = (double)(P(j,1)-rsize/2);
Eout(outcont.rows()+j,0) = outcont.rows()+j;
Eout(j+outcont.rows(),1) = outcont.rows()+j+1;
}

Eout(outcont.rows()+P.rows()-1,1) = Eout(outcont.rows(),0);

//manually set boundary
Vout(outcont.rows()+P.rows(),0) = -csize/2;
Vout(outcont.rows()+P.rows(),1) = -rsize/2;
Vout(outcont.rows()+P.rows()+1,0) = csize/2;
Vout(outcont.rows()+P.rows()+1,1) = -rsize/2;
Vout(outcont.rows()+P.rows()+2,0) = csize/2;
Vout(outcont.rows()+P.rows()+2,1) = rsize/2;
Vout(outcont.rows()+P.rows()+3,0) = -csize/2;
Vout(outcont.rows()+P.rows()+3,1) = rsize/2;


Eout(outcont.rows()+P.rows(),0) = P.rows()+100;
Eout(outcont.rows()+P.rows(),1) = P.rows()+100+1;
Eout(outcont.rows()+P.rows()+1,0) = P.rows()+100+1;
Eout(outcont.rows()+P.rows()+1,1) = P.rows()+100+2;
Eout(outcont.rows()+P.rows()+2,0) = P.rows()+100+2;
Eout(outcont.rows()+P.rows()+2,1) = P.rows()+100+3;
Eout(outcont.rows()+P.rows()+3,0) = P.rows()+100+3;
Eout(outcont.rows()+P.rows()+3,1) = P.rows()+100;

Hout.resize(1,2);
Hout << 0,0;
//igl::triangle::triangulate(Vout,Eout,Hout,"a1q",V2,F2);
std::cout<<"End of Snake"<<std::endl;
   // viewer.data.clear();
   // viewer.data.set_mesh(V2,F2);
   // viewer.core.show_texture = true;

    // Use the image as a texture
  }

*/
 viewer2.data.clear();
    viewer2.data.set_mesh(outV,F);
    viewer2.core.show_texture = false;
    //viewer2.core.invert_normal = true;
viewer2.launch();
	//viewer2.core.show_texture = false;
      //  viewer2.callback_key_down = &key_down;
//	viewer2.launch();

  return false;
}

}


int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
 ///Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R1,G1,B1,A1,I1;
  //  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> temp1;
    // Read the PNG
 //   temp1 = igl::png::readPNG("myim2.png",R1,G1,B1,A1,temp1);

//double rsize1 =(double)temp1.rows();
//double csize1 =(double)temp1.cols();
/*
std::cout<<"r = "<<rsize1<<" c = "<<csize1<<std::endl;
Eigen::MatrixXd Vout1;
Eigen::MatrixXi Eout1;
Eigen::MatrixXd Hout1;
Vout1.resize(5,2);
Eout1.resize(5,2);

//Vout1(1,0) = -csize1/2;
//Vout1(1,1) = -rsize1/2;
Vout1(2,0) = csize1/2;
Vout1(2,1) = -rsize1/2;
Vout1(3,0) = csize1/2;
Vout1(3,1) = rsize1/2;
Vout1(4,0) = -csize1/2;
Vout1(4,1) = rsize1/2;
Vout1(0,0) = 0;
Vout1(0,1) = 0;

Eout1(1,0) = 1;
Eout1(1,1) = 2;
Eout1(2,0) = 2;
Eout1(2,1) = 3;
Eout1(3,0) = 3;
Eout1(3,1) = 4;
Eout1(4,0) = 4;
Eout1(4,1) = 1;
Eout1(0,0) = 0;
Eout1(0,1) = 0;



Hout1.resize(1,2);
Hout1 << 0,0;
  // Create the boundary of a square
 // V.resize(12,2);
 // E.resize(12,2);
 // H.resize(1,2);

//xy.resize(10,2);
//xyedge.resize(10,2);

// V << 100,100, 200,100, 200,200, 100, 200,
  //     50,50, 250,50, 250,250, 50, 250,
  //     0,0, 300,0, 300,300, 0, 300;

/*  V << -1,-1, 1,-1, 1,1, -1, 1,
       -2,-2, 2,-2, 2,2, -2, 2,
       -3,-3, 3,-3, 3,3, -3, 3;
*/
//V = V-150;
//  E << 0,1, 1,2, 2,3, 3,0,
  //     4,5, 5,6, 6,7, 7,4,
  //     8,9, 9,10, 10,11, 11,8;

 // H << 0,1;

 /*
Eigen::MatrixXd V22;
Eigen::MatrixXi F22;

  // Triangulate the interior
        igl::triangle::triangulate(Vout1,Eout1,Hout1,"a1q",V22,F22);
	igl::viewer::Viewer viewer;
//igl::triangle::triangulate(V,E,H,"a5q",V2,F2);
  // Plot the generated mesh
  */

        
	viewer2.data.clear();
  viewer2.callback_init = [&](igl::viewer::Viewer& viewer){
    viewer.ngui->addVariable("ALPHA : ", alpha);
    viewer.ngui->addVariable("BETA : ", beta);
    viewer.ngui->addVariable("ITERATIONS : ", iterations);
    viewer.screen->performLayout();
    return false;
  };
  
	//viewer2.data.set_mesh(V22,F22);
//	viewer2.core.align_camera_center(V22,F22);
	viewer2.core.show_texture = false;
  viewer2.callback_key_down = &key_down;
	viewer2.launch();

// Wait for Key? 
  
}
#endif
