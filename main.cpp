#include <iostream>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/triangle/triangulate.h>
#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "interpcont.h"
#include "extForce.h"
#include "imDev.h"
#include "snake3D.h"

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
outcont = snake3D(P, //initial contour
temp, //input gray-scale image
1, //time step, default 1
400, //# of iteration, default 100
100, //# of pts to interpolate contours 
3, // sigma to calculate img derivative default 10
0, //attraction to lines, < 0 to black line; > 0 to white line; default 0.04 wl
2, //attraction to edge, default 2.0 wedge
0, //attraction to end points, default 0.01 wterm
3, //sigma to calculate gradient of edge energy image (give image force), default 20
0.1, //membrane energy, default 0.2 alpha
0.1, //thin plate energy, default 0.2 beta
-0.1, //baloon force, default 0.1 delta
4, //weight of external img force, default 2 kappa 
// the following is used for GVF snake
0.2, //tradeoff between real edge vectors and noise vectors, default 0.2
100, //GVF iteration, default 0
1, //sigma used to calculate laplacian in GVF, default 1
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


  return false;
}




int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
 Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R1,G1,B1,A1,I1;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> temp1;
    // Read the PNG
    temp1 = igl::png::readPNG("myim2.png",R1,G1,B1,A1,temp1);

double rsize1 =(double)temp1.rows();
double csize1 =(double)temp1.cols();
std::cout<<"r = "<<rsize1<<" c = "<<csize1<<std::endl;
Eigen::MatrixXd Vout1;
Eigen::MatrixXi Eout1;
Eigen::MatrixXd Hout1;
Vout1.resize(5,2);
Eout1.resize(5,2);

Vout1(1,0) = -csize1/2;
Vout1(1,1) = -rsize1/2;
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
Eigen::MatrixXd V22;
Eigen::MatrixXi F22;

  // Triangulate the interior
        igl::triangle::triangulate(Vout1,Eout1,Hout1,"a1q",V22,F22);
	igl::viewer::Viewer viewer;
//igl::triangle::triangulate(V,E,H,"a5q",V2,F2);
  // Plot the generated mesh
	igl::viewer::Viewer viewer2;
        
	viewer2.data.clear();
	viewer2.data.set_mesh(V22,F22);
	viewer2.core.align_camera_center(V22,F22);
	viewer2.core.show_texture = false;
        viewer2.callback_key_down = &key_down;
	viewer2.launch();

// Wait for Key? 
	
	
        
}
