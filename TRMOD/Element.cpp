#include "Element.h"
#include <string>

// ------------------------------------------------------------------------------------------------------------------ 
CElement::CElement(void)
// ------------------------------------------------------------------------------------------------------------------ 
// --- constuctor initialize all fixed (precalculated constant) values 
// ------------------------------------------------------------------------------------------------------------------ 
{
   int i,j,k,l;
   int num1,num2;

	// precalculated shape function values for all faces and gaussian points
	// stored in the order
	// g[ 8*numgp + 32*numface]
	// f1: xi  =  1.0, gp1: eta=+gp, psi=+gp, gp2: eta=+gp, psi=-gp, gp3: eta=-gp, psi=+gp, gp4: eta=-gp, psi=-gp 
	// f2: xi  = -1.0, gp1: eta=+gp, psi=+gp, gp2: eta=+gp, psi=-gp, gp3: eta=-gp, psi=+gp, gp4: eta=-gp, psi=-gp  
	// f3: eta =  1.0, gp1: xi =+gp, psi=+gp, gp2: xi =+gp, psi=-gp, gp3: xi =-gp, psi=+gp, gp4: xi =-gp, psi=-gp  
	// f4: eta = -1.0, gp1: xi =+gp, psi=+gp, gp2: xi =+gp, psi=-gp, gp3: xi =-gp, psi=+gp, gp4: xi =-gp, psi=-gp  
	// f5: psi =  1.0, gp1: xi =+gp, eta=+gp, gp2: xi =+gp, eta=-gp, gp3: xi =-gp, eta=+gp, gp4: xi =-gp, eta=-gp  
	// f6: psi = -1.0, gp1: xi =+gp, eta=+gp, gp2: xi =+gp, eta=-gp, gp3: xi =-gp, eta=+gp, gp4: xi =-gp, eta=-gp  
	const double g[192] = {
		0.00000000000000, 0.04465819873852, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.62200846792815, 0.00000000000000,
		0.00000000000000, 0.16666666666667, 0.62200846792815, 0.00000000000000, 0.00000000000000, 0.04465819873852, 0.16666666666667, 0.00000000000000,
		0.00000000000000, 0.16666666666667, 0.04465819873852, 0.00000000000000, 0.00000000000000, 0.62200846792815, 0.16666666666667, 0.00000000000000,
		0.00000000000000, 0.62200846792815, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.04465819873852, 0.00000000000000,
		0.04465819873852, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.62200846792815,
		0.16666666666667, 0.00000000000000, 0.00000000000000, 0.62200846792815, 0.04465819873852, 0.00000000000000, 0.00000000000000, 0.16666666666667,
		0.16666666666667, 0.00000000000000, 0.00000000000000, 0.04465819873852, 0.62200846792815, 0.00000000000000, 0.00000000000000, 0.16666666666667,
		0.62200846792815, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.04465819873852,
		0.00000000000000, 0.00000000000000, 0.16666666666667, 0.04465819873852, 0.00000000000000, 0.00000000000000, 0.62200846792815, 0.16666666666667,
		0.00000000000000, 0.00000000000000, 0.62200846792815, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.04465819873852,
		0.00000000000000, 0.00000000000000, 0.04465819873852, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.62200846792815,
		0.00000000000000, 0.00000000000000, 0.16666666666667, 0.62200846792815, 0.00000000000000, 0.00000000000000, 0.04465819873852, 0.16666666666667,
		0.04465819873852, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.62200846792815, 0.00000000000000, 0.00000000000000,
		0.16666666666667, 0.62200846792815, 0.00000000000000, 0.00000000000000, 0.04465819873852, 0.16666666666667, 0.00000000000000, 0.00000000000000,
		0.16666666666667, 0.04465819873852, 0.00000000000000, 0.00000000000000, 0.62200846792815, 0.16666666666667, 0.00000000000000, 0.00000000000000,
		0.62200846792815, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.04465819873852, 0.00000000000000, 0.00000000000000,
		0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.04465819873852, 0.16666666666667, 0.62200846792815, 0.16666666666667,
		0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.62200846792815, 0.16666666666667, 0.04465819873852,
		0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.16666666666667, 0.04465819873852, 0.16666666666667, 0.62200846792815,
		0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.62200846792815, 0.16666666666667, 0.04465819873852, 0.16666666666667,
		0.04465819873852, 0.16666666666667, 0.62200846792815, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000,
		0.16666666666667, 0.62200846792815, 0.16666666666667, 0.04465819873852, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000,
		0.16666666666667, 0.04465819873852, 0.16666666666667, 0.62200846792815, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000,
		0.62200846792815, 0.16666666666667, 0.04465819873852, 0.16666666666667, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000 };


   // precalculated shape function derivate values for all faces and gaussian points
   // only for computing the tangential directions of the area, two derivates per face
   // G1 derivate to xi, G2 derivate to eta, G3 derivate to psi
   // stored in the order
   // g[ 8*numgp + 64*numface]
   // f1(G2,G3): xi  =  1.0, gp1: eta=+gp, psi=+gp, gp2: eta=+gp, psi=-gp, gp3: eta=-gp, psi=+gp, gp4: eta=-gp, psi=-gp 
   // f2(G2,G3): xi  = -1.0, gp1: eta=+gp, psi=+gp, gp2: eta=+gp, psi=-gp, gp3: eta=-gp, psi=+gp, gp4: eta=-gp, psi=-gp  
   // f3(G1,G3): eta =  1.0, gp1: xi =+gp, psi=+gp, gp2: xi =+gp, psi=-gp, gp3: xi =-gp, psi=+gp, gp4: xi =-gp, psi=-gp  
   // f4(G1,G3): eta = -1.0, gp1: xi =+gp, psi=+gp, gp2: xi =+gp, psi=-gp, gp3: xi =-gp, psi=+gp, gp4: xi =-gp, psi=-gp  
   // f5(G1,G2): psi =  1.0, gp1: xi =+gp, eta=+gp, gp2: xi =+gp, eta=-gp, gp3: xi =-gp, eta=+gp, gp4: xi =-gp, eta=-gp  
   // f6(G1,G2): psi = -1.0, gp1: xi =+gp, eta=+gp, gp2: xi =+gp, eta=-gp, gp3: xi =-gp, eta=+gp, gp4: xi =-gp, eta=-gp 
	const double gd[384] = {
		 0.00000000000000, -0.10566243270259,  0.10566243270259,  0.00000000000000,  0.00000000000000, -0.39433756729741, 0.39433756729741,  0.00000000000000,
		 0.00000000000000, -0.39433756729741,  0.39433756729741,  0.00000000000000,  0.00000000000000, -0.10566243270259, 0.10566243270259,  0.00000000000000,
		 0.00000000000000, -0.10566243270259,  0.10566243270259,  0.00000000000000,  0.00000000000000, -0.39433756729741, 0.39433756729741,  0.00000000000000,
		 0.00000000000000, -0.39433756729741,  0.39433756729741,  0.00000000000000,  0.00000000000000, -0.10566243270259, 0.10566243270259,  0.00000000000000,
		 0.00000000000000, -0.10566243270259, -0.39433756729741,  0.00000000000000,  0.00000000000000,  0.10566243270259, 0.39433756729741,  0.00000000000000,
		 0.00000000000000, -0.10566243270259, -0.39433756729741,  0.00000000000000,  0.00000000000000,  0.10566243270259, 0.39433756729741,  0.00000000000000,
		 0.00000000000000, -0.39433756729741, -0.10566243270259,  0.00000000000000,  0.00000000000000,  0.39433756729741, 0.10566243270259,  0.00000000000000,
		 0.00000000000000, -0.39433756729741, -0.10566243270259,  0.00000000000000,  0.00000000000000,  0.39433756729741, 0.10566243270259,  0.00000000000000,
		-0.10566243270259,  0.00000000000000,  0.00000000000000,  0.10566243270259, -0.39433756729741,  0.00000000000000, 0.00000000000000,  0.39433756729741,
		-0.39433756729741,  0.00000000000000,  0.00000000000000,  0.39433756729741, -0.10566243270259,  0.00000000000000, 0.00000000000000,  0.10566243270259,
		-0.10566243270259,  0.00000000000000,  0.00000000000000,  0.10566243270259, -0.39433756729741,  0.00000000000000, 0.00000000000000,  0.39433756729741,
		-0.39433756729741,  0.00000000000000,  0.00000000000000,  0.39433756729741, -0.10566243270259,  0.00000000000000, 0.00000000000000,  0.10566243270259,
		-0.10566243270259,  0.00000000000000,  0.00000000000000, -0.39433756729741,  0.10566243270259,  0.00000000000000, 0.00000000000000,  0.39433756729741,
		-0.10566243270259,  0.00000000000000,  0.00000000000000, -0.39433756729741,  0.10566243270259,  0.00000000000000, 0.00000000000000,  0.39433756729741,
		-0.39433756729741,  0.00000000000000,  0.00000000000000, -0.10566243270259,  0.39433756729741,  0.00000000000000, 0.00000000000000,  0.10566243270259,
		-0.39433756729741,  0.00000000000000,  0.00000000000000, -0.10566243270259,  0.39433756729741,  0.00000000000000, 0.00000000000000,  0.10566243270259,
		 0.00000000000000,  0.00000000000000,  0.10566243270259, -0.10566243270259,  0.00000000000000,  0.00000000000000, 0.39433756729741, -0.39433756729741,
		 0.00000000000000,  0.00000000000000,  0.39433756729741, -0.39433756729741,  0.00000000000000,  0.00000000000000, 0.10566243270259, -0.10566243270259,
		 0.00000000000000,  0.00000000000000,  0.10566243270259, -0.10566243270259,  0.00000000000000,  0.00000000000000, 0.39433756729741, -0.39433756729741,
		 0.00000000000000,  0.00000000000000,  0.39433756729741, -0.39433756729741,  0.00000000000000,  0.00000000000000, 0.10566243270259, -0.10566243270259,
		 0.00000000000000,  0.00000000000000, -0.39433756729741, -0.10566243270259,  0.00000000000000,  0.00000000000000, 0.39433756729741,  0.10566243270259,
		 0.00000000000000,  0.00000000000000, -0.39433756729741, -0.10566243270259,  0.00000000000000,  0.00000000000000, 0.39433756729741,  0.10566243270259,
		 0.00000000000000,  0.00000000000000, -0.10566243270259, -0.39433756729741,  0.00000000000000,  0.00000000000000, 0.10566243270259,  0.39433756729741,
		 0.00000000000000,  0.00000000000000, -0.10566243270259, -0.39433756729741,  0.00000000000000,  0.00000000000000, 0.10566243270259,  0.39433756729741,
		-0.10566243270259,  0.10566243270259,  0.00000000000000,  0.00000000000000, -0.39433756729741,  0.39433756729741, 0.00000000000000,  0.00000000000000,
		-0.39433756729741,  0.39433756729741,  0.00000000000000,  0.00000000000000, -0.10566243270259,  0.10566243270259, 0.00000000000000,  0.00000000000000,
		-0.10566243270259,  0.10566243270259,  0.00000000000000,  0.00000000000000, -0.39433756729741,  0.39433756729741, 0.00000000000000,  0.00000000000000,
		-0.39433756729741,  0.39433756729741,  0.00000000000000,  0.00000000000000, -0.10566243270259,  0.10566243270259, 0.00000000000000,  0.00000000000000,
		-0.10566243270259, -0.39433756729741,  0.00000000000000,  0.00000000000000,  0.10566243270259,  0.39433756729741, 0.00000000000000,  0.00000000000000,
		-0.10566243270259, -0.39433756729741,  0.00000000000000,  0.00000000000000,  0.10566243270259,  0.39433756729741, 0.00000000000000,  0.00000000000000,
		-0.39433756729741, -0.10566243270259,  0.00000000000000,  0.00000000000000,  0.39433756729741,  0.10566243270259, 0.00000000000000,  0.00000000000000,
		-0.39433756729741, -0.10566243270259,  0.00000000000000,  0.00000000000000,  0.39433756729741,  0.10566243270259, 0.00000000000000,  0.00000000000000,
		 0.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000, -0.10566243270259,  0.10566243270259, 0.39433756729741, -0.39433756729741,
		 0.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000, -0.39433756729741,  0.39433756729741, 0.10566243270259, -0.10566243270259,
		 0.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000, -0.10566243270259,  0.10566243270259, 0.39433756729741, -0.39433756729741,
		 0.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000, -0.39433756729741,  0.39433756729741, 0.10566243270259, -0.10566243270259,
		 0.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000, -0.10566243270259, -0.39433756729741, 0.39433756729741,  0.10566243270259,
		 0.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000, -0.10566243270259, -0.39433756729741, 0.39433756729741,  0.10566243270259,
		 0.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000, -0.39433756729741, -0.10566243270259, 0.10566243270259,  0.39433756729741,
		 0.00000000000000,  0.00000000000000,  0.00000000000000,  0.00000000000000, -0.39433756729741, -0.10566243270259, 0.10566243270259,  0.39433756729741,
		-0.10566243270259,  0.10566243270259,  0.39433756729741, -0.39433756729741,  0.00000000000000,  0.00000000000000, 0.00000000000000,  0.00000000000000,
		-0.39433756729741,  0.39433756729741,  0.10566243270259, -0.10566243270259,  0.00000000000000,  0.00000000000000, 0.00000000000000,  0.00000000000000,
		-0.10566243270259,  0.10566243270259,  0.39433756729741, -0.39433756729741,  0.00000000000000,  0.00000000000000, 0.00000000000000,  0.00000000000000,
		-0.39433756729741,  0.39433756729741,  0.10566243270259, -0.10566243270259,  0.00000000000000,  0.00000000000000, 0.00000000000000,  0.00000000000000,
		-0.10566243270259, -0.39433756729741,  0.39433756729741,  0.10566243270259,  0.00000000000000,  0.00000000000000, 0.00000000000000,  0.00000000000000,
		-0.10566243270259, -0.39433756729741,  0.39433756729741,  0.10566243270259,  0.00000000000000,  0.00000000000000, 0.00000000000000,  0.00000000000000,
		-0.39433756729741, -0.10566243270259,  0.10566243270259,  0.39433756729741,  0.00000000000000,  0.00000000000000, 0.00000000000000,  0.00000000000000,
		-0.39433756729741, -0.10566243270259,  0.10566243270259,  0.39433756729741,  0.00000000000000,  0.00000000000000, 0.00000000000000,  0.00000000000000};

   /* ---------------------------------------------------------------------
   * set shape derivates at gaussain points 1 to 8 for matrix material
   * index operation: gg1[8*i+j] -> i = number of gaussian point,
   *                                j = number of shape function
   * ---------------------------------------------------------------------*/
   const double gg1[64] = {
      -0.02232909936926,0.02232909936926,0.083333333333333,-0.083333333333333,
      -0.083333333333333,0.083333333333333,0.31100423396407,-0.31100423396407,
      -0.083333333333333,0.083333333333333,0.31100423396407,-0.31100423396407,
      -0.02232909936926,0.02232909936926,0.083333333333333,-0.083333333333333,
      -0.083333333333333,0.083333333333333,0.02232909936926,-0.02232909936926,
      -0.31100423396407,0.31100423396407,0.083333333333333,-0.083333333333333,
      -0.31100423396407,0.31100423396407,0.083333333333333,-0.083333333333333,
      -0.083333333333333,0.083333333333333,0.02232909936926,-0.02232909936926,
      -0.02232909936926,0.02232909936926,0.083333333333333,-0.083333333333333,
      -0.083333333333333,0.083333333333333,0.31100423396407,-0.31100423396407,
      -0.083333333333333,0.083333333333333,0.31100423396407,-0.31100423396407,
      -0.02232909936926,0.02232909936926,0.083333333333333,-0.083333333333333,
      -0.083333333333333,0.083333333333333,0.02232909936926,-0.02232909936926,
      -0.31100423396407,0.31100423396407,0.083333333333333,-0.083333333333333,
      -0.31100423396407,0.31100423396407,0.083333333333333,-0.083333333333333,
      -0.083333333333333,0.083333333333333,0.02232909936926,-0.02232909936926};

   const double gg2[64] = {
      -0.02232909936926,-0.083333333333333,0.083333333333333,0.02232909936926,
      -0.083333333333333,-0.31100423396407,0.31100423396407,0.083333333333333,
      -0.083333333333333,-0.31100423396407,0.31100423396407,0.083333333333333,
      -0.02232909936926,-0.083333333333333,0.083333333333333,0.02232909936926,
      -0.02232909936926,-0.083333333333333,0.083333333333333,0.02232909936926,
      -0.083333333333333,-0.31100423396407,0.31100423396407,0.083333333333333,
      -0.083333333333333,-0.31100423396407,0.31100423396407,0.083333333333333,
      -0.02232909936926,-0.083333333333333,0.083333333333333,0.02232909936926,
      -0.083333333333333,-0.02232909936926,0.02232909936926,0.083333333333333,
      -0.31100423396407,-0.083333333333333,0.083333333333333,0.31100423396407,
      -0.31100423396407,-0.083333333333333,0.083333333333333,0.31100423396407,
      -0.083333333333333,-0.02232909936926,0.02232909936926,0.083333333333333,
      -0.083333333333333,-0.02232909936926,0.02232909936926,0.083333333333333,
      -0.31100423396407,-0.083333333333333,0.083333333333333,0.31100423396407,
      -0.31100423396407,-0.083333333333333,0.083333333333333,0.31100423396407,
      -0.083333333333333,-0.02232909936926,0.02232909936926,0.083333333333333};

   const double gg3[64] = {
      -0.02232909936926,-0.083333333333333,-0.31100423396407,-0.083333333333333,
      0.02232909936926,0.083333333333333,0.31100423396407,0.083333333333333,
      -0.02232909936926,-0.083333333333333,-0.31100423396407,-0.083333333333333,
      0.02232909936926,0.083333333333333,0.31100423396407,0.083333333333333,
      -0.083333333333333,-0.31100423396407,-0.083333333333333,-0.02232909936926,
      0.083333333333333,0.31100423396407,0.083333333333333,0.02232909936926,
      -0.083333333333333,-0.31100423396407,-0.083333333333333,-0.02232909936926,
      0.083333333333333,0.31100423396407,0.083333333333333,0.02232909936926,
      -0.083333333333333,-0.02232909936926,-0.083333333333333,-0.31100423396407,
      0.083333333333333,0.02232909936926,0.083333333333333,0.31100423396407,
      -0.083333333333333,-0.02232909936926,-0.083333333333333,-0.31100423396407,
      0.083333333333333,0.02232909936926,0.083333333333333,0.31100423396407,
      -0.31100423396407,-0.083333333333333,-0.02232909936926,-0.083333333333333,
      0.31100423396407,0.083333333333333,0.02232909936926,0.083333333333333,
      -0.31100423396407,-0.083333333333333,-0.02232909936926,-0.083333333333333,
      0.31100423396407,0.083333333333333,0.02232909936926,0.083333333333333};

   /* ---------------------------------------------------------------------
   * set shape function values at gaussain points
   * index operation: gg1[8*i+j] -> i = number of gaussian point,
   *                                j = number of shape function
   * ---------------------------------------------------------------------*/
   const double ga[64]={
      0.0094373878376559,0.035220810900865,0.1314458557658,0.035220810900865,
      0.035220810900865,0.1314458557658,0.49056261216234,0.1314458557658,
      0.035220810900865,0.1314458557658,0.49056261216234,0.1314458557658,
      0.0094373878376559,0.035220810900865,0.1314458557658,0.035220810900865,
      0.035220810900865,0.1314458557658,0.035220810900865,0.0094373878376559,
      0.1314458557658,0.49056261216234,0.1314458557658,0.035220810900865,
      0.1314458557658,0.49056261216234,0.1314458557658,0.035220810900865,
      0.035220810900865,0.1314458557658,0.035220810900865,0.0094373878376559,
      0.035220810900865,0.0094373878376559,0.035220810900865,0.1314458557658,
      0.1314458557658,0.035220810900865,0.1314458557658,0.49056261216234,
      0.1314458557658,0.035220810900865,0.1314458557658,0.49056261216234,
      0.035220810900865,0.0094373878376559,0.035220810900865,0.1314458557658,
      0.1314458557658,0.035220810900865,0.0094373878376559,0.035220810900865,
      0.49056261216234,0.1314458557658,0.035220810900865,0.1314458557658,
      0.49056261216234,0.1314458557658,0.035220810900865,0.1314458557658,
      0.1314458557658,0.035220810900865,0.0094373878376559,0.035220810900865};


   // --- copy constants to arrays
	memcpy(GF , g , 192 * sizeof(double));
	memcpy(GFD, gd, 384 * sizeof(double));
   memcpy(g1md,gg1, 64 * sizeof(double));
   memcpy(g2md,gg2, 64 * sizeof(double));
   memcpy(g3md,gg3, 64 * sizeof(double));
   memcpy(GA , ga,  64 * sizeof(double));

   // --- generate GF times GT'
   // --- loop over 6 faces
   for (i=0; i<6; i++)
   {
      // --- loop over gaussian points
      for (j=0; j<4; j++)
      {
         num1 = i*256+j*64;
         num2 = i*32 +j*8;
         for (k=0; k<8; k++)
         {
            for (l=0; l<8; l++)
            {
               GFTGF[num1 + k*8 + l] = GF[num2 + k]*GF[num2 + l];
            }
         }
      }
   }
}
// ------------------------------------------------------------------------------------------------------------------ 
void CElement::Init(double *x,double *y,double *z)
// ------------------------------------------------------------------------------------------------------------------ 
// --- store the geometry
// ------------------------------------------------------------------------------------------------------------------ 
{
   // --- store x,y,z coordinate pointers
   xnode = x;
   ynode = y;
   znode = z;
}
// ------------------------------------------------------------------------------------------------------------------ 
void CElement::Eval_Jacob_Iso_dV(const int numgp)
// ------------------------------------------------------------------------------------------------------------------ 
// --- compute jacobian at point xi,eta,psi and store it as matrix for isopar
// ---      | 0 1 2 |
// ---  J = | 3 4 5 |
// ---      | 6 7 8 |
// ------------------------------------------------------------------------------------------------------------------
{
   J[0] = compute_CATtimesB(xnode,g1md,numgp);
   J[1] = compute_CATtimesB(ynode,g1md,numgp);
   J[2] = compute_CATtimesB(znode,g1md,numgp);

   J[3] = compute_CATtimesB(xnode,g2md,numgp);
   J[4] = compute_CATtimesB(ynode,g2md,numgp);
   J[5] = compute_CATtimesB(znode,g2md,numgp);

   J[6] = compute_CATtimesB(xnode,g3md,numgp);
   J[7] = compute_CATtimesB(ynode,g3md,numgp);
   J[8] = compute_CATtimesB(znode,g3md,numgp);

   compute_Ainvers(J,Ji,&Jdet);
}
// ------------------------------------------------------------------------------------------------------------------ 
void CElement::Eval_Jacob_Iso_dA(const int face, const int numgp)
// ------------------------------------------------------------------------------------------------------------------ 
// --- function calculate the surface area for a gaussian point on a certain face
// --- face 1: xi=+1.0, face 2: xi=-1.0, face 3: eta=+1.0, face 4: eta=-1.0, face 5: psi=+1.0, face 6: psi=-1.0
// --- order of gaussian points look to the class constructer
// ------------------------------------------------------------------------------------------------------------------ 
{
	 int    i,num;
	 double t1[3];
	 double t2[3];
	 double  n[3];

    // --- init vectors
	 memset(t1, 0, 3 * sizeof(double));
	 memset(t2, 0, 3 * sizeof(double));

	 // --- comp tangential vectors, using structured storing of percalculated values of GFD
    // --- f1: t1 use G2, t2 use G3; f2: t1 use G2, t2 use G3; f3: t1 use G1, t2 use G3
    // --- f4: t1 use G1, t2 use G3; f5: t1 use G1, t2 use G2; f6: t1 use G1, t2 use G2
    num = numgp*8 + face*64;
	 for (i=0; i<8; i++)
	 {
       
		 t1[0] += xnode[i] * GFD[num + i];
		 t1[1] += ynode[i] * GFD[num + i];
		 t1[2] += znode[i] * GFD[num + i];

		 t2[0] += xnode[i] * GFD[num + i + 32];
		 t2[1] += ynode[i] * GFD[num + i + 32];
		 t2[2] += znode[i] * GFD[num + i + 32];
	 }

    // --- get normal vector 
	 n[0] = t1[2]*t2[1] - t1[1]*t2[2];
	 n[1] = t1[0]*t2[2] - t1[2]*t2[0];
	 n[2] = t1[1]*t2[0] - t1[0]*t2[1];

    // --- comp area of gp and return value
	 dA = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
 }
// ------------------------------------------------------------------------------------------------------------------ 
void CElement::Comp_Flux(double flux, double *flux_vector, int face)
// ------------------------------------------------------------------------------------------------------------------ 
// --- function calculate the flux_vector (load vector of flux) for one face denoted by the number "face"
// --- face 1: xi=+1.0, face 2: xi=-1.0, face 3: eta=+1.0, face 4: eta=-1.0, face 5: psi=+1.0, face 6: psi=-1.0
// ------------------------------------------------------------------------------------------------------------------ 
{
	 int    i,j,num;
    double factor;
	memset(flux_vector, 0, 8 * sizeof(double));
	 // --- face zero based, no face return function
	 if (face == 0) return;
	 face = face - 1;

	 // --- loop over 4 gaussian points 
	 for (i=0; i<4; i++)
	 {
		 // --- get dA and add constant flux value
		 Eval_Jacob_Iso_dA(face,i);
		 factor = dA*flux;
		 // --- add dA times g to the fluxvector, GF contains the shape function values
       num = face*32 + i*8;
       for (j=0; j<8; j++) flux_vector[j] += factor*GF[num + j];
	 }
}
// ------------------------------------------------------------------------------------------------------------------ 
void CElement::Comp_Convection(double h,double Tfluid,double *con_matrix,double *con_vector,int face)
// ------------------------------------------------------------------------------------------------------------------ 
// --- function calculate the con_vector (load vector of convection) for one face denoted by the number "face"
// --- function calculate the con_matrix (temp matrix of convection) for one face denoted by the number "face"
// --- face 1: xi=+1.0, face 2: xi=-1.0, face 3: eta=+1.0, face 4: eta=-1.0, face 5: psi=+1.0, face 6: psi=-1.0
// ------------------------------------------------------------------------------------------------------------------ 
{
   int    i,j,num1,num2;
   double factor1,factor2;

   // --- face zero based, no face return function
   if (face == 0) return;
   face = face - 1;

   // --- loop over 4 gaussian points 
   for (i=0; i<4; i++)
   {
      // --- get dA and add constant flux value
      Eval_Jacob_Iso_dA(face,i);

      // --- comp constant values 
      factor1 = dA*h*Tfluid;
      factor2 = dA*h;

      // --- sum up the values and generate load vector and convection matrix
      num1 = face*32  + i*8;
      num2 = face*256 + i*64;
      for (j = 0; j<8;  j++) con_vector[j] += GF   [num1+j]*factor1;
      for (j = 0; j<64; j++) con_matrix[j] += GFTGF[num2+j]*factor2;
   }
}
// ------------------------------------------------------------------------------------------------------------------ 
void CElement::Comp_Convection_Diagonal(double h,double Tfluid,double *con_matrix,double *con_vector,int face)
// ------------------------------------------------------------------------------------------------------------------ 
// --- function calculate the con_vector (load vector of convection) for one face denoted by the number "face"
// --- function calculate the con_matrix (temp matrix of convection in diagonalized way) for one face denoted by the number "face"
// --- face 1: xi=+1.0, face 2: xi=-1.0, face 3: eta=+1.0, face 4: eta=-1.0, face 5: psi=+1.0, face 6: psi=-1.0
// ------------------------------------------------------------------------------------------------------------------ 
{
   int    i,j,k,num1,num2;
   double factor1,factor2;
   memset(con_vector, 0, 8 * sizeof(double));
   memset(con_matrix, 0, 64 * sizeof(double));

   // --- face zero based, no face return function
   if(face == 0) return;
   face = face - 1;

   // --- loop over 4 gaussian points 
   for(i=0; i<4; i++)
   {
      // --- get dA and add constant flux value
      Eval_Jacob_Iso_dA(face,i);

      // --- comp constant values 
      factor1 = dA*h*Tfluid;
      factor2 = dA*h;
      // --- sum up the values and generate load vector and convection matrix
      num1 = face*32  + i*8;
      num2 = face*256 + i*64;
      for (j=0; j<8; j++)
      {
         con_vector[j] += GF[num1+j]*factor1;
         for (k=0; k<8; k++)
			 con_matrix[j*9] += GFTGF[num2+j*8+k]*factor2;
      }
   }
}
// ------------------------------------------------------------------------------------------------------------------ 
void CElement::Comp_HeatGeneration(double Q, double *heat_vector, int element) 
// ------------------------------------------------------------------------------------------------------------------ 
// --- function calculate the heat_vector (load vector of heat) for whole element
// ------------------------------------------------------------------------------------------------------------------ 
{
   int    i,j;
   double factor;

   // --- loop over 8 gaussian points (volume)
   for (i=0; i<8; i++)
   {
      // --- get Jacobian, inverted Jacobian and dV
      Eval_Jacob_Iso_dV(i);

      // --- add constant heat value to Jdet
      factor = Q*Jdet*element;

      // --- add values to heat vector
      for (j=0; j<8; j++) heat_vector[j] += GA[8*i+j]*factor;
   }
}
// ------------------------------------------------------------------------------------------------------------------
void CElement::Comp_Conductivity(double Conductivity,double *con_matrix)
// ------------------------------------------------------------------------------------------------------------------ 
// --- function calculate the conductivity_matrix for whole element
// --- only for same conductivity to each three directions
// --- conductivity matrix is stored in the order S[i][j] of i*8 + j 
// ------------------------------------------------------------------------------------------------------------------
{
   int    i,j,k;
   double factor;
   // --- loop over 8 gaussian points (volume)
   for (i=0; i<8; i++)
   {
      // --- get Jacobian, inverted Jacobian and dV
      Eval_Jacob_Iso_dV(i);

      // --- comp shape function derivates with inv. Jacobian
      Eval_GM(G1,1,i);
      Eval_GM(G2,2,i);
      Eval_GM(G3,3,i);

      // --- B*D*B' (D only diagonal elements, same conductivity for all directions)
      for (j=0; j<8; j++)
      {
         for (k=0; k<8; k++)
         {
            // --- sparse computation because of diagonal filled matrix, B.D.B' * Jdet
            con_matrix[j*8 + k] += G1[j]*G1[k]*Conductivity*Jdet; // G1.D(1,1).G1' * Jdet
            con_matrix[j*8 + k] += G2[j]*G2[k]*Conductivity*Jdet; // G2.D(2,2).G2' * Jdet
            con_matrix[j*8 + k] += G3[j]*G3[k]*Conductivity*Jdet; // G3.D(3,3).G3' * Jdet
         }
      }
   }
}
// ------------------------------------------------------------------------------------------------------------------ 
void CElement::Eval_GM(double *G,const int i,const int j)
// ------------------------------------------------------------------------------------------------------------------ 
//   multiply invers of Jacobian with shape function gj derivate for Gi
//   i can be 1,2,3 for x,y or z and j denotes the gaussian point
// ------------------------------------------------------------------------------------------------------------------ 
{
   double a11,a12,a13;

   switch(i)
   {
   case 1: a11 = Ji[0]; a12 = Ji[1]; a13 = Ji[2]; break;
   case 2: a11 = Ji[3]; a12 = Ji[4]; a13 = Ji[5]; break;
   case 3: a11 = Ji[6]; a12 = Ji[7]; a13 = Ji[8]; break;
   }

   switch(j)
   {
   case 0:
      G[0] = a11*g1md[0] + a12*g2md[0] + a13*g3md[0];
      G[1] = a11*g1md[1] + a12*g2md[1] + a13*g3md[1];
      G[2] = a11*g1md[2] + a12*g2md[2] + a13*g3md[2];
      G[3] = a11*g1md[3] + a12*g2md[3] + a13*g3md[3];
      G[4] = a11*g1md[4] + a12*g2md[4] + a13*g3md[4];
      G[5] = a11*g1md[5] + a12*g2md[5] + a13*g3md[5];
      G[6] = a11*g1md[6] + a12*g2md[6] + a13*g3md[6];
      G[7] = a11*g1md[7] + a12*g2md[7] + a13*g3md[7];
      break;
   case 1:
      G[0] = a11*g1md[8] + a12*g2md[8] + a13*g3md[8];
      G[1] = a11*g1md[9] + a12*g2md[9] + a13*g3md[9];
      G[2] = a11*g1md[10] + a12*g2md[10] + a13*g3md[10];
      G[3] = a11*g1md[11] + a12*g2md[11] + a13*g3md[11];
      G[4] = a11*g1md[12] + a12*g2md[12] + a13*g3md[12];
      G[5] = a11*g1md[13] + a12*g2md[13] + a13*g3md[13];
      G[6] = a11*g1md[14] + a12*g2md[14] + a13*g3md[14];
      G[7] = a11*g1md[15] + a12*g2md[15] + a13*g3md[15];
      break;
   case 2:
      G[0] = a11*g1md[16] + a12*g2md[16] + a13*g3md[16];
      G[1] = a11*g1md[17] + a12*g2md[17] + a13*g3md[17];
      G[2] = a11*g1md[18] + a12*g2md[18] + a13*g3md[18];
      G[3] = a11*g1md[19] + a12*g2md[19] + a13*g3md[19];
      G[4] = a11*g1md[20] + a12*g2md[20] + a13*g3md[20];
      G[5] = a11*g1md[21] + a12*g2md[21] + a13*g3md[21];
      G[6] = a11*g1md[22] + a12*g2md[22] + a13*g3md[22];
      G[7] = a11*g1md[23] + a12*g2md[23] + a13*g3md[23];
      break;
   case 3:
      G[0] = a11*g1md[24] + a12*g2md[24] + a13*g3md[24];
      G[1] = a11*g1md[25] + a12*g2md[25] + a13*g3md[25];
      G[2] = a11*g1md[26] + a12*g2md[26] + a13*g3md[26];
      G[3] = a11*g1md[27] + a12*g2md[27] + a13*g3md[27];
      G[4] = a11*g1md[28] + a12*g2md[28] + a13*g3md[28];
      G[5] = a11*g1md[29] + a12*g2md[29] + a13*g3md[29];
      G[6] = a11*g1md[30] + a12*g2md[30] + a13*g3md[30];
      G[7] = a11*g1md[31] + a12*g2md[31] + a13*g3md[31];
      break;
   case 4:
      G[0] = a11*g1md[32] + a12*g2md[32] + a13*g3md[32];
      G[1] = a11*g1md[33] + a12*g2md[33] + a13*g3md[33];
      G[2] = a11*g1md[34] + a12*g2md[34] + a13*g3md[34];
      G[3] = a11*g1md[35] + a12*g2md[35] + a13*g3md[35];
      G[4] = a11*g1md[36] + a12*g2md[36] + a13*g3md[36];
      G[5] = a11*g1md[37] + a12*g2md[37] + a13*g3md[37];
      G[6] = a11*g1md[38] + a12*g2md[38] + a13*g3md[38];
      G[7] = a11*g1md[39] + a12*g2md[39] + a13*g3md[39];
      break;
   case 5:
      G[0] = a11*g1md[40] + a12*g2md[40] + a13*g3md[40];
      G[1] = a11*g1md[41] + a12*g2md[41] + a13*g3md[41];
      G[2] = a11*g1md[42] + a12*g2md[42] + a13*g3md[42];
      G[3] = a11*g1md[43] + a12*g2md[43] + a13*g3md[43];
      G[4] = a11*g1md[44] + a12*g2md[44] + a13*g3md[44];
      G[5] = a11*g1md[45] + a12*g2md[45] + a13*g3md[45];
      G[6] = a11*g1md[46] + a12*g2md[46] + a13*g3md[46];
      G[7] = a11*g1md[47] + a12*g2md[47] + a13*g3md[47];
      break;
   case 6:
      G[0] = a11*g1md[48] + a12*g2md[48] + a13*g3md[48];
      G[1] = a11*g1md[49] + a12*g2md[49] + a13*g3md[49];
      G[2] = a11*g1md[50] + a12*g2md[50] + a13*g3md[50];
      G[3] = a11*g1md[51] + a12*g2md[51] + a13*g3md[51];
      G[4] = a11*g1md[52] + a12*g2md[52] + a13*g3md[52];
      G[5] = a11*g1md[53] + a12*g2md[53] + a13*g3md[53];
      G[6] = a11*g1md[54] + a12*g2md[54] + a13*g3md[54];
      G[7] = a11*g1md[55] + a12*g2md[55] + a13*g3md[55];
      break;
   case 7:
      G[0] = a11*g1md[56] + a12*g2md[56] + a13*g3md[56];
      G[1] = a11*g1md[57] + a12*g2md[57] + a13*g3md[57];
      G[2] = a11*g1md[58] + a12*g2md[58] + a13*g3md[58];
      G[3] = a11*g1md[59] + a12*g2md[59] + a13*g3md[59];
      G[4] = a11*g1md[60] + a12*g2md[60] + a13*g3md[60];
      G[5] = a11*g1md[61] + a12*g2md[61] + a13*g3md[61];
      G[6] = a11*g1md[62] + a12*g2md[62] + a13*g3md[62];
      G[7] = a11*g1md[63] + a12*g2md[63] + a13*g3md[63];
      break;
   }
}