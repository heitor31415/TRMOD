#include <iostream>
#include "Math8.h"

using namespace std;

#ifndef ELEMENT_H
#define ELEMENT_H

class CElement : public CMath8
{

public:

	CElement(void);

   void Init                       (double *x,double *y,double *z);
   void Comp_Flux                  (double flux,double *flux_vector,int face);
   void Comp_Convection            (double h,double Tfluid,double *con_matrix,double *con_vector,int face);
   void Comp_Convection_Diagonal   (double h,double Tfluid,double *con_matrix,double *con_vector,int face);
   void Comp_HeatGeneration        (double Q,double *heat_vector,int element);
   void Comp_Conductivity          (double Conductivity,double *con_matrix);

private:
   void   Eval_Jacob_Iso_dA   (const int face,const int numgp);
   void   Eval_Jacob_Iso_dV   (const int numgp);
   void   Eval_GM             (double *G,const int i,const int j);


protected:

   double  J[9],Ji[9];
   double  G1[8],G2[8],G3[8];
   double  g1md[64],g2md[64],g3md[64];
	double  GF[192],GFD[384];
   double  GFTGF[1536];
   double  GA[64];
   double  *xnode,*ynode,*znode;
   double  Jdet,dA;
};

#endif