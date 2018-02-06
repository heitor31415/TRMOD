#pragma once
#include<math.h>

class CMath8 
{
public:
    CMath8(void);
   ~CMath8(void);

// --- functions used for algebraic computation with precalculated shape functions etc.

   void compute_G1      (double *G, double *Ji, double *g1, double *g2, double *g3);
   void compute_G2      (double *G, double *Ji, double *g1, double *g2, double *g3);
   void compute_G3      (double *G, double *Ji, double *g1, double *g2, double *g3);
   void compute_Ainvers (double *A, double *Ai, double *Adet);    
   void compute_B (double *G, double *Ji, const double *gd1, const double *gd2, const double *gd3);
   void compute_T (double *B,double *BT);
   void compute_sp1(double *Bt,int *d,double *sp1);
// --- inline functions used for algebraic computation

   inline double compute_mag (double *A)  // value of norm of a vector
   {
      return (sqrt (A[0]*A[0] + A[1]*A[1] + A[2]*A[2]));
   }

   inline double compute_norm (double *A)  // norm of a vector applied to the vector
   {
      double val;
      val  = sqrt (A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
      A[0] /= val;
      A[1] /= val;
      A[2] /= val;
      return (val);
   }


   inline void compute_norm (double *A, double *B)  // norm of a vector applied to the vector
   {
      double val;
      val  = sqrt (A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
      B[0] = A[0]/val;
      B[1] = A[1]/val;
      B[2] = A[2]/val;
   }

   inline void compute_norm (double *a1 , double *a2, double *a3)  // norm of a vector applied to the vector
   {
      double val  = sqrt ((*a1)*(*a1) + (*a2)*(*a2) + (*a3)*(*a3));
      *a1 /= val;
      *a2 /= val;
      *a3 /= val;
   }

   inline void compute_GIJ (double *GIJ, double *GI, double *GJ) //   dyade multiply GI and GJ in row storage order
   {
      for (int i=0; i<8; i++) for (int j=0; j<8; j++) GIJ[8*i+j] = GI[i]*GJ[j];
   }

   inline double compute_CAtB (double *A, double *B, double a1) //   scalar multiply At and B with a factor 
   {
      double val=0.0;
      for (int i=0; i<8; i++) val += A[i]*B[i];
      return (val*a1);
   }
   inline double compute_CAtB (double *A, double *B) //   scalar multiply At and B without a factor 
   {
      double val=0.0;
      for (int i=0; i<8; i++) val += A[i]*B[i];
      return val;
   }

   inline void compute_CAcrossB (double *C, double *A, double *B) //   cross product of A and B 
   {
      C[0] = A[1]*B[2]-A[2]*B[1];
      C[1] = A[2]*B[0]-A[0]*B[2];
      C[2] = A[0]*B[1]-A[1]*B[0];
   }

   inline double compute_CAtB1 (double *A, double *B) //   multiply At and B   
   {
      double val=0.0;
      for (int i=0; i<8; i++) val += A[i]*B[i];
      return (val);
   }

   inline void compute_CAAt (double *C, double *A, double a1, double a2, double a3) 
   {
      for (int i=0; i<8; i++) for (int j=0; j<8; j++) 
      {
         int k = 8*i+j;
         C[k] = a3*(a1*A[k] + a2*A[8*j+i]);
      }
   }

   inline void compute_DABC (double *D, double *A, double *B, double *C, double a1, double a2, double a3, double a4) 
   {
      for (int i=0; i<64; i++)  D[i] = a4*(a1*A[i] + a2*B[i] + a3*C[i]);
   }

   inline void compute_CAB1 (double *C, double *A, double *B, double a1) 
   {
      for (int i=0; i<8; i++) C[i] = a1*A[i] + B[i];
   }

   inline void compute_CAB2 (double *C, double *A, double *B, double a1, double a2) 
   {
      for (int i=0; i<8; i++) C[i] = a2*(a1*A[i] + B[i]);  
   }

   inline void compute_CAtimesBT (double *C, double *A, double *B) 
   {
      for (int i=0; i<8; i++) for (int j=0; j<8; j++) C[8*i+j] = A[i]*B[j];
   }

   inline double compute_CATtimesB (double *A, double *B, const int num) 
   {
      double c = 0.0; for (int i=0; i<8; i++) c += A[i]*B[8*num+i]; return (c);
   }

   inline double compute_CATtimesBS (double *A, double *B, const int num, const int dir) 
   {
		int i,k;

		double c = 0.0; 
		for (i=0; i<8; i++) 
		{
			k=24*num+8*dir+i; 
			c += A[i]*B[k];
		} 
		return (c);
   }
  inline void compute_sum(double *s,double *LS)   // summation of 8*8 matrix
 {

	 for(int i=0;i<64;i++)
         LS[i] = LS[i] + s[i];	 
 }

  inline void compute_sumVec(double *f,double *Lf)   // summation of 8*1 vector
 {

	 for(int i=0;i<8;i++)
         Lf[i] = Lf[i] + f[i];	 
 }

  inline double ArrToMat(double *arr, int m, int n)             //convert array to 2d array 8*8 matrix
  {
	  return arr[m * 8 + n];
	 /*double mat[8][8];
     for(int i=0;i<8;i++)
		 for(int j=0;j<8;j++)
			 mat[i][j] = arr[(i*8)+j];

	       return mat[m][n];
    */
  }

	inline void compute_ShapeFunctions (const double xi, const double eta, const double psi, double *g)
	{
		g[0] = (1.0-xi)*(1.0-eta)*(1.0-psi)/8.0;
		g[1] = (1.0+xi)*(1.0-eta)*(1.0-psi)/8.0;
		g[2] = (1.0+xi)*(1.0+eta)*(1.0-psi)/8.0;
		g[3] = (1.0-xi)*(1.0+eta)*(1.0-psi)/8.0;
		g[4] = (1.0-xi)*(1.0-eta)*(1.0+psi)/8.0;
		g[5] = (1.0+xi)*(1.0-eta)*(1.0+psi)/8.0;
		g[6] = (1.0+xi)*(1.0+eta)*(1.0+psi)/8.0;
		g[7] = (1.0-xi)*(1.0+eta)*(1.0+psi)/8.0;
	}
	inline void get_NodalPosition (double *xi, double *eta, const int node)
	{
		switch (node)
		{
			case 0 :  {*xi =-1.0; *eta =-1.0;  break;}
			case 1 :  {*xi = 1.0; *eta =-1.0;  break;}
			case 2 :  {*xi = 1.0; *eta = 1.0;  break;}
			case 3 :  {*xi =-1.0; *eta = 1.0;  break;}
		}
	}
};


