
#include "Math8.h"


CMath8::CMath8(void)
{
}

CMath8::~CMath8(void)
{
}

void CMath8::compute_G1 (double *G, double *Ji, double *g1md, double *g2md, double *g3md)
//---------------------------------------------------------------------
//   multiply Jacobian with shape function derivate for G1 
//---------------------------------------------------------------------
{
   double a11,a12,a13;

   a11 = Ji[0];
   a12 = Ji[1];
   a13 = Ji[2];

   G[0] = a11*g1md[0] + a12*g2md[0] + a13*g3md[0];
   G[1] = a11*g1md[1] + a12*g2md[1] + a13*g3md[1];
   G[2] = a11*g1md[2] + a12*g2md[2] + a13*g3md[2];
   G[3] = a11*g1md[3] + a12*g2md[3] + a13*g3md[3];
   G[4] = a11*g1md[4] + a12*g2md[4] + a13*g3md[4];
   G[5] = a11*g1md[5] + a12*g2md[5] + a13*g3md[5];
   G[6] = a11*g1md[6] + a12*g2md[6] + a13*g3md[6];
   G[7] = a11*g1md[7] + a12*g2md[7] + a13*g3md[7];
}

void CMath8::compute_G2 (double *G, double *Ji, double *g1md, double *g2md, double *g3md)
//---------------------------------------------------------------------
//   multiply Jacobian with shape function derivate for G2
//---------------------------------------------------------------------
{
   double a11,a12,a13;

   a11 = Ji[3];
   a12 = Ji[4];
   a13 = Ji[5];

   G[0] = a11*g1md[0] + a12*g2md[0] + a13*g3md[0];
   G[1] = a11*g1md[1] + a12*g2md[1] + a13*g3md[1];
   G[2] = a11*g1md[2] + a12*g2md[2] + a13*g3md[2];
   G[3] = a11*g1md[3] + a12*g2md[3] + a13*g3md[3];
   G[4] = a11*g1md[4] + a12*g2md[4] + a13*g3md[4];
   G[5] = a11*g1md[5] + a12*g2md[5] + a13*g3md[5];
   G[6] = a11*g1md[6] + a12*g2md[6] + a13*g3md[6];
   G[7] = a11*g1md[7] + a12*g2md[7] + a13*g3md[7];
}

void CMath8::compute_G3 (double *G, double *Ji, double *g1md, double *g2md, double *g3md)
//---------------------------------------------------------------------
//   multiply Jacobian with shape function derivate for G3
//---------------------------------------------------------------------
{
   double a11,a12,a13;

   a11 = Ji[6];
   a12 = Ji[7];
   a13 = Ji[8];

   G[0] = a11*g1md[0] + a12*g2md[0] + a13*g3md[0];
   G[1] = a11*g1md[1] + a12*g2md[1] + a13*g3md[1];
   G[2] = a11*g1md[2] + a12*g2md[2] + a13*g3md[2];
   G[3] = a11*g1md[3] + a12*g2md[3] + a13*g3md[3];
   G[4] = a11*g1md[4] + a12*g2md[4] + a13*g3md[4];
   G[5] = a11*g1md[5] + a12*g2md[5] + a13*g3md[5];
   G[6] = a11*g1md[6] + a12*g2md[6] + a13*g3md[6];
   G[7] = a11*g1md[7] + a12*g2md[7] + a13*g3md[7];
}

void CMath8::compute_Ainvers (double *A, double *Ai, double *Adet)
/* ---------------------------------------------------------------------
* compute invers of matrix A stored in linear array 
*      | 0 1 2 |
*  A = | 3 4 5 |
*      | 6 7 8 |
* ---------------------------------------------------------------------*/
{      
      double t1,t2,t3;

      t1 = A[4]*A[8]-A[5]*A[7];
      t2 = A[5]*A[6]-A[3]*A[8];
      t3 = A[3]*A[7]-A[4]*A[6];

      *Adet = A[0]*t1+A[1]*t2+A[2]*t3;

      Ai[0] =  t1/(*Adet);
      Ai[1] = (A[2]*A[7]-A[1]*A[8])/(*Adet);
      Ai[2] = (A[1]*A[5]-A[2]*A[4])/(*Adet);
      Ai[3] =  t2/(*Adet);
      Ai[4] = (A[0]*A[8]-A[2]*A[6])/(*Adet);
      Ai[5] = (A[2]*A[3]-A[0]*A[5])/(*Adet);
      Ai[6] =  t3/(*Adet);
      Ai[7] = (A[1]*A[6]-A[0]*A[7])/(*Adet);
      Ai[8] = (A[0]*A[4]-A[1]*A[3])/(*Adet);

      for (int i=0; i<9; i++)  if ( fabs(Ai[i]) < 1.0e-12 ) Ai[i] = 0.0;
}

void CMath8::compute_B (double *G, double *Ji, const double *gd1, const double *gd2, const double *gd3) // compute B
{
 
	 double a11,a12,a13,a21,a22,a23,a31,a32,a33;

   a11 = Ji[0];
   a12 = Ji[1];
   a13 = Ji[2];
   a21 = Ji[3];
   a22 = Ji[4];
   a23 = Ji[5];
   a31 = Ji[6];
   a32 = Ji[7];
   a33 = Ji[8];

   G[0] = a11*gd1[0] + a12*gd2[0] + a13*gd3[0];
   G[1] = a11*gd1[1] + a12*gd2[1] + a13*gd3[1];
   G[2] = a11*gd1[2] + a12*gd2[2] + a13*gd3[2];
   G[3] = a11*gd1[3] + a12*gd2[3] + a13*gd3[3];
   G[4] = a11*gd1[4] + a12*gd2[4] + a13*gd3[4];
   G[5] = a11*gd1[5] + a12*gd2[5] + a13*gd3[5];
   G[6] = a11*gd1[6] + a12*gd2[6] + a13*gd3[6];
   G[7] = a11*gd1[7] + a12*gd2[7] + a13*gd3[7];
   G[8] = a21*gd1[0] + a22*gd2[0] + a23*gd3[0];
   G[9] = a21*gd1[1] + a22*gd2[1] + a23*gd3[1];
   G[10] = a21*gd1[2] + a22*gd2[2] + a23*gd3[2];
   G[11] = a21*gd1[3] + a22*gd2[3] + a23*gd3[3];
   G[12] = a21*gd1[4] + a22*gd2[4] + a23*gd3[4];
   G[13] = a21*gd1[5] + a22*gd2[5] + a23*gd3[5];
   G[14] = a21*gd1[6] + a22*gd2[6] + a23*gd3[6];
   G[15] = a21*gd1[7] + a22*gd2[7] + a23*gd3[7];
   G[16] = a31*gd1[0] + a32*gd2[0] + a33*gd3[0];
   G[17] = a31*gd1[1] + a32*gd2[1] + a33*gd3[1];
   G[18] = a31*gd1[2] + a32*gd2[2] + a33*gd3[2];
   G[19] = a31*gd1[3] + a32*gd2[3] + a33*gd3[3];
   G[20] = a31*gd1[4] + a32*gd2[4] + a33*gd3[4];
   G[21] = a31*gd1[5] + a32*gd2[5] + a33*gd3[5];
   G[22] = a31*gd1[6] + a32*gd2[6] + a33*gd3[6];
   G[23] = a31*gd1[7] + a32*gd2[7] + a33*gd3[7];

 }

 void CMath8::compute_T (double *B,double *BT)   //compute Bt
 { 
 
   BT[0] = B[0];
   BT[1] = B[8];
   BT[2] = B[16];
   BT[3] = B[1];
   BT[4] = B[9];
   BT[5] = B[17];
   BT[6] = B[2];
   BT[7] = B[10];
   BT[8] = B[18];
   BT[9] = B[3];
   BT[10] = B[11];
   BT[11] = B[19];
   BT[12] = B[4];
   BT[13] = B[12];
   BT[14] = B[20];
   BT[15] = B[5];
   BT[16] = B[13];
   BT[17] = B[21];
   BT[18] = B[6];
   BT[19] = B[14];
   BT[20] = B[22];
   BT[21] = B[7];
   BT[22] = B[15];
   BT[23] = B[23];

 }

 void CMath8::compute_sp1(double *Bt,int *d,double *sp1) //comput Bt * Dmat
 {
 
  int d11,d12,d13,d21,d22,d23,d31,d32,d33;

   d11 = d[0];
   d12 = d[1];
   d13 = d[2];
   d21 = d[3];
   d22 = d[4];
   d23 = d[5];
   d31 = d[6];
   d32 = d[7];
   d33 = d[8];

  sp1[0] = d11*Bt[0] + d21*Bt[1] + d31*Bt[2];
  sp1[1] = d12*Bt[0] + d22*Bt[1] + d32*Bt[2];
  sp1[2] = d13*Bt[0] + d23*Bt[1] + d33*Bt[2];

  sp1[3] = d11*Bt[3] + d21*Bt[4] + d31*Bt[5];
  sp1[4] = d12*Bt[3] + d22*Bt[4] + d32*Bt[5];
  sp1[5] = d13*Bt[3] + d23*Bt[4] + d33*Bt[5];

  sp1[6] = d11*Bt[6] + d21*Bt[7] + d31*Bt[8];
  sp1[7] = d12*Bt[6] + d22*Bt[7] + d32*Bt[8];
  sp1[8] = d13*Bt[6] + d23*Bt[7] + d33*Bt[8];

  sp1[9] = d11*Bt[9] + d21*Bt[10] + d31*Bt[11];
  sp1[10] = d12*Bt[9] + d22*Bt[10] + d32*Bt[11];
  sp1[11] = d13*Bt[9] + d23*Bt[10] + d33*Bt[11];

  sp1[12] = d11*Bt[12] + d21*Bt[13] + d31*Bt[14];
  sp1[13] = d12*Bt[12] + d22*Bt[13] + d32*Bt[14];
  sp1[14] = d13*Bt[12] + d23*Bt[13] + d33*Bt[14];

  sp1[15] = d11*Bt[15] + d21*Bt[16] + d31*Bt[17];
  sp1[16] = d12*Bt[15] + d22*Bt[16] + d32*Bt[17];
  sp1[17] = d13*Bt[15] + d23*Bt[16] + d33*Bt[17];

  sp1[18] = d11*Bt[18] + d21*Bt[19] + d31*Bt[20];
  sp1[19] = d12*Bt[18] + d22*Bt[19] + d32*Bt[20];
  sp1[20] = d13*Bt[18] + d23*Bt[19] + d33*Bt[20];

  sp1[21] = d11*Bt[21] + d21*Bt[22] + d31*Bt[23];
  sp1[22] = d12*Bt[21] + d22*Bt[22] + d32*Bt[23];
  sp1[23] = d13*Bt[21] + d23*Bt[22] + d33*Bt[23];

  /*cout<<endl;
  cout<<sp1[0]<<" "<<sp1[1]<<" "<<sp1[2]<<endl;
  cout<<sp1[3]<<" "<<sp1[4]<<" "<<sp1[5]<<endl;
  cout<<sp1[6]<<" "<<sp1[7]<<" "<<sp1[8]<<endl;
  cout<<sp1[9]<<" "<<sp1[10]<<" "<<sp1[11]<<endl;
  cout<<sp1[12]<<" "<<sp1[13]<<" "<<sp1[14]<<endl;
  cout<<sp1[15]<<" "<<sp1[16]<<" "<<sp1[17]<<endl;
  cout<<sp1[18]<<" "<<sp1[19]<<" "<<sp1[20]<<endl;
  cout<<sp1[21]<<" "<<sp1[22]<<" "<<sp1[23]<<endl;*/

 }