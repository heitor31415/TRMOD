
#include "MatrixStuff_Modal.h"

CMatrixStuff_Modal::CMatrixStuff_Modal(void)
{
}

CMatrixStuff_Modal::~CMatrixStuff_Modal(void)
{
}
//-----------------------------------------------------------------------------------------------------
double** CMatrixStuff_Modal::Matrix_Allocate(int ir, int ic, double initfactor, bool *status)
//-----------------------------------------------------------------------------------------------------
// --- allocate space for irow times icolumn matrix and return pointer
//-----------------------------------------------------------------------------------------------------
{
   int    i,j;
   double** p;

   *status = true;
   p = ( double ** ) malloc (ir*sizeof(double *));

   if ( p == NULL ) *status = false;

   for (i=0; i<ir; i++ )
   {
      p[i] = (double*) malloc (ic*sizeof(double));

      if ( p[i] == NULL ) *status = false;
   }

   for ( i=0; i<ir; i++ )
   {
      for ( j=0; j<ic; j++ )
      {
         p[i][j] = (double) initfactor;
      }
   }

   return (p);
}
//-----------------------------------------------------------------------------------------------------
int** CMatrixStuff_Modal::Matrix_Allocate_Int(int ir, int ic, int initfactor, bool *status)
//-----------------------------------------------------------------------------------------------------
// --- allocate space for irow times icolumn matrix and return pointer
//-----------------------------------------------------------------------------------------------------
{
	int    i, j;
	int** p;

	*status = true;
	p = (int **)malloc(ir*sizeof(int *));

	if (p == NULL) *status = false;

	for (i = 0; i<ir; i++)
	{
		p[i] = (int *)malloc(ic*sizeof(int));

		if (p[i] == NULL) *status = false;
	}

	for (i = 0; i<ir; i++)
	{
		for (j = 0; j<ic; j++)
		{
			p[i][j] = (int)initfactor;
		}
	}

	return (p);
}
//-----------------------------------------------------------------------------------------------------
double* CMatrixStuff_Modal::Vector_Allocate(int ir, double initfactor, bool *status)
//-----------------------------------------------------------------------------------------------------
// --- allocate space for irow times icolumn matrix and return pointer
//-----------------------------------------------------------------------------------------------------
{
	int    i;
	double* p;

	*status = true;
	p = (double *)malloc(ir*sizeof(double));

	if (p == NULL) *status = false;

	for (i = 0; i<ir; i++)
	{
			p[i] = (double)initfactor;
	}

	return (p);
}
//-----------------------------------------------------------------------------------------------------
int* CMatrixStuff_Modal::Vector_Allocate_Int(int ir, int initfactor, bool *status)
//-----------------------------------------------------------------------------------------------------
// --- allocate space for irow times icolumn matrix and return pointer
//-----------------------------------------------------------------------------------------------------
{
	int    i;
	int* p;

	*status = true;
	p = (int *)malloc(ir*sizeof(int));

	if (p == NULL) *status = false;

	for (i = 0; i<ir; i++)
	{
		p[i] = (int)initfactor;
	}

	return (p);
}
//-----------------------------------------------------------------------------------------------------
void CMatrixStuff_Modal::Matrix_SetValues(int ir, int ic, double** matrix, double* lineararray)
// ----------------------------------------------------------------------------------------------------
// --- store values of an array into matrix
// ----------------------------------------------------------------------------------------------------
{
   int i,j;

   for (i=0; i<ir; i++) for (j=0; j<ic; j++) matrix[i][j] = lineararray[j*ir + i];
}
// ----------------------------------------------------------------------------------------------------
void CMatrixStuff_Modal::Matrix_SetValues_T(int ir, int ic, double** matrix, double* lineararray)
// ----------------------------------------------------------------------------------------------------
// --- store values of an array into matrix(Transposed)
// ----------------------------------------------------------------------------------------------------
{
   int i,j;

   for (i=0; i<ir; i++) for (j=0; j<ic; j++) matrix[j][i] = lineararray[j*ir + i];
}
// ----------------------------------------------------------------------------------------------------
void CMatrixStuff_Modal::Matrix_Deallocate(int ir, double **p)
// ----------------------------------------------------------------------------------------------------
// --- deallocate space for irow times icolumn matrix and return pointer
// ----------------------------------------------------------------------------------------------------
{
   for (int i=0; i<ir; i++)
   {
      free(p[i]);
   }
   free(p);
}
// ----------------------------------------------------------------------------------------------------
void CMatrixStuff_Modal::Matrix_Times_Vector(int ir, int ic, double* result, double** matrix, double* vector)
// ----------------------------------------------------------------------------------------------------
// --- deallocate space for irow times icolumn matrix and return pointer
// ----------------------------------------------------------------------------------------------------
{
   int    i,j;
   double sum;

   for (i=0; i<ir; i++)
   {
      sum = 0.0;

      // --- sum up delta values
      for (j=0; j<ic; j++) sum += matrix[i][j]*vector[j];

      result[i] = sum;
   }
}
// ----------------------------------------------------------------------------------------------------
void CMatrixStuff_Modal::M_Times_M_Times_M(int ir1,int ir2,int ir3, double** result, double **mat1, double** mat2, double **mat3)
// ----------------------------------------------------------------------------------------------------
// --- deallocate space for irow times icolumn matrix and return pointer
// ----------------------------------------------------------------------------------------------------
{
   int    i,j,k;
   bool   status;
   double sum;
   double **HMat;

   // --- allocate help matrix
   HMat = Matrix_Allocate(ir2,ir3,0.0,&status);

   // --- compute first matrix HMat = mat2*mat3
   for (i=0; i<ir2; i++)
   {
      for (j=0; j<ir3; j++)
      {
         sum = 0.0;
         for (k=0; k<ir2; k++)
         {
            sum += mat2[i][k]*mat3[k][j];
         }
         HMat[i][j] = sum;
      }
   }
   
   // --- compute second matrix result = mat1*HMat
   for (i=0; i<ir1; i++)
   {
      for (j=0; j<ir3; j++)
      {
         sum = 0.0;
         for (k=0; k<ir2; k++)
         {
            sum += mat1[i][k]*HMat[k][j];
         }
         result[i][j] = sum;
      }
   }

   Matrix_Deallocate(ir2,HMat);
}