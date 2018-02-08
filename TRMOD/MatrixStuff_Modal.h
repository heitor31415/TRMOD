#pragma once

#include <malloc.h>
#include <stdio.h>

class CMatrixStuff_Modal
{
public:
   CMatrixStuff_Modal(void);
   ~CMatrixStuff_Modal(void);

protected:
   void     Matrix_SetValues   (int ir, int ic, double** matrix, double* lineararray);
   void     Matrix_SetValues_T (int ir, int ic, double** matrix, double* lineararray);
   
public:
   void     Matrix_Times_Vector	(int ir, int ic, double* result, double** matrix, double* vector);
   void     M_Times_M_Times_M	(int ir1,int ir2,int ir3, double** result, double **mat1, double** mat2, double **mat3);
   void     Matrix_Deallocate	(int ir, double **p);
   double** Matrix_Allocate		(int ir, int ic, double initfactor, bool *status);
   double*	Vector_Allocate		(int ir, double initfacetor, bool *status);
   int**	Matrix_Allocate_Int	(int ir, int ic, int initfactor, bool *status);
   int*		Vector_Allocate_Int	(int ir, int initfactor, bool *status);
};
