#include "gsl/gsl_math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_spmatrix.h"
#include "gsl/gsl_splinalg.h"

int ItSolver(gsl_spmatrix *A, double *Load, double *Temp, int size)
{
	int    i, j, result=0;

	//gsl_spmatrix *A = gsl_spmatrix_alloc(size, size); /* triplet format */
	gsl_spmatrix *C;
	gsl_vector *b = gsl_vector_alloc(size); /* right hand side vector */
	gsl_vector *x = gsl_vector_alloc(size); /* solution vector */

	for (i = 0; i<size; i++)
	{
		//for (j = 0; j<size; j++)
		{
			//gsl_spmatrix_set(A, i, j, AMat[i][j]);
		}
		gsl_vector_set(b, i, Load[i]);
	}


	/* convert to compressed column format */
	C = gsl_spmatrix_compcol(A);

	const double tol = 1.0e-12; /* solution relative tolerance */
	const size_t max_iter = 1000; /* maximum iterations */
	const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
	gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, size, 0);
	size_t iter = 0;
	double residual;
	int status;

	/* initial guess x = 0 */
	gsl_vector_set_zero(x);

	/* solve the system C x = b */
	do
	{
		status = gsl_splinalg_itersolve_iterate(C, b, tol, x, work);
		/* print out residual norm ||C*x - b|| */
		residual = gsl_splinalg_itersolve_normr(work);
		//fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);
		if (status == GSL_SUCCESS)
			fprintf(stderr, "Converged\n");
	} while (status == GSL_CONTINUE && ++iter < max_iter);

	/* output solution */
	for (i = 0; i<size; i++)
	{
		Temp[i] = gsl_vector_get(x, i);
	}

	gsl_splinalg_itersolve_free(work);
	//gsl_spmatrix_free(A);
	gsl_spmatrix_free(C);
	gsl_vector_free(x);
	gsl_vector_free(b);


	return result;
}
//
// --- function converts MARCs face numbering to the used Face Numbering
//
int faceNumberConverter(int MarcFaceNumber)
{
	const int marcFace[6] = { 4, 1, 3, 2, 6, 5 };
	return marcFace[MarcFaceNumber];
}

int condenseGauss(int node, double **AMat, double *Load, int size)
{
	int ii, i, j;
	if (Load[node] != 0)
	{
		printf("External node selected");
		return 1;
	}
	for (i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i != j)
				AMat[i][j] -= AMat[i][node] * AMat[node][j] / AMat[node][node];

}

void BinaryOutputTemperature(FILE *fp, const double* temperature, const int counter)
/* ---------------------------------------------------------------------
* generate the binary output of the accelerations at the index
* ---------------------------------------------------------------------*/
{
	char s[512];
	int  value;

	memset(s, ' ', 512 * sizeof(char));
	sprintf(s, "TemperatureInformation\n");
	value = counter;

	fwrite(s, sizeof(char), 512, fp);
	fwrite(&value, sizeof(int), 1, fp);
	fwrite(temperature, sizeof(double), counter, fp);

}

void BinaryOutput(FILE *fp, const double* temperature, const int counter)
/* ---------------------------------------------------------------------
* generate the binary output of the Loadcase LoadcaseName
* ---------------------------------------------------------------------*/
{
	char s[512];
	double SimulationTime = 1.0;
	int NumberOfOutputEntries = 1;
	int BinaryDataEntryCount = 1;
	int BinaryEntryCount = 1;

	memset(s, ' ', 512 * sizeof(char));

	// --- write local header
	sprintf(s, "Dynamics\n");
	fwrite(s, sizeof(char), 512, fp);
	fwrite(&BinaryEntryCount, sizeof(int), 1, fp);
	fwrite(&BinaryDataEntryCount, sizeof(int), 1, fp);
	fwrite(&SimulationTime, sizeof(double), 1, fp);
	fwrite(&NumberOfOutputEntries, sizeof(int), 1, fp);

	// --- write values
	BinaryOutputTemperature(fp, temperature, counter);            // --- write displacement data
}