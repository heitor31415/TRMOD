#include "Tmesh.h"
#include "Math8.h"

CTmesh::CTmesh(unsigned int nNodes, unsigned int nElements, unsigned int nElCS)
// ------------------------------------------------------------------------------------------------------------------ 
// --- constuctor allocates the memory to store the coarse mesh (nodes that will be exchanged) 
// ------------------------------------------------------------------------------------------------------------------ 
{
	// Sizing the class
	nDoF = nNodes;
	nEl = nElements;
	nElCrossSec = nElCS;
	tailElements = nElements;
	tailNodes = nNodes;


	CMatrixStuff_Modal Matrix; // for alocating vectors and matrices (int or double)
	
	/* Allocating memory to store node's coordinates*/
	xn = Matrix.Vector_Allocate(nDoF, 0.0, &allocStatus); // nDoF is the number of nodes/Deg. of Freedom(nDoF)
	yn = Matrix.Vector_Allocate(nDoF, 0.0, &allocStatus);
	zn = Matrix.Vector_Allocate(nDoF, 0.0, &allocStatus);

	/* Allocating memory to store element's global nodes*/
	cMat = Matrix.Matrix_Allocate_Int(nEl, 9, 0, &allocStatus); // nEl is the number of Elements
}
// ------------------------------------------------------------------------------------------------------------------ 
// --- addLayer inserts a layer with layerMat conductivity and 'LayerThickness' on the middle of element 'element') 
// ------------------------------------------------------------------------------------------------------------------ 

void CTmesh::addLayer(unsigned int element, double layerMat, double layerThic)
{
	int oldElement[8];					//Connectivity matrix of the old Element
	double x[8], y[8], z[8];			//Layer nodes position
	for (int k = 0; k < 8; k++)
		oldElement[k] = cMat[element][k];
	//memcpy(oldElement, cMat[element], 8);


	/* Computing coordinates of layer`s nodes
	for (int i = 0; i < 4; i++)			// Loop over each  edge
	{
		double edge[3], norm;
		edge[0] = xn[oldElement[i] - 1] - xn[oldElement[i + 4] - 1]; // Calculating edge x component
		edge[1] = yn[oldElement[i] - 1] - yn[oldElement[i + 4] - 1];
		edge[2] = zn[oldElement[i] - 1] - zn[oldElement[i + 4] - 1];
		norm = compute_norm(edge);									// computing edge norm
		x[i] = (xn[oldElement[i] - 1] + xn[oldElement[i + 4] - 1]) / 2 - (layerThic / 2)*edge[1] / norm; // Layer nodes coordinates
		x[i + 4] = (xn[oldElement[i] - 1] + xn[oldElement[i + 4] - 1]) / 2 + (layerThic / 2)*edge[1] / norm;
		y[i] = (yn[oldElement[i] - 1] + yn[oldElement[i + 4] - 1]) / 2 - (layerThic / 2)*edge[2] / norm;
		y[i + 4] = (yn[oldElement[i] - 1] + yn[oldElement[i + 4] - 1]) / 2 + (layerThic / 2)*edge[2] / norm;
		z[i] = (zn[oldElement[i] - 1] + zn[oldElement[i + 4] - 1]) / 2 - (layerThic / 2)*edge[3] / norm;
		z[i + 4] = (zn[oldElement[i] - 1] + zn[oldElement[i + 4] - 1]) / 2 + (layerThic / 2)*edge[3] / norm;
	}
	*/
	for (int i = 0; i < 8; i++)			// Loop over each  edge
	{
		double edge[3], norm;
		//	using variables to represent the nodes would improve code maintenability
		if (i % 2 == 0)
		{
			edge[0] = xn[oldElement[i] - 1] - xn[oldElement[i + 3] - 1]; // Calculating edge x component
			edge[1] = yn[oldElement[i] - 1] - yn[oldElement[i + 3] - 1];
			edge[2] = zn[oldElement[i] - 1] - zn[oldElement[i + 3] - 1];
			norm = compute_norm(edge);
			x[i] = (xn[oldElement[i] - 1] + xn[oldElement[i + 3] - 1]) / 2 - (layerThic / 2)*edge[1] / norm; // Layer nodes coordinates
			x[i + 3] = (xn[oldElement[i] - 1] + xn[oldElement[i + 3] - 1]) / 2 + (layerThic / 2)*edge[1] / norm;
			y[i] = (yn[oldElement[i] - 1] + yn[oldElement[i + 3] - 1]) / 2 - (layerThic / 2)*edge[2] / norm;
			y[i + 3] = (yn[oldElement[i] - 1] + yn[oldElement[i + 3] - 1]) / 2 + (layerThic / 2)*edge[2] / norm;
			z[i] = (zn[oldElement[i] - 1] + zn[oldElement[i + 3] - 1]) / 2 - (layerThic / 2)*edge[3] / norm;
			z[i + 3] = (zn[oldElement[i] - 1] + zn[oldElement[i + 3] - 1]) / 2 + (layerThic / 2)*edge[3] / norm;
		}
		else
		{
			edge[0] = xn[oldElement[i] - 1] - xn[oldElement[i + 1] - 1]; // Calculating edge x component
			edge[1] = yn[oldElement[i] - 1] - yn[oldElement[i + 1] - 1];
			edge[2] = zn[oldElement[i] - 1] - zn[oldElement[i + 1] - 1];
			norm = compute_norm(edge);
			x[i] = (xn[oldElement[i] - 1] + xn[oldElement[i + 1] - 1]) / 2 - (layerThic / 2)*edge[1] / norm; // Layer nodes coordinates
			x[i + 1] = (xn[oldElement[i] - 1] + xn[oldElement[i + 1] - 1]) / 2 + (layerThic / 2)*edge[1] / norm;
			y[i] = (yn[oldElement[i] - 1] + yn[oldElement[i + 1] - 1]) / 2 - (layerThic / 2)*edge[2] / norm;
			y[i + 1] = (yn[oldElement[i] - 1] + yn[oldElement[i + 1] - 1]) / 2 + (layerThic / 2)*edge[2] / norm;
			z[i] = (zn[oldElement[i] - 1] + zn[oldElement[i + 1] - 1]) / 2 - (layerThic / 2)*edge[3] / norm;
			z[i + 1] = (zn[oldElement[i] - 1] + zn[oldElement[i + 1] - 1]) / 2 + (layerThic / 2)*edge[3] / norm;
			i += 2;
		}
	}
#define GIGA_BYTE (1024 * 1024 * 1024)

	unsigned long long mallocSize = 0, numGigaBytes = 0;
	void *mallocMemory = NULL;

	do
	{
		mallocSize += GIGA_BYTE;
		numGigaBytes = mallocSize / GIGA_BYTE;
		mallocMemory = malloc(mallocSize);
		if (mallocMemory)
		{
			printf("Dynamically allocated %llu GBs\n", numGigaBytes);
			free(mallocMemory);
		}
		else
		{
			printf("Failed to allocate %llu GBs\n", numGigaBytes);
			break;
		}
	} while (true);
	system("PAUSE");
	/* Checking if there is sufficient space on current matrices*/
	if (tailElements <= nEl+2)				 // Global connect. matrix check
		expandConnectMatrix();
	if (tailNodes <= nDoF+8)				 // Global nodes coordiantes vector check
		expandNodesArray();

	for (int i = 0; i < 8; i++)				//Loop over local nodes
		cMat[nEl][7] = nDoF + i+1;			// Global numbering on the layer

	/*Updating Conectivity matrix*/
	for (int i = 0; i < 4; i++)
	{
		cMat[element][i + 4] = cMat[nEl][i];				// Updating element

		cMat[nEl + 1][i] = cMat[nEl][i+4];					// New element, with same material of 'element'
		cMat[nEl + 1][i+4] = oldElement[i + 4];	

	}
	nEl += 2;									// Updating number of elements
	for (int i = 0; i < 8; i++)
	{
		xn[nDoF+i] = x[i];
		yn[nDoF + i] = y[i];
		zn[nDoF + i] = z[i];
		nDoF++;									// Updating number of nodes
	}
	
}

void CTmesh::expandConnectMatrix()
{
	tailElements += 500;						// Optimized for Amortized complexity
	int i;
	int **newMatrix = (int**)realloc(cMat, tailElements * sizeof(int*));
	for (i = 0; i < tailElements; i++)
		newMatrix[i] = (int*)realloc(cMat[i], 9* sizeof(int));
	if (newMatrix != NULL)
		cMat = newMatrix;
	else
	{
		printf("Realloc error");
		system("PAUSE");
	}
}

void CTmesh::expandNodesArray()
{
	tailNodes *= 2;							// Optimized for Amortized complexity
	xn=(double*)realloc(xn, tailNodes*sizeof(double));
	yn=(double*)realloc(yn, tailNodes*sizeof(double));
	zn=(double*)realloc(zn, tailNodes*sizeof(double));
}