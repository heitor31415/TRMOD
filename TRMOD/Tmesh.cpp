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
	tailElements = 2*nElements;
	tailNodes = 2*nNodes;


	CMatrixStuff_Modal Matrix; // for alocating vectors and matrices (int or double)
	
	/* Allocating memory to store node's coordinates*/
	xn = Matrix.Vector_Allocate(tailNodes, 0.0, &allocStatus); // nDoF is the number of nodes/Deg. of Freedom(nDoF)
	yn = Matrix.Vector_Allocate(tailNodes, 0.0, &allocStatus);
	zn = Matrix.Vector_Allocate(tailNodes, 0.0, &allocStatus);

	/* Allocating memory to store element's global nodes*/
	cMat = Matrix.Matrix_Allocate_Int(tailElements, 9, 0, &allocStatus); // nEl is the number of Elements
}
// ------------------------------------------------------------------------------------------------------------------ 
// --- addLayer inserts a layer with layerMat conductivity and 'LayerThickness' on the middle of element 'element+1') 
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
			x[i] = (xn[oldElement[i] - 1] + xn[oldElement[i + 3] - 1]) / 2 - (layerThic / 2)*edge[0] / norm; // Layer nodes coordinates
			x[i + 3] = (xn[oldElement[i] - 1] + xn[oldElement[i + 3] - 1]) / 2 + (layerThic / 2)*edge[0] / norm;
			y[i] = (yn[oldElement[i] - 1] + yn[oldElement[i + 3] - 1]) / 2 - (layerThic / 2)*edge[1] / norm;
			y[i + 3] = (yn[oldElement[i] - 1] + yn[oldElement[i + 3] - 1]) / 2 + (layerThic / 2)*edge[1] / norm;
			z[i] = (zn[oldElement[i] - 1] + zn[oldElement[i + 3] - 1]) / 2 - (layerThic / 2)*edge[2] / norm;
			z[i + 3] = (zn[oldElement[i] - 1] + zn[oldElement[i + 3] - 1]) / 2 + (layerThic / 2)*edge[2] / norm;
		}
		else
		{
			edge[0] = xn[oldElement[i] - 1] - xn[oldElement[i + 1] - 1]; // Calculating edge x component
			edge[1] = yn[oldElement[i] - 1] - yn[oldElement[i + 1] - 1];
			edge[2] = zn[oldElement[i] - 1] - zn[oldElement[i + 1] - 1];
			norm = compute_norm(edge);
			x[i] = (xn[oldElement[i] - 1] + xn[oldElement[i + 1] - 1]) / 2 - (layerThic / 2)*edge[0] / norm; // Layer nodes coordinates
			x[i + 1] = (xn[oldElement[i] - 1] + xn[oldElement[i + 1] - 1]) / 2 + (layerThic / 2)*edge[0] / norm;
			y[i] = (yn[oldElement[i] - 1] + yn[oldElement[i + 1] - 1]) / 2 - (layerThic / 2)*edge[1] / norm;
			y[i + 1] = (yn[oldElement[i] - 1] + yn[oldElement[i + 1] - 1]) / 2 + (layerThic / 2)*edge[1] / norm;
			z[i] = (zn[oldElement[i] - 1] + zn[oldElement[i + 1] - 1]) / 2 - (layerThic / 2)*edge[2] / norm;
			z[i + 1] = (zn[oldElement[i] - 1] + zn[oldElement[i + 1] - 1]) / 2 + (layerThic / 2)*edge[2] / norm;
			i += 2;
		}
	}

	/* Checking if there is sufficient space on current matrices*/
	if (tailElements <= nEl+2)				 // Global connect. matrix check
		expandConnectMatrix();
	if (tailNodes <= nDoF+8)				 // Global nodes coordiantes vector check
		expandNodesArray();

	for (int i = 0; i < 8; i++)				//Loop over local nodes
		cMat[nEl][i] = nDoF + i+1;			// Global numbering on the layer

	cMat[element][2] = cMat[nEl][1];				// Updating element
	cMat[element][3] = cMat[nEl][0];
	cMat[element][6] = cMat[nEl][5];
	cMat[element][7] = cMat[nEl][4];

	cMat[nEl + 1][0] = cMat[nEl][3];
	cMat[nEl + 1][1] = cMat[nEl][2];
	cMat[nEl + 1][4] = cMat[nEl][7];
	cMat[nEl + 1][5] = cMat[nEl][6];
	cMat[nEl + 1][2] = oldElement[2];
	cMat[nEl + 1][3] = oldElement[3];
	cMat[nEl + 1][6] = oldElement[6];
	cMat[nEl + 1][7] = oldElement[7];
	
	nEl += 2;									// Updating number of elements
	for (int i = 0; i < 8; i++)
	{
		xn[nDoF+i] = x[i];
		yn[nDoF + i] = y[i];
		zn[nDoF + i] = z[i];
	}
	nDoF+=8;									// Updating number of nodes
	
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
	printf("ELEMENTS REALLOC\n");
}

void CTmesh::expandNodesArray()
{
	tailNodes *= 2;							// Optimized for Amortized complexity
	xn=(double*)realloc(xn, tailNodes*sizeof(double));
	yn=(double*)realloc(yn, tailNodes*sizeof(double));
	zn=(double*)realloc(zn, tailNodes*sizeof(double));
	printf("NODES REALLOC\n");
}
void CTmesh::exportNodes(void)
{
	FILE *NODESEXP;

	NODESEXP = fopen("nodes.out", "w+");
	for (int i = 0; i < nDoF; i++)
		fprintf(NODESEXP, "%+12.6e %+12.6e %+12.6e\n", xn[i],yn[i],zn[i]);

	fclose(NODESEXP);
}
void CTmesh::exportElements(void)
{
	FILE *ELEMEXP;

	ELEMEXP = fopen("elements.out", "w+");
	for (int i = 0; i < nEl; i++)
	{ 
		for (int j = 0; j < 8;j++)
			fprintf(ELEMEXP, "%d ", cMat[i][j]);
		fprintf(ELEMEXP, "\n");
	}
	fclose(ELEMEXP);
}