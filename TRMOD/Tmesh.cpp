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
	int oldElement[8];
	double x[8], y[8], z[8];			//New element nodes position
	memcpy(oldElement, cMat[element], 8);

	if (tailElements == nEl)
		expandConnectMatrix(nEl);
	if (tailNodes == nDoF)
		expandNodesArray(nDoF);
	for (int i = 0; i < 4; i++)			// Loop over each local node
	{
		double edge[3];
		edge[1] = xn[oldElement[i] - 1] - xn[oldElement[i + 4] - 1];
		edge[2] = yn[oldElement[i] - 1] - yn[oldElement[i + 4] - 1];
		edge[3] = zn[oldElement[i] - 1] - zn[oldElement[i + 4] - 1];
		double norm;
		norm = compute_norm(edge);
		
		x[i] = xn[oldElement[i] - 1] + layerThic*edge[1];
		x[i+4]=
	}
	
	nEl += 2;							// 1 element turns in 3 elements
	nDoF += 8;							// 1 new element = 8 nodes

	
}

void CTmesh::expandConnectMatrix(unsigned int size)
{
	tailElements *= 2;						// Optimized for Amortized complexity
	for (int c = 0; c < 9; c++)
		realloc(cMat[c], tailElements);
}

void CTmesh::expandNodesArray(unsigned int size)
{
	tailNodes *= 2;							// Optimized for Amortized complexity
	realloc(xn, tailNodes);
	realloc(yn, tailNodes);
	realloc(zn, tailNodes);
}