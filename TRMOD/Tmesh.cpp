#include "Tmesh.h"
CTmesh::CTmesh(unsigned int nNodes, unsigned int nElements, unsigned int nElCS)
// ------------------------------------------------------------------------------------------------------------------ 
// --- constuctor allocates the memory to store the coarse mesh (nodes that will be exchanged) 
// ------------------------------------------------------------------------------------------------------------------ 
{
	// Sizing the class
	nDoF = nNodes;
	nEl = nElements;
	nElCrossSec = nElCS;
	tail = nElements;


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
	memcpy(oldElement, cMat[element], 8);
	//if
	
}

void CTmesh::expandConnectMatrix(unsigned int size)
{
	tail *= 2;						// Amortized complexity
	for (int c = 0; c < 9; c++)
		realloc(cMat[c], tail);
}