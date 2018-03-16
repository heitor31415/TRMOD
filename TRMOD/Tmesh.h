#include <iostream>
#include "Math8.h"
#include "MatrixStuff_Modal.h"

using namespace std;
class CTmesh : public CMath8
{
public:
	CTmesh(unsigned int nNodes, unsigned int nElements, unsigned int nElCS);
	//~CTmesh(void);

	void addLayer(unsigned int element, double layerMat, double layerThic);
	void exportNodes(void);
	void exportElements(void);
	void removeDuplicatedNodes(int initialNode);


	unsigned int nDoF, nEl, nElCrossSec;	 // Mesh sizing
	double *xn, *yn, *zn;					 // Nodes global coordinates
	int **cMat;								 // Connectivity matrix
	bool allocStatus;						// Allocation Status

private:
	unsigned int tailNodes, tailElements;
	void expandConnectMatrix();
	void expandNodesArray();
};