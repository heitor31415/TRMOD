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

	unsigned int nDoF, nEl, nElCrossSec;
	double *xn, *yn, *zn;			   // Nodes global coordinates
	int **cMat;						  // Connectivity 
	bool allocStatus;

private:
	int tailNodes, tailElements;
	void expandConnectMatrix(unsigned int size);
	void expandNodesArray(unsigned int size);
};