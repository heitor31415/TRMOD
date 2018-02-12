#include <iostream>
#include "Element.h"
#include "misc.h"
#include "MatrixStuff_Modal.h"
#include <string>
#include "Tmesh.h"

using namespace std;
int main(){
	/*  MESH SIZING AND CONNECTIVITY MATRIX */
	int nDoF, nEl, nElCrossSec; // Number of degrees of Freedom(= num.Nodes in thermal analysis) and number of elements

	/* NODES POSITIONS AND GENERAL VARIABLES */
	double matConductivity = 35, ambTemp; // Material property and Ambient Temperature (used in convection)

	/* NATURAL/NEUMANN BOUNDARY CONDITIONS VARIABLES */
	double *heatGenList, **fluxList, **convList; // Vectors that stores B.C.s information for each element

	/* ESSENTIAL/DIRICHLET BOUNDARY CONDITIONS VARIABLES */
	int *fixedNodesList, *freeNodesList; // Lists with nodes with fixed temp. (Dirichlet BC) and 'free' nodes repect. 
	double *tempFixedNodes; // List with the fixed nodes' temperatures.
	int nFixedNodes, nFreeNodes; // Number of nodes with natural/dirichlet Boundary conditions and those without it, respect.

	/* STIFFNESS MATRIXES AND FORCE VECTORS */
	double **SG, **SSol; // SG is the Global Stiffness matrix and SSol is SG after applied the Boundary Conditions
	double *FG, *FSol; // FG is the Global Force Vector
	double *tempSol, *temp; // tempSol vector contains unkown temperatures and temp contain all nodal temperatures 
	// Therefore SG is nDoFxnDof and SSol is nFreeNodesxnFreeNodes.


	char s[256]; // We will use this 'string' to read the input
	bool allocStatus;

	CMatrixStuff_Modal Matrix; // for alocating vectors and matrices (int or double)


	/* ####################################################################################
	########################## READING .dat INPUT FILE #################################
	#################################################################################### */
	FILE *INPUT;
	//INPUT = fopen("modifiedDATAFlux.dat", "r"); // "fileName.dat"
	//INPUT = fopen("modifiedDATAFlux.dat", "r");
	INPUT = fopen("TireNEG.dat", "r");
	fgets(s, 256, INPUT);// >>Text line= MESH SIZING

	/* Read the sizing values, which are Number of Elements and number of nodes. */
	fgets(s, 256, INPUT);
	sscanf(s, "%d %d %d", &nDoF, &nEl, &nElCrossSec);	// Number of elements & nodes on the model, resp.
	CTmesh cMesh(nDoF, nEl, nElCrossSec);				// coarse/common mesh


	/* Store nodes' global coordinates >>(NODES HAVE TO BE SORTED)<<*/
	int globalNodeNumber, elementType, elementNumber;	// Not used in the current version of the code
	fgets(s, 256, INPUT); // >>Text line= GLOBAL NODES COORDINATES

	int layer, fixedNodesCounter = 0;					// Layer represents which layer(s) the node belongs
	bool fixedNodesStart = true;						// Used to identify the first fixed node
	/* Nodes coordinates input */
	for (int i = 0; i < nDoF; i++){ //Loop over all nodes
		fgets(s, 256, INPUT);
		sscanf(s, "%d %lf %lf %lf %d", &globalNodeNumber, &cMesh.xn[i], &cMesh.yn[i], &cMesh.zn[i], &layer); // Scan each node coordinate
		if (layer < 0){
			if (fixedNodesStart){ // just happens one time and allocate 
				fixedNodesStart = false;
				nFreeNodes = i;
				nFixedNodes = nDoF - nFreeNodes;
				fixedNodesList = Matrix.Vector_Allocate_Int(nFixedNodes, 0, &allocStatus); // List with 'fixed' nodes
				tempFixedNodes = Matrix.Vector_Allocate(nFixedNodes, 0.0, &allocStatus); // List with each fixed node temperature
				freeNodesList = Matrix.Vector_Allocate_Int(nFreeNodes, 0, &allocStatus); // List with free nodes
			}
			fixedNodesList[fixedNodesCounter] = globalNodeNumber;
			tempFixedNodes[fixedNodesCounter] = 120.0;						// A map can be used to fixed nodes*
			fixedNodesCounter++;
		}
	}

	fgets(s, 256, INPUT); // >>Text Line= CONNECTIVITY MATRIX
	int elemPosition;
	double corFactor; // not used in current version
	for (int i = 0; i < nEl; i++){ // Loop over all the elements
		fgets(s, 256, INPUT);
		sscanf(s, "%d %d %d %d %d %d %d %d %d %d %lf", &elementNumber, &cMesh.cMat[i][0], &cMesh.cMat[i][1], &cMesh.cMat[i][2], &cMesh.cMat[i][3], &cMesh.cMat[i][4], &cMesh.cMat[i][5], &cMesh.cMat[i][6], &cMesh.cMat[i][7], &cMesh.cMat[i][8], &corFactor);
	}

	/*LAYER INPUT*/
	int oldnEl = nEl;
	for (int i = 0; i < oldnEl; i++)
	{
		if (cMesh.cMat[i][8] == 303) //inner element
		{
			cMesh.addLayer(i, 80, 0.1);
		}
	}

	double intTemp = 180.0, outTemp = 30.0;
	/* ####### READING BOUNDARY CONDITIONS #########*/
	heatGenList = Matrix.Vector_Allocate(nEl, 0.0, &allocStatus); // Stores the Heat Generation of each Element
	fluxList = Matrix.Matrix_Allocate(nEl, 6, 0.0, &allocStatus); //     "      Heat flux in each Face of each Element (HEX8)
	convList = Matrix.Matrix_Allocate(nEl, 6, 0.0, &allocStatus); //     "      convection coef. in each Face (HEX8)
	for (int i = 0; i < nEl; i++)
	{
		if (cMesh.cMat[i][8] == 301) //inner element
		{
			convList[i][0] = 30.0;
		}
		else if (cMesh.cMat[i][8] == 300 + nElCrossSec) // outter element
			convList[i][2] = 50.0;
	}

	int nBC, bElem; // Counter: Number of boundary conditions (nBC) and index: boundary Element (bElem) 
	double BCvalue; // Variable that store the BC value, it is used when the BC value and its index are read at the same time.


	// Free nodes *OBS: The nodes number have to be complete and sorted, so if there is
	// 100 nodes, the nodes have to be numbered from 1 to 100.
	int j = 0;
	for (int i = 0; i < nFreeNodes; i++){
		while ((fixedNodesList[j] <= i + 1 + j) && j < nFixedNodes) // This loop eliminates the nodes with Dirichlet BC
			j++;
		freeNodesList[i] = i + 1 + j;
	}
	fclose(INPUT);
	/* #################################################################################################  */
	printf("READ FINISHED\n");
	printf("nDoF: %d, freeDoF: %d, fixedDoF: %d\n", nDoF, nFreeNodes, nFixedNodes);
	system("PAUSE");

	// open binary file for storing the temperatures
	FILE *BINARYRESULT;
	BINARYRESULT = fopen("FEM-Tyre_Job1.tbi", "wb+");

	int EntryNumbers[3], Values_n2[3];
	double Values_n3[2];
	int BinaryEntryCount;
	int BinaryEntryTestCount;
	char JobTitle[512];

	// --- reset Count for BinaryEntries 
	BinaryEntryCount = 1;
	BinaryEntryTestCount = 1;

	// --- Lenght of all Entries
	EntryNumbers[0] = 1;    //n1
	EntryNumbers[1] = 3;    //n2
	EntryNumbers[2] = 1;    //n3

	// --- Values of n2
	Values_n2[0] = nFreeNodes * 3;
	Values_n2[1] = nFixedNodes * 3;
	Values_n2[2] = 100000;

	// --- Values of n3
	Values_n3[0] = (180.0) / (100.0);

	// --- write Entries into binary File
	fwrite(EntryNumbers, sizeof(int), 3, BINARYRESULT);          // lenght
	fwrite(JobTitle, sizeof(char), 512, BINARYRESULT);          // n1
	fwrite(Values_n2, sizeof(int), EntryNumbers[1], BINARYRESULT);          // n2
	fwrite(Values_n3, sizeof(double), EntryNumbers[2], BINARYRESULT);          // n3

	/*Allocating matrices and vectors to assembly, solve and output the problem.*/
	//SG = Matrix.Matrix_Allocate(nDoF, nDoF, 0.0, &allocStatus); // Full Stiffness matrix
	//gsl_spmatrix *SGt = gsl_spmatrix_alloc(nDoF, nDoF);
	gsl_spmatrix *SGt = gsl_spmatrix_alloc_nzmax(nDoF, nDoF, 64 * nEl, GSL_SPMATRIX_TRIPLET);
	gsl_spmatrix_set_zero(SGt);
	if (allocStatus) printf("SG Allocated succefully\n"); system("PAUSE");
	//SSol = Matrix.Matrix_Allocate(nFreeNodes, nFreeNodes, 0.0, &allocStatus); // Just solving free Nodes
	gsl_spmatrix *SSolt = gsl_spmatrix_alloc_nzmax(nFreeNodes, nFreeNodes, 64 * nEl, GSL_SPMATRIX_TRIPLET);
	if (allocStatus) printf("SSol Allocated succefully\n");
	system("PAUSE");

	FG = Matrix.Vector_Allocate(nDoF, 0.0, &allocStatus);
	FSol = Matrix.Vector_Allocate(nFreeNodes, 0.0, &allocStatus);


	temp = Matrix.Vector_Allocate(nDoF, 0.0, &allocStatus);
	tempSol = Matrix.Vector_Allocate(nFreeNodes, 0.0, &allocStatus);

	/* Initializing classes objects */
	CMath8 math; //Math8 class init
	CElement elem; //Element class init

	/* ####################### ASSEMBLY (HEX8) ##########################*/
	/* Declaring and allocating local vectors and matrices*/
	double x[8], y[8], z[8]; // Vectors with nodes coordinates of 'ElementCount'

	double localS[64], localSConvection[64], heatGenFlocal[8], convFlocal[8], fluxFlocal[8], Fsum[8], Ssum[64];
	/* ASSEMBLY VARIABLES
	LocalS: Local Stiffness matrix IN VECTOR NOTATION M[m][n]=V[8*m+n]
	LocalSConvection: Local Stiffness matrix due to convection
	heatGenFlocal:  Local Heat generation Vector
	convFlocal: Local Convection Force vector
	fluxFlocal:  Local Heat Flux foce vector
	Fsum: Local Force Vector with contributions of all BCs on the element.
	Ssum:  Local Stif. matrix with contribution of conduction & convection       */
	FILE *CONV;

	CONV = fopen("conv.out", "w+");

	for (int elementcount = 1; elementcount <= nEl; elementcount++) //loop over all elements
	{
		for (int h = 0; h < 8; h++) // Loop over each node
		{
			// Get global coordinates of every node in element 'elementcount'
			x[h] = cMesh.xn[cMesh.cMat[elementcount - 1][h] - 1];
			y[h] = cMesh.yn[cMesh.cMat[elementcount - 1][h] - 1];
			z[h] = cMesh.zn[cMesh.cMat[elementcount - 1][h] - 1];

		}
		elem.Init(x, y, z); // Initializing object elem, this avoid passing pointer to x, y, z every call

		memset(heatGenFlocal, 0, 8 * sizeof(double));
		memset(Fsum, 0, 8 * sizeof(double));
		memset(Ssum, 0, 64 * sizeof(double));
		memset(localS, 0, 64 * sizeof(double));

		elem.Comp_Conductivity(matConductivity, localS); // Compute the Local S (Local Stiffness matrix)

		if (heatGenList[elementcount - 1] != 0.0) // Check if there is Heat generation on the element
			elem.Comp_HeatGeneration(heatGenList[elementcount - 1], heatGenFlocal, 1); // Compute Heat gen.

		for (int w = 0; w < 6; w++) // Loop over the faces to calc. flux BCs
			if (fluxList[elementcount - 1][w] != 0.0){
				elem.Comp_Flux(fluxList[elementcount - 1][w], fluxFlocal, faceNumberConverter(w));
				printf("Flux on element %d of %f\n", elementcount, fluxList[elementcount - 1][w]);
				math.compute_sumVec(fluxFlocal, Fsum);
			}


		for (int w = 0; w < 6; w++) // Loop over the faces to calculate the Convection
			if (convList[elementcount - 1][w] != 0){
				if (cMesh.cMat[elementcount - 1][8] == 301) //inner element
				{
					ambTemp = intTemp;
					//fprintf(CONV, "%d:0\n", elementcount);
				}
				else if (cMesh.cMat[elementcount - 1][8] == 300 + nElCrossSec) // outter element
				{
					ambTemp = outTemp;
					fprintf(CONV, "%d:2\n", elementcount);
				}
				elem.Comp_Convection_Diagonal(convList[elementcount - 1][w], ambTemp, localSConvection, convFlocal, faceNumberConverter(w));
				math.compute_sumVec(convFlocal, Fsum);
				math.compute_sum(localSConvection, Ssum);
				//printf("Convecttion on element %d of %f\n", elementcount, convList[elementcount - 1][w]);
			}

		for (int m = 0; m < 8; m++)
		{
			for (int n = 0; n < 8; n++)
			{
				double v = gsl_spmatrix_get(SGt, cMesh.cMat[elementcount - 1][m] - 1, cMesh.cMat[elementcount - 1][n] - 1);
				//SG[cMesh.cMat[elementcount - 1][m] - 1][cMesh.cMat[elementcount - 1][n] - 1] += math.ArrToMat(localS, m, n) + math.ArrToMat(Ssum, m, n);
				gsl_spmatrix_set(SGt, cMesh.cMat[elementcount - 1][m] - 1, cMesh.cMat[elementcount - 1][n] - 1, gsl_spmatrix_get(SGt, cMesh.cMat[elementcount - 1][m] - 1, cMesh.cMat[elementcount - 1][n] - 1) + math.ArrToMat(localS, m, n) + math.ArrToMat(Ssum, m, n));
				v = gsl_spmatrix_get(SGt, cMesh.cMat[elementcount - 1][m] - 1, cMesh.cMat[elementcount - 1][n] - 1);
			}
			FG[cMesh.cMat[elementcount - 1][m] - 1] += Fsum[m] + heatGenFlocal[m];
		}
	}
	fclose(CONV);
	printf("ASSEMBLY FINISHED, APPLYING BC\n");
	system("PAUSE");
	/*Applying boundary conditions on the force matrix*/
	for (int fNode = 0; fNode < nFreeNodes; fNode++) // Loop over the free nodes (fNode)
		for (int bNode = 0; bNode < nFixedNodes; bNode++) // Loop over the nodes with Diri. BC (bNode stands for Bounded node)
		{
			//FG[freeNodesList[fNode] - 1] -= tempFixedNodes[bNode] * SG[freeNodesList[fNode] - 1][fixedNodesList[bNode] - 1];
			FG[freeNodesList[fNode] - 1] -= tempFixedNodes[bNode] * gsl_spmatrix_get(SGt, freeNodesList[fNode] - 1, fixedNodesList[bNode] - 1);
		}

	/* Taking SSol out of SG, because smaller the matrix faster the solver*/
	for (int i = 0; i < nFreeNodes; i++){
		for (int j = 0; j < nFreeNodes; j++)
			gsl_spmatrix_set(SSolt, i, j, gsl_spmatrix_get(SGt, freeNodesList[i] - 1, freeNodesList[j] - 1));
		//SSol[i][j] = SG[freeNodesList[i] - 1][freeNodesList[j] - 1];

		FSol[i] = FG[freeNodesList[i] - 1];
	}
	printf("BC APPLIED, INIT SOLVER\n");
	system("PAUSE");
	ItSolver(SSolt, FSol, tempSol, nFreeNodes); // Calling the solver, the result will be stored in tempSol vector

	double *bigtemp;
	bigtemp = Matrix.Vector_Allocate((nFixedNodes + nFreeNodes) * 3, 0.0, &allocStatus);

	/* Gathering all nodes temperatures in one vector */
	for (int i = 0; i < nFixedNodes; i++)
	{
		temp[fixedNodesList[i] - 1] = tempFixedNodes[i];
		bigtemp[(fixedNodesList[i] - 1) * 3] = tempFixedNodes[i];
		bigtemp[(fixedNodesList[i] - 1) * 3 + 1] = tempFixedNodes[i];
		bigtemp[(fixedNodesList[i] - 1) * 3 + 2] = tempFixedNodes[i];
	}

	for (int i = 0; i < nFreeNodes; i++)
	{
		temp[freeNodesList[i] - 1] = tempSol[i];
		bigtemp[(freeNodesList[i] - 1) * 3] = tempSol[i];
		bigtemp[(freeNodesList[i] - 1) * 3 + 1] = tempSol[i];
		bigtemp[(freeNodesList[i] - 1) * 3 + 2] = tempSol[i];
	}

	BinaryOutput(BINARYRESULT, bigtemp, (nFixedNodes + nFreeNodes) * 3);

	/* Writing Result file*/
	FILE *OUT;

	OUT = fopen("test.out", "w+");

	for (int i = 0; i<nDoF; i++)
		fprintf(OUT, "%+12.6e\n", temp[i]);

	fclose(OUT);

	return 0;
}