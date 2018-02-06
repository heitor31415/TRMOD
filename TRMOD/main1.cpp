#include <iostream>
#include "Element.h"
#include "misc.h"
#include "MatrixStuff_Modal.h"
#include <string>
using namespace std;
int main(){
	/*  MESH SIZING AND CONNECTIVITY MATRIX */
	int nDoF, nEl, nElCrossSec; // Number of degrees of Freedom(= num.Nodes in thermal analysis) and number of elements
	int **connectMatrix; // Connectivity matrix, which contains elements global nodes numbers.

	/* NODES POSITIONS AND GENERAL VARIABLES */
	double *xn, *yn, *zn; // Global nodes' coordinates
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
	INPUT = fopen("TireLayer.dat", "r");// "fileName.dat"
	fgets(s, 256, INPUT);// >>Text line= MESH SIZING

	/* Read the sizing values, which are Number of Elements and number of nodes. */
	fgets(s, 256, INPUT);
	sscanf(s, "%d %d %d", &nDoF, &nEl, &nElCrossSec); // Number of elements & nodes on the model, resp.

	/* Allocating memory to store node's coordinates*/
	xn = Matrix.Vector_Allocate(nDoF, 0.0, &allocStatus); // nDoF is the number of nodes/Deg. of Freedom(nDoF)
	yn = Matrix.Vector_Allocate(nDoF, 0.0, &allocStatus);
	zn = Matrix.Vector_Allocate(nDoF, 0.0, &allocStatus);

	/* Store nodes' global coordinates >>(NODES HAVE TO BE SORTED)<<*/
	int globalNodeNumber, elementType, elementNumber; // Not used in the current version of the code
	fgets(s, 256, INPUT); // >>Text line= GLOBAL NODES COORDINATES

	int layer, fixedNodesCounter = 0; // Layer not used in the current version
	bool fixedNodesStart = true;
	/* Nodes coordinates input */
	for (int i = 0; i < nDoF; i++){ //Loop over all nodes
		fgets(s, 256, INPUT);
		sscanf(s, "%d %lf %lf %lf %d", &globalNodeNumber, &xn[i], &yn[i], &zn[i], &layer); // Scan each node coordinate
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
			tempFixedNodes[fixedNodesCounter] = 120.0;
			fixedNodesCounter++;
		}
	}

	/* Allocating and scanning the Connectivity matrix (Matrix with element's nodes)*/
	connectMatrix = Matrix.Matrix_Allocate_Int(nEl, 9, 0, &allocStatus); //8 = Number of nodes per element (HEX8)

	fgets(s, 256, INPUT); // >>Text Line= CONNECTIVITY MATRIX
	int elemPosition;
	double corFactor; // not used in current version
	for (int i = 0; i < nEl; i++){ // Loop over all the elements
		fgets(s, 256, INPUT);
		sscanf(s, "%d %d %d %d %d %d %d %d %d %d %lf", &elementNumber, &connectMatrix[i][0], &connectMatrix[i][1], &connectMatrix[i][2], &connectMatrix[i][3], &connectMatrix[i][4], &connectMatrix[i][5], &connectMatrix[i][6], &connectMatrix[i][7], &connectMatrix[i][8], &corFactor);
	}


	double intTemp = 180.0, outTemp = 30.0;
	/* ####### READING BOUNDARY CONDITIONS #########*/
	heatGenList = Matrix.Vector_Allocate(nEl, 0.0, &allocStatus); // Stores the Heat Generation of each Element
	fluxList = Matrix.Matrix_Allocate(nEl, 6, 0.0, &allocStatus); //     "      Heat flux in each Face of each Element (HEX8)
	convList = Matrix.Matrix_Allocate(nEl, 6, 0.0, &allocStatus); //     "      convection coef. in each Face (HEX8)
	for (int i = 0; i < nEl; i++)
	{
		if (connectMatrix[i][8] == 301) //inner element
		{
			convList[i][0] = 30.0;
		}
		else if (connectMatrix[i][8] == 300 + nElCrossSec) // outter element
			convList[i][2] = 50.0;
	}

	int nBC, bElem; // Counter: Number of boundary conditions (nBC) and index: boundary Element (bElem) 
	double BCvalue; // Variable that store the BC value, it is used when the BC value and its index are read at the same time.

	/* Reading Heat Generation Boundary conditions
	fgets(s, 256, INPUT); // TXT = heat generation list

	fgets(s, 256, INPUT);
	sscanf(s, "%d", &nBC); // Number of elements with heat generation
	// Faces with flux
	fgets(s, 256, INPUT); //TXT= choosing elements with heat gen.
	for (int i = 0; i < nBC; i++){ // Loop over number of elements with Heat generation
	fgets(s, 256, INPUT);
	sscanf(s, "%d %lf", &bElem, &BCvalue);
	heatGenList[bElem - 1] = BCvalue; // 0 indexed list.
	}

	/* Reading flux Boundary conditions*/
	/*
	fgets(s, 256, INPUT); // >>Text line= Elements with flux
	fgets(s, 256, INPUT);
	sscanf(s, "%d", &nBC); // Number of ELEMENTS with heat flux BC
	fgets(s, 256, INPUT); // >>Text line= choosing the faces to put the flux on
	for (int i = 0; i < nBC; i++){ //Loop over all ELEMENTS with heat flux BC
	fgets(s, 256, INPUT);
	sscanf(s, "%d", &bElem); // Element with flux BC
	bElem--; // 0 indexed vector
	fgets(s, 256, INPUT); // Each column contains the normal heat flux value. 0 if there is no heat flux on the face.
	sscanf(s, "%lf %lf %lf %lf %lf %lf", &fluxList[bElem][0], &fluxList[bElem][1], &fluxList[bElem][2], &fluxList[bElem][3], &fluxList[bElem][4], &fluxList[bElem][5]);
	}
	fgets(s, 256, INPUT); // >>Text line= Elements with flux
	fgets(s, 256, INPUT);
	sscanf(s, "%d", &nBC); // Number of ELEMENTS with heat flux BC
	/* Reading Convection(Film) Boundary conditions*/
	/*
	fgets(s, 256, INPUT);//>>Text line= choosing the faces to put convection on
	for (int i = 0; i < nBC; i++){ //Loop over all ELEMENTS with convection BC
	fgets(s, 256, INPUT);
	sscanf(s, "%d", &bElem);// Element with convection(film) BC
	bElem--; // 0 indexed vector
	fgets(s, 256, INPUT); // Each column must contain the convection(film) coeficient (h). ambTemp will be used in all faces
	sscanf(s, "%lf %lf %lf %lf %lf %lf", &convList[bElem][0], &convList[bElem][1], &convList[bElem][2], &convList[bElem][3], &convList[bElem][4], &convList[bElem][5]);
	}
	fgets(s, 256, INPUT); // TXT=Num of Nodes with Dirichlet BC
	fgets(s, 256, INPUT);
	sscanf(s, "%d", &nFixedNodes);

	/* Reading Dirichlet boundary conditions*/
	/*
	// List of nodes with Dirichlet BC
	fgets(s, 256, INPUT); // >>Text line: List of nodes with Dirichlet BC
	for (int i = 0; i < nFixedNodes; i++){
	fgets(s, 256, INPUT);
	sscanf(s, "%d %lf", &fixedNodesList[i], &tempFixedNodes[i]);
	}*/

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
			x[h] = xn[connectMatrix[elementcount - 1][h] - 1];
			y[h] = yn[connectMatrix[elementcount - 1][h] - 1];
			z[h] = zn[connectMatrix[elementcount - 1][h] - 1];

		}
		elem.Init(x, y, z); // Initializing object elem, this avoid passing pointer to x, y, z every call

		memset(heatGenFlocal, 0, 8 * sizeof(double));
		memset(Fsum, 0, 8 * sizeof(double));
		memset(Ssum, 0, 64 * sizeof(double));
		memset(localS, 0, 64 * sizeof(double));

		if (connectMatrix[elementcount - 1][8] == 302)
		{ 
			elem.Comp_Conductivity(90, localS); // Compute the Local S (Local Stiffness matrix)
			fprintf(CONV, "%d\n", elementcount);
		}
		else
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
				if (connectMatrix[elementcount - 1][8] == 301) //inner element
				{
					ambTemp = intTemp;
					//fprintf(CONV, "%d:0\n", elementcount);
				}
				else if (connectMatrix[elementcount - 1][8] == 300 + nElCrossSec) // outter element
				{
					ambTemp = outTemp;
					//fprintf(CONV, "%d:2\n", elementcount);
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
				double v = gsl_spmatrix_get(SGt, connectMatrix[elementcount - 1][m] - 1, connectMatrix[elementcount - 1][n] - 1);
				//SG[connectMatrix[elementcount - 1][m] - 1][connectMatrix[elementcount - 1][n] - 1] += math.ArrToMat(localS, m, n) + math.ArrToMat(Ssum, m, n);
				gsl_spmatrix_set(SGt, connectMatrix[elementcount - 1][m] - 1, connectMatrix[elementcount - 1][n] - 1, gsl_spmatrix_get(SGt, connectMatrix[elementcount - 1][m] - 1, connectMatrix[elementcount - 1][n] - 1) + math.ArrToMat(localS, m, n) + math.ArrToMat(Ssum, m, n));
				v = gsl_spmatrix_get(SGt, connectMatrix[elementcount - 1][m] - 1, connectMatrix[elementcount - 1][n] - 1);
			}
			FG[connectMatrix[elementcount - 1][m] - 1] += Fsum[m] + heatGenFlocal[m];
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