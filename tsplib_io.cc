	/* To use TSPLIB, set the environ variable "TSPLIB_DIR"
	   for example, "$setenv TSPLIB_DIR $HOME/ga/TSPLIB"
	 */

const double EPS = 1e-9;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include "macro.h"
#include "tsplib_io.h"

	/*************************** Debug Flag ******************************/
#define affirm(x)
#define trace

	//#define affirm(x)           assert(x)
	//#define trace               printf


	/*************************** Macros **********************************/
#define GRAPH_INFO_SIZE      (sizeof(gGraphInfo)/sizeof(GRAPH_INFO))
#define MAX_COORD_SYSTEM     20

	/* tsplib structure */
typedef struct
{
	char       szGraphName[64];
	char       szType[64];
	char       szComment[256];
	int        nDimension;
	int        nCoordDim;
	char       szGraphType[64];
	char       szEdgeType[64];
	char       szEdgeWeightType[64];
	char       szEdgeWeightFormat[64];
	char       szEdgeDataFormat[64];
	char       szNodeType[64];
	char       szNodeCoordType[64];
	char       szDisplayDataType[64];
} TSP_FILE_INFO;

typedef struct
{
	double     pt[MAX_COORD_SYSTEM];
} POINT;

typedef struct
{
	char*      graph_name;
	int        opt_cost;
} GRAPH_INFO;

	/********************* External Global Variables ********************/
unsigned int    gNumCity;
int             gType;          /* is Symmetric or Asysmmetric ? */
int             gNumNN;         /* is # of nearest neighbors */
double*            gDistMat;       /* distance matrix */
int*            gNNI;           /* is nearest neighbors arrrays */

	/********************* Internal Global Varialbes ********************/
POINT*          gNodeCoords;    /* When NODE_COORD_SECTION exists, used */
TSP_FILE_INFO   gtfi;
	/* function pointer return distance of (i, j) */
static int      (*internalCalcOnline)(int i, int j);             
static char     gGraphName[64]; /* graph name (without file extension) */
int             gOptimumCost;   /* opt. tour cost or lower bound if not, 0 */
int*            gOrgInfo;       /* index : remapped city's id
								   value : original city's id */

	/************************ Internal Function Declaration ***************/
	// Error printing and exiting
void Error(const char *format, const char *str);

	// Neighbor List related.
void constructNormal(int* NNI);
void constructQuadNeighbors(int* NNI);
void sortNeighbors(int* neighbors, double* ndist, int size, int numNN);

	// Reading File related.
void GetTspFileInfo(FILE *fd, TSP_FILE_INFO *tfi, char *nextl);
void ReadDataInNodeCoord(FILE *fd, TSP_FILE_INFO *tfi);
void ReadDataInEdgeWeight(FILE *fd, TSP_FILE_INFO *tfi);

	// Calculating Distances related
double GetEuc2DDist(int c1, int c2);
int GetEuc3DDist(int c1, int c2);
int GetEucNDDist(int c1, int c2);
int GetCeil2DDist(int c1, int c2);
int GetGeoDist(int c1, int c2);
int GetAttDist(int c1, int c2);
int GetFullMatDist(int c1, int c2);
int GetHalfMatDist(int c1, int c2);

	// Optimal Cost related
void SetOptimumCost(char* graphname);

	// distance cache related
void initDistanceCache(int num_city, int num_nn);

	/************************** Fuction Definition ***************/
	/****************** Reading File related ********************/
void ReadTspFile(char* graphname)
{
	FILE *fd;
	char dummy[256], graph[256], *tsplib_dir;
	
	trace("Entering ReadTspFile\n");
	tsplib_dir = getenv("TSPLIB_DIR");

	strcpy(graph, graphname);  // try to find out in the current directory.
	fd = fopen(graph, "r");
	if( !fd && tsplib_dir) // try to find out in the specified directory.
	{
		sprintf(graph, "%s/%s", tsplib_dir, graphname);
		fd = fopen(graph, "r");
	}
	if( !fd && tsplib_dir) // try to find out in the second spec. dir.
	{
		sprintf(graph, "%s/MYLIB/%s", tsplib_dir, graphname);
		fd = fopen(graph, "r");
	}
	if( !fd) Error(" %s data file not found..\n", graph);

	GetTspFileInfo(fd, &gtfi, dummy);

	gNumCity = gtfi.nDimension;
/*
	if( strcmp(gtfi.szEdgeWeightType, "EXPLICIT")!=0 &&
		strcmp(gtfi.szEdgeWeightType, "EUC_2D")!=0 &&
		strcmp(gtfi.szEdgeWeightType, "EUC_3D")!=0 &&
		strcmp(gtfi.szEdgeWeightType, "EUC_ND")!=0 &&
		strcmp(gtfi.szEdgeWeightType, "CEIL_2D")!=0 &&
		strcmp(gtfi.szEdgeWeightType, "GEO")!=0 &&
		strcmp(gtfi.szEdgeWeightType, "ATT")!=0)
		Error(" not supported EDGE_WEIGHT_TYPE : %s\n", gtfi.szEdgeWeightType);

	if( strcmp(gtfi.szEdgeWeightType, "EXPLICIT")==0 &&
		strcmp(gtfi.szEdgeWeightFormat, "FULL_MATRIX")!=0)
		Error(" not supported EDGE_WEIGHT_FORMAT with EXPLICIT type: %s\n",
			   gtfi.szEdgeWeightFormat);

	if( strncmp(dummy, "NODE_COORD_SECTION", 18) == 0)
		ReadDataInNodeCoord(fd, &gtfi);
	else if( strncmp(dummy, "EDGE_WEIGHT_SECTION", 19) == 0)
		ReadDataInEdgeWeight(fd, &gtfi);
	else
		Error(" unimplemented section ... %s\n", dummy);

	fclose(fd);
	gType = (strcmp(gtfi.szType, "ATSP")==0 ? TYPE_ATSP : TYPE_STSP);
*/
	ReadDataInNodeCoord(fd, &gtfi);
	gType = TYPE_STSP;
	trace("Quitting ReadTspFile\n");
}

int ReadTourFile(char* file, int *tour, int size)
{
	trace("Entering ReadTourFile\n");
	
	FILE*   fd;
	char    buf[1024];
	int     i, n;

	fd = fopen(file, "r");
	if( !fd)
	{
		printf(" %s data file not found..\n", file);
		return FALSE;
	}
	if( size < gNumCity)
	{
		printf(" tour array's size must be bigger.\n");
		return FALSE;
	}

	char *temp = fgets(buf, 1023, fd);
	while( strncmp(buf, "TOUR_SECTION", 12) != 0)
		char *temp = fgets(buf, 1023, fd);
	
	for( i=0; i<gNumCity; i++)
	{
		fscanf(fd, "%d", &tour[i]);
		tour[i]--;
	}
	fclose(fd);
	trace("Quitting ReadTourFile\n");
	return TRUE;
}

int ReadOptTourFile(int *opt_tour, int size)
{
	trace("Entering ReadOptTourFile\n");

	FILE*   fd;
	char    dummy[256], graph[256], *tsplib_dir;
	int     res;

	CopyGraphName(dummy);
	tsplib_dir = getenv("TSPLIB_DIR");

	strcat(dummy, ".opt.tour");
	strcpy(graph, dummy);
	fd = fopen(graph, "r");
	if( !fd && tsplib_dir)
	{
		sprintf(graph, "%s/%s", tsplib_dir, dummy);
		fd = fopen(graph, "r");
	}
	if( !fd && tsplib_dir)
	{
		sprintf(graph, "%s/MYLIB/%s", tsplib_dir, dummy);
		fd = fopen(graph, "r");
	}

	if( !fd)
	{
		printf(" %s data file not found..\n", graph);
		return 0;
	}
	fclose(fd);

	res = ReadTourFile(graph, opt_tour, size);
	trace("Quitting ReadOptTourFile\n");
	return res;
}

void ReadDataInNodeCoord(FILE* fd, TSP_FILE_INFO* tfi)
{
	int i, j, r, c, id, n;
	double (*ComputeDist)(int c1, int c2);

	trace("Entering ReadDataInNodeCoord\n");
	n = tfi->nDimension;
	/*
	if( strcmp(tfi->szEdgeWeightType, "EUC_2D")==0)
		ComputeDist = GetEuc2DDist;
	else if( strcmp(tfi->szEdgeWeightType, "EUC_3D")==0)
	{
		ComputeDist = GetEuc3DDist;
		tfi->nCoordDim = 3;
	}
	else if( strcmp(tfi->szEdgeWeightType, "EUC_ND")==0)
		ComputeDist = GetEucNDDist;
	else if( strcmp(tfi->szEdgeWeightType, "CEIL_2D")==0)
		ComputeDist = GetCeil2DDist;
	else if( strcmp(tfi->szEdgeWeightType, "GEO")==0)
		ComputeDist = GetGeoDist;
	else if( strcmp(tfi->szEdgeWeightType, "ATT")==0)
		ComputeDist = GetAttDist;
	else
		Error("Error : In ReadDataInNodeCoord()%s\n", "");
	if( tfi->nCoordDim<=0) tfi->nCoordDim = 2;
	if( tfi->nCoordDim>MAX_COORD_SYSTEM)
		Error("Error : Too Large Coord. Dimension%s\n", "");
*/
	gNodeCoords = new POINT[n];
	tfi->nCoordDim = 2;
	ComputeDist = GetEuc2DDist;
	for( i=0; i<n; i++)
	{
//        fscanf(fd, "%d", &id);
		for( j=0; j<tfi->nCoordDim; j++)
			fscanf(fd, "%lf", &gNodeCoords[i].pt[j]);
	}
	
	gDistMat = new double[(n*(n-1))/2];
	for( r=0; r<n-1; r++)
		for( c=r+1; c<n; c++)
			gDistMat[c*(c-1)/2 + r] = ComputeDist(r, c);

	trace("Quitting ReadDataInNodeCoord\n");
}

void ReadDataInEdgeWeight(FILE* fd, TSP_FILE_INFO* tfi)
{
}

	/************* Printing, Error, Util Function Definition *************/ 
void PrintTspInfo()
{
	TSP_FILE_INFO* tfi = &gtfi;

	trace("Entering PrintTspInfo\n");
	puts("================ TSP File INFO =================");
	printf("Graph Name         : %s\n", tfi->szGraphName);
	printf("Type               : %s\n", tfi->szType);
	puts(tfi->szComment);
	printf("Dimension          : %d\n", tfi->nDimension);
	if( tfi->szGraphType[0]) 
	printf("Graph Type         : %s\n", tfi->szGraphType);
	if( tfi->szEdgeType[0])  
		printf("Edge Type          : %s\n", tfi->szEdgeType);
	if( tfi->szEdgeWeightType[0]) 
		printf("Edge Weight Type   : %s\n", tfi->szEdgeWeightType);
	if( tfi->szEdgeWeightFormat[0])
		printf("Edge Weight Format : %s\n", tfi->szEdgeWeightFormat);
	if( tfi->szEdgeDataFormat[0])
		printf("Edge Data Format   : %s\n", tfi->szEdgeDataFormat);
	if( tfi->szNodeType[0])
		printf("Node Type          : %s\n", tfi->szNodeType);
	if( tfi->szNodeCoordType[0])
		printf("Node Coord. Type   : %s\n", tfi->szNodeCoordType);
	if( tfi->szDisplayDataType[0])
		printf("Display Data Type  : %s\n", tfi->szDisplayDataType);
	puts("================================================");
	trace("Quitting PrintTspInfo\n");
}

void Error(const char* format, const char* str)
{
	printf(format, str);
	exit(1);
} 

void PrintCoords(FILE* fd, int* tour, int size)
{
	int i, j, t;
	trace("Entering PrintCoords\n");
	if( gNodeCoords)
	{
		for( i=0; i<size; i++)
		{
			t = (gOrgInfo) ? gOrgInfo[tour[i]] : tour[i];
			for( j=0; j<gtfi.nCoordDim; j++)
				fprintf(fd, "%f ", gNodeCoords[t].pt[j]);
			fprintf(fd, "\n");
		}
	}
	trace("Quitting PrintCoords\n");
}

int GetCoords(int city, double* pt, int size)
{
	int j;
	if( gtfi.nCoordDim>size) return -gtfi.nCoordDim;
	for( j=0; j<gtfi.nCoordDim; j++)
		pt[j] = gNodeCoords[city].pt[j];
	return gtfi.nCoordDim;
}

void PrintTour(FILE* fd, int* tour, int size,
			   char* name, char* comment)
{
	trace("Entering PrintTour\n");
	int i, t;
	char   graph[256];

	CopyGraphName(graph);

	if( name) fprintf(fd, "NAME : %s\n", name);
	if( comment) fprintf(fd, "COMMENT : %s\n", comment);
	else         fprintf(fd, "COMMENT : %s\n", graph);
	fprintf(fd, "TYPE : TOUR\n");
	fprintf(fd, "DIMENSION : %d\n", size);
	fprintf(fd, "TOUR_SECTION\n");

	for( i=0; i<size; i++)
	{
		t = (gOrgInfo) ? gOrgInfo[tour[i]] : tour[i];
		fprintf(fd, "%d ", t+1);
		if( i%16 == 15 ) fprintf(fd, "\n");
	}
	fprintf(fd, "\n");
	fprintf(fd, "-1\nEOF\n");
	trace("Quitting PrintTour\n");
}

void GetTspFileInfo(FILE* fd, TSP_FILE_INFO* tfi, char* nextl)
{
//	char buf[1024], keyword[64], value[128], dummy[128];
//	int n;
//
//	trace("Entering GetTspFileInfo\n");
//	memset(tfi, 0, sizeof(*tfi));
//	
//	while( 1)
//	{
//		fgets(buf, 1023, fd);
//		if( buf[strlen(buf)-1] == '\n')
//			buf[strlen(buf)-1] = 0;
//		n = sscanf(buf, "%s%s%s", keyword, dummy, value);
//		if( n == 1) /* This is not of the specification part. */
//		{
//			strcpy(nextl, buf);
//			/*if( nextl[strlen(nextl)-1] == '\n')
//				nextl[strlen(nextl)-1] = 0;*/
//			trace("Quitting GetTspFileInfo\n");
//			return;
//		}
//		else if( n==2) /* keyword contains <keyword> following by ":". */
//		{              /* so, dummy contains <value>. fix!! */
//			keyword[strlen(keyword)-1] = 0; /* remove ":" */
//			strcpy(value, dummy);
//		}
//		/*if( value[strlen(value)-1] == '\n')
//			value[strlen(value)-1] = 0;*/
//		
//		if( strcmp(keyword, "NAME")==0) strcpy(tfi->szGraphName, value);
//		else if( strcmp(keyword, "TYPE")==0) strcpy(tfi->szType, value);
//		else if( strcmp(keyword, "COMMENT")==0) strcpy(tfi->szComment, buf);
//		else if( strcmp(keyword, "COMMENT:")==0) strcpy(tfi->szComment, buf);
//		else if( strcmp(keyword, "DIMENSION")==0) tfi->nDimension = atoi(value);
//		else if( strcmp(keyword, "COORD_DIM")==0) tfi->nCoordDim = atoi(value);
//		else if( strcmp(keyword, "GRAPH_TYPE")==0) strcpy(tfi->szGraphType, value);
//		else if( strcmp(keyword, "EDGE_TYPE")==0) strcpy(tfi->szEdgeType, value);
//		else if( strcmp(keyword, "EDGE_WEIGHT_TYPE")==0) strcpy(tfi->szEdgeWeightType, value);
//		else if( strcmp(keyword, "EDGE_WEIGHT_FORMAT")==0) strcpy(tfi->szEdgeWeightFormat, value);
//		else if( strcmp(keyword, "EDGE_DATA_FORMAT")==0) strcpy(tfi->szEdgeDataFormat, value);
//		else if( strcmp(keyword, "NODE_TYPE")==0) strcpy(tfi->szNodeType, value);
//		else if( strcmp(keyword, "NODE_COORD_TYPE")==0) strcpy(tfi->szNodeCoordType, value);
//		else if( strcmp(keyword, "DISPLAY_DATA_TYPE")==0) strcpy(tfi->szDisplayDataType, value);
//	}
	fscanf(fd,"%d",&(tfi->nDimension));
//	tfi->nDimension++;
}

void CopyGraphName(char* dest)
{
	strcpy(dest, gGraphName);
}

	/*************** Calculating Distances Function Definition **************/
double GetEuc2DDist(int c1, int c2)
{
	double xd, yd;
	xd = gNodeCoords[c1].pt[0] - gNodeCoords[c2].pt[0]; 
	yd = gNodeCoords[c1].pt[1] - gNodeCoords[c2].pt[1];
    return sqrt(xd*xd + yd*yd);
}

int GetEuc3DDist(int c1, int c2)
{
	double xd, yd, zd;
	xd = gNodeCoords[c1].pt[0] - gNodeCoords[c2].pt[0]; 
	yd = gNodeCoords[c1].pt[1] - gNodeCoords[c2].pt[1];
	zd = gNodeCoords[c1].pt[2] - gNodeCoords[c2].pt[2];
    return nint(sqrt(xd*xd + yd*yd + zd*zd));
}

int GetEucNDDist(int c1, int c2)
{
	int j;
	double d, sd = 0;
	for( j=0; j<gtfi.nCoordDim; j++)
	{
		d = gNodeCoords[c1].pt[j] - gNodeCoords[c2].pt[j];
		sd += (d*d);
	}
    return nint(sqrt(sd));
}

int GetCeil2DDist(int c1, int c2)
{
	double xd, yd, t1, t2;
	xd = gNodeCoords[c1].pt[0] - gNodeCoords[c2].pt[0]; 
	yd = gNodeCoords[c1].pt[1] - gNodeCoords[c2].pt[1];
    return cint(sqrt(xd*xd + yd*yd));
}

int GetGeoDist(int c1, int c2)
{
	double deg, min, latitude1, longitude1, q1, q2, q3;
	double latitude2, longitude2;
	int t;
	
	deg = t = gNodeCoords[c1].pt[0];
	min = gNodeCoords[c1].pt[0] - deg;
	latitude1 = 3.141592654*(deg + 5.0*min/3.0)/180.0;
	deg = t = gNodeCoords[c1].pt[1];
	min = gNodeCoords[c1].pt[1] - deg;
	longitude1 = 3.141592654*(deg + 5.0*min/3.0)/180.0;
	
	deg = t = gNodeCoords[c2].pt[0];
	min = gNodeCoords[c2].pt[0] - deg;
	latitude2 = 3.141592654*(deg + 5.0*min/3.0)/180.0;
	deg = t = gNodeCoords[c2].pt[1];
	min = gNodeCoords[c2].pt[1] - deg;
	longitude2 = 3.141592654*(deg + 5.0*min/3.0)/180.0;

	q1 = cos(longitude1 - longitude2);
	q2 = cos(latitude1 - latitude2);
	q3 = cos(latitude1 + latitude2);

	return (int)(6378.388*acos(0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0);
}
	
int GetAttDist(int c1, int c2)
{
	static double xd, yd, t1, t2, rij;
	int tij;

	xd = gNodeCoords[c1].pt[0] - gNodeCoords[c2].pt[0]; 
	yd = gNodeCoords[c1].pt[1] - gNodeCoords[c2].pt[1];
	rij = sqrt((xd*xd + yd*yd)/10.0);
	tij = nint(rij);
	if( tij<rij)
		return tij + 1;
	else
		return tij;
}

	/******************* Online Calculating Distances ***********************/
	/* online calculation is allowed only in case of EUC-?D instance. */
#define ANOTHER_CACHES          4
struct DISTANCE_CACHE
{
	IntPair*    nn;   // nearest cities. i1 -> city number, i2 -> distance.
	int         cur;  // current index of another cached cities.
		// another cached cities. i1-> city number, i2 -> distance
	IntPair     another[ANOTHER_CACHES];
};
DISTANCE_CACHE* gCache;

void initDistanceCache(int num_city, int num_nn)
{
	affirm( gType == TYPE_STSP);

	int             i, j;
	IntPair         *home, *nns;

	gCache = new DISTANCE_CACHE[num_city];
	home   = new IntPair[num_city*num_nn];
	for( i=0; i<num_city; i++)
	{
		gCache[i].cur = 0;
		for( j=0; j<ANOTHER_CACHES; j++)
		{
			gCache[i].another[j].i1 = -1;
			gCache[i].another[j].i2 = -1;
		}

		gCache[i].nn  = &home[num_nn*i];
		nns           = gCache[i].nn;
		for( j=0; j<num_nn; j++)
		{
			nns[j].i1 = nni(i, j);  // save neighbor city.
			nns[j].i2 = internalCalcOnline(i, nns[j].i1); // save distances.
		}
	}
}

	//int gFailCount, gTotalCount;
int calcOnline(int c1, int c2)
{
	int         i;
	IntPair     *nns;

	if( gOrgInfo)
	{
		c1 = gOrgInfo[c1]; c2 = gOrgInfo[c2];
	}

		// gTotalCount++;
		// first trial : search among the nearest neighbor cities.
	nns = gCache[c1].nn;
	for( i=0; i<gNumNN; i++)
		if( nns[i].i1 == c2)
			return nns[i].i2;

		// second trial : search among another cities.
	nns = gCache[c1].another;
	for( i=0; i<ANOTHER_CACHES; i++)
		if( nns[i].i1 == c2)
			return nns[i].i2;

		//gFailCount++;
		// fail!!  save this city and distance.
	nns[gCache[c1].cur].i1 = c2;
	i                      = internalCalcOnline(c1, c2); // calc. distance.
	nns[gCache[c1].cur].i2 = i;

		// update current index.
	gCache[c1].cur++;
	if( gCache[c1].cur == ANOTHER_CACHES) gCache[c1].cur = 0;

	return i;
}

	/***************************************************************/
	/********************** remap implementation *******************/
	/* map is remapped city -> org city. */
void RemapCities(int* map)
{
	int i, j, c1, c2, mi, mx, rcity;
	int *temp, *rmap;
	
	gOrgInfo = new int[gNumCity];
	rmap = new int[gNumCity];
	memcpy(gOrgInfo, map, sizeof(int)*gNumCity);
	
	for( i=0; i<gNumCity; i++)
		rmap[map[i]] = i;
	
	temp = new int[(gNumCity*(gNumCity-1))/2];
	memcpy(temp, gDistMat, sizeof(int)*(gNumCity*(gNumCity-1))/2);
	for( i=0; i<gNumCity-1; i++)
	{
		for( j=i+1; j<gNumCity; j++)
		{
			c1 = min(rmap[i], rmap[j]);
			c2 = max(rmap[i], rmap[j]);
			mi = min(i, j); mx = max(i, j);
			
			gDistMat[c2*(c2-1)/2 + c1] = temp[j*(j-1)/2 + i];
		}
	}
	delete[] temp;
	
	if( gNNI)
	{
		int     rcity;
		
		temp = new int[gNumCity*gNumNN];
		
		memcpy(temp, gNNI, sizeof(int)*gNumCity*gNumNN);
		for( i=0; i<gNumCity; i++)
		{
			rcity = rmap[i];
			for( j=0; j<gNumNN; j++)
				gNNI[rcity*gNumNN + j] = rmap[temp[i*gNumNN + j]];
		}
		delete[] temp;
	}
	delete[] rmap;
}

int GetOrgCity(int city)
{
	return ((gOrgInfo)? gOrgInfo[city] : city);
}
			
	/***********************************************************/
	/********************* optimum cost ***********************/
void SetOptimumCost(char* gn)
{
}

int GetOptimumCost()
{
	return gOptimumCost;
}

	/***********************************************************/
	/*************** construct Nearest Neighbors array *********/
void ConstructNN(int numNN, int is_quadrant)
{
	trace("Entering constructNN\n");
	gNumNN = numNN;		
	if( is_quadrant &&
		(strcmp(gtfi.szEdgeWeightType, "EUC_2D")==0 ||
		strcmp(gtfi.szEdgeWeightType, "CEIL_2D")==0 ||
		strcmp(gtfi.szEdgeWeightType, "ATT")==0))
	{

		if( numNN%20) gNumNN = (numNN/20 + 1)*20;
		if( gNumNN >= gNumCity) gNumNN -= 20;
		gNNI = new int[gNumCity*gNumNN];
		constructQuadNeighbors(gNNI);
		
	}
	else
	{
		if( gNumNN >= gNumCity) gNumNN = gNumCity - 1;
		gNNI = new int[gNumCity*gNumNN];
		constructNormal(gNNI);
	}

		/* in the case of online calculation, setup a cache. */
#if     (DISTANCE_CALC_METHOD==2)
    initDistanceCache(gNumCity, gNumNN);
#endif

	trace("Quitting constructNN\n");
}

void constructNormal(int* NNI)
{
	trace("Entering constructNormal\n");
	
	int             i, j, count;
	int             *neighbors;
	double *ndist;
	
	neighbors = new int[gNumCity-1];
	ndist     = new double[gNumCity-1];
	
	for( i=0; i<gNumCity; i++)
	{
		for( j=0, count=0; j<gNumCity; j++)
			if( i != j)
			{
				neighbors[count] = j;

#if     (DISTANCE_CALC_METHOD==2)
				ndist[count++] = internalCalcOnline(i, j);
#else
				ndist[count++] = dist(i, j);
#endif
			}
		
			// sort neighbors to "numNN"th.
		sortNeighbors(neighbors, ndist, gNumCity-1, gNumNN);
		
		for( j=0; j<gNumNN; j++)
			NNI[i*gNumNN + j] = neighbors[j];
	}
	delete[] neighbors;
	delete[] ndist;
	trace("Quitting constructNormal\n");
}

void constructQuadNeighbors(int* NNI)
{
	int             i, j, k, t, count, dir;
	int             *nns[4],  size[4], alloc[4];
	double *ndist[4];
	trace("Entering constructQuadNeighbors\n");

	for( i=0; i<4; i++)
	{
		nns[i]   = new int[gNumCity-1];
		ndist[i] = new double[gNumCity-1];
	}
	
	for( i=0; i<gNumCity; i++)
	{
		for( j=0; j<4; j++)
		{
			alloc[j] = gNumNN/4;
			size[j] = 0;
		}
		for( j=0; j<gNumCity; j++)
			if( i != j)
			{
				dir  = (gNodeCoords[i].pt[0] < gNodeCoords[j].pt[0]) ? 0 : 1;
				dir += (gNodeCoords[i].pt[1] < gNodeCoords[j].pt[1]) ? 0 : 2;
				nns[dir][size[dir]] = j;
#if     (DISTANCE_CALC_METHOD==2)
				ndist[dir][size[dir]++] = internalCalcOnline(i, j);
#else
				ndist[dir][size[dir]++] = dist(i, j);
#endif
			}
			// If size[?] < alloc[?], adjust size[*] and alloc[*].
		for( j=0; j<4; j++)
		{
			if( size[j] < alloc[j])
				while( size[j] != alloc[j])
				{
					alloc[j]--;
					dir = alloc[j]%4;
					while( size[dir] <= alloc[dir]) dir = (dir+1)%4;
					alloc[dir]++;
				}
		}
		affirm( alloc[0]+alloc[1]+alloc[2]+alloc[3] == gNumNN);
		
			// sort neighbors to "alloc[i]"th. and combine...
		t = 0;
		for( j=0; j<4; j++)
		{
			sortNeighbors(nns[j], ndist[j], size[j], alloc[j]);
			for( k=0; k<alloc[j]; k++)
			{
				nns[0][t]     = nns[j][k];
				ndist[0][t++] = ndist[j][k];
			}
		}
		affirm( t== gNumNN);
		
		sortNeighbors(nns[0], ndist[0], gNumNN, gNumNN);
		for( j=0; j<gNumNN; j++)
			NNI[i*gNumNN + j] = nns[0][j];
	}
	
	for( i=0; i<4; i++)
	{
		delete[] nns[i]; delete[] ndist[i];
	}
	trace("Quitting constructQuadNeighbors\n");
}
void sortNeighbors(int* neighbors, double* ndist, int size, int numNN)
{
	int mini, j, k, t;
	int  left, right, prev_right;
	double pivot, doublet;
		//int count=0; float swap_ratio;
	
	affirm(size >= numNN);
	
	if( size == 0 || size == 1) return; // A sort is not need.
	
		// quick sort variants 
	pivot = 0; left = 0; prev_right = right = size;
	while( right > numNN)
	{
		// save right 
		prev_right = right--;
		
		// calc. pivot value
		pivot = (ndist[0]+ndist[right/4]+ndist[right/2]+
				 ndist[(right*3)/4]+ndist[right])/5;
		
		while( left < right)
		{
			while( ndist[left] < pivot - EPS) left++;
			while( ndist[right] >= pivot - EPS) right--;
			
			if( left < right)
			{
				t = neighbors[left]; neighbors[left] = neighbors[right];
				neighbors[right] = t;
				doublet = ndist[left]; ndist[left] = ndist[right];
				ndist[right] = doublet;
					//count++;
			}
		}
		left = 0; right++;
	}
	// end of quick sort
	
	// we have the interest only between 0 and prev_right.
	size = prev_right;
	
		/* selection sort..... */
	for( j=0; j<numNN; j++)
	{
		mini = j;
		for( k=j+1; k<size; k++)
			if( ndist[k] < ndist[mini] - EPS)
				mini = k;
		t = neighbors[mini]; neighbors[mini] = neighbors[j];
		neighbors[j] = t;
		doublet = ndist[mini]; ndist[mini] = ndist[j];
		ndist[j] = doublet;
			//count++;
	}
		//printf("prev_right = %4d, ratio= %f\n",
		//size, (swap_ratio=(float)count/(float)gNumCity));
}

	/*************************** Test Driver ******************************/
#ifdef TEST_DRIVER
/* Test Driver routine
   compiling : gcc tsplib_io.cc -DTEST_DRIVER -lm
   usage : a.out tsp-file
   ex) a.out lin105.tsp 
*/
int main(int argc, char* argv[])
{
	int i, cost, t;
	int* opt_tour, *neighbors;
	FILE* fd1, *fd2;
	char opt_file[256];

	if( argc <2) return 0;
	
	ReadTspFile(argv[1]);
	PrintTspInfo();


	// To test nearest neighborhood lists
	ConstructNN(20);

		// To test quadrant-neighbor lists..
		//ConstructNN(40, 1);
	neighbors = new int[gNumNN];

	printf("# of city = %d\n", gNumCity);
	printf("Known optimum cost = %d\n", GetOptimumCost());
	
	opt_tour = new int[gNumCity];
	if( ReadOptTourFile(opt_tour, gNumCity))
	{
			// Print optimal tour's cost
		printf(" Getting optimal tour cost.\n");
		cost = 0;
		for( i=0; i<gNumCity-1; i++)
			cost += dist(opt_tour[i], opt_tour[i+1]);
		cost += dist(opt_tour[gNumCity-1], opt_tour[0]);
		printf(" optimal cost is %d\n", cost);
		
			// To test nearest neighborhood lists
		puts("Entering TEST : neighborhood list");
		t = 200;
		printf("%f %f\n", gNodeCoords[t].pt[0], gNodeCoords[t].pt[1]);
		for( i=0; i<gNumNN; i++)
			neighbors[i] = nni(t, i);
		fd1 = fopen("nn_test.xy", "w"); fd2 = fopen("nn_test.nn_xy", "w");
		PrintCoords(fd1, opt_tour, gNumCity);
		PrintCoords(fd2, neighbors, gNumNN);
		fclose(fd1); fclose(fd2);
		
			// To test remapping cities.
		int temp_cost, city;
		
		RemapCities(opt_tour);
			// Now, optimal tour sequence is "0123456789.....".
		temp_cost = 0;
		for( i=0; i<gNumCity-1; i++)
			temp_cost += dist(i, i+1);
		temp_cost += dist(gNumCity-1, 0);
		printf("After remapped, optimal cost is %d\n", temp_cost);
		
			// After remapped, neighbor list test
		for( i=0; i<gNumCity; i++)
			if( GetOrgCity(i) == gNumCity/2)
				city = i;
		for( i=0; i<gNumNN; i++)
			neighbors[i] = nni(city, i);
		fd1 = fopen("remap_test.nn_xy", "w");
		PrintCoords(fd1, neighbors, gNumNN);
		fclose(fd1);
	}
	
	puts("end...................");
	return 0;
}
#endif // TEST_DRIVER
		

		


































