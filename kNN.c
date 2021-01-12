#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include <time.h>
#include "lib.h"
#include "VPT.h"
#include "timediff.h"

char MatName[128];

// Definition of the kNN result struct
typedef struct knnresult{
  int    * nidx;    //!< Indices (0-based) of nearest neighbors [m-by-k]
  double * ndist;   //!< Distance of nearest neighbors          [m-by-k]
  int      m;       //!< Number of query points                 [scalar]
  int      k;       //!< Number of nearest neighbors            [scalar]
} knnresult;

void ebeMultiply(const int n, const double *a, const double *x, double *y)
{
    

    blasint k = 0; // Just the diagonal; 0 super-diagonal bands
    static const double alpha = 1.0;
    blasint lda = 1;
    blasint incx = 1;
    static const double beta = 0.0;
    blasint incy = 1;

	cblas_dsbmv(CblasRowMajor,CblasLower,n,k,alpha,a,lda,x,incx,beta,y,incy);

}


void distanceMatrix(double * D, const double * X, const double * Y, const int n, const int m, const int d){

	double * X2 = malloc(n*d*sizeof(double));
	double * Y2 = malloc(m*d*sizeof(double));

	double *  e = malloc(d  *sizeof(double));

	double * XX = malloc(n  *sizeof(double));
	double * YY = malloc(m  *sizeof(double));


	ebeMultiply(n*d,X,X,X2);
	ebeMultiply(m*d,Y,Y,Y2);

	/* for(int i=0; i<n*d; i++){
		X2[i]=pow(X[i],2);
	}
	for(int i=0; i<m*d; i++){
		Y2[i]=pow(Y[i],2);
	} */

	for(int i=0;i<d;i++) e[i]=1;

	cblas_dgemv(CblasRowMajor,CblasNoTrans,n,d,1,X2,d,e,1,0,XX,1);
	cblas_dgemv(CblasRowMajor,CblasNoTrans,m,d,1,Y2,d,e,1,0,YY,1);
	
	free(X2);
	free(Y2);
	free(e);



	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,n,m,d,-2,X,d,Y,d,0,D,m);

	for(int i=0; i<n; i++)
		for(int j=0; j<m; j++){
			D[i*m + j] += XX[i] + YY[j];
			D[i*m + j] = sqrt(D[i*m + j]);
		}
	
	free(XX);
	free(YY);
	



}


/*!
	Initial kNN algorithm. 
 	Compute k nearest neighbors of each point in X [n-by-d]

  \param  X      Corpus data points              [n-by-d]
  \param  Y      Query data points               [m-by-d]
  \param  n      Number of corpus points         [scalar]
  \param  m      Number of query points          [scalar]
  \param  d      Number of dimensions            [scalar]
  \param  k      Number of neighbors             [scalar]

  \return  The kNN result
*/
knnresult kNN(const double * X, const double * Y, const int n, const int m, const int d, const int k){

    knnresult Results;
	Results.k=k;
	Results.m=m;
	Results.ndist=(double *)malloc(m*k*sizeof(double));
	Results.nidx =(int    *)malloc(m*k*sizeof(int));

	double D[m*n];

	int indices[n];

	for(int i=0; i<n; i++){
		indices[i]=i;
	}

	distanceMatrix(D,Y,X,m,n,d);

	printf("\n");
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++)
			printf("%3.0lf ",D[i*m+j]);
		printf("\n");
	}
	printf("\n");

	for(int Yindex=0; Yindex<m; Yindex++){
		int indicesTemp[n];
		memcpy(indicesTemp,indices,n*sizeof(int));

		kthSmallest(D+Yindex*n,indicesTemp,0,n-1,k);

		memcpy(Results.nidx +Yindex*k,indicesTemp,k*sizeof(int));
		memcpy(Results.ndist+Yindex*k,D+Yindex*n, k*sizeof(double));

	}

	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++)
			printf("%3.0lf ",D[i*m+j]);
		printf("\n");
	}
	printf("\n"); 

    return Results;
}

/*!
	kNN algorithm, splits query points into p blocks. indStart and indEnd are the start and end indices of the corpus set.
 	Compute k nearest neighbors of each point in X [n-by-d]

  \param  X      Corpus data points              [n-by-d]
  \param  Y      Query data points               [m-by-d]
  \param  n      Number of corpus points         [scalar]
  \param  m      Number of query points          [scalar]
  \param  d      Number of dimensions            [scalar]
  \param  k      Number of neighbors             [scalar]

  \return  The kNN result
*/
knnresult kNNblocked(const double * X, const double * Y, const int n, const int m, const int d, const int k, const int p, const int indStart, const int indEnd){

	//allocating space for the results
    knnresult Results;
	Results.k=k;
	Results.m=m;
	Results.ndist=(double *)malloc(m*k*sizeof(double));
	Results.nidx =(int    *)malloc(m*k*sizeof(int));

	//initialising results to max values
	for(int i=0; i<m*k; i++){
		Results.ndist[i]=999;
		Results.nidx[i]=999;
	}
	
	//computing block sizes
	int chunk=floor(m/p);
	int extra=m%p;
	
	double *D=(double *)malloc(n*(chunk+1)*sizeof(double));

	//---------------------------------------------------------------------------------------------------------------------
	// computing block parameters, scounts[] being the size of each block and displs[] their starting index in Y.
	int scounts[p];
	int displs[p+1];
	
	displs[0]=0;
	displs[p]=m;

	for(int i=0; i<p; i++){
		scounts[i]=chunk;
	}
	for(int i=0; i<extra; i++){
		scounts[i]+=1;
	}
	for(int i=1; i<p; i++){
		displs[i]=displs[i-1]+scounts[i-1];
	}
	//---------------------------------------------------------------------------------------------------------------------

	// creating the index array of the corpus set ONCE and then copying it to a temp array every time it's needed
	int indices[n];
	for(int i= 0,j=indStart; j<indEnd; i++, j++){
		indices[i]=j;
	}

	int indicesTemp[n];

	//---------------------------------------------------------------------------------------------------------------------
	// main algorithm:
	// if k<n then we compute the distance matrix, find the k nearest for every point in the current block, copy the k nearest into resutls and repeat for each block
	if(k<n){
		for(int times=0; times<p; times++){

			distanceMatrix(D,Y+displs[times]*d,X,scounts[times],n,d);

			for(int Yindex=displs[times], counter=0; Yindex<displs[times+1]; Yindex++, counter++){
				memcpy(indicesTemp,indices,n*sizeof(int));

				kthSmallest(D+counter*n,indicesTemp,0,n-1,k);

				memcpy(Results.nidx +Yindex*k, indicesTemp, k*sizeof(int));
				memcpy(Results.ndist+Yindex*k, D+counter*n, k*sizeof(double));
			}

		}
	}else{ 		// if k>=n then we dont need to compute the k nearest, it is all of them so we just copy them
		for(int times=0; times<p; times++){

			distanceMatrix(D,Y+displs[times]*d,X,scounts[times],n,d);

			for(int Yindex=displs[times], counter=0; Yindex<displs[times+1]; Yindex++, counter++){
				memcpy(indicesTemp,indices,n*sizeof(int));

				memcpy(Results.nidx +Yindex*k, indicesTemp, n*sizeof(int));
				memcpy(Results.ndist+Yindex*k, D+counter*n, n*sizeof(double));
			}

		}
	}
	free(D);

    return Results;
}

// make list of pknnresult structs(queries)
Pknnresult** make_list_local_queries(const double *X, const int n, const int d, const int k){
	Pknnresult **queryList=(Pknnresult**)malloc(n*sizeof(Pknnresult*));
	for(int i=0; i<n; i++){
		queryList[i] = create_pknnresult(d,k);
		queryList[i] = initialize_pknnresult(queryList[i], X+i*d , d, k);
	}
	return queryList;
}

// distributed kNN algorithm, every process should run this. nblocks is the number of blocks each process must split its data into for kNNblocked
knnresult distrAllkNN(double *X, const int n, const int d, const int k, const int nblocks){

	struct timespec tSearchS,tSearchE,tWaitS,tWaitE,tSearchS2,tSearchE2;
	double totalWait=0,totalSearch=0,totalSearch2=0;

	// basic mpi IO variable initialization
	int numtasks,rank,prev,next;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	prev = rank-1;
   	next = rank+1;
   	if (rank == 0)  prev = numtasks - 1;
   	if (rank == (numtasks - 1))  next = 0;

	// compute size of data for each process
	int p=numtasks;
	int baseChunk = floor(n/p), extra = n%p;


	// allocate memory for process data (static, send and receive)
	double *xp=(double *)malloc((baseChunk+1)*d*sizeof(double));
	double *yp=(double *)malloc((baseChunk+1)*d*sizeof(double));
	double *zp=(double *)malloc((baseChunk+1)*d*sizeof(double));

	// initialise extra values in xp, may be overwriten if needed
	for(int i=0; i<d; i++){
		xp[(baseChunk+1)*d -1 -i]=999;
	}


	//-----------------------------------------------------------------------------------------------------------------------------------------------
	// compute size and displacement for each process's data, only process 0 needs this for scatter but everyone needs it for later computations.
	int  *displs=(int*)malloc((p+1)*sizeof(int));
	int *scounts=(int*)malloc(p*sizeof(int));

	displs[0]=0;
	for(int i=0; i<p; i++){
		scounts[i]=baseChunk*d;
	}
	for(int i=0; i<extra; i++){
		scounts[i]+=d;
	}
	for(int i=1; i<p+1; i++){
		displs[i]=displs[i-1]+scounts[i-1];
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------------


	// process 0 scatters the data and the rest receive
	if(rank==0){


		MPI_Scatterv(X,scounts,displs,MPI_DOUBLE,xp,(baseChunk+1)*d,MPI_DOUBLE,0,MPI_COMM_WORLD);

	}else{
		MPI_Scatterv(NULL,NULL,NULL,  MPI_DOUBLE,xp,(baseChunk+1)*d,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}


	// first time each process runs compares the data with its own
	memcpy(yp,xp,(baseChunk+1)*d*sizeof(double));


	// mpi variables for sending and receiving
	MPI_Request reqs[2];
	MPI_Status stats[2];

	// final result of each process space allocation
	knnresult pkNN;
	pkNN.k=k;
	pkNN.m=(baseChunk+1);
	pkNN.ndist=(double *)malloc((baseChunk+1)*k*sizeof(double));
	pkNN.nidx= (int *)   malloc((baseChunk+1)*k*sizeof(int));
	
	// value initialization
	for(int i=0; i<(baseChunk+1)*k; i++){
		pkNN.nidx[i]=999;
		pkNN.ndist[i]=999;
	}

	// temporary result, needed for computation each time a process compares new data. Space gets allocated for it inside kNNblocked.
	knnresult pkNNtemp;
	pkNNtemp.k=k;
	pkNNtemp.m=(baseChunk+1);

	// temporary arrays for k nearest neighbors of final and temp in order for them to get compared		
	double 	*tempD=(double *)malloc(2*k*sizeof(double));
	int 	*tempIDX=(int *) malloc(2*k*sizeof(int));

	//! \param counter contains the id of whose data is contained in yp[]
	int counter=rank;

	//! \param xsize contains the block size of the static data of each process
	int xsize;

	//! \param ysize contans the block size of yp[], essentially the xsize of the data we have received each time
	int ysize;

	// if our rank is less than the extra points we have xsize should be 1 bigger to split the points as equally as possible
	if(rank<extra){
		xsize=(baseChunk+1);
	}else{
		xsize=baseChunk;
	}

	//printf("extra=%d, xsize=%d \n",extra,xsize);

	// start and end of indices of our data
	int indStart;
	int indEnd;
	
	
	//-----------------------------------------------------------------------------------------------------------------------------------------------
	// main algorithm:
	// each process compares data with each other process by sending the data in a ring and keeping the kNN each time 
	for(int i=0; i<p; i++){

		// non blocking send and receives so we can have IO while we compute
		MPI_Isend(yp,(baseChunk+1)*d,MPI_DOUBLE,next,1,MPI_COMM_WORLD,&reqs[0]);
		MPI_Irecv(zp,(baseChunk+1)*d,MPI_DOUBLE,prev,1,MPI_COMM_WORLD,&reqs[1]);

		// do the same as in xsize
		if(counter<extra){
			ysize=(baseChunk+1);
		}else{
			ysize=baseChunk;
		}
		
		// compute start and end indices
		indStart=displs[counter]/d;
		indEnd=displs[counter+1]/d;

		clock_gettime(CLOCK_MONOTONIC, &tSearchS);
		// compute the kNN of our data
		pkNNtemp = kNNblocked(yp,xp,ysize,xsize,d,k,nblocks,indStart,indEnd);
		
		//pkNNtemp=kNN(yp,xp,ysize,xsize,d,k); <- non blocked kNN, maybe faster but needs more memory
		clock_gettime(CLOCK_MONOTONIC, &tSearchE);

		totalSearch += timeConv(diff(tSearchS,tSearchE));

		// counter computation and loop arround when needed
		counter--;
		if(counter==-1) counter=p-1;

	

		clock_gettime(CLOCK_MONOTONIC, &tSearchS2);
		// comparison of procesess's final kNN and the temporary kNN, copy both into temporary arrays and find the k smalles of the bunch
		for(int i=0; i<xsize; i++){
			memcpy(&tempD[0],   &pkNN.ndist[i*k],     k*sizeof(double));
			memcpy(&tempD[k],   &pkNNtemp.ndist[i*k], k*sizeof(double));
			memcpy(&tempIDX[0], &pkNN.nidx[i*k],      k*sizeof(int));
			memcpy(&tempIDX[k], &pkNNtemp.nidx[i*k],  k*sizeof(int));

			kthSmallest(tempD, tempIDX, 0, 2*k-1, k);

			memcpy(&pkNN.ndist[i*k], tempD,   k*sizeof(double));
			memcpy(&pkNN.nidx[i*k],  tempIDX, k*sizeof(int));
		}
		clock_gettime(CLOCK_MONOTONIC, &tSearchE2);

		totalSearch2 += timeConv(diff(tSearchS2,tSearchE2));


		clock_gettime(CLOCK_MONOTONIC, &tWaitS);
		// make sure we have finished our data transfer and switch yp[] with zp[]
		MPI_Waitall(2,reqs,stats);
		clock_gettime(CLOCK_MONOTONIC, &tWaitE);

		totalWait += timeConv(diff(tWaitS,tWaitE));
		double *temp=yp;
		yp=zp;
		zp=temp;	
	}
	free(pkNNtemp.ndist);
	free(pkNNtemp.nidx);
	free(tempD);
	free(tempIDX);

	// FINAL FINAL result computed by process 0 from the data of all processes
	knnresult pkNNFinal;
	int rcounts[p];
	int disps[p];
	if(rank==0){
		
		// memory allocation
		pkNNFinal.m=n;
		pkNNFinal.k=k;
		pkNNFinal.nidx=(int *)malloc(n*k*sizeof(int));
		pkNNFinal.ndist=(double *)malloc(n*k*sizeof(double));

		
		disps[0]=0;
		for(int i=0; i<p; i++){
				rcounts[i]=baseChunk*k;
		}
		for(int i=0; i<extra; i++){
				rcounts[i]+=k;
		}
		for(int i=1; i<p; i++){
			disps[i]=disps[i-1]+rcounts[i-1];
		}

	}
	// gather the data in a single pkNN struct in process 0 and return
	MPI_Gatherv(pkNN.ndist,xsize*k,MPI_DOUBLE,pkNNFinal.ndist,rcounts,disps,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gatherv(pkNN.nidx, xsize*k,MPI_INT,   pkNNFinal.nidx, rcounts,disps,MPI_INT,   0,MPI_COMM_WORLD);
	free(pkNN.ndist);
	free(pkNN.nidx);

	
	printf("%s,DIST,0,%lf,%lf,%lf,%lf,%lf,%d,%d,%d,%d,%d\n",MatName,totalSearch,totalSearch2,totalSearch+totalSearch2,totalWait,totalSearch+totalSearch2+totalWait,n,d,k,numtasks,rank);

	return pkNNFinal;
}

// distr serial vpt
knnresult *distrAllkNN_serial_vpt(const double *X, const int n, const int d, const int k, const int B){

	struct timespec tStart1,tStart2,tEnd1,tEnd2,tSearchS,tSearchE,tWaitS,tWaitE;
	double totalWait=0,totalSearch=0;

	// basic mpi IO variable initialization
	int numtasks,rank,prev,next;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	prev = rank-1;
   	next = rank+1;
   	if (rank == 0)  prev = numtasks - 1;
   	if (rank == (numtasks - 1))  next = 0;

	// compute size of data for each process
	const int p = numtasks;
	const int baseChunk = floor(n/p), extra = n%p;


	// allocate memory for process data (static, send and receive)
	double *xp = malloc((baseChunk+1)*d*sizeof(double));

	// initialise extra values in xp, may be overwriten if needed
	for(int i=0; i<d; i++){
		xp[(baseChunk+1)*d -1 -i]=999;
	}


	//-----------------------------------------------------------------------------------------------------------------------------------------------
	// compute size and displacement for each process's data, only process 0 needs this for scatter but everyone needs it for later computations.
	int  *displs = malloc((p+1)*sizeof(int));
	int *scounts = malloc( p*   sizeof(int));

	displs[0]=0;
	for(int i=0; i<p; i++){
		scounts[i]=baseChunk*d;
	}
	for(int i=0; i<extra; i++){
		scounts[i]+=d;
	}
	for(int i=1; i<p+1; i++){
		displs[i]=displs[i-1]+scounts[i-1];
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------------


	// process 0 scatters the data and the rest receive
	if(rank==0){
		MPI_Scatterv(X,scounts,displs,MPI_DOUBLE,xp,(baseChunk+1)*d,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}else{
		MPI_Scatterv(NULL,NULL,NULL,  MPI_DOUBLE,xp,(baseChunk+1)*d,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	int counter=rank;
	int xsize;
	int ysize;

	// if our rank is less than the extra points we have xsize should be 1 bigger to split the points as equally as possible
	if(rank<extra){
		xsize=(baseChunk+1);
	}else{
		xsize=baseChunk;
	}

	int indStart = displs[rank  ]/d;
	int indEnd   = displs[rank+1]/d;

	int *tempInd = malloc(xsize*  sizeof(int));

	for(int i=0, j=indStart; j<indEnd; j++, i++){
		tempInd[i]=j;
	}

	// make the local queries
	Pknnresult** localQueryList = make_list_local_queries(xp, xsize, d, k);

	// make the tree which we will later send to other nodes
	VPtree *localTree;
	localTree =    initialize_vpt(localTree, d, (baseChunk+1));


	clock_gettime(CLOCK_MONOTONIC, &tStart1);
	make_serial_bucket_vpt_select(localTree, xp, tempInd, 0, d, xsize, B);
	clock_gettime(CLOCK_MONOTONIC, &tEnd1);

	free(tempInd);

	free(xp);
	//printVPtree_serial(localTree,0,0,5);

	// reserve space for the tree to be received
	VPtree *recTree;
	recTree = initialize_vpt(recTree, d, (baseChunk+1));

	// mpi variables for sending and receiving
	MPI_Request reqs[10];
	MPI_Status stats[10];

	// tags for sending and receiving different parts of the tree
	const int cordTag=1, indTag=2, LptrTag=3, RptrTag=4, muTag=5;
	
	clock_gettime(CLOCK_MONOTONIC, &tStart2);
	//-----------------------------------------------------------------------------------------------------------------------------------------------
	// main algorithm:
	// each process compares data with each other process by sending the data in a ring and keeping the kNN each time 
	for(int i=0; i<p; i++){

		if(counter<extra){
			ysize=(baseChunk+1);
		}else{
			ysize=baseChunk;
		}

		MPI_Isend(localTree->cordsP,   ysize*d,         MPI_DOUBLE, next, cordTag, MPI_COMM_WORLD, &reqs[0]);
		MPI_Irecv(  recTree->cordsP,   (baseChunk+1)*d, MPI_DOUBLE, prev, cordTag, MPI_COMM_WORLD, &reqs[1]);

		MPI_Isend(localTree->indP,     ysize,           MPI_INT,    next, indTag,  MPI_COMM_WORLD, &reqs[2]);
		MPI_Irecv(  recTree->indP,     (baseChunk+1),   MPI_INT,    prev, indTag,  MPI_COMM_WORLD, &reqs[3]);

		MPI_Isend(localTree->mu,       ysize,           MPI_DOUBLE, next, muTag,   MPI_COMM_WORLD, &reqs[4]);
		MPI_Irecv(  recTree->mu,       (baseChunk+1),   MPI_DOUBLE, prev, muTag,   MPI_COMM_WORLD, &reqs[5]);

		MPI_Isend(localTree->leftptr,  ysize,           MPI_INT,    next, LptrTag, MPI_COMM_WORLD, &reqs[6]);
		MPI_Irecv(  recTree->leftptr,  (baseChunk+1),   MPI_INT,    prev, LptrTag, MPI_COMM_WORLD, &reqs[7]);

		MPI_Isend(localTree->rightptr, ysize,           MPI_INT,    next, RptrTag, MPI_COMM_WORLD, &reqs[8]);
		MPI_Irecv(  recTree->rightptr, (baseChunk+1),   MPI_INT,    prev, RptrTag, MPI_COMM_WORLD, &reqs[9]);

		clock_gettime(CLOCK_MONOTONIC, &tSearchS);
		// update our local query list
		for(int j=0; j<xsize; j++){
			localQueryList[j] = search_serial_bucket_vpt(localTree, 0, localQueryList[j], d, k);
		}
		clock_gettime(CLOCK_MONOTONIC, &tSearchE);

		totalSearch += timeConv(diff(tSearchS,tSearchE));

		counter--;
		if(counter==-1) counter=p-1;

		/* printVPtree_serial(localTree,0,0,5);

		printf("\n---------------\np=%d:\n",p);
		for(int j=0; j<xsize; j++){
			for(int i=0; i<k+1; i++){
				printf("[%d: %.3lf]",localQueryList[j]->idx[i],localQueryList[j]->dist[i]);
			}
			printf("\n");
		} */
		clock_gettime(CLOCK_MONOTONIC, &tWaitS);
		// make sure we have finished our data transfer and switch yp[] with zp[]
		MPI_Waitall(10,reqs,stats);
		clock_gettime(CLOCK_MONOTONIC, &tWaitE);

		totalWait += timeConv(diff(tWaitS,tWaitE));

		VPtree *temp = recTree;
		     recTree = localTree;
		   localTree = temp;	
	}
	clock_gettime(CLOCK_MONOTONIC, &tEnd2);

	//printVPtree_serial(localTree,0,0,5);
	free_searial_vpt(localTree);
	free_searial_vpt(recTree);
	free(displs);
	free(scounts);


	// final result of each process space allocation
	knnresult *pkNNlocal=malloc(sizeof(knnresult));
	pkNNlocal->k=k;
	pkNNlocal->m=xsize;
	pkNNlocal->ndist=(double *)malloc((baseChunk+1)*k*sizeof(double));
	pkNNlocal->nidx= (int *)   malloc((baseChunk+1)*k*sizeof(int));

	for(int i=0; i<xsize; i++){
		memcpy(pkNNlocal->ndist + i*k, localQueryList[i]->dist, k*sizeof(double));
		memcpy(pkNNlocal->nidx  + i*k, localQueryList[i]->idx,  k*sizeof(int));
	}

	for(int i=0; i<xsize; i++){
		free_pknnresult(localQueryList[i]);
	}
	free(localQueryList);


	// FINAL FINAL result gathered in process 0 from the data of all processes
	knnresult *pkNNFinal=malloc(sizeof(knnresult));

	int rcounts[p];
	int disps[p];
	if(rank==0){
		
		// memory allocation
		pkNNFinal->m=n;
		pkNNFinal->k=k;
		pkNNFinal->nidx =(int *)   malloc(n*k*sizeof(int));
		pkNNFinal->ndist=(double *)malloc(n*k*sizeof(double));

		
		disps[0]=0;
		for(int i=0; i<p; i++){
				rcounts[i]=baseChunk*k;
		}
		for(int i=0; i<extra; i++){
				rcounts[i]+=k;
		}
		for(int i=1; i<p; i++){
			disps[i]=disps[i-1]+rcounts[i-1];
		}

	}
	// gather the data in a single pkNN struct in process 0 and return
	MPI_Gatherv(pkNNlocal->ndist,xsize*k,MPI_DOUBLE,pkNNFinal->ndist,rcounts,disps,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gatherv(pkNNlocal->nidx, xsize*k,MPI_INT,   pkNNFinal->nidx, rcounts,disps,MPI_INT,   0,MPI_COMM_WORLD);
	free(pkNNlocal->ndist);
	free(pkNNlocal->nidx);
	free(pkNNlocal);

	double tBuild;

	tBuild = timeConv(diff(tStart1,tEnd1));

	printf("%s,VPT,%lf,%lf,0,0,%lf,%lf,%d,%d,%d,%d,%d\n",MatName,tBuild,totalSearch,totalWait,tBuild+totalSearch+totalWait,n,d,k,numtasks,rank);
	//printf("%s,VPT,%lf,%lf,0,0,%lf,%lf,%d,%d,%d,%d,%d,%d\n",MatName,tBuild,totalSearch,totalWait,tBuild+totalSearch+totalWait,n,d,k,numtasks,rank,BB);

	return pkNNFinal;
}

//testing main
int main1(int argc, char *argv[])
{

	int numtasks, rank;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	//srand(time(NULL));

	//getchar();

	const int n=30000/* 20000*sqrt(numtasks) */,d=20,k=10;
	const int B=0.02*n/numtasks;
	
	const int vpt=atoi(argv[1]);

	int nblocks;
	if(vpt==0) nblocks=atoi(argv[2]);

	knnresult *result;
	knnresult simple;

	struct timespec tStart,tEnd;

	if(rank==0){
		if(vpt==1){
			printf("n=%d, B=%d\n",n,B);
		}else{
			printf("n=%d, nblocks=%d\n",n,nblocks);
		}

		double *X = malloc(n*d*sizeof(double));
		for(int i=0; i<n; i++){
			for(int j=0; j<d; j++){
				X[i*d+j]=(double)rand()/RAND_MAX;
			}
		}

		clock_gettime(CLOCK_MONOTONIC, &tStart);

		if(vpt==1){
			result = distrAllkNN_serial_vpt(X,n,d,k,B);
		}else{
			simple = distrAllkNN(X,n,d,k,nblocks);
		}

		
		clock_gettime(CLOCK_MONOTONIC, &tEnd);

		free(X);
	}else{
		clock_gettime(CLOCK_MONOTONIC, &tStart);
		
		if(vpt==1){
			result = distrAllkNN_serial_vpt(NULL,n,d,k,B);
			free(result);
		}
		else
			simple = distrAllkNN(NULL,n,d,k,nblocks);

		clock_gettime(CLOCK_MONOTONIC, &tEnd);

	}


	if(rank==0){
		if(vpt==1){
			/* printf("\n-----VPT-----\n");
			for(int i=0; i<n; i++){
				for(int j=0;j<k;j++)
					printf("[%2d: %3.5lf]",result->nidx[i*k+j],result->ndist[i*k+j]);
				printf("\n");
			}
			printf("\n-----VPTend-----\n"); */

			free(result->ndist);
			free(result->nidx);
			free(result);
		}else
		{
			/* printf("\nhi-----\n");
			for(int i=0; i<n; i++){
				for(int j=0;j<k;j++)
					printf("[%2d: %3.20lf]",simple.nidx[i*k+j],simple.ndist[i*k+j]);
				printf("\n");
			}
			printf("\nbye-----\n"); */
			free(simple.ndist);
			free(simple.nidx);
		}
		
		
	}

	struct timespec tResult = diff(tStart,tEnd);

	printf("node: %d took %lf seconds\n",rank,timeConv(tResult));

	MPI_Finalize();
    return 0;
}

//main main 
int main(int argc, char *argv[])
{

	int numtasks, rank;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	srand(time(NULL));

	//getchar();
	if(argc!=8){
		printf("Usage: ./kNN matFile (corel->1, tv->2, PID->3, feature->4) n d k (0 for distance/1 for VPT) (nblocks/B*0.0001 of n)\n");
		exit(1);
	}

	strcpy(MatName,argv[1]);

	const int numbered=atoi(argv[2]);

	const int n=atoi(argv[3]),d=atoi(argv[4]),k=atoi(argv[5]);
	const int B=atoi(argv[7]);
	
	const int vpt=atoi(argv[6]);

	int nblocks;
	if(vpt==0) nblocks=atoi(argv[7]);

	knnresult *result;
	knnresult simple;


	if(rank==0){
		/* if(vpt==1){
			printf("n=%d, d=%d, k=%d, B=%d\n",n,d,k,B);
		}else{
			printf("n=%d, d=%d, k=%d, nblocks=%d\n",n,d,k,nblocks);
		} */

		double *X = malloc(n*d*sizeof(double));
		FILE *matFile=fopen(argv[1],"r");

		if (matFile == NULL)
		{
			printf("Couldnt open file\n");
			exit(1);
		}

		double num,temp;
		int i=0,j=0;
		if(numbered==1){
			while ( fscanf(matFile,"%lf",&num) !=EOF )
			{
				if(i%(d+1)!=0){
					X[j]=num;
					j++;
				}
				i++;
			}
		}else if(numbered==2){
			for(i=0; i<n; i++){
				fscanf(matFile,"%lf",&num);
				for(j=0; j<d; j++){
					if(fscanf(matFile," %lf:%lf",&temp,&num)==EOF) break;
					X[i*d+j]=num;
        		}
				fscanf(matFile,"%*[^\n]\n");
			}
		}else if(numbered==3){
			fscanf(matFile,"%*[^\n]\n");
			while ( fscanf(matFile," %lf",&num) !=EOF )
			{
				X[j]=num;
				j++;
				
			}
		}else{
			for(int skip=0;skip<4;skip++){
				fscanf(matFile,"%*[^\n]\n");
			}
			for(i=0; i<n; i++){
				fscanf(matFile,"%lf",&num);
				for(j=0; j<d; j++){
					if(fscanf(matFile,",%lf",&num)==EOF) break;
					X[i*d+j]=num;
        		}
				fscanf(matFile,"%*[^\n]\n");
			}
		}
    
		
		fclose(matFile);

		//clock_gettime(CLOCK_MONOTONIC, &tStart);

		if(vpt==1){
			result = distrAllkNN_serial_vpt(X,n,d,k,B);
		}else{
			simple = distrAllkNN(X,n,d,k,nblocks);
		}

		
		//clock_gettime(CLOCK_MONOTONIC, &tEnd);

		free(X);
	}else{
		//clock_gettime(CLOCK_MONOTONIC, &tStart);
		
		if(vpt==1){
			result = distrAllkNN_serial_vpt(NULL,n,d,k,B);
			free(result);
		}
		else
			simple = distrAllkNN(NULL,n,d,k,nblocks);

		//clock_gettime(CLOCK_MONOTONIC, &tEnd);

	}


	if(rank==0){
		if(vpt==1){
			/* printf("\n-----VPT-----\n");
			for(int i=0; i<n; i++){
				for(int j=0;j<k;j++)
					printf("[%2d: %3.5lf]",result->nidx[i*k+j],result->ndist[i*k+j]);
				printf("\n");
			}
			printf("\n-----VPTend-----\n"); */

			free(result->ndist);
			free(result->nidx);
			free(result);
		}else
		{
			/* printf("\nhi-----\n");
			for(int i=0; i<n; i++){
				for(int j=0;j<k;j++)
					printf("[%2d: %3.20lf]",simple.nidx[i*k+j],simple.ndist[i*k+j]);
				printf("\n");
			}
			printf("\nbye-----\n"); */
			free(simple.ndist);
			free(simple.nidx);
		}
		
		
	}


	//printf("node: %d took %lf seconds\n",rank,timeConv(tResult));

	MPI_Finalize();
    return 0;
}

//testing main
int main3(int argc, char *argv[])
{

	int numtasks, rank;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	strcpy(MatName,"None");
	srand(time(NULL));

	//getchar();

	const int n=atoi(argv[2]),d=atoi(argv[3]),k=atoi(argv[4]);
	const int B=0.0002*n/numtasks;
	
	const int vpt=atoi(argv[1]);

	int nblocks;
	if(vpt==0) nblocks=atoi(argv[5]);

	knnresult *result;
	knnresult simple;

	struct timespec tStart,tEnd;

	if(rank==0){
		

		double *X = malloc(n*d*sizeof(double));
		for(int i=0; i<n; i++){
			for(int j=0; j<d; j++){
				X[i*d+j]=(double)rand()/RAND_MAX;
			}
		}

		

		if(vpt==1){
			result = distrAllkNN_serial_vpt(X,n,d,k,B);
		}else{
			simple = distrAllkNN(X,n,d,k,nblocks);
		}

		
		

		free(X);
	}else{
		
		
		if(vpt==1){
			result = distrAllkNN_serial_vpt(NULL,n,d,k,B);
			free(result);
		}
		else
			simple = distrAllkNN(NULL,n,d,k,nblocks);

		

	}


	if(rank==0){
		if(vpt==1){
			/* printf("\n-----VPT-----\n");
			for(int i=0; i<n; i++){
				for(int j=0;j<k;j++)
					printf("[%2d: %3.5lf]",result->nidx[i*k+j],result->ndist[i*k+j]);
				printf("\n");
			}
			printf("\n-----VPTend-----\n"); */

			free(result->ndist);
			free(result->nidx);
			free(result);
		}else
		{
			/* printf("\nhi-----\n");
			for(int i=0; i<n; i++){
				for(int j=0;j<k;j++)
					printf("[%2d: %3.20lf]",simple.nidx[i*k+j],simple.ndist[i*k+j]);
				printf("\n");
			}
			printf("\nbye-----\n"); */
			free(simple.ndist);
			free(simple.nidx);
		}
		
		
	}

	struct timespec tResult = diff(tStart,tEnd);

	//printf("node: %d took %lf seconds\n",rank,timeConv(tResult));

	MPI_Finalize();
    return 0;
}
