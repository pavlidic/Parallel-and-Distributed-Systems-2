#include "VPT.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
//#include <cblas.h>
#include "lib.h"
#include "insert.h"

// struct for the nodes of the VP tree



// calculaates distance between 2 d-dimensional points
double distCalc(const double *x, const double *y, const int d){
    double dist=0;
    for(int i=0; i<d; i++){
        dist+=pow(x[i]-y[i],2);
    }
    return sqrt(dist);
}

/* double distCalc(const double *x, const double *y, double *z,const int d){
    memcpy(z,y,d*sizeof(double));

    cblas_daxpy(d,-1,x,1,z,1);
    
    return cblas_dnrm2(d,z,1);
} */




// returns which element of X has the best spread for a randomly chosen sample
int select_vp(const double *X, const int d, const int n){
    double temp;
    double spread=0;
    double best_spread=0;
    int    best_spread_int=0;
    int    rand_ind;
    

    const int sample_size=(int)floor(log2(n))+1;

    double *P    = malloc(sample_size*d*sizeof(double));
    int    *Pind = malloc(sample_size*  sizeof(int));
    double *Q    = malloc(sample_size*d*sizeof(double));
    double *dist = malloc(sample_size*  sizeof(double));

    for(int i=0; i<sample_size; i++){

        rand_ind=rand()%n;
        Pind[i]=rand_ind;
        memcpy(P+i*d, X+rand_ind*d, d*sizeof(double));

    }

    for(int i=0; i<sample_size; i++){

        for(int j=0; j<sample_size; j++){
            dist[j]=distCalc(P+i*d, X+(rand()%n)*d, d);
        }

        kthSmallestNoInd(dist, 0, sample_size-1, sample_size/2+1 );

        const double mu = dist[sample_size/2+1];

        for(int i=0; i<sample_size; i++){
            temp=dist[i]-mu;
            spread+=temp*temp;
        }
        spread/=sample_size;

        if(spread > best_spread){
            best_spread = spread;
            best_spread_int = Pind[i];
        }

        spread=0;
    }

    free(P);
    free(Q);
    free(Pind);
    free(dist);

    return best_spread_int;

}

// ALL VP TREE MAKERS ARE DESTRICTIVE (for speed and simplicity), parse a COPY of X if you want to retain the element's order
// makes a vp tree and return pointer to the head
VPnode* make_vp_tree(double *X, int *ind, const int d, const int n){
    if(n==1){
        VPnode *end=(VPnode*)malloc(sizeof(VPnode));

        end->isBucket=0;
        end->mu=0;
        end->indP=ind[0];
        end->P=(double*)malloc(d*sizeof(double));
        memcpy(end->P,X,d*sizeof(double));

        end->right=NULL;
        end->left=NULL;

        return end;

    }else if(n==0){
        return NULL;
    }

    double *LX=(double*)malloc(n*d*sizeof(double));
    double *RX=(double*)malloc(n*d*sizeof(double));
    int  *Lind=(int*)   malloc(n*sizeof(int));
    int  *Rind=(int*)   malloc(n*sizeof(int));
    
    double *dist=(double*)malloc(n*sizeof(double));
    double *distTemp=(double*)malloc(n*sizeof(double));
    dist[0]=0;
    double temp,tempSq;

    VPnode *node = (VPnode*)malloc(sizeof(VPnode));
    
    node->isBucket=0;
    node->indP=ind[0];
    node->P=(double*)malloc(d*sizeof(double));
    memcpy(node->P,X,d*sizeof(double));


    for(int i=1; i<n; i++){
        tempSq=0;
        for(int j=0; j<d; j++){
            temp=(X[0*d+j] - X[i*d + j]);
            tempSq+=temp*temp;
        }
        dist[i]=sqrt(tempSq);
    }

    memcpy(distTemp,dist,n*sizeof(double));

    kthSmallestNoInd(distTemp,1,n-1,(n-1)/2);
    node->mu=distTemp[(n-1)/2+1];
    free(distTemp);

    int Ln=0, Rn=0;

    for(int i=0; i<n-1; i++){
        if(dist[i+1] < node->mu){
            memcpy(&LX[Ln*d], &X[(i+1)*d], d*sizeof(double));
            Lind[Ln]=ind[i+1];
            Ln++;
        }else{
            memcpy(&RX[Rn*d], &X[(i+1)*d], d*sizeof(double));
            Rind[Rn]=ind[i+1];
            Rn++;
        }
    }
    free(dist);

    LX=(double*)realloc(LX,   Ln*d*sizeof(double));
    Lind= (int*)realloc(Lind, Ln*  sizeof(int));
    RX=(double*)realloc(RX,   Rn*d*sizeof(double));
    Rind= (int*)realloc(Rind, Rn*  sizeof(int));

    node->left=  make_vp_tree(LX,Lind,d,Ln);
    free(LX);
    free(Lind);
    node->right= make_vp_tree(RX,Rind,d,Rn);
    free(RX);
    free(Rind);

    return node;
}

// ALL VP TREE MAKERS ARE DESTRICTIVE (for speed and simplicity), parse a COPY of X if you want to retain the element's order
// makes a vp tree with a statistically good selection for VP and return pointer to the head
VPnode* make_vp_tree_select(double *X, int *ind, const int d, const int n){
    if(n==1){
        VPnode *end=(VPnode*)malloc(sizeof(VPnode));

        end->isBucket=0;
        end->mu=0;
        end->indP=ind[0];
        end->P=(double*)malloc(d*sizeof(double));
        memcpy(end->P,X,d*sizeof(double));

        end->right=NULL;
        end->left=NULL;

        return end;

    }else if(n==0){
        return NULL;
    }

    double *LX=(double*)malloc(n*d*sizeof(double));
    double *RX=(double*)malloc(n*d*sizeof(double));
    int  *Lind=(int*)   malloc(n*sizeof(int));
    int  *Rind=(int*)   malloc(n*sizeof(int));
    
    double *dist=(double*)malloc(n*sizeof(double));
    double *distTemp=(double*)malloc(n*sizeof(double));
    dist[0]=0;
    double temp,tempSq;

    VPnode *node = (VPnode*)malloc(sizeof(VPnode));
    

    //----------------------------------------------------------------
    // using select vp to find a suitable vp and swap its data
    // with the first element for ease of use later
    int vp=select_vp(X,d,n);
    int tempInd;
    double *tempX=(double*)malloc(d*sizeof(double));

    // move X[vp]
    memcpy(tempX,  X+vp*d, d*sizeof(double));
    memcpy(X+vp*d, X+0*d,  d*sizeof(double));
    memcpy(X+0*d,  tempX,  d*sizeof(double));
    free(tempX);
    // move its index
    tempInd=ind[vp];
    ind[vp]=ind[0];
    ind[0]=tempInd;

    //----------------------------------------------------------------


    node->isBucket=0;
    node->indP=ind[0];
    node->P=(double*)malloc(d*sizeof(double));
    memcpy(node->P,X,d*sizeof(double));


    for(int i=1; i<n; i++){
        tempSq=0;
        for(int j=0; j<d; j++){
            temp=(X[0*d+j] - X[i*d + j]);
            tempSq+=temp*temp;
        }
        dist[i]=sqrt(tempSq);
    }

    memcpy(distTemp,dist,n*sizeof(double));

    kthSmallestNoInd(distTemp,1,n-1,(n-1)/2);
    node->mu=distTemp[(n-1)/2+1];
    free(distTemp);

    int Ln=0, Rn=0;

    for(int i=0; i<n-1; i++){
        if(dist[i+1] < node->mu){
            memcpy(&LX[Ln*d], &X[(i+1)*d], d*sizeof(double));
            Lind[Ln]=ind[i+1];
            Ln++;
        }else{
            memcpy(&RX[Rn*d], &X[(i+1)*d], d*sizeof(double));
            Rind[Rn]=ind[i+1];
            Rn++;
        }
    }
    free(dist);

    LX=(double*)realloc(LX,   Ln*d*sizeof(double));
    Lind= (int*)realloc(Lind, Ln*  sizeof(int));
    RX=(double*)realloc(RX,   Rn*d*sizeof(double));
    Rind= (int*)realloc(Rind, Rn*  sizeof(int));

    node->left=  make_vp_tree_select(LX,Lind,d,Ln);
    free(LX);
    free(Lind);
    node->right= make_vp_tree_select(RX,Rind,d,Rn);
    free(RX);
    free(Rind);

    return node;
}

// ALL VP TREE MAKERS ARE DESTRICTIVE (for speed and simplicity), parse a COPY of X if you want to retain the element's order
// makes a vp tree with buckets of size B, returns pointer to the head
VPnode* make_vp_bucket_tree(double *X, int *ind, const int d, const int n, const int B){

    if(n<=B){
        VPnode *bucket=(VPnode*)malloc(sizeof(VPnode));
        bucket->bucketX=(double*)malloc(n*d*sizeof(double));
        bucket->bucketIndX=(int*)malloc(n*  sizeof(int));

        bucket->right=NULL;
        bucket->left=NULL;

        bucket->isBucket=1;

        bucket->bucketN=n;
        memcpy(bucket->bucketX,      X, n*d*sizeof(double));
        memcpy(bucket->bucketIndX, ind, n*  sizeof(int));

        return bucket;

    }else if(n==0){
        return NULL;
    }

    double *LX=(double*)malloc(n*d*sizeof(double));
    double *RX=(double*)malloc(n*d*sizeof(double));
    int  *Lind=(int*)   malloc(n*sizeof(int));
    int  *Rind=(int*)   malloc(n*sizeof(int));
    
    double *dist=(double*)malloc(n*sizeof(double));
    double *distTemp=(double*)malloc(n*sizeof(double));
    dist[0]=0;
    double temp,tempSq;

    VPnode *node = (VPnode*)malloc(sizeof(VPnode));
    
    node->isBucket=0;
    node->indP=ind[0];
    node->P=(double*)malloc(d*sizeof(double));
    memcpy(node->P,X,d*sizeof(double));


    for(int i=1; i<n; i++){
        tempSq=0;
        for(int j=0; j<d; j++){
            temp=(X[0*d+j] - X[i*d + j]);
            tempSq+=temp*temp;
        }
        dist[i]=sqrt(tempSq);
    }

    memcpy(distTemp,dist,n*sizeof(double));

    kthSmallestNoInd(distTemp,1,n-1,(n-1)/2);
    node->mu=distTemp[(n-1)/2+1];
    free(distTemp);

    int Ln=0, Rn=0;

    for(int i=0; i<n-1; i++){
        if(dist[i+1] < node->mu){
            memcpy(&LX[Ln*d], &X[(i+1)*d], d*sizeof(double));
            Lind[Ln]=ind[i+1];
            Ln++;
        }else{
            memcpy(&RX[Rn*d], &X[(i+1)*d], d*sizeof(double));
            Rind[Rn]=ind[i+1];
            Rn++;
        }
    }
    free(dist);

    LX=(double*)realloc(LX,   Ln*d*sizeof(double));
    Lind= (int*)realloc(Lind, Ln*  sizeof(int));
    RX=(double*)realloc(RX,   Rn*d*sizeof(double));
    Rind= (int*)realloc(Rind, Rn*  sizeof(int));

    node->left=  make_vp_bucket_tree(LX,Lind,d,Ln,B);
    free(LX);
    free(Lind);
    node->right= make_vp_bucket_tree(RX,Rind,d,Rn,B);
    free(RX);
    free(Rind);

    return node;
}

// ALL VP TREE MAKERS ARE DESTRICTIVE (for speed and simplicity), parse a COPY of X if you want to retain the element's order
// makes a vp tree with buckets of size B with a statistically good selection for VP, returns pointer to the head
VPnode* make_vp_bucket_tree_select(double *X, int *ind, const int d, const int n, const int B){

    if(n<=B){
        VPnode *bucket=(VPnode*)malloc(sizeof(VPnode));
        bucket->bucketX=(double*)malloc(n*d*sizeof(double));
        bucket->bucketIndX=(int*)malloc(n*  sizeof(int));

        bucket->right=NULL;
        bucket->left=NULL;

        bucket->isBucket=1;

        bucket->bucketN=n;
        memcpy(bucket->bucketX,      X, n*d*sizeof(double));
        memcpy(bucket->bucketIndX, ind, n*  sizeof(int));

        return bucket;

    }else if(n==0){
        return NULL;
    }

    double *LX=(double*)malloc(n*d*sizeof(double));
    double *RX=(double*)malloc(n*d*sizeof(double));
    int  *Lind=(int*)   malloc(n*sizeof(int));
    int  *Rind=(int*)   malloc(n*sizeof(int));
    
    double *dist=(double*)malloc(n*sizeof(double));
    double *distTemp=(double*)malloc(n*sizeof(double));
    dist[0]=0;
    double temp,tempSq;

    VPnode *node = (VPnode*)malloc(sizeof(VPnode));
    

    //----------------------------------------------------------------
    // using select vp to find a suitable vp and swap its data
    // with the first element for ease of use later
    int vp=select_vp(X,d,n);
    int tempInd;
    double *tempX=(double*)malloc(d*sizeof(double));

    // move X[vp]
    memcpy(tempX,  X+vp*d, d*sizeof(double));
    memcpy(X+vp*d, X+0*d,  d*sizeof(double));
    memcpy(X+0*d,  tempX,  d*sizeof(double));
    free(tempX);
    // move its index
    tempInd=ind[vp];
    ind[vp]=ind[0];
    ind[0]=tempInd;

    //----------------------------------------------------------------

    node->isBucket=0;
    node->indP=ind[0];
    node->P=(double*)malloc(d*sizeof(double));
    memcpy(node->P,X,d*sizeof(double));


    for(int i=1; i<n; i++){
        tempSq=0;
        for(int j=0; j<d; j++){
            temp=(X[0*d+j] - X[i*d + j]);
            tempSq+=temp*temp;
        }
        dist[i]=sqrt(tempSq);
    }

    memcpy(distTemp,dist,n*sizeof(double));

    kthSmallestNoInd(distTemp,1,n-1,(n-1)/2);
    node->mu=distTemp[(n-1)/2+1];
    free(distTemp);

    int Ln=0, Rn=0;

    for(int i=0; i<n-1; i++){
        if(dist[i+1] < node->mu){
            memcpy(&LX[Ln*d], &X[(i+1)*d], d*sizeof(double));
            Lind[Ln]=ind[i+1];
            Ln++;
        }else{
            memcpy(&RX[Rn*d], &X[(i+1)*d], d*sizeof(double));
            Rind[Rn]=ind[i+1];
            Rn++;
        }
    }
    free(dist);

    LX=(double*)realloc(LX,   Ln*d*sizeof(double));
    Lind= (int*)realloc(Lind, Ln*  sizeof(int));
    RX=(double*)realloc(RX,   Rn*d*sizeof(double));
    Rind= (int*)realloc(Rind, Rn*  sizeof(int));

    node->left=  make_vp_bucket_tree_select(LX,Lind,d,Ln,B);
    free(LX);
    free(Lind);
    node->right= make_vp_bucket_tree_select(RX,Rind,d,Rn,B);
    free(RX);
    free(Rind);

    return node;
}

// frees a vp tree
void freeVP(VPnode *head){
    if(head==NULL) return;

    if(head->left!=NULL){
        freeVP(head->left);
    }
    if(head->right!=NULL){
        freeVP(head->right);
    }

    if(head->isBucket==1){
        free(head->bucketX);
        free(head->bucketIndX);
    }else{
        free(head->P);
    }
    free(head);
}

// inserts a distance into the k smallest of a querry and updates tau
void pknn_insert(VPnode *node, Pknnresult *query, const double dist, const int d, const int k){
    query->dist[k] = dist;
    query->idx[k]  = node->indP;
    binInsert(query->dist, query->idx, k+1);

    query->tau=query->dist[k-1];
    //if(query->tauInd < k-1) query->tauInd++;
}

// searches a normal vp tree and updates querry's kNN
Pknnresult* search_vp_tree(VPnode* node, Pknnresult *query, const int d, const int k){


    if(node== NULL){
        return query;
    }else{
        double dist = distCalc(node->P,query->cord,d);
        if(dist < query->dist[k-1]){
            pknn_insert(node,query,dist,d,k);
        }

        if(dist < node->mu){
            query=search_vp_tree(node->left, query, d, k);
            if((dist + query->tau) >= node->mu){
                query=search_vp_tree(node->right, query, d, k);
            }
        }else{
            query=search_vp_tree(node->right, query, d, k);
            if((dist - query->tau) <= node->mu){
                query=search_vp_tree(node->left, query, d, k);
            }
        }

        return query;

    }

}

// searches a bucket vp tree and updates querry's kNN
Pknnresult* search_vp_bucket_tree(VPnode* node, Pknnresult *query, const int d, const int k){


    if(node== NULL){
        return query;
    }else if(node->isBucket==0){
        double dist = distCalc(node->P,query->cord,d);
        if(dist < query->dist[k-1]){
            pknn_insert(node,query,dist,d,k);
        }

        if(dist < node->mu){
            query=search_vp_bucket_tree(node->left, query, d, k);
            if((dist + query->tau) >= node->mu){
                query=search_vp_bucket_tree(node->right, query, d, k);
            }
        }else{
            query=search_vp_bucket_tree(node->right, query, d, k);
            if((dist - query->tau) <= node->mu){
                query=search_vp_bucket_tree(node->left, query, d, k);
            }
        }

        return query;

    }else{
        double dist;

        for(int i=0; i< node->bucketN; i++){

            dist=distCalc(query->cord, node->bucketX+i*d, d);

            if(dist < query->dist[k-1]){
                query->dist[k]=dist;
                query->idx[k]=node->bucketIndX[i];
                binInsert(query->dist, query->idx, k+1);

                query->tau=query->dist[k-1];
                //if(query->tauInd < k-1) query->tauInd++;

            }
        }

        return query;

    }

}

// initializes a pkNN result needed to form a query point;
Pknnresult* create_pknnresult(const int d, const int k){
    Pknnresult *pknn=(Pknnresult*)malloc(sizeof(Pknnresult));
    pknn->cord=(double*)malloc(d*    sizeof(double));
    pknn->dist=(double*)malloc((k+1)*sizeof(double));
    pknn->idx= (int*)   malloc((k+1)*sizeof(int));
    return pknn;
}

Pknnresult* initialize_pknnresult(Pknnresult* pknn, const double *cord, const int d, const int k){
    //pknn->tauInd=0;
    pknn->tau=999;

    for(int i=0; i<k+1; i++){
        pknn->dist[i]=999;
        pknn->idx[i]=999;
    }
    memcpy(pknn->cord,cord,d*sizeof(double));
    return pknn;
}

// frees a pkNN result
void free_pknnresult(Pknnresult *node){
    free(node->dist);
    free(node->idx);
    free(node->cord);
    free(node);
}

// prints VP tree structure, space should be 0, count >5
void printVPtree(VPnode *root, int space, const int COUNT) 
{ 
    // Base case 
    if (root == NULL) 
        return; 
  
    // Increase distance between levels 
    space += COUNT; 
  
    // Process right child first 
    printVPtree(root->right, space, COUNT); 
  
    // Print current node after space 
    // count 
    printf("\n"); 

    for(int i = COUNT; i < space; i++) 
        printf(" "); 

    if(root->isBucket==0)
        printf("%d\n", root->indP); 
    else if(root->isBucket==1){
        for(int i=0; i< root->bucketN; i++){
            printf("[%d] ", root->bucketIndX[i]);
        }
        printf("\n");
    }
    // Process left child 
    printVPtree(root->left, space, COUNT); 
} 


VPtree *initialize_vpt(VPtree *tree, const int d, const int n){
    tree=(VPtree*)malloc(sizeof(VPtree));
    tree->indP    =(int*) malloc(n*sizeof(double));
    tree->mu    =(double*)malloc(n*sizeof(double));
    tree->cordsP=(double*)malloc(n*d*sizeof(double));
    tree->leftptr= (int*) malloc(n*sizeof(int));
    tree->rightptr=(int*) malloc(n*sizeof(int));

    for(int i=0; i<n; i++){
        tree->mu[i]=0;
        tree->leftptr[i]=0;
    }
    memcpy(tree->rightptr, tree->leftptr, n*sizeof(int));

    return tree;
}

void make_serial_vpt(VPtree *tree, double *X, int *ind, const int VPptr, const int d, const int n){

    if(n==1){
        tree->leftptr[VPptr]=-1;
        tree->rightptr[VPptr]=-1;
        //tree->isBucket[VPptr]=0;
        tree->mu[VPptr]=0;
        tree->indP[VPptr]=ind[0];
        memcpy(&tree->cordsP[VPptr*d], &X[0], d*sizeof(double));
        return;
    }else if(n==0){
        return;
    }
    double *LX=(double*)malloc(n*d*sizeof(double));
    double *RX=(double*)malloc(n*d*sizeof(double));
    int  *Lind=(int*)   malloc(n*sizeof(int));
    int  *Rind=(int*)   malloc(n*sizeof(int));

    //tree->isBucket[VPptr] = 0;
    tree->indP    [VPptr] = ind[0];
    memcpy(&tree->cordsP[VPptr*d], &X[0], d*sizeof(double));

    double *dist=    (double*)malloc(n*sizeof(double));
    double *distTemp=(double*)malloc(n*sizeof(double));
    dist[0]=0;
    double tempSq,temp;

    for(int i=1; i<n; i++){
        tempSq=0;
        for(int j=0; j<d; j++){
            temp=(X[0*d+j] - X[i*d + j]);
            tempSq+=temp*temp;
        }
        dist[i]=sqrt(tempSq);
    }


    memcpy(distTemp,dist,n*sizeof(double));

    kthSmallestNoInd(distTemp, 1, n-1, (n-1)/2);
    tree->mu[VPptr] = distTemp[(n-1)/2+1];
    free(distTemp);


    int Ln=0, Rn=0;

    for(int i=0; i<n-1; i++){
        if(dist[i+1] < tree->mu[VPptr]){
            memcpy(&LX[Ln*d], &X[(i+1)*d], d*sizeof(double));
            Lind[Ln]=ind[i+1];
            Ln++;
        }else{
            memcpy(&RX[Rn*d], &X[(i+1)*d], d*sizeof(double));
            Rind[Rn]=ind[i+1];
            Rn++;
        }
    }
    free(dist);

    LX=(double*)realloc(LX,   Ln*d*sizeof(double));
    Lind= (int*)realloc(Lind, Ln*  sizeof(int));
    RX=(double*)realloc(RX,   Rn*d*sizeof(double));
    Rind= (int*)realloc(Rind, Rn*  sizeof(int));

    if(Ln==0){
        tree->leftptr [VPptr]=-1;
    }else{
        tree->leftptr [VPptr]=VPptr+1;
    }
    
    if(Rn==0){
        tree->rightptr [VPptr]=-1;
    }else{
        tree->rightptr [VPptr]=VPptr+1+Ln;
    }
    
    make_serial_vpt(tree, LX, Lind, VPptr+1,    d, Ln);
    free(LX);
    free(Lind);
    make_serial_vpt(tree, RX, Rind, VPptr+1+Ln, d, Rn);
    free(RX);
    free(Rind);

    return;
}

void make_serial_vpt_select(VPtree *tree, double *X, int *ind, const int VPptr, const int d, const int n){

    if(n==1){
        tree->leftptr[VPptr]=-1;
        tree->rightptr[VPptr]=-1;
        //tree->isBucket[VPptr]=0;
        tree->mu[VPptr]=0;
        tree->indP[VPptr]=ind[0];
        memcpy(&tree->cordsP[VPptr*d], &X[0], d*sizeof(double));
        return;
    }else if(n==0){
        return;
    }
    double *LX=(double*)malloc(n*d*sizeof(double));
    double *RX=(double*)malloc(n*d*sizeof(double));
    int  *Lind=(int*)   malloc(n*sizeof(int));
    int  *Rind=(int*)   malloc(n*sizeof(int));

    //----------------------------------------------------------------
    // using select vp to find a suitable vp and swap its data
    // with the first element for ease of use later
    int vp=select_vp(X,d,n);
    int tempInd;
    double *tempX=(double*)malloc(d*sizeof(double));

    // move X[vp]
    memcpy(tempX,  X+vp*d, d*sizeof(double));
    memcpy(X+vp*d, X+0*d,  d*sizeof(double));
    memcpy(X+0*d,  tempX,  d*sizeof(double));
    free(tempX);
    // move its index
    tempInd=ind[vp];
    ind[vp]=ind[0];
    ind[0]=tempInd;

    //----------------------------------------------------------------

    //tree->isBucket[VPptr] = 0;
    tree->indP    [VPptr] = ind[0];
    memcpy(&tree->cordsP[VPptr*d], &X[0], d*sizeof(double));

    double *dist=    (double*)malloc(n*sizeof(double));
    double *distTemp=(double*)malloc(n*sizeof(double));
    dist[0]=0;
    double tempSq,temp;

    for(int i=1; i<n; i++){
        tempSq=0;
        for(int j=0; j<d; j++){
            temp=(X[0*d+j] - X[i*d + j]);
            tempSq+=temp*temp;
        }
        dist[i]=sqrt(tempSq);
    }


    memcpy(distTemp,dist,n*sizeof(double));

    kthSmallestNoInd(distTemp, 1, n-1, (n-1)/2);
    tree->mu[VPptr] = distTemp[(n-1)/2+1];
    free(distTemp);


    int Ln=0, Rn=0;

    for(int i=0; i<n-1; i++){
        if(dist[i+1] < tree->mu[VPptr]){
            memcpy(&LX[Ln*d], &X[(i+1)*d], d*sizeof(double));
            Lind[Ln]=ind[i+1];
            Ln++;
        }else{
            memcpy(&RX[Rn*d], &X[(i+1)*d], d*sizeof(double));
            Rind[Rn]=ind[i+1];
            Rn++;
        }
    }
    free(dist);

    LX=(double*)realloc(LX,   Ln*d*sizeof(double));
    Lind= (int*)realloc(Lind, Ln*  sizeof(int));
    RX=(double*)realloc(RX,   Rn*d*sizeof(double));
    Rind= (int*)realloc(Rind, Rn*  sizeof(int));

    if(Ln==0){
        tree->leftptr [VPptr]=-1;
    }else{
        tree->leftptr [VPptr]=VPptr+1;
    }
    
    if(Rn==0){
        tree->rightptr [VPptr]=-1;
    }else{
        tree->rightptr [VPptr]=VPptr+1+Ln;
    }
    
    make_serial_vpt_select(tree, LX, Lind, VPptr+1,    d, Ln);
    free(LX);
    free(Lind);
    make_serial_vpt_select(tree, RX, Rind, VPptr+1+Ln, d, Rn);
    free(RX);
    free(Rind);

    return;
}

void make_serial_bucket_vpt(VPtree *tree, double *X, int *ind, const int VPptr, const int d, const int n, const int B){

    if(n==0){
        return;
    }else if(n<=B){
        memcpy(&tree->cordsP[VPptr*d],   X, n*d*sizeof(double));
        memcpy(&tree->indP  [VPptr], ind, n*  sizeof(int));
        
        tree->leftptr[VPptr]=-2;
        tree->rightptr[VPptr]=n;//we only need n here so we save it in rightptr to save space since rightptr isnt needed here
        tree->mu[VPptr]=0; 
        memcpy(&tree->cordsP[VPptr*d], &X[0], n*d*sizeof(double));

        return;
    }

    double *LX=(double*)malloc(n*d*sizeof(double));
    double *RX=(double*)malloc(n*d*sizeof(double));
    int  *Lind=(int*)   malloc(n*sizeof(int));
    int  *Rind=(int*)   malloc(n*sizeof(int));

    //tree->isBucket[VPptr] = 0;
    tree->indP    [VPptr] = ind[0];
    memcpy(&tree->cordsP[VPptr*d], &X[0], d*sizeof(double));

    double *dist=    (double*)malloc(n*sizeof(double));
    double *distTemp=(double*)malloc(n*sizeof(double));
    dist[0]=0;
    double tempSq,temp;

    for(int i=1; i<n; i++){
        tempSq=0;
        for(int j=0; j<d; j++){
            temp=(X[0*d+j] - X[i*d + j]);
            tempSq+=temp*temp;
        }
        dist[i]=sqrt(tempSq);
    }


    memcpy(distTemp,dist,n*sizeof(double));

    kthSmallestNoInd(distTemp, 1, n-1, (n-1)/2);
    tree->mu[VPptr] = distTemp[(n-1)/2+1];
    free(distTemp);


    int Ln=0, Rn=0;

    for(int i=0; i<n-1; i++){
        if(dist[i+1] < tree->mu[VPptr]){
            memcpy(&LX[Ln*d], &X[(i+1)*d], d*sizeof(double));
            Lind[Ln]=ind[i+1];
            Ln++;
        }else{
            memcpy(&RX[Rn*d], &X[(i+1)*d], d*sizeof(double));
            Rind[Rn]=ind[i+1];
            Rn++;
        }
    }
    free(dist);

    LX=(double*)realloc(LX,   Ln*d*sizeof(double));
    Lind= (int*)realloc(Lind, Ln*  sizeof(int));
    RX=(double*)realloc(RX,   Rn*d*sizeof(double));
    Rind= (int*)realloc(Rind, Rn*  sizeof(int));

    if(Ln==0){
        tree->leftptr [VPptr]=-1;
    }else{
        tree->leftptr [VPptr]=VPptr+1;
    }
    
    if(Rn==0){
        tree->rightptr [VPptr]=-1;
    }else{
        tree->rightptr [VPptr]=VPptr+1+Ln;
    }
    
    make_serial_bucket_vpt(tree, LX, Lind, VPptr+1,    d, Ln, B);
    free(LX);
    free(Lind);
    make_serial_bucket_vpt(tree, RX, Rind, VPptr+1+Ln, d, Rn, B);
    free(RX);
    free(Rind);

    return;
}

void make_serial_bucket_vpt_select(VPtree *tree, double *X, int *ind, const int VPptr, const int d, const int n, const int B){

    if(n==0){
        return;
    }else if(n<=B){
        memcpy(&tree->cordsP[VPptr*d],   X, n*d*sizeof(double));
        memcpy(&tree->indP  [VPptr], ind, n*  sizeof(int));
        
        tree->leftptr[VPptr]=-2;
        tree->rightptr[VPptr]=n;//we only need n here so we save it in rightptr to save space since rightptr isnt needed here
        tree->mu[VPptr]=0; 
        memcpy(&tree->cordsP[VPptr*d], &X[0], n*d*sizeof(double));

        return;
    }

    double *LX = malloc(n*d*sizeof(double));
    double *RX = malloc(n*d*sizeof(double));
    int  *Lind = malloc(n*  sizeof(int));
    int  *Rind = malloc(n*  sizeof(int));

    //----------------------------------------------------------------
    // using select vp to find a suitable vp and swap its data
    // with the first element for ease of use later
    int vp = select_vp(X,d,n);
    int tempInd;
    double *tempX = malloc(d*sizeof(double));

    // move X[vp]
    memcpy(tempX,  X+vp*d, d*sizeof(double));
    memcpy(X+vp*d, X+0*d,  d*sizeof(double));
    memcpy(X+0*d,  tempX,  d*sizeof(double));
    free(tempX);
    // move its index
    tempInd=ind[vp];
    ind[vp]=ind[0];
    ind[0]=tempInd;

    //----------------------------------------------------------------


    tree->indP    [VPptr] = ind[0];
    memcpy(&tree->cordsP[VPptr*d], &X[0], d*sizeof(double));

    double *dist=    (double*)malloc(n*sizeof(double));
    double *distTemp=(double*)malloc(n*sizeof(double));
    dist[0]=0;
    double tempSq,temp;

    for(int i=1; i<n; i++){
        tempSq=0;
        for(int j=0; j<d; j++){
            temp=(X[0*d+j] - X[i*d + j]);
            tempSq+=temp*temp;
        }
        dist[i]=sqrt(tempSq);
    }


    memcpy(distTemp,dist,n*sizeof(double));

    kthSmallestNoInd(distTemp, 1, n-1, (n-1)/2);
    tree->mu[VPptr] = distTemp[(n-1)/2+1];
    free(distTemp);


    int Ln=0, Rn=0;

    for(int i=0; i<n-1; i++){
        if(dist[i+1] < tree->mu[VPptr]){
            memcpy(&LX[Ln*d], &X[(i+1)*d], d*sizeof(double));
            Lind[Ln]=ind[i+1];
            Ln++;
        }else{
            memcpy(&RX[Rn*d], &X[(i+1)*d], d*sizeof(double));
            Rind[Rn]=ind[i+1];
            Rn++;
        }
    }
    free(dist);

    LX=(double*)realloc(LX,   Ln*d*sizeof(double));
    Lind= (int*)realloc(Lind, Ln*  sizeof(int));
    RX=(double*)realloc(RX,   Rn*d*sizeof(double));
    Rind= (int*)realloc(Rind, Rn*  sizeof(int));

    if(Ln==0){
        tree->leftptr [VPptr]=-1;
    }else{
        tree->leftptr [VPptr]=VPptr+1;
    }
    
    if(Rn==0){
        tree->rightptr [VPptr]=-1;
    }else{
        tree->rightptr [VPptr]=VPptr+1+Ln;
    }
    
    make_serial_bucket_vpt_select(tree, LX, Lind, VPptr+1,    d, Ln, B);
    free(LX);
    free(Lind);
    make_serial_bucket_vpt_select(tree, RX, Rind, VPptr+1+Ln, d, Rn, B);
    free(RX);
    free(Rind);

    return;
}

void free_searial_vpt(VPtree *tree){
    free(tree->cordsP);
    free(tree->indP);
    free(tree->leftptr);
    free(tree->rightptr);
    free(tree->mu);
    free(tree);
}

void pknn_insert_serial(const int indP, Pknnresult *query, const double dist, const int d, const int k){
    query->dist[k] = dist;
    query-> idx[k] = indP;
    binInsert(query->dist, query->idx, k+1);

    query->tau = query->dist[k-1];
    //if(query->tauInd < k-1) query->tauInd++;
}

Pknnresult* search_serial_vpt(VPtree *tree, const int VPptr, Pknnresult *query, const int d, const int k){


    if(VPptr==-1){
        return query;
    }else{
        double dist = distCalc(&tree->cordsP[VPptr*d], query->cord, d);
        if(dist < query->dist[k-1]){
            pknn_insert_serial(tree->indP[VPptr], query, dist, d, k);
        }

        if(dist < tree->mu[VPptr]){
            query=search_serial_vpt(tree, tree->leftptr[VPptr], query, d, k);
            if((dist + query->tau) > tree->mu[VPptr]){
                query=search_serial_vpt(tree, tree->rightptr[VPptr], query, d, k);
            }
        }else{
            query=search_serial_vpt(tree, tree->rightptr[VPptr], query, d, k);
            if((dist - query->tau) < tree->mu[VPptr]){
                query=search_serial_vpt(tree, tree->leftptr[VPptr], query, d, k);
            }
        }

        return query;

    }

}

Pknnresult* search_serial_bucket_vpt(VPtree *tree, const int VPptr, Pknnresult *query, const int d, const int k){

// leftptr==-2 indicates bucket, rightptr stores n in this occasion to save space
    if(VPptr==-1){
        return query;
    }else if(tree->leftptr[VPptr]==-2){
        double dist;

        for(int i=0; i< tree->rightptr[VPptr]; i++){

            dist = distCalc(query->cord, &tree->cordsP[(VPptr + i) * d], d);

            if(dist < query->dist[k-1]){
                query->dist[k] = dist;
                query-> idx[k] = tree->indP[VPptr + i];
                binInsert(query->dist, query->idx, k+1);

                query->tau = query->dist[k-1];
            }
        }

        return query;
    }else{
        double dist = distCalc(&tree->cordsP[VPptr*d], query->cord, d);
        if(dist < query->dist[k-1]){
            pknn_insert_serial(tree->indP[VPptr], query, dist, d, k);
        }

        if(dist < tree->mu[VPptr]){
            query = search_serial_bucket_vpt(tree, tree->leftptr[VPptr], query, d, k);
            if((dist + query->tau) > tree->mu[VPptr]){
                query = search_serial_bucket_vpt(tree, tree->rightptr[VPptr], query, d, k);
            }
        }else{
            query = search_serial_bucket_vpt(tree, tree->rightptr[VPptr], query, d, k);
            if((dist - query->tau) < tree->mu[VPptr]){
                query = search_serial_bucket_vpt(tree, tree->leftptr[VPptr], query, d, k);
            }
        }

        return query;

    }

}

void printVPtree_serial(VPtree *tree, const int VPptr, int space, const int COUNT) 
{ 
    // Base case 
    if (VPptr == -1) 
        return; 
  
    // Increase distance between levels 
    
    space += COUNT; 
    // Process right child first 
    if(tree->leftptr[VPptr]==-2){
        for(int i = COUNT; i < space; i++) 
            printf(" ");
        for(int i=0; i< tree->rightptr[VPptr]; i++){
            printf("[%d] ", tree->indP[VPptr+i]);
        }
        printf("\n");
        return;
    }

    
    printVPtree_serial(tree,tree->rightptr[VPptr],space,COUNT);
  
    for(int i = COUNT; i < space; i++) 
        printf(" ");
    // Print current node after space 
    // count 
    printf("\n"); 

    for(int i = COUNT; i < space; i++) 
        printf(" "); 

    printf("%d\n", tree->indP[VPptr]); 


    printVPtree_serial(tree,tree->leftptr[VPptr], space, COUNT); 
} 

/* int main(int argc, char const *argv[])
{
    //srand(time(NULL));

    const int n=9,d=2;
    double X[]={1,1,2,2,3,3,4,4,5,5,7,6,7,7,8,8,9,9};
    int ind[]={1,2,3,4,5,6,7,8,9};

    VPnode *head=make_vp_tree(X,ind,d,n);
    printf("\nhi\n");

    printVPtree(head,0,5);
    printf("\nhi\n");
    freeVP(head);
    //-------------------------------------------------------

    double y[]={7,6};
    int k=5;
    Pknnresult *test=create_pknnresult(d,k);

    test=initialize_pknnresult(test,y,d,k);



    test=search_vp_tree(head,test,d,k);
    
    printf("\n");
    for(int i=0; i<k+1; i++){
        printf("[%d: %.3lf]",test->idx[i],test->dist[i]);
    }
    printf("\n");  
    free_pknnresult(test);  
    freeVP(head);

    //-------------------------------------------------------

    double y[]={9,9};
    const int k=6;

    VPnode *headB=make_vp_bucket_tree_select(X,ind,d,n,3);

    printVPtree(headB,0,5);

    Pknnresult *testB=create_pknnresult(d,k);
    testB=initialize_pknnresult(testB,y,d,k);

    testB=search_vp_bucket_tree(headB,testB,d,k);
    printf("\n");  
    for(int i=0; i<k+1; i++){
        printf("[%d: %.3lf]",testB->idx[i],testB->dist[i]);
    }
    printf("\n"); 
    free_pknnresult(testB);
    freeVP(headB);


    return 0;
} */
