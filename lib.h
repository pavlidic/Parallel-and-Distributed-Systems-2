#ifndef LIB__
#define LIB__

void swapD(double *a,double *b);

void swapI(int *a,int *b);

int partition(double arr[],int indices[], const int l, const int r);

int partitionNoInd(double arr[], const int l, const int r);

int kthSmallest(double arr[],int indices[], const int l, const int r, const int k);

int kthSmallestNoInd(double arr[], const int l, const int r, const int k);

#endif