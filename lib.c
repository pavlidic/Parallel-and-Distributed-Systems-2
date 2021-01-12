#include <limits.h>
#include "lib.h"

void swapD(double *a,double *b){
    const double temp=*a;
    *a=*b;
    *b=temp;
}

void swapI(int *a,int *b){
    const int temp=*a;
    *a=*b;
    *b=temp;
}

int partition(double arr[],int indices[], const int l, const int r) 
{ 
    const double x = arr[r];
    int i = l; 

    for (int j = l; j < r; j++) { 
        if (arr[j] <= x) { 
            swapD(&arr[i], &arr[j]);
            swapI(&indices[i], &indices[j]);  
            i++; 
        } 
    } 
    swapD(&arr[i], &arr[r]);
    swapI(&indices[i], &indices[r]);
    return i; 
} 

int partitionNoInd(double arr[], const int l, const int r) 
{ 
    const double x = arr[r];
    int i = l; 
    for (int j = l; j < r; j++) { //was j<=r-1
        if (arr[j] <= x) { 
            swapD(&arr[i], &arr[j]);
            i++; 
        } 
    } 
    swapD(&arr[i], &arr[r]);
    return i; 
} 

int kthSmallest(double arr[],int indices[], const int l, const int r, const int k) 
{ 
    // If k is smaller than number of  
    // elements in array 
    if (k > 0 && k <= r - l + 1) { 
  
        // Partition the array around last  
        // element and get position of pivot  
        // element in sorted array 
        const int index = partition(arr,indices, l, r); 
  
        // If position is same as k 
        if (index - l == k - 1) 
            return arr[index]; 
  
        // If position is more, recur  
        // for left subarray 
        if (index - l > k - 1)  
            return kthSmallest(arr,indices, l, index - 1, k); 
  
        // Else recur for right subarray 
        return kthSmallest(arr,indices, index + 1, r,  
                            k - index + l - 1); 
    } 
  
    // If k is more than number of  
    // elements in array 
    return INT_MAX; 
} 

int kthSmallestNoInd(double arr[], const int l, const int r, const int k) 
{ 
    // If k is smaller than number of  
    // elements in array 
    if (k > 0 && k <= r - l + 1) { 
  
        // Partition the array around last  
        // element and get position of pivot  
        // element in sorted array 
        const int index = partitionNoInd(arr, l, r); 
  
        // If position is same as k 
        if (index - l == k - 1) 
            return arr[index]; 
  
        // If position is more, recur  
        // for left subarray 
        if (index - l > k - 1)  
            return kthSmallestNoInd(arr, l, index - 1, k); 
  
        // Else recur for right subarray 
        return kthSmallestNoInd(arr, index + 1, r,  
                            k - index + l - 1); 
    } 
  
    // If k is more than number of  
    // elements in array 
    return INT_MAX; 
} 


