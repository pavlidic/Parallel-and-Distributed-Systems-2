
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// A binary search based function
// to find the position
// where item should be inserted 
// in a[low..high]
int binarySearch(const double a[], const double item, 
				const int low, const int high)
{
	if (high <= low)
		return (item > a[low]) ? 
				(low + 1) : low;

	const int mid = (low + high) / 2;

	if (item == a[mid])
		return mid + 1;

	if (item > a[mid])
		return binarySearch(a, item, 
							mid + 1, high);
	return binarySearch(a, item, low, 
						mid - 1);
}


void binInsert(double a[], int ind[], const int n)
{
	int j;
    double selected_element;
    int selected_index;
	
    j=n-1;
    selected_element = a[j];
    selected_index = ind[j];

    // find location where selected sould be inseretd
    const int loc = binarySearch(a, selected_element, 0, j-1);

    // Move all elements after location to create space
    while (j > loc) 
    {
        a[j] = a[j-1];
        ind[j] = ind[j-1];
        j--;
    }
    a[j] = selected_element;
    ind[j]= selected_index;

} 
void binInsert1(double a[], int ind[], const int n)
{
	int j;
    double selected_element;
    int selected_index;
	
    j=n-1;
    selected_element = a[j];
    selected_index = ind[j];

    // find location where selected sould be inseretd
    const int loc = binarySearch(a, selected_element, 0, j-1);

    // Move all elements after location to create space
    memcpy(a+loc+1,a+loc,(n-loc-1)*sizeof(double));
    memcpy(ind+loc+1,ind+loc,(n-loc-1)*sizeof(int));
    
    a[loc]=selected_element;
    ind[loc]=selected_index;

} 


// Driver Code
int insertTesting(int argc, char const *argv[])
{

    double A[]={1.2,2.2,3.6,3.7,3.9,4,5};
    int ind[]={1,3,6,4,2,7,9};

	int n = sizeof(A) / sizeof(A[0]), i;
    
    A[n-1]= atof(argv[1]);
    
    for(int i=0; i<n; i++){
        printf("[%d: %.0lf]",ind[i],A[i]);
    }
    printf("\n");
	binInsert(A,ind,n);

    for(int i=0; i<n; i++){
        printf("[%d: %.0lf]",ind[i],A[i]);
    }
    printf("\n");


	return 0;
}
