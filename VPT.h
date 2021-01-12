#ifndef VPT__
#define VPT__
#include <stdbool.h>

typedef struct VPnode{

    int indP;
    double *P;
    double mu;

    bool isBucket;
    int bucketN;
    double *bucketX;
    int *bucketIndX;

    struct VPnode *left,*right;

} VPnode;

typedef struct Pknnresult{
  double   tau;
  double * cord;
  int    * idx;
  double * dist;    
} Pknnresult;

typedef struct VPtree{
    int    *indP;
    double *cordsP;
    double *mu;
    int    *leftptr, *rightptr;
} VPtree;

double distCalc(const double *x,const double *y, const int d);

int select_vp(const double *X, const int d, const int n);

// node stuff here

VPnode* make_vp_tree(double *X, int *ind, const int d, const int n);

VPnode* make_vp_tree_select(double *X, int *ind, const int d, const int n);

VPnode* make_vp_bucket_tree(double *X, int *ind, const int d, const int n, const int B);

VPnode* make_vp_bucket_tree_select(double *X, int *ind, const int d, const int n, const int B);

void freeVP(VPnode *head);

// query stuff here

void pknn_insert(VPnode *node, Pknnresult *query, const double dist, const int d, const int k);

Pknnresult* search_vp_tree(VPnode* node, Pknnresult *query, const int d, const int k);

Pknnresult* search_vp_bucket_tree(VPnode* node, Pknnresult *query, const int d, const int k);

Pknnresult* create_pknnresult(const int d, const int k);

Pknnresult* initialize_pknnresult(Pknnresult* pknn, const double *cord, const int d, const int k);

void free_pknnresult(Pknnresult *node);

void printVPtree(VPnode *root, int space, const int COUNT);

// serial stuff here

VPtree *initialize_vpt(VPtree *tree, const int d, const int n);

void make_serial_vpt(VPtree *tree, double *X, int *ind, const int VPptr, const int d, const int n);

void make_serial_vpt_select(VPtree *tree, double *X, int *ind, const int VPptr, const int d, const int n);

void make_serial_bucket_vpt(VPtree *tree, double *X, int *ind, const int VPptr, const int d, const int n, const int B);

void make_serial_bucket_vpt_select(VPtree *tree, double *X, int *ind, const int VPptr, const int d, const int n, const int B);

void free_searial_vpt(VPtree *tree);

void pknn_insert_serial(int indP, Pknnresult *query, const double dist, const int d, const int k);

Pknnresult* search_serial_vpt(VPtree *tree, const int VPptr, Pknnresult *query, const int d, const int k);

Pknnresult* search_serial_bucket_vpt(VPtree *tree, const int VPptr, Pknnresult *query, const int d, const int k);

void printVPtree_serial(VPtree *tree, const int VPptr, int space, const int COUNT);

#endif