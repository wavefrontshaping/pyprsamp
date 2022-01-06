#include <string.h>
#include <math.h>
#include <complex.h>
#include "gsl_sf_bessel.h"
//#include <time.h>
#include <stdlib.h>
#include <stdio.h>
//#include <stdbool.h>
#include <Python.h>


//Python functions
PyObject* calibration(PyObject* self, PyObject* args);
PyObject* calibrationFromFile(PyObject* self, PyObject* args);
PyMODINIT_FUNC initprsamp(void);

// SOLVERS
double cgamp (
        int n, int m, int strIndY, int strInda,double *y, int *ir, int *jc,
        double delta, double *pr_prmts,
        int t_max, double eps, double damp, double vnf, int disp, complex double *a
    );

// PRIORS
void prior( complex double r, double sig, double *p_prms, complex double *a, double *c, double vnf);

// CHANNELS
void channel( double y, complex double o, double v, double delta,
        complex double *g, double *dg, double vnf );


// COMMON
void sort_rand( int n, int *seq );
double drand(double dMin, double dMax);
double gaussrand();
void readDimensions(FILE *fp, int *m, int *n, int *p, int *nnz);
void readFromText(FILE *fp, int n, int nnz, int *ir, int *jc, double *p_prms, int *t_max, double *eps, double *damp, double *vnf, double *delta, int *disp);
void readFromTextY(FILE *fp, int m, int p, double *y);
void writeRow(FILE *fp, int n, complex double *a);

static inline double max( double a, double b ) { return a > b ? a : b; }
static inline double min( double a, double b ) { return a < b ? a : b; }
static inline int rand_int( int k ) { return rand() / (RAND_MAX / k + 1); }

#define send_data_tag 2001
#define return_data_tag 2002

