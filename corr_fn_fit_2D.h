#ifndef CORR_FN_FIT_2D_H
#define CORR_FN_FIT_2D_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#include <time.h>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <complex>
#include <limits>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

using namespace std;

typedef struct
{
   double qx, qy, qz;
   double Corr_fn, error;
}corr_fcn_data;

struct data
{
	size_t n;
	double * y;
	double * sigma;
};

vector<double>* corr_fn_fit_2D(vector<double>* passed_corr_fn_ptr, vector<double>* q1_pts_ptr, vector<double>* q2_pts_ptr);
int print_state (size_t iter, gsl_multifit_fdfsolver * s);
int my_gaussian_f (const gsl_vector * x, void *params, gsl_vector * f);
int my_gaussian_df (const gsl_vector * x, void *params, gsl_matrix * J);
int my_gaussian_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);
int get_filelength(string filepath);
void read_corrfn(int length, corr_fcn_data* corrfn_ptr);
void read_corrfn(corr_fcn_data* corrfn_ptr, vector<double>* passed_corr_fn_ptr, vector<double>* q1_pts_ptr, vector<double>* q2_pts_ptr);

#endif
