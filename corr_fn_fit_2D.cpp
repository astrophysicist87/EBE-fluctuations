//Program: EBE_fluctuations_1.0
//Author: Christopher Plumberg
//Date: January 27, 2014

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

#include "corr_fn_fit_2D.h"

size_t n;
string filename;
string path = "/n/home00/plumberg.1/research/ToChris_HBT_package/src_HBT_class/results_MC/";  //or results_MC
const size_t p = 3;
corr_fcn_data * corrfn_ptr;

vector<double>* corr_fn_fit_2D(vector<double>* passed_corr_fn_ptr, vector<double>* q1_pts_ptr, vector<double>* q2_pts_ptr)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

n=(*passed_corr_fn_ptr).size();
cout << "n = " << n << endl;
corrfn_ptr = new corr_fcn_data[n];
read_corrfn(corrfn_ptr, passed_corr_fn_ptr, q1_pts_ptr, q2_pts_ptr);


  int status;
  size_t i, j, k, iter = 0;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  double y[n], sigma[n];

  struct data d = { n, y, sigma};
  
  gsl_multifit_function_fdf f;

  double x_init[p] = { 1.0, 1.0, 1.0 };

  gsl_vector_view x = gsl_vector_view_array (x_init, p);

  const gsl_rng_type * type;
  gsl_rng * r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &my_gaussian_f;
  f.df = &my_gaussian_df;
  f.fdf = &my_gaussian_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  /* This is the data to be fitted */

int idx = 0;
  for (idx = 0; idx < n; idx++)
    {
      y[idx] = (corrfn_ptr[idx]).Corr_fn;
      sigma[idx] = 0.0001;
	idx++;
    }

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  //print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      //printf ("status = %s\n", gsl_strerror (status));

      //print_state (iter, s);

//cout << "jacobian!" << endl;
 // gsl_matrix_fprintf (stdout, s->J, "%g");
//cout << "end jacobian!" << endl;

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
      //status = gsl_multifit_test_delta (s->dx, s->x, 1e-1, 1e-1);
    }
  while (status == GSL_CONTINUE && iter < 500);

print_state (iter, s);

//cout << "status = " << status << endl;

//cout << "s->J: size1 = " << (*(s->J)).size1 << " and size2 = " << (*(s->J)).size2 << endl;
//cout << "covar: size1 = " << (*covar).size1 << " and size2 = " << (*covar).size2 << endl;

  gsl_multifit_covar (s->J, 0.0, covar);

  //gsl_matrix_fprintf (stdout, covar, "%g");

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

//output << FIT(0) << "\t" << ERR(0) << "\t" << FIT(1) << "\t" << ERR(1) << "\t" << FIT(2) << "\t" << ERR(2) << "\t" << FIT(3) << "\t" << ERR(3) << endl;
//cout << FIT(0) << "\t" << ERR(0) << "\t" << FIT(1) << "\t" << ERR(1) << "\t" << FIT(2) << "\t" << ERR(2) << endl;
//cout << "FITS: " << FIT(0) << "\t" << FIT(1) << "\t" << FIT(2) << endl;
//cout << "ERRS: " << ERR(0) << "\t" << ERR(1) << "\t" << ERR(2) << endl;

vector<double>* finalHBTresult_ptr;
finalHBTresult_ptr = new vector<double> (2*p);
for (int idx1 = 0; idx1 < p; idx1++)
{
		(*finalHBTresult_ptr)[2*idx1] = FIT(idx1);
		(*finalHBTresult_ptr)[2*idx1+1] = ERR(idx1);
}

  gsl_multifit_fdfsolver_free (s);

return (finalHBTresult_ptr);
}

int
my_gaussian_f (const gsl_vector * x, void *params, 
        gsl_vector * f)
{
  size_t n = ((struct data *)params)->n;
  double *y = ((struct data *)params)->y;
  double *sigma = ((struct data *) params)->sigma;

  double R2o = gsl_vector_get (x, 0);
  double R2s = gsl_vector_get (x, 1);
  double R2os = gsl_vector_get (x, 2);

int idx = 0;
  for (idx = 0; idx < n; idx++)
    {
	double xp = corrfn_ptr[idx].qx;
	double yp = corrfn_ptr[idx].qy;
      /* Model Yi = exp (-R2o* t * t + R2s) */
      //double t = i;
	//double x = i;
	//double x = q_interval * t + (corrfn_ptr[0]).qval;
	//double t = (corrfn_ptr[i]).qval;
      double Yi = exp (- R2o*xp*xp - R2s*yp*yp - 2.*R2os*xp*yp);
      gsl_vector_set (f, idx, (Yi - y[idx])/sigma[idx]);
	idx++;
    }

  return GSL_SUCCESS;
}

int
my_gaussian_df (const gsl_vector * x, void *params, 
         gsl_matrix * J)
{
  size_t n = ((struct data *)params)->n;
  double *sigma = ((struct data *) params)->sigma;

  double R2o = gsl_vector_get (x, 0);
  double R2s = gsl_vector_get (x, 1);
  double R2os = gsl_vector_get (x, 2);

int idx = 0;
  for (idx = 0; idx < n; idx++)
    {
	double xp = corrfn_ptr[idx].qx;
	double yp = corrfn_ptr[idx].qy;
	//double z = qz_lower + qz_interval * k;
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      //double t = i;
	//double x = q_interval * t + (corrfn_ptr[0]).qval;
	//double x = i;
	//double t = (corrfn_ptr[i]).qval;
      double s = sigma[idx];
      double e = exp(- R2o*xp*xp - R2s*yp*yp - 2.*R2os*xp*yp);
      gsl_matrix_set (J, idx, 0, -xp*xp*e/s); 
      gsl_matrix_set (J, idx, 1, -yp*yp*e/s);
      gsl_matrix_set (J, idx, 2, -2.*xp*yp*e/s);
idx++;
    }
  return GSL_SUCCESS;
}

int
my_gaussian_fdf (const gsl_vector * x, void *params,
          gsl_vector * f, gsl_matrix * J)
{
  my_gaussian_f (x, params, f);
  my_gaussian_df (x, params, J);

  return GSL_SUCCESS;
}

int
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          //gsl_vector_get (s->x, 3),
          gsl_blas_dnrm2 (s->f));
}

void read_corrfn(int length, corr_fcn_data* corrfn_ptr)
{ 
  cout<<"read in correlation function data..." << endl;
  ostringstream corrfn_stream;
  //corrfn_stream << "../EBE_fluctuations/EBE-fluctuations-results/EBE-fluctuations-01-30-2014-1/corr_fns_qslice_0.00.dat";
  corrfn_stream << filename;
  ifstream corrfn(corrfn_stream.str().c_str());
  for(int i=0; i<length; i++)
  {
     corrfn >> corrfn_ptr[i].qx;
     corrfn >> corrfn_ptr[i].qy;
     corrfn >> corrfn_ptr[i].qz;
     corrfn >> corrfn_ptr[i].Corr_fn;
     corrfn >> corrfn_ptr[i].error;
//cout << "read in: " << endl << corrfn_ptr[i].qval << "\t" << corrfn_ptr[i].Cbar << "\t" << corrfn_ptr[i].Cavgd << "\t" << endl;
  }
  corrfn.close();
  cout<<"done"<<endl;
  return;
}

void read_corrfn(corr_fcn_data* corrfn_ptr, vector<double>* passed_corr_fn_ptr, vector<double>* q1_pts_ptr, vector<double>* q2_pts_ptr)
{ 
  for(int i=0; i<n; i++)
  {
     //corrfn_ptr[i].qx = (*q1_pts_ptr)[i];
     //corrfn_ptr[i].qy = (*q2_pts_ptr)[i];
     corrfn_ptr[i].qx = (*q1_pts_ptr)[i] * cos( (*q2_pts_ptr)[i] );  //ONLY WORKS IN 2D POLAR COORDS!!!
     corrfn_ptr[i].qy = (*q1_pts_ptr)[i] * sin( (*q2_pts_ptr)[i] );
     corrfn_ptr[i].Corr_fn = (*passed_corr_fn_ptr)[i]-1.;
//corrfn_ptr[i].Corr_fn = exp( -2.*corrfn_ptr[i].qx*corrfn_ptr[i].qx - 0.5*corrfn_ptr[i].qy*corrfn_ptr[i].qy - 3.*corrfn_ptr[i].qx*corrfn_ptr[i].qy );
     corrfn_ptr[i].error = 1.e-1;
//cout << "read in (in read_corrfn): " << (*q1_pts_ptr)[i] << "               " << (*q2_pts_ptr)[i] << "               " << corrfn_ptr[i].qx  << "               " << corrfn_ptr[i].qy << endl;
//cout << "read in (n = " << n << "): " << "\t" << corrfn_ptr[i].qx << "\t" << corrfn_ptr[i].qy << "\t" << corrfn_ptr[i].Corr_fn << "\t" << endl;
  }
  return;
}

int get_filelength(string filepath)
{
   int length=0; 
   char line[512];
   ostringstream filepath_stream;
   filepath_stream << filepath;
   ifstream infile(filepath_stream.str().c_str());
   //Determine the length of the file
   while (!infile.eof ())
   {
      infile.getline(line, 512);
      length++;
   }
   length = length-1;
   infile.close();
   return(length);
}

//End of file
