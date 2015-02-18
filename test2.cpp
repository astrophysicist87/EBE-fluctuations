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

typedef struct
{
   double qo, qs, ql;
   double Corr_fn, error;
}corr_fcn_data;

struct data
{
	size_t n;
	double * y;
	double * sigma;
};

  size_t n;
const size_t npts=20;
size_t n1=npts, n2=npts, n3=npts;
string filename;
string path = "/n/home00/plumberg.1/research/ToChris_HBT_package/src_HBT_class/results_MC/";  //or results_MC
//string path = "";  //or results_MC
  const size_t p = 4;
	double q_interval;
corr_fcn_data * corrfn_ptr;
double qo_lower, qo_upper, qs_lower, qs_upper, ql_lower, ql_upper;
double qo_interval, qs_interval, ql_interval;

int print_state (size_t iter, gsl_multifit_fdfsolver * s);
int my_gaussian_f (const gsl_vector * x, void *params, gsl_vector * f);
int my_gaussian_df (const gsl_vector * x, void *params, gsl_matrix * J);
int my_gaussian_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);
int print_state (size_t iter, gsl_multifit_fdfsolver * s);
int get_filelength(string filepath);
void read_corrfn(int length, corr_fcn_data* corrfn_ptr);

int
main (int argc, char *argv[])
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

string output_filename = "HBT_radii_fit.dat";
string outputfilepath = path + output_filename;
ofstream output (outputfilepath.data(), ios::app);

string file = string(argv[1]);

filename = path + file;

n=get_filelength(filename);
//cout << "n = " << n << endl;
corrfn_ptr = new corr_fcn_data[n];
read_corrfn(n, corrfn_ptr);

//n=n1*n2*n3;
qo_lower = (corrfn_ptr)[0].qo;
qo_upper = (corrfn_ptr)[n-1].qo;
qs_lower = (corrfn_ptr)[0].qs;
qs_upper = (corrfn_ptr)[n-1].qs;
ql_lower = (corrfn_ptr)[0].ql;
ql_upper = (corrfn_ptr)[n-1].ql;
qo_interval = (qo_upper - qo_lower)/double(n1-1);
qs_interval = (qs_upper - qs_lower)/double(n2-1);
ql_interval = (ql_upper - ql_lower)/double(n3-1);

  int status;
  size_t i, j, k, iter = 0;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  double y[n], sigma[n];

  struct data d = { n, y, sigma};
  
  gsl_multifit_function_fdf f;

  double x_init[p] = { 1.0, 1.0, 1.0, 1.0 };

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

//read_corrfn(n, corrfn_ptr);
//q_interval = ((corrfn_ptr[n-1]).qval - (corrfn_ptr[0]).qval) / double(n-1);

  /* This is the data to be fitted */

int idx = 0;
  for (i = 0; i < n1; i++)
    {
  for (j = 0; j < n2; j++)
    {
  for (k = 0; k < n3; k++)
    {
	double xo = qo_lower + qo_interval * i;
	double xs = qs_lower + qs_interval * j;
	double xl = ql_lower + ql_interval * k;
      //double t = i;
	//double t = (corrfn_ptr[i]).qval;
      //y[i] = exp (-0.1 * t * t + 0.5) 
      //           + gsl_ran_gaussian(r, 0.1);
      //y[i] = exp (- 10.*x*x*x*x -0.1 * x * x + 0.5);
      y[idx] = (corrfn_ptr[idx]).Corr_fn;
	//y[idx] = exp(-4.*xo*xo - xs*xs - 2.*xl*xl) + gsl_ran_gaussian(r, 0.1);
	//y[i] = exp(-0.1*t*t+0.5)+0.2;
      sigma[idx] = 0.0001;
      //sigma[idx] = (corrfn_ptr[idx]).error;
      //printf("data: %d %g %g %g\n", idx, y[idx], sigma[idx]);
	idx++;
    }
    }
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

      print_state (iter, s);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 5000);
//cout << "status = " << (status == GSL_CONTINUE) << endl;

//cout << "s->J: size1 = " << (*(s->J)).size1 << " and size2 = " << (*(s->J)).size2 << endl;
//cout << "covar: size1 = " << (*covar).size1 << " and size2 = " << (*covar).size2 << endl;

  gsl_multifit_covar (s->J, 0.0, covar);

  gsl_matrix_fprintf (stdout, covar, "%g");

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

output << FIT(0) << "\t" << ERR(0) << "\t" << FIT(1) << "\t" << ERR(1) << "\t" << FIT(2) << "\t" << ERR(2) << "\t" << FIT(3) << "\t" << ERR(3) << endl;

  //printf("R2o	= %.5f +/- %.5f\n", FIT(0), ERR(0));
  //printf("R2s	= %.5f +/- %.5f\n", FIT(1), ERR(1));
  //printf("R2l	= %.5f +/- %.5f\n", FIT(2), ERR(2));
  //printf("R2os	= %.5f +/- %.5f\n", FIT(3), ERR(3));

  //printf ("status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);
output.close();
  return 0;
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
  double R2l = gsl_vector_get (x, 2);
  double R2os = gsl_vector_get (x, 3);

  size_t i, j, k;

int idx = 0;
  for (i = 0; i < n1; i++)
    {
  for (j = 0; j < n2; j++)
    {
  for (k = 0; k < n3; k++)
    {
	double xo = qo_lower + qo_interval * i;
	double xs = qs_lower + qs_interval * j;
	double xl = ql_lower + ql_interval * k;
      /* Model Yi = exp (-R2o* t * t + R2s) */
      //double t = i;
	//double x = i;
	//double x = q_interval * t + (corrfn_ptr[0]).qval;
	//double t = (corrfn_ptr[i]).qval;
      double Yi = exp (- R2o*xo*xo - R2s*xs*xs - R2l*xl*xl - R2os*xo*xs);
      gsl_vector_set (f, idx, (Yi - y[idx])/sigma[idx]);
	idx++;
    }
    }
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
  double R2l = gsl_vector_get (x, 2);
  double R2os = gsl_vector_get (x, 3);

  size_t i, j, k;

int idx = 0;
  for (i = 0; i < n1; i++)
    {
  for (j = 0; j < n2; j++)
    {
  for (k = 0; k < n3; k++)
    {
	double xo = qo_lower + qo_interval * i;
	double xs = qs_lower + qs_interval * j;
	double xl = ql_lower + ql_interval * k;
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      //double t = i;
	//double x = q_interval * t + (corrfn_ptr[0]).qval;
	//double x = i;
	//double t = (corrfn_ptr[i]).qval;
      double s = sigma[idx];
      double e = exp (- R2o*xo*xo - R2s*xs*xs - R2l*xl*xl - R2os*xo*xs);
      gsl_matrix_set (J, idx, 0, -xo*xo*e/s); 
      gsl_matrix_set (J, idx, 1, -xs*xs*e/s);
      gsl_matrix_set (J, idx, 2, -xl*xl*e/s);
      gsl_matrix_set (J, idx, 3, -xo*xs*e/s);
idx++;
    }
    }
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
  printf ("iter: %3u x = % 15.8f % 15.8f "
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
     corrfn >> corrfn_ptr[i].qo;
     corrfn >> corrfn_ptr[i].qs;
     corrfn >> corrfn_ptr[i].ql;
     corrfn >> corrfn_ptr[i].Corr_fn;
     corrfn >> corrfn_ptr[i].error;
//cout << "read in: " << endl << corrfn_ptr[i].qval << "\t" << corrfn_ptr[i].Cbar << "\t" << corrfn_ptr[i].Cavgd << "\t" << endl;
  }
  corrfn.close();
  cout<<"done"<<endl;
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
