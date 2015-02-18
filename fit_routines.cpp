#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

using namespace std;

#include "declarations_2.3.h"
#include "fit_routines.h"

void Fit_Correlationfunction2D()
{
  const size_t data_length = nq_points1*nq_points2;  // # of points
  const size_t n_para = 3;	// # of parameters
				// n_para = 4 if you want to include lambda

  // allocate space for a covariance matrix of size p by p
  gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);

  // allocate and setup for generating gaussian distibuted random numbers
  gsl_rng_env_setup ();
  const gsl_rng_type *type = gsl_rng_default;
  gsl_rng *rng_ptr = gsl_rng_alloc (type);

  //set up test data
  struct Correlationfunction2D_data Correlfun2D_data;
  Correlfun2D_data.data_length = data_length;
  Correlfun2D_data.q_o = new double [data_length];
  Correlfun2D_data.q_s = new double [data_length];
  Correlfun2D_data.y = new double [data_length];
  Correlfun2D_data.sigma = new double [data_length];

  int idx = 0;
  for(int i=0; i<nq_points1; i++)
  {
    for(int j=0; j<nq_points2; j++)
    {
         Correlfun2D_data.q_o[idx] = q_out[i];
         Correlfun2D_data.q_s[idx] = q_side[j];
         // This sets up the data to be fitted, with gaussian noise added
         // Correlfun2D_data.y[idx] = 1.0*exp( - 0.81*q_out[i]*q_out[i] - 1.21*q_side[j]*q_side[j] - 4.0*q_long[k]*q_long[k] - 0.25*q_out[i]*q_side[j]) + gsl_ran_gaussian(rng_ptr, error);
         Correlfun2D_data.y[idx] = Correl_2D[i][j];
         Correlfun2D_data.sigma[idx] = Correl_2D_err[i][j];
         //Correlfun2D_data.sigma[idx] = 1e-2;
         idx++;
    }
  }

  double para_init[n_para] = { 1.0, 1.0, 1.0 };  // initial guesses of parameters

  gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);
  
  // set up the function to be fit 
  gsl_multifit_function_fdf target_func;
  target_func.f = &Fittarget_correlfun2D_f;        // the function of residuals
  target_func.df = &Fittarget_correlfun2D_df;      // the gradient of this function
  target_func.fdf = &Fittarget_correlfun2D_fdf;    // combined function and gradient
  target_func.n = data_length;              // number of points in the data set
  target_func.p = n_para;              // number of parameters in the fit function
  target_func.params = &Correlfun2D_data;  // structure with the data and error bars

  const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver *solver_ptr 
       = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
  gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);

  size_t iteration = 0;         // initialize iteration counter
  print_fit_state_2D (iteration, solver_ptr);
  int status;  		// return value from gsl function calls (e.g., error)
  do
  {
      iteration++;
      
      // perform a single iteration of the fitting routine
      status = gsl_multifit_fdfsolver_iterate (solver_ptr);

      // print out the status of the fit
      cout << "status = " << gsl_strerror (status) << endl;

      // customized routine to print out current parameters
      print_fit_state_2D (iteration, solver_ptr);

      if (status)    // check for a nonzero status code
      {
          break;  // this should only happen if an error code is returned 
      }

      // test for convergence with an absolute and relative error (see manual)
      status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, 
                                        fit_tolerance, fit_tolerance);
  }
  while (status == GSL_CONTINUE && iteration < fit_max_iterations);

  // calculate the covariance matrix of the best-fit parameters
  gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

  // print out the covariance matrix using the gsl function (not elegant!)
  cout << endl << "Covariance matrix: " << endl;
  gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

  cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
  cout.precision (5);		                // # of digits in doubles

  int width = 7;		// setw width for output
  cout << endl << "Best fit results:" << endl;
  cout << "Ro = " << setw (width) << get_fit_results (0, solver_ptr)
    << " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

  cout << "Rs      = " << setw (width) << get_fit_results (1, solver_ptr)
    << " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;
  
  cout << "Ros      = " << setw (width) << get_fit_results (2, solver_ptr)
    << " +/- " << setw (width) << get_fit_err (2, covariance_ptr) << endl;
    
  cout << "status = " << gsl_strerror (status) << endl;
  cout << "--------------------------------------------------------------------" << endl;

  double chi = gsl_blas_dnrm2(solver_ptr->f);
  double dof = data_length - n_para;
  double c = GSL_MAX_DBL(1, chi/sqrt(dof));

  R_out_Correl = fabs(get_fit_results(0, solver_ptr))*hbarC;
  R_side_Correl = fabs(get_fit_results(1, solver_ptr))*hbarC;
  R_os_Correl = fabs(get_fit_results(2, solver_ptr))*hbarC;
  R_out_Correl_err = c*get_fit_err(0, covariance_ptr)*hbarC;
  R_side_Correl_err = c*get_fit_err(1, covariance_ptr)*hbarC;
  R_os_Correl_err = c*get_fit_err(2, covariance_ptr)*hbarC;

  cout << "final results: " << endl;
  cout << scientific << setw(10) << setprecision(5) 
       << "chisq/dof = " << chi*chi/dof << endl;
  cout << scientific << setw(10) << setprecision(5);
  cout << " R_out = " << R_out_Correl << " +/- " << R_out_Correl_err << endl;
  cout << " R_side = " << R_side_Correl << " +/- " << R_side_Correl_err << endl;
  cout << " R_os = " << R_os_Correl << " +/- " << R_os_Correl_err << endl;

  //clean up
  gsl_matrix_free (covariance_ptr);
  gsl_rng_free (rng_ptr);

  delete[] Correlfun2D_data.q_o;
  delete[] Correlfun2D_data.q_s;
  delete[] Correlfun2D_data.y;
  delete[] Correlfun2D_data.sigma;

  gsl_multifit_fdfsolver_free (solver_ptr);  // free up the solver

  return;
}

//*********************************************************************
// 2D case
//*********************************************************************
//  Simple function to print results of each iteration in nice format
int print_fit_state_2D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
  cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
  cout.precision (5);		// digits in doubles

  int width = 15;		// setw width for output
  cout << scientific
    << "iteration " << iteration << ": "
    << "  x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
    << setw (width) << gsl_vector_get (solver_ptr->x, 1)
    << setw (width) << gsl_vector_get (solver_ptr->x, 2)
    << "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
    << endl << endl;

  return 0;
}
//*********************************************************************
//  Function returning the residuals for each point; that is, the 
//  difference of the fit function using the current parameters
//  and the data to be fit.
int Fittarget_correlfun2D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
  size_t n = ((struct Correlationfunction2D_data *) params_ptr)->data_length;
  double *q_o = ((struct Correlationfunction2D_data *) params_ptr)->q_o;
  double *q_s = ((struct Correlationfunction2D_data *) params_ptr)->q_s;
  double *y = ((struct Correlationfunction2D_data *) params_ptr)->y;
  double *sigma = ((struct Correlationfunction2D_data *) params_ptr)->sigma;

  //fit parameters
  double R_o = gsl_vector_get (xvec_ptr, 0);
  double R_s = gsl_vector_get (xvec_ptr, 1);
  double R_os = gsl_vector_get (xvec_ptr, 2);

  size_t i;

  for (i = 0; i < n; i++)
  {
      double Yi = exp(- q_s[i]*q_s[i]*R_s*R_s
                   - q_o[i]*q_o[i]*R_o*R_o - q_o[i]*q_s[i]*R_os*R_os);
      gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
  }

  return GSL_SUCCESS;
}

//*********************************************************************
//  Function returning the Jacobian of the residual function
int Fittarget_correlfun2D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
  size_t n = ((struct Correlationfunction2D_data *) params_ptr)->data_length;
  double *q_o = ((struct Correlationfunction2D_data *) params_ptr)->q_o;
  double *q_s = ((struct Correlationfunction2D_data *) params_ptr)->q_s;
  double *sigma = ((struct Correlationfunction2D_data *) params_ptr)->sigma;

  //fit parameters
  double R_o = gsl_vector_get (xvec_ptr, 0);
  double R_s = gsl_vector_get (xvec_ptr, 1);
  double R_os = gsl_vector_get (xvec_ptr, 2);

  size_t i;

  for (i = 0; i < n; i++)
  {
      // Jacobian matrix J(i,j) = dfi / dxj, 
      // where fi = (Yi - yi)/sigma[i],      
      //       Yi = A * exp(-lambda * i) + b 
      // and the xj are the parameters (A,lambda,b) 
      double sig = sigma[i];

      //derivatives
      double common_elemt = exp(- q_s[i]*q_s[i]*R_s*R_s
                   - q_o[i]*q_o[i]*R_o*R_o - q_o[i]*q_s[i]*R_os*R_os);
      
      gsl_matrix_set (Jacobian_ptr, i, 0, - q_o[i]*q_o[i]*2.0*R_o*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 1, - q_s[i]*q_s[i]*2.0*R_s*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 2, - q_o[i]*q_s[i]*2.0*R_os*common_elemt/sig);
  }
  return GSL_SUCCESS;
}
//*********************************************************************
//  Function combining the residual function and its Jacobian
int Fittarget_correlfun2D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
  Fittarget_correlfun2D_f(xvec_ptr, params_ptr, f_ptr);
  Fittarget_correlfun2D_df(xvec_ptr, params_ptr, Jacobian_ptr);

  return GSL_SUCCESS;
}

//*********************************************************************
//  Function to return the i'th best-fit parameter
inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr)
{
  return gsl_vector_get (solver_ptr->x, i);
}

//*********************************************************************
//  Function to retrieve the square root of the diagonal elements of
//   the covariance matrix.
inline double get_fit_err (int i, gsl_matrix * covariance_ptr)
{
  return sqrt (gsl_matrix_get (covariance_ptr, i, i));
}

//End of file
