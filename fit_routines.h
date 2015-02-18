#ifndef FIT_ROUTINES_H
#define FIT_ROUTINES_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "declarations_2.3.h"
//#include "readindata.h"
//#include "parameters.h"
//#include "Arsenal.h"

using namespace std;

struct Correlationfunction2D_data
{
  size_t data_length;
  double *q_o;
  double *q_s;
  double *y;
  double *sigma;
};

int Fittarget_correlfun2D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun2D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun2D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);

      //particle information 
      string particle_name;
      double particle_mass;
      int particle_id;     //particle id
      double particle_sign;   //+/- 1 for Fermi/Bose statistic for baryon/meson
      double particle_gspin;  //particle degeneracy 
     
      //pair momentum
      double K_T, K_phi, K_y;
     
      //Emission function
      double* Emissionfunction_Data;
      double* Emissionfunction_t;
      double* Emissionfunction_x;
      double* Emissionfunction_y;
      double* Emissionfunction_z;
      double* Emissionfunction_Data_CDF;
      int Emissionfunction_length;
      int FO_length;

      double* eta_s;
      double* eta_s_weight;

      double spectra;

      double* q_out;
      double* q_side;

      //store correlation functions
      double* Correl_2D_out;
      double* Correl_2D_out_err;
      double* Correl_2D_side;
      double* Correl_2D_side_err;
      double** Correl_2D;
      double** Correl_2D_err;



      //HBT radii calculated from fitting correlation functions
      double R_out_Correl;
      double R_side_Correl;
      double R_os_Correl;
      double R_out_Correl_err;
      double R_side_Correl_err;
      double R_os_Correl_err;

      double get_Rout_Correl() {return(R_out_Correl);};
      double get_Rout_Correl_err() {return(R_out_Correl_err);};
      double get_Rside_Correl() {return(R_side_Correl);};
      double get_Rside_Correl_err() {return(R_side_Correl_err);};
      double get_Ros_Correl() {return(R_os_Correl);};
      double get_Ros_Correl_err() {return(R_os_Correl_err);};

      void Fit_Correlationfunction2D();
      int print_fit_state_2D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
      inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr);
      inline double get_fit_err (int i, gsl_matrix * covariance_ptr);

#endif
