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

using namespace std;

#include "declarations.h"
//#include "generate_vector_1D.h"
//#include "generate_vector_2D.h"
//#include "filename_new.h"
//#include "outfilename_1.0.h"
//#include "qng_1d_vec.h"
#include "qng_2d_vec.h"
#include "integration_routines.h"
#include "gauss.h"
#include "fourier_kernels.h"

//function prototype(s)
void initialize_parameters();
void set_gaussian_points();
void fill_pos_coords();
void fill_ens_avg_vector();
double S_function (vector<double>*, void *);
bool is_multiple(int x, int y);
double pdf(double r, void * params_ptr);
double get_normalization(double (*function_original) (double, double, double, void *), void * params_ptr);
int index(int i, int j, int k, int l);
void fill_fourier_transform_ens_avg();
void fill_fourier_transform_model_function(int i, void * params_ptr);
void fill_ens_avg_FTmf_vector();
double ens_avg_FT2_integrand(double fluct, void * params_ptr);
double ens_avg_integrand(double fluct, void * params_ptr);
double check_normalization();
double check_qng();
double test_integrand(double x, double y, void * params_ptr);
void new_points_and_weights(int mod_index);

int
main (int argc, char *argv[])
{
//grid_ptr = new vector<position> (order*order*order*order);
xi_ptr = new vector<double> (order);
wi_ptr = new vector<double> (order);
xi_0pinf_ptr = new vector<double> (order);
wi_0pinf_ptr = new vector<double> (order);
xi_minfpinf_ptr = new vector<double> (order);
wi_minfpinf_ptr = new vector<double> (order);
grid_ptr = new vector<position> (int(pow(order,4)));
ens_avg_ptr = new vector<double> (int(pow(order,dim)));
FT_ens_avg_ptr = new vector<complex<double> > (nq_points);
FT_model_function_ptr = new vector<complex<double> > (nq_points);
FT_ens_avg2_ptr = new vector<double> (nq_points);
FT_model_function2_ens_avg_ptr = new vector<double> (nq_points);
Cbar_ptr = new vector<double> (nq_points);
C_ens_avg_ptr = new vector<double> (nq_points);

/*
//take arguments from command-line
if (argc == 1) { //if no arguments included in command-line, run with defaults
	eps_2_bar_cmd = eps_2_bar_default;
}
else if (is_multiple(argc , 2)) {
	//if incorrect number of arguments, end program and return (1)
	cerr << "Incorrect number of arguments: expected (function_call) (int param_key) (double param_val)" << endl
	     << "param_key: 1 - eps_2_bar" << endl;

	return (1);
}
else {
	eps_2_bar_cmd = eps_2_bar_default;

int arg_index = 0;
do {
	switch (atoi(argv[arg_index+1]))
	{
		case 1:
			eps_2_bar_cmd = atof(argv[arg_index+2]);
			break;
	}

//	cout << "Using eps_2_bar = " << eps_2_bar_cmd << endl << endl;
arg_index += 2;  //advance to the next two input arguments
  } while (arg_index < argc-1);
}
*/
//  string output_filename = outputfilename();  //trying a possible speed-up for code
  time_t start, end, start2, end2;
	time(&start);

	//set parameters in model function and pdf function
	initialize_parameters();
cout << "Checkpoint 1" << endl;

if (coords >= 2) { cerr << "Only Cartesian coordinates presently supported!  Ending run..." << endl;
				return (1);
			}

fluct_lower = my_pdf_params.mean - my_pdf_params.sigma;
fluct_upper = my_pdf_params.mean + my_pdf_params.sigma;

	//set the points needed to perform the integrations
	set_gaussian_points();
cout << "Checkpoint 2" << endl;

	//define position vectors over coordinates corresponding to Gaussian points
	//needed for integrations
	//fill_pos_coords(x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, t_lower, t_upper);
	fill_pos_coords();
cout << "Checkpoint 3" << endl;

	//define a vector over spatial coordinates for doing remaining integrals
	fill_ens_avg_vector();
cout << "Checkpoint 4" << endl;

cout << "normalization is " << check_normalization() << endl;

//if (1) return (0);

	//fourier transform ensemble average over q
	//fill vector over q values
	fill_fourier_transform_ens_avg();
cout << "Checkpoint 5" << endl;

	//compute normalization of ensemble average
	double ens_avg_norm = real((*FT_ens_avg_ptr).at(0));
cout << "ens_avg_norm = " << ens_avg_norm << endl;
cout << "Checkpoint 6" << endl;

	//computes pointers to hold | |^2 vectors
	for (int i = 0; i <= nq_points-1; i++)
	{
		(*FT_ens_avg2_ptr).at(i) = real((*FT_ens_avg_ptr).at(i) * conj((*FT_ens_avg_ptr).at(i)));
	}

/******************************************************************************/

	//ensemble-average the FT of the model-function
	//need to treat entire | |^2 as integrand of ensemble-average
	//--> need to compute for each value of fluctuating quantity in question
	fill_ens_avg_FTmf_vector();
cout << "Checkpoint 7" << endl;

//finally, compute correlation functions
for (int i = 0; i <= nq_points-1; i++)
{
	(*Cbar_ptr).at(i) = 1 + (*FT_ens_avg2_ptr).at(i) / (ens_avg_norm*ens_avg_norm);
	(*C_ens_avg_ptr).at(i) = 1 + (*FT_model_function2_ens_avg_ptr).at(i) / (ens_avg_norm*ens_avg_norm);
	cout << q_lower + i * q_interval << "\t" << (*Cbar_ptr).at(i) << "\t" << (*C_ens_avg_ptr).at(i) << endl;
}
cout << "Checkpoint 8" << endl;

	time(&end);
	return (0);
}

/************************************************************************/

void initialize_parameters()
{
	//Define default parameters for probability density function here
	my_pdf_params.sigma = 4.;
	my_pdf_params.mean = 5.;

	//Define default parameters for EBE emission function here
	my_model_params.sigma = 5.;
	my_model_params.norm = 1.;
	my_model_params.mean = 0.;
}

/************************************************************************/

void set_gaussian_points()
{
	//calculate zeroes of legendre polynomials (abscissae) and corresponding weights
	const double midpoint = my_pdf_params.mean;
	double xpts[order];
	double wts[order];
	gauss(order,3,-1,1,xpts,wts);

	for(int i = 0; i < order; i++) {
		//finite
		(*xi_ptr)[i] = xpts[i];
		(*wi_ptr)[i] = wts[i];
//if (i==0) cout << (*xi_ptr)[i] << "\t" << (*wi_ptr)[i] << endl;

		//half-infinite
		(*xi_0pinf_ptr)[i] = midpoint * (1. + xpts[i]) / (1. - xpts[i]);
		(*wi_0pinf_ptr)[i] = 2. * midpoint * wts[i] / ( (1. - xpts[i]) * (1. - xpts[i]) );
//if (i==0) cout << (*xi_0pinf_ptr)[i] << "\t" << (*wi_0pinf_ptr)[i] << endl;

		//full-infinite
		(*xi_minfpinf_ptr)[i] = midpoint * xpts[i] / (1. - xpts[i] * xpts[i]);
		(*wi_minfpinf_ptr)[i] = midpoint * wts[i] * (1. + xpts[i] * xpts[i]) / ( (1. - xpts[i] * xpts[i]) * (1. - xpts[i] * xpts[i]) );
//if (i==0) cout << (*xi_minfpinf_ptr)[i] << "\t" << (*wi_minfpinf_ptr)[i] << endl;
	}

	cout<< "Computed Gaussian abscissae and corresponding weights..." << endl;
}

void fill_pos_coords()
{
	vector<double>* local_xpts_ptr1;
	vector<double>* local_wts_ptr1;
	vector<double>* local_xpts_ptr2;
	vector<double>* local_wts_ptr2;
	vector<double>* local_xpts_ptr3;
	vector<double>* local_wts_ptr3;
	vector<double>* local_xpts_ptr4;
	vector<double>* local_wts_ptr4;

/*
	double a1 = lower_limits.r;  //all of these things will change with coords as well
	double b1 = upper_limits.r;
	double a2 = lower_limits.phi;
	double b2 = upper_limits.phi;
	double a3 = lower_limits.theta;
	double b3 = upper_limits.theta;
	double a4 = lower_limits.t;
	double b4 = upper_limits.t;
*/
/*
	double a1 = lower_limits.x;  //all of these things will change with coords as well
	double b1 = upper_limits.x;
	double a2 = lower_limits.y;
	double b2 = upper_limits.y;
	double a3 = lower_limits.z;
	double b3 = upper_limits.z;
	double a4 = lower_limits.t;
	double b4 = upper_limits.t;
*/
	double half_length1, center1, half_length2, center2, half_length3, center3, half_length4, center4;
switch (coords)
{
	case 0:
		local_xpts_ptr1 = xi_minfpinf_ptr;
		local_wts_ptr1 = wi_minfpinf_ptr;
		local_xpts_ptr2 = xi_minfpinf_ptr;
		local_wts_ptr2 = wi_minfpinf_ptr;
		local_xpts_ptr3 = xi_minfpinf_ptr;
		local_wts_ptr3 = wi_minfpinf_ptr;
		local_xpts_ptr4 = xi_minfpinf_ptr;
		local_wts_ptr4 = wi_minfpinf_ptr;
		half_length1 = 1.;
		half_length2 = 1.;
		half_length3 = 1.;
		half_length4 = 1.;
		center1 = 0.;
		center2 = 0.;
		center3 = 0.;
		center4 = 0.;
		break;
	case 1:
		local_xpts_ptr1 = xi_0pinf_ptr;
		local_wts_ptr1 = wi_0pinf_ptr;
		local_xpts_ptr2 = xi_ptr;
		local_wts_ptr2 = wi_ptr;
		local_xpts_ptr3 = xi_ptr;
		local_wts_ptr3 = wi_ptr;
		local_xpts_ptr4 = xi_minfpinf_ptr;
		local_wts_ptr4 = wi_minfpinf_ptr;
		half_length1 = 1.;
		half_length2 = 0.5 * (PI - (-PI));
		half_length3 = 0.5 * (PI - 0.);
		half_length4 = 1.;
		center1 = 0.;
		center2 = 0.5 * (PI + (-PI));
		center3 = 0.5 * (PI + 0.);
		center4 = 0.;
		break;
	case 2:
		local_xpts_ptr1 = xi_0pinf_ptr;
		local_wts_ptr1 = wi_0pinf_ptr;
		local_xpts_ptr2 = xi_ptr;
		local_wts_ptr2 = wi_ptr;
		local_xpts_ptr3 = xi_minfpinf_ptr;
		local_wts_ptr3 = wi_minfpinf_ptr;
		local_xpts_ptr4 = xi_minfpinf_ptr;
		local_wts_ptr4 = wi_minfpinf_ptr;
		half_length1 = 1.;
		half_length2 = 0.5 * (PI - (-PI));
		half_length3 = 1.;
		half_length4 = 1.;
		center1 = 0.;
		center2 = 0.5 * (PI + (-PI));
		center3 = 0.;
		center4 = 0.;
		break;
}

	for (int i = 0; i <= order-1; i++) {
		double abscissa1 = half_length1 * (*local_xpts_ptr1)[i];
	for (int j = 0; j <= order-1; j++) {
		double abscissa2 = half_length2 * (*local_xpts_ptr2)[j];
	for (int k = 0; k <= order-1; k++) {
		double abscissa3 = half_length3 * (*local_xpts_ptr3)[k];
	for (int l = 0; l <= order-1; l++) {
		double abscissa4 = half_length4 * (*local_xpts_ptr4)[l];	//coords = 0, 1, 2
		((*grid_ptr).at(index(i,j,k,l))).coord1 = center1 + abscissa1;	//x, r, r
		((*grid_ptr).at(index(i,j,k,l))).coord2 = center2 + abscissa2;	//y, phi, phi
		((*grid_ptr).at(index(i,j,k,l))).coord3 = center3 + abscissa3;	//z, theta, z
		((*grid_ptr).at(index(i,j,k,l))).coord4 = center4 + abscissa4;	//t, t, t
	}
	}
	}
	}
}

/*
double S_function (double r, void * params_ptr)
{
//assume 1D Gaussian, instantaneous freeze-out for simplicity --> t-dependence suppressed for timebeing
//for now, assume only fluctuations of normalization
	model_params params = * (struct model_params *) params_ptr;

	double R = params.sigma;
	double N = params.norm;
	double r0 = params.mean;

	const double result = ( N / (sqrt(2.*PI)*R) ) * exp(-(r-r0)*(r-r0)/(2.*R*R));

	return (result);
	//return (1.);
	//return (r);
	//return (exp(-r*r));
	//return (N);
}
*/

double S_function (double r, double phi, void * params_ptr)
{
//assume 2D Gaussian, instantaneous freeze-out for simplicity --> t-dependence suppressed for timebeing
//for now, assume only fluctuations of normalization
	model_params params = * (struct model_params *) params_ptr;
	//double r = (*x_ptr).at(0);  //of course, this will change depending on coords
	//double phi = (*x_ptr).at(1);

	double R = params.sigma;
	double N = params.norm;
	double r0 = params.mean;

	const double result = ( N / (2.*PI*R*R) ) * exp(-(r-r0)*(r-r0)/(2.*R*R));

	return (result);
}

/*
double S_function (double x, double y, void * params_ptr)
{
//assume 2D Gaussian, instantaneous freeze-out for simplicity --> t-dependence suppressed for timebeing
//for now, assume only fluctuations of normalization
	model_params params = * (struct model_params *) params_ptr;
	//double r = (*x_ptr).at(0);  //of course, this will change depending on coords
	//double phi = (*x_ptr).at(1);

	double R = params.sigma;
	double N = params.norm;
	double r0 = params.mean;

	const double result = ( N / (2.*PI*R*R) ) * exp(-(x*x+y*y)/(2.*R*R));

	return (result);
	//return (1.);
}
*/

double pdf(double fluct, void * params_ptr)
{
//currently set to uniform distribution --> generalize eventually
	pdf_params params = * (struct pdf_params *) params_ptr;
	double width = params.sigma;
	double center = params.mean;

	const double result = ( abs(fluct-center) > width ) ? 0. : 1/(2.*width);

	return (result);
}

void fill_ens_avg_vector()
{
	vector<double>* local_xpts_ptr1;
	vector<double>* local_wts_ptr1;
	vector<double>* local_xpts_ptr2;
	vector<double>* local_wts_ptr2;
	double half_length1, center1, half_length2, center2, half_length3, center3, half_length4, center4;

switch (coords)
{

	case 0:
		local_xpts_ptr1 = xi_minfpinf_ptr;
		local_wts_ptr1 = wi_minfpinf_ptr;
		local_xpts_ptr2 = xi_minfpinf_ptr;
		local_wts_ptr2 = wi_minfpinf_ptr;
		half_length1 = 1.;
		half_length2 = 1.;
		center1 = 0.;
		center2 = 0.;
		break;
	case 1:
		local_xpts_ptr1 = xi_0pinf_ptr;
		local_wts_ptr1 = wi_0pinf_ptr;
		local_xpts_ptr2 = xi_ptr;
		local_wts_ptr2 = wi_ptr;
		half_length1 = 1.;
		half_length2 = 0.5 * (PI - (-PI));
		center1 = 0.;
		center2 = 0.5 * (PI + (-PI));
		break;
	case 2:
		local_xpts_ptr1 = xi_0pinf_ptr;
		local_wts_ptr1 = wi_0pinf_ptr;
		local_xpts_ptr2 = xi_ptr;
		local_wts_ptr2 = wi_ptr;
		half_length1 = 1.;
		half_length2 = 0.5 * (PI - (-PI));
		center1 = 0.;
		center2 = 0.5 * (PI + (-PI));
		break;
}

	double abscissa1, abscissa2;
	position my_pos;
	for (int i = 0; i <= order-1; i++)
	{
		abscissa1 = half_length1 * (*local_xpts_ptr1).at(i);
		my_pos.coord1 = center1 + abscissa1;
	for (int j = 0; j <= order-1; j++)
	{
		abscissa2 = half_length2 * (*local_xpts_ptr2).at(j);
		my_pos.coord2 = center2 + abscissa2;

//cout << "fluct_lower = " << fluct_lower << " and fluct_upper = " << fluct_upper << endl;
		(*ens_avg_ptr).at(index(i,j)) = integrate_1D(&ens_avg_integrand, &my_pos, fluct_lower, fluct_upper, 0);  //0 indicates finite integration interval
		//if (j==41) cout << sqrt(my_pos.coord1*my_pos.coord1 + my_pos.coord2*my_pos.coord2) << "\t" << (*ens_avg_ptr).at(index(i,j)) << endl;
	}
	}
}

double ens_avg_integrand(double fluct, void * params_ptr)
{
	position my_pos = * (struct position *) params_ptr;
	model_params my_model_params_copy = my_model_params;
	my_model_params_copy.sigma = fluct;

	double S_val = S_function(my_pos.coord1, my_pos.coord2, &my_model_params_copy);
	double prob = pdf(fluct, &my_pdf_params);

	return (prob*S_val);
	//return (prob*exp(-my_pos.coord1*my_pos.coord1-my_pos.coord2*my_pos.coord2));
}

double check_normalization()
{
	vector<double>* temp_ptr;
	temp_ptr = new vector<double> (int(pow(order,dim)));

	for (int j = 0; j <= order-1; j++)
	{
		for (int k = 0; k <= order-1; k++)
		{
			(*temp_ptr).at(index(j,k)) = (*ens_avg_ptr).at(index(j,k));
		}
	}

	//TWO full-infinite integration intervals
	double result = qng_2d_vec(temp_ptr, 0., 0., -PI, PI, 1, 0);  //this changes with coords!!
	//cout << x_lower << "\t" << x_upper << "\t" << y_lower << "\t" << y_upper << endl;
	//cout << lower_limit_vec.x << "\t" << upper_limit_vec.x << "\t" << lower_limit_vec.y << "\t" << upper_limit_vec.y << endl;

	delete temp_ptr;

	return (result);
}

double test_integrand(double x, void * params_ptr)
{
	return (exp(-x*x));
}

double check_qng()
{
	double result = integrate_1D(&test_integrand, &my_pos, 0, 0, 2);  //full-infinite integration interval
	cout << x_lower << "\t" << x_upper << "\t" << y_lower << "\t" << y_upper << endl;
	//cout << lower_limit_vec.x << "\t" << upper_limit_vec.x << "\t" << lower_limit_vec.y << "\t" << upper_limit_vec.y << endl;

	return (result);
}

void fill_fourier_transform_ens_avg()
{
	//cout << "size = " << int(pow(order,dim)) << endl;
	for (int i = 0;  i <= nq_points - 1; i++)
	{
		vector<double>* temp_re_ptr;
		temp_re_ptr = new vector<double> (int(pow(order,dim)));
		vector<double>* temp_im_ptr;
		temp_im_ptr = new vector<double> (int(pow(order,dim)));
		vector<double>* q_ptr;
		q_ptr = new vector<double> (dim);
		vector<double>* x_ptr;
		x_ptr = new vector<double> (dim);
		(*q_ptr).at(0) = q_lower + double(i)*q_interval;
		(*q_ptr).at(1) = 0.;

		for (int j = 0; j <= order-1; j++)
		{
		for (int k = 0; k <= order-1; k++)
		{
			//(*x_ptr).at(0) = (*grid_ptr).at(index(j,k,0,0)).r;
			//(*x_ptr).at(1) = (*grid_ptr).at(index(j,k,0,0)).phi;
			(*x_ptr).at(0) = (*grid_ptr).at(index(j,k,0,0)).coord1;
			(*x_ptr).at(1) = (*grid_ptr).at(index(j,k,0,0)).coord2;  //these indices change with dim!!
			(*temp_re_ptr).at(index(j,k)) = fourier_kernel_real(q_ptr, x_ptr, dim)*(*ens_avg_ptr).at(index(j,k));
			(*temp_im_ptr).at(index(j,k)) = fourier_kernel_imag(q_ptr, x_ptr, dim)*(*ens_avg_ptr).at(index(j,k));
		}
		}
		//TWO full-infinite integration intervals
		(*FT_ens_avg_ptr).at(i) = qng_2d_vec(temp_re_ptr, 0., 0., -PI, PI, 1, 0) + complex_i*qng_2d_vec(temp_im_ptr, 0., 0., -PI, PI, 1, 0);
		//if (i==0) cout << (*FT_ens_avg_ptr).at(i) << endl;

		delete temp_re_ptr;
		delete temp_im_ptr;
		delete q_ptr;
		delete x_ptr;
	}
}

void fill_ens_avg_FTmf_vector()
{
	for ( int i = 0;  i <= nq_points-1; i++)
	{
		my_pos.q_index = i;
		vector<double>* temp_ptr;	//vector to hold all possible
						//values of fluctuating quantity
		//temp_ptr = new vector<double> (order);

		//generate_vector(&ens_avg_FT2_integrand, &my_pos, temp_ptr, fluct_lower, fluct_upper);
		//(*FT_model_function2_ens_avg_ptr).at(i) = qng_1d_vec(temp_ptr, fluct_lower, fluct_upper);
		(*FT_model_function2_ens_avg_ptr).at(i) = integrate_1D(&ens_avg_FT2_integrand, &my_pos, fluct_lower, fluct_upper, 0);  //finite integration interval

		//delete temp_ptr;
	}
}

double ens_avg_FT2_integrand(double fluct, void * params_ptr)
{
	position my_pos = * (struct position *) params_ptr;
	int i = my_pos.q_index;
	model_params my_model_params_copy = my_model_params;
	my_model_params_copy.sigma = fluct;  //this will change in general!!

	fill_fourier_transform_model_function(i, &my_model_params_copy);
	//fourier transform of given EBE emission function is stored as
	//function of q in FT_model_function_ptr

	double prob = pdf(fluct, &my_pdf_params);
	//probability of given value of fluctuation

	const double result = prob * real((*FT_model_function_ptr).at(i) * conj((*FT_model_function_ptr).at(i)));

	return (result);
}

void fill_fourier_transform_model_function(int i, void * params_ptr)
{
		vector<double>* temp_re_ptr;
		temp_re_ptr = new vector<double> (int(pow(order,dim)));
		vector<double>* temp_im_ptr;
		temp_im_ptr = new vector<double> (int(pow(order,dim)));
		vector<double>* q_ptr;
		q_ptr = new vector<double> (dim);
		vector<double>* x_ptr;
		x_ptr = new vector<double> (dim);
		(*q_ptr).at(0) = q_lower + i*q_interval;
		(*q_ptr).at(1) = 0.;

		for (int j = 0; j <= order-1; j++)
		{
		for (int k = 0; k <= order-1; k++)
		{
			(*x_ptr).at(0) = (*grid_ptr).at(index(j,k,0,0)).coord1;
			(*x_ptr).at(1) = (*grid_ptr).at(index(j,k,0,0)).coord2;
			//(*x_ptr).at(0) = (*grid_ptr).at(index(j,k,0,0)).r;
			//(*x_ptr).at(1) = (*grid_ptr).at(index(j,k,0,0)).phi;
			(*temp_re_ptr).at(index(j,k)) = fourier_kernel_real(q_ptr, x_ptr, dim)*S_function((*x_ptr).at(0), (*x_ptr).at(1), params_ptr);
			(*temp_im_ptr).at(index(j,k)) = fourier_kernel_imag(q_ptr, x_ptr, dim)*S_function((*x_ptr).at(0), (*x_ptr).at(1), params_ptr);
		}
		}
		//TWO full-infinite integration intervals
		(*FT_model_function_ptr).at(i) = qng_2d_vec(temp_re_ptr, 0., 0., -PI, PI, 1, 0) + complex_i*qng_2d_vec(temp_im_ptr, 0., 0., -PI, PI, 1, 0);

		delete temp_re_ptr;
		delete temp_im_ptr;
		delete q_ptr;
		delete x_ptr;
}

bool is_multiple(int x, int y)
{
	//x is number which might be a multiple of y
	return (double(x/y) == double(x)/double(y));
}

//End of file
