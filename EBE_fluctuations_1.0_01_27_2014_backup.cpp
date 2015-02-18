//Program: EBE_fluctuations_1.0
//Author: Christopher Plumberg
//Date: January 24, 2014

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
#include "generate_vector_1D.h"
//#include "filename_new.h"
//#include "outfilename_1.0.h"
#include "qng_1d_vec.h"
#include "gauss.h"
#include "fourier_kernels.h"

//function prototype(s)
void initialize_parameters();
void set_gaussian_points();
//void fill_pos_coords(double,double,double,double,double,double,double,double);
void fill_pos_coords(void * params_ptr, void * params_ptr);
void fill_ens_avg_vector();
double S_function (double r, void *);
bool is_multiple(int x, int y);
double pdf(double r, void * params_ptr);
double get_normalization(double (*function_original) (double, double, double, void *), void * params_ptr);
int index(int i, int j, int k, int l);
void fill_fourier_transform_ens_avg(int dim);
void fill_fourier_transform_model_function(int i, int dim, void * params_ptr);
void fill_ens_avg_FTmf_vector();
double ens_avg_FT2_integrand(double fluct, void * params_ptr);
double ens_avg_integrand(double fluct, void * params_ptr);

int
main (int argc, char *argv[])
{
grid_ptr = new vector<position> (order*order*order*order);
ens_avg_ptr = new vector<double> (order);
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

if (coords != 0) { cerr << "Only Cartesian coordinates presently supported!  Ending run..." << endl;
				return (1);
			}

lower_limit_vec.x = 0.;
upper_limit_vec.x = 20.;
fluct_lower = my_pdf_params.mean - my_pdf_params.sigma;
fluct_upper = my_pdf_params.mean + my_pdf_params.sigma;

	//set the points needed to perform the integrations
	set_gaussian_points();
cout << "Checkpoint 2" << endl;

	//define position vectors over coordinates corresponding to Gaussian points
	//needed for integrations
	//fill_pos_coords(x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, t_lower, t_upper);
	fill_pos_coords(&lower_limit_vec, &upper_limit_vec);
cout << "Checkpoint 3" << endl;

	//define a vector over spatial coordinates for doing remaining integrals
	fill_ens_avg_vector();
cout << "Checkpoint 4" << endl;

	//fourier transform ensemble average over q
	//fill vector over q values
	fill_fourier_transform_ens_avg(1);
cout << "Checkpoint 5" << endl;

	//compute normalization of ensemble average
	double ens_avg_norm = real((*FT_ens_avg_ptr).at(0));
//cout << "ens_avg_norm = " << ens_avg_norm << endl;
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

	lower_limit_vec.r = 0.;
	lower_limit_vec.theta = 0.;
	lower_limit_vec.phi = -PI;
	lower_limit_vec.eta = -limit;
	lower_limit_vec.tau = 0.;
	lower_limit_vec.x = -4.*limit;
	lower_limit_vec.y = -4.*limit;
	lower_limit_vec.z = -4.*limit;
	lower_limit_vec.t = -limit;

	upper_limit_vec.r = 4.*limit;
	upper_limit_vec.theta = PI;
	upper_limit_vec.phi = PI;
	upper_limit_vec.eta = limit;
	upper_limit_vec.tau = 4.*limit;
	upper_limit_vec.x = 4.*limit;
	upper_limit_vec.y = 4.*limit;
	upper_limit_vec.z = 4.*limit;
	upper_limit_vec.t = limit;

/*
double r_lower=0., r_upper=4.*limit, phi_lower=-PI, phi_upper=PI;
double eta_lower=-limit, eta_upper=limit, tau_lower=0., tau_upper=4.*limit;
double t_lower=-limit, t_upper=limit, x_lower=-4.*limit, x_upper=4.*limit;
double y_lower=-limit, y_upper=limit, z_lower=-4.*limit, z_upper=4.*limit;
double theta_lower=0., theta_upper=PI;
*/

	phi_interval = (phi_upper-phi_lower)/double(order);
}

/************************************************************************/

void set_gaussian_points()
{
	//calculate zeroes of legendre polynomials (abscissae) and corresponding weights
	double xpts[order];
	double wts[order];
	gauss(order,3,-1,1,xpts,wts);
	//cout << "i \t xi[i] \t wi[i]" << endl;
	for(int i = 0; i < half_order; i++) {
		xi[i] = xpts[i+half_order+1];
		wi[i] = wts[i+half_order+1];
		//cout << i << "\t" << xi[i] << "\t" << wi[i] << endl;
	}
	wi[half_order] = wts[half_order];

	//cout<< "Computed Gaussian abscissae and corresponding weights..." << endl;
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


double S_function (double r, void * params_ptr)
{
//assume 2D Gaussian, instantaneous freeze-out for simplicity --> t-dependence suppressed for timebeing
//for now, assume only fluctuations of normalization
	model_params params = * (struct model_params *) params_ptr;

	double R = params.sigma;
	double N = params.norm;
	double r0 = params.mean;

	const double result = ( N / (2.*PI*R*R) ) * exp(-(r-r0)*(r-r0)/(2.*R*R));

	return (result);
}


double pdf(double fluct, void * params_ptr)
{
//currently set to uniform distribution --> generalize eventually
	pdf_params params = * (struct pdf_params *) params_ptr;
	double width = params.sigma;
	double center = params.mean;

	const double result = ( abs(fluct-center) > width ) ? 0. : 1/(2.*width);

	return (result);
	//return (1.);
}

void fill_pos_coords(void * params_ptr1, void * params_ptr2)
{
	position lower_limits = * (struct position *) params_ptr1;
	position upper_limits = * (struct position *) params_ptr2;
	double a1 = lower_limits.r;
	double b1 = upper_limits.r;
	double a2 = lower_limits.phi;
	double b2 = upper_limits.phi;
	double a3 = lower_limits.theta;
	double b3 = upper_limits.theta;
	double a4 = lower_limits.t;
	double b4 = upper_limits.t;


	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	const double half_length2 = 0.5 * (b2 - a2);
	const double center2 = 0.5 * (a2 + b2);
	const double half_length3 = 0.5 * (b3 - a3);
	const double center3 = 0.5 * (a3 + b3);
	const double half_length4 = 0.5 * (b4 - a4);
	const double center4 = 0.5 * (a4 + b4);

	for (int i = 0; i <= order-1; i++) {
		double abscissa1;
		if (i < half_order) {abscissa1 = -half_length1 * xi[i];}
		else if (i == half_order) {abscissa1 = 0.;}
		else {abscissa1 = half_length1 * xi[order-1-i];}
	for (int j = 0; j <= order-1; j++) {
		double abscissa2;
		if (j < half_order) {abscissa2 = -half_length2 * xi[j];}
		else if (j == half_order) {abscissa2 = 0.;}
		else {abscissa2 = half_length2 * xi[order-1-j];}
	for (int k = 0; k <= order-1; k++) {
		double abscissa3;
		if (k < half_order) {abscissa3 = -half_length3 * xi[k];}
		else if (k == half_order) {abscissa3 = 0.;}
		else {abscissa3 = half_length3 * xi[order-1-k];}
	for (int l = 0; l <= order-1; l++) {
		double abscissa4;
		if (l < half_order) {abscissa4 = -half_length4 * xi[l];}
		else if (l == half_order) {abscissa4 = 0.;}
		else {abscissa4 = half_length4 * xi[order-1-l];}
		switch (coords)
		{
			case 0: 
				((*grid_ptr).at(index(i,j,k,l))).x = center1 + abscissa1;
				((*grid_ptr).at(index(i,j,k,l))).y = center2 + abscissa2;
				((*grid_ptr).at(index(i,j,k,l))).z = center3 + abscissa3;
				((*grid_ptr).at(index(i,j,k,l))).t = center4 + abscissa4;
				break;
			case 1: 
				((*grid_ptr).at(index(i,j,k,l))).r = center1 + abscissa1;
				((*grid_ptr).at(index(i,j,k,l))).phi = center2 + abscissa2;
				((*grid_ptr).at(index(i,j,k,l))).theta = center3 + abscissa3;
				((*grid_ptr).at(index(i,j,k,l))).t = center4 + abscissa4;
				break;
			case 2: 
				((*grid_ptr).at(index(i,j,k,l))).r = center1 + abscissa1;
				((*grid_ptr).at(index(i,j,k,l))).phi = center2 + abscissa2;
				((*grid_ptr).at(index(i,j,k,l))).z = center3 + abscissa3;
				((*grid_ptr).at(index(i,j,k,l))).t = center4 + abscissa4;
				break;
		}
	}
	}
	}
	}
}

void fill_ens_avg_vector()
{
	const double x_half_length = 0.5 * (x_upper - x_lower);
	const double x_center = 0.5 * (x_lower + x_upper);
	double abscissa;
	position my_pos;
	for (int i = 0; i <= order-1; i++)
	{
		if (i < half_order) {abscissa = -x_half_length * xi[i];}
		else if (i == half_order) {abscissa = 0.;}
		else {abscissa = x_half_length * xi[order-1-i];}
		my_pos.x = x_center + abscissa;

		vector<double>* temp_ptr;
		temp_ptr = new vector<double> (order);

		generate_vector(&ens_avg_integrand, &my_pos, temp_ptr, fluct_lower, fluct_upper);
		(*ens_avg_ptr).at(i) = qng_1d_vec(temp_ptr, fluct_lower, fluct_upper);

		delete temp_ptr;
	}
}

double ens_avg_integrand(double fluct, void * params_ptr)
{
	position my_pos = * (struct position *) params_ptr;
	double x = my_pos.x;
	model_params my_model_params_copy = my_model_params;
	my_model_params_copy.sigma = fluct;

	double S_val = S_function(x, &my_model_params_copy);
	double prob = pdf(fluct, &my_pdf_params);

	return (prob*S_val);
}

void fill_fourier_transform_ens_avg(int dim)
{
	for (int i = 0;  i <= nq_points - 1; i++)
	{
		vector<double>* temp_re_ptr;
		temp_re_ptr = new vector<double> (order);
		vector<double>* temp_im_ptr;
		temp_im_ptr = new vector<double> (order);
		vector<double>* q_ptr;
		q_ptr = new vector<double> (dim);
		vector<double>* x_ptr;
		x_ptr = new vector<double> (dim);
		(*q_ptr).at(0) = q_lower + i*q_interval;
		(*q_ptr).at(1) = 0.;

		for (int j = 0; j <= order-1; j++)
		{
		for (int k = 0; k <= order-1; k++)
			(*x_ptr).at(0) = (*grid_ptr).at(index(j,0,k,0)).r;
			(*x_ptr).at(1) = (*grid_ptr).at(index(j,0,k,0)).phi;
			(*temp_re_ptr).at(j) = fourier_kernel_real(q_ptr,x_ptr,2)*(*ens_avg_ptr).at(j);
			(*temp_im_ptr).at(j) = fourier_kernel_imag(q_ptr,x_ptr,2)*(*ens_avg_ptr).at(j);
		}
		}
		(*FT_ens_avg_ptr).at(i) = qng_1d_vec(temp_re_ptr, x_lower, x_upper) + complex_i*qng_1d_vec(temp_im_ptr, x_lower, x_upper);

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
		temp_ptr = new vector<double> (order);

		generate_vector(&ens_avg_FT2_integrand, &my_pos, temp_ptr, fluct_lower, fluct_upper);
		(*FT_model_function2_ens_avg_ptr).at(i) = qng_1d_vec(temp_ptr, fluct_lower, fluct_upper);

		delete temp_ptr;
	}
}

double ens_avg_FT2_integrand(double fluct, void * params_ptr)
{
	position my_pos = * (struct position *) params_ptr;
	double x = my_pos.x;
	int i = my_pos.q_index;
	model_params my_model_params_copy = my_model_params;
	my_model_params_copy.sigma = fluct;  //this will change in general!!

	fill_fourier_transform_model_function(i, 1, &my_model_params_copy);
	//fourier transform of given EBE emission function is stored as
	//function of q in FT_model_function_ptr

	double prob = pdf(fluct, &my_pdf_params);
	//probability of given value of fluctuation

	const double result = prob * real((*FT_model_function_ptr).at(i) * conj((*FT_model_function_ptr).at(i)));

	return (result);
}

void fill_fourier_transform_model_function(int i, int dim, void * params_ptr)
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

		for (int j = 0; j <= order-1; j++)
		{
			(*x_ptr).at(0) = (*grid_ptr).at(index(j,0,k,0)).r;
			(*x_ptr).at(1) = (*grid_ptr).at(index(j,0,k,0)).phi;
			(*temp_re_ptr).at(j) = fourier_kernel_real(q_ptr,x_ptr,2)*S_function((*x_ptr).at(0), params_ptr);
			(*temp_im_ptr).at(j) = fourier_kernel_imag(q_ptr,x_ptr,2)*S_function((*x_ptr).at(0), params_ptr);
		}
		(*FT_model_function_ptr).at(i) = qng_1d_vec(temp_re_ptr, x_lower, x_upper) + complex_i*qng_1d_vec(temp_im_ptr, x_lower, x_upper);

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

int index(int i, int j, int k, int l)
{
	return (order*order*order*i+order*order*j+order*k+l);
}

//End of file
