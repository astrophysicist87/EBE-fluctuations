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
#include "outfilename_1.0.h"
#include "qng_1d_vec.h"
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
vector<double> cartesian_to_polar(double, double);
vector<double> polar_to_cartesian(double, double);
void fill_q_vector();
void get_ens_avg_norm();
vector<double> get_corr_fn_moments(int k1, int k2);
double do_integrations(vector<double>* params_ptr);
vector<double> get_corr_fn_moments();
int local_index(int i, int j);

vector<int> convert_base(int number, int base)
{
	const int length = 1 + int(log(number)/log(base));
	vector<int> result(length);
	for (int dim = 1; dim <= length; dim++) result[dim-1] = (number%int(pow(base,dim)))/int(pow(base,dim-1));

	return (result);
}

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
xiq_0pinf_ptr = new vector<double> (order);
wiq_0pinf_ptr = new vector<double> (order);
xiq_minfpinf_ptr = new vector<double> (order);
wiq_minfpinf_ptr = new vector<double> (order);
grid_ptr = new vector<position> (int(pow(order,dim)));
ens_avg_ptr = new vector<double> (int(pow(order,dim)));
if (program_function == 0) {
	Cbar_ptr = new vector<double> (int(pow(nq_points,dim)));
	C_ens_avg_ptr = new vector<double> (int(pow(nq_points,dim)));
	q_vec_ptr = new vector<qposition> (int(pow(nq_points,dim)));
	FT_ens_avg_ptr = new vector<complex<double> > (int(pow(nq_points,dim)));
	FT_model_function_ptr = new vector<complex<double> > (int(pow(nq_points,dim)));
	FT_ens_avg2_ptr = new vector<double> (int(pow(nq_points,dim)));
	FT_model_function2_ens_avg_ptr = new vector<double> (int(pow(nq_points,dim)));
	}
else if (program_function == 1 || program_function == 2) {
	Cbar_ptr = new vector<double> (int(pow(order,dim)));
	C_ens_avg_ptr = new vector<double> (int(pow(order,dim)));
	q_vec_ptr = new vector<qposition> (int(pow(order,dim)));
	FT_ens_avg_ptr = new vector<complex<double> > (int(pow(order,dim)));
	FT_model_function_ptr = new vector<complex<double> > (int(pow(order,dim)));
	FT_ens_avg2_ptr = new vector<double> (int(pow(order,dim)));
	FT_model_function2_ens_avg_ptr = new vector<double> (int(pow(order,dim)));
	}

/*int number, base;
cout << "Enter a number in base 10: " << endl;
cin >> number;
cout << "Enter a base to convert it to: " << endl;
cin >> base;
vector<int> new_number = convert_base(number,base);
int length = new_number.size();
cout << "The result is ";
for (int i = 0; i < length; i++) cout << new_number[length-i-1] << " ";
cout << endl;
if (1) return (0);*/

//take arguments from command-line
if (argc == 1) { //if no arguments included in command-line, run with defaults
	sigma_cmd = sigma_default;
	mean_cmd = mean_default;
	R_cmd = R_default;
	N_cmd = N_default;
	r0_cmd = r0_default;
	eps_n_bar_cmd = eps_n_bar_default;
	harmonic_n_cmd = harmonic_n_default;
	psi_n_bar_cmd = psi_n_bar_default;
}
else if (is_multiple(argc , 2)) {
	//if incorrect number of arguments, end program and return (1)
	cerr << "Incorrect number of arguments: expected (function_call) (int param_key) (double param_val)" << endl
	     << "param_key: 1 - sigma, 2 - mean, 3 - R, 4 - N, 5 - r0" << endl
	     << "6 - eps_n_bar, 7 - harmonic_n, 8 - psi_n_bar" << endl;

	return (1);
}
else {
	sigma_cmd = sigma_default;
	mean_cmd = mean_default;
	R_cmd = R_default;
	N_cmd = N_default;
	r0_cmd = r0_default;
	eps_n_bar_cmd = eps_n_bar_default;
	harmonic_n_cmd = harmonic_n_default;
	psi_n_bar_cmd = psi_n_bar_default;

int arg_index = 0;
do {
	switch (atoi(argv[arg_index+1]))
	{
		case 1:
			sigma_cmd = atof(argv[arg_index+2]);
			break;
		case 2:
			mean_cmd = atof(argv[arg_index+2]);
			break;
		case 3:
			R_cmd = atof(argv[arg_index+2]);
			break;
		case 4:
			N_cmd = atof(argv[arg_index+2]);
			break;
		case 5:
			r0_cmd = atof(argv[arg_index+2]);
			break;
		case 6:
			eps_n_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 7:
			harmonic_n_cmd = atof(argv[arg_index+2]);
			break;
		case 8:
			psi_n_bar_cmd = atof(argv[arg_index+2]);
			break;
	}
arg_index += 2;  //advance to the next two input arguments
  } while (arg_index < argc-1);
}
cout << "sigma = " << sigma_cmd << endl
	<< "mean = " << mean_cmd << endl
	<< "R = " << R_cmd << endl
	<< "N = " << N_cmd << endl
	<< "r0 = " << r0_cmd << endl
	<< "eps_n_bar = " << eps_n_bar_cmd << endl
	<< "harmonic_n = " << harmonic_n_cmd << endl
	<< "psi_n_bar = " << psi_n_bar_cmd << endl;

	string output_filename = "corr_fns";
	string paramsfile = output_filename + "_params.dat";
	string corrfnsfile = output_filename + ".dat";
	time_t start, end, start2, end2;
	time(&start);

	//set parameters in model function and pdf function
	initialize_parameters();

	if ( !fexists( paramsfile.data() ) ) print_params_to_output(paramsfile, &my_pdf_params, &my_model_params);
	ofstream corrfnsoutput (corrfnsfile.data());

cout << "Checkpoint 1" << endl;

if (coords >= 2) { cerr << "Only Cartesian and polar coordinates presently supported!  Ending run..." << endl;
				return (1);
			}

	//set the points needed to perform the integrations
	set_gaussian_points();
cout << "Checkpoint 2" << endl;

	//set all values of q needed for desired computations
	fill_q_vector();

cout << "Checkpoint 2b" << endl;

	//define position vectors over coordinates corresponding to Gaussian points needed for integrations
	fill_pos_coords();
cout << "Checkpoint 3" << endl;

	//define a vector over spatial coordinates for doing remaining integrals
	fill_ens_avg_vector();
cout << "Checkpoint 4" << endl;

cout << "normalization is " << check_normalization() << endl;

//if (1) return (0);

	//fourier transform ensemble average over q fill vector over q values
	fill_fourier_transform_ens_avg();
cout << "Checkpoint 5" << endl;

	//compute normalization of ensemble average
	if (ens_avg_norm == 0.) get_ens_avg_norm();
cout << "ens_avg_norm = " << ens_avg_norm << endl;
cout << "Checkpoint 6" << endl;

	//computes pointers to hold | |^2 vectors
	for (int i = 0; i <= (*q_vec_ptr).size()-1; i++)
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
for (int i = 0; i <= (*q_vec_ptr).size()-1; i++)
{
	(*Cbar_ptr).at(i) = 1 + (*FT_ens_avg2_ptr).at(i) / (ens_avg_norm*ens_avg_norm);
	(*C_ens_avg_ptr).at(i) = 1 + (*FT_model_function2_ens_avg_ptr).at(i) / (ens_avg_norm*ens_avg_norm);
	//cout << q_lower + i * q_interval << "\t" << (*Cbar_ptr).at(i) << "\t" << (*C_ens_avg_ptr).at(i) << endl;
	//corrfnsoutput << q_lower + i * q_interval << "\t" << (*Cbar_ptr).at(i) << "\t" << (*C_ens_avg_ptr).at(i) << endl;
	//cout << (*q_vec_ptr).at(i).qcoord1 << "\t" << (*q_vec_ptr).at(i).qcoord2 << "\t" << (*Cbar_ptr).at(i) << "\t" << (*C_ens_avg_ptr).at(i) << endl;
	corrfnsoutput << (*q_vec_ptr).at(i).qcoord1 << "\t" << (*q_vec_ptr).at(i).qcoord2 << "\t" << (*Cbar_ptr).at(i) << "\t" << (*C_ens_avg_ptr).at(i) << endl;
}
cout << "Checkpoint 8" << endl;

if (program_function == 1 || program_function == 2 || program_function == 0)
{
	vector<double> result1 = get_corr_fn_moments(1, 1);
	vector<double> result2 = get_corr_fn_moments(1, 2);	// = get_corr_fn_moments(2, 1) by symmetry
	vector<double> result3 = get_corr_fn_moments(2, 2);

	cout << "RESULTS:" << endl;
	cout << result1.at(0) << "\t" << result1.at(1) << endl;
	cout << result2.at(0) << "\t" << result2.at(1) << endl;
	cout << result3.at(0) << "\t" << result3.at(1) << endl;

/*vector<double> result = get_corr_fn_moments();
	cout << "RESULTS:" << endl;
	cout << result.at(0) << "\t" << result.at(1) << endl;*/
}

	time(&end);
	return (0);
}

/************************************************************************/

void initialize_parameters()
{
	//Define default parameters for probability density function here
	my_pdf_params.sigma = sigma_cmd;
	my_pdf_params.mean = mean_cmd;

	//Define default parameters for EBE emission function here
	my_model_params.R = R_cmd;
	my_model_params.N = N_cmd;
	my_model_params.r0 = r0_cmd;
	my_model_params.eps_n_bar = eps_n_bar_cmd;
	my_model_params.harmonic_n = harmonic_n_cmd;
	my_model_params.psi_n_bar = psi_n_bar_cmd;

	//integration limits for ensemble-averaging (assumes uniform distribution)
	fluct_lower = my_pdf_params.mean - my_pdf_params.sigma;
	fluct_upper = my_pdf_params.mean + my_pdf_params.sigma;
}

/************************************************************************/

void set_gaussian_points()
{
	//calculate zeroes of legendre polynomials (abscissae) and corresponding weights
	const double midpoint = my_model_params.R;
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

//do it for q-integrations too
	const double midpointq = 1./midpoint;
		//half-infinite
		(*xiq_0pinf_ptr)[i] = midpointq * (1. + xpts[i]) / (1. - xpts[i]);
		(*wiq_0pinf_ptr)[i] = 2. * midpointq * wts[i] / ( (1. - xpts[i]) * (1. - xpts[i]) );
//if (i==0) cout << (*xi_0pinf_ptr)[i] << "\t" << (*wi_0pinf_ptr)[i] << endl;

		//full-infinite
		(*xiq_minfpinf_ptr)[i] = midpointq * xpts[i] / (1. - xpts[i] * xpts[i]);
		(*wiq_minfpinf_ptr)[i] = midpointq * wts[i] * (1. + xpts[i] * xpts[i]) / ( (1. - xpts[i] * xpts[i]) * (1. - xpts[i] * xpts[i]) );
//if (i==0) cout << (*xi_minfpinf_ptr)[i] << "\t" << (*wi_minfpinf_ptr)[i] << endl;
	}

	//cout<< "Computed Gaussian abscissae and corresponding weights..." << endl;
}

void fill_q_vector()
{
	vector<double>* temp_q_pts_ptr1;
	vector<double>* temp_q_pts_ptr2;
	double half_length1 = 1., half_length2 = 1., q1_lower, q1_upper, q2_lower, q2_upper, center1, center2;
	double q1_interval, q2_interval;

	switch (coords)
	{
		case 0:
			//temp_q_pts_ptr1 = xiq_minfpinf_ptr;  //testing...
			temp_q_pts_ptr1 = xi_ptr;
			half_length1 = 0.5 * (q_upper - q_lower);
			center1 = 0.5 * (q_upper + q_lower);
			half_length2 = 0.5 * (q_upper - q_lower);
			center2 = 0.5 * (q_upper + q_lower);
			//temp_q_pts_ptr2 = xiq_minfpinf_ptr;	//infinite limits don't work so well...
			temp_q_pts_ptr2 = xi_ptr;			//using finite limits instead...
			q1_lower = q_lower;
			q2_lower = q_lower;
			q1_upper = q_upper;
			q2_upper = q_upper;
			break;
		case 1:
			//temp_q_pts_ptr1 = xiq_0pinf_ptr;	//infinite limits don't work so well...
			temp_q_pts_ptr1 = xi_ptr;			//using finite limits instead...
			temp_q_pts_ptr2 = xi_ptr;
			half_length1 = 0.5 * (q_upper - q_lower);
			center1 = 0.5 * (q_upper + q_lower);
			half_length2 = PI;
			q1_lower = q_lower;
			q2_lower = -PI;
			q1_upper = q_upper;
			q2_upper = PI;
			break;
		case 2:
			temp_q_pts_ptr1 = xiq_0pinf_ptr;
			temp_q_pts_ptr2 = xi_ptr;
			half_length2 = PI;
			q1_lower = q_lower;
			q2_lower = -PI;
			q1_upper = q_upper;
			q2_upper = PI;
			break;
	}



	if (program_function == 0) {
		q1_interval = (q1_upper - q1_lower) / double(nq_points-1);
		q2_interval = (q2_upper - q2_lower) / double(nq_points-1);
	for (int i=0; i<int(pow(nq_points,dim)); i++)
		{
			(*q_vec_ptr).at(i).qcoord1 = q1_lower + double(i % nq_points)*q1_interval;
			if (dim >= 2) (*q_vec_ptr).at(i).qcoord2 = q2_lower + double(i / nq_points)*q2_interval;
		}
	}
	else if (program_function == 1 || program_function == 2) {
	for (int i=0; i<int(pow(order,dim)); i++)
		{
			(*q_vec_ptr).at(i).qcoord1 = center1 + half_length1 * (*temp_q_pts_ptr1).at(i % order);
			if (dim >= 2) (*q_vec_ptr).at(i).qcoord2 = half_length2 * (*temp_q_pts_ptr2).at(i / order);
		}
	}
}

void fill_pos_coords()
{
	vector<double>* local_xpts_ptr1;
	vector<double>* local_wts_ptr1;
	vector<double>* local_xpts_ptr2;
	vector<double>* local_wts_ptr2;

	double half_length1, center1, half_length2, center2;
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

	for (int i = 0; i <= order-1; i++) {
		double abscissa1 = half_length1 * (*local_xpts_ptr1)[i];
	for (int j = 0; j <= order-1; j++) {
		double abscissa2 = half_length2 * (*local_xpts_ptr2)[j];
		((*grid_ptr).at(local_index(i,j))).coord1 = center1 + abscissa1;	//x, r, r
		((*grid_ptr).at(local_index(i,j))).coord2 = center2 + abscissa2;	//y, phi, phi
//cout << ((*grid_ptr).at(local_index(i,j))).coord1 << endl;
	}
	}
}

/*
double S_function (vector<double>* x_ptr, void * params_ptr)
{
//assume 1D Gaussian, instantaneous freeze-out for simplicity --> t-dependence suppressed for timebeing
//for now, assume only fluctuations of normalization
	model_params params = * (struct model_params *) params_ptr;

	double R = params.R;
	double N = params.N;
	double r0 = params.r0;

	double r = (*x_ptr).at(0);

	//N.B. - double dummy is not used for anything!!

	const double result = ( N / (sqrt(2.*PI)*R) ) * exp(-(r-r0)*(r-r0)/(2.*R*R));
//cout << "The result is " << result << endl;
//cout << "Specifically, N = " << N << ", R = " << R << ", r0 = " << r0 << " and r = " << r << endl;

	return (result);
	//return (1.);
	//return (r);
	//return (exp(-r*r));
	//return (N);
}*/


double S_function (vector<double>* x_ptr, void * params_ptr)
{
//assume 2D Gaussian, instantaneous freeze-out for simplicity --> t-dependence suppressed for timebeing
//for now, assume only fluctuations of normalization
	model_params params = * (struct model_params *) params_ptr;

	double R = params.R;
	double N = params.N;
	double r0 = params.r0;
	double eps_n_bar = params.eps_n_bar;
	double harmonic_n = params.harmonic_n;
	double psi_n_bar = params.psi_n_bar;

	double r = (*x_ptr).at(0);
	double phi = (*x_ptr).at(1);

	const double result = sqrt(1-4.*eps_n_bar*eps_n_bar)*( N / (2.*PI*R*R) ) * exp(-(r*r)/(2.*R*R)*(1.+2.*eps_n_bar*cos(harmonic_n*(phi-psi_n_bar))));

	return (result);
	//return (exp(-r*r));
}

/*
double S_function (double x, double y, void * params_ptr)
{
//assume 2D Gaussian, instantaneous freeze-out for simplicity --> t-dependence suppressed for timebeing
//for now, assume only fluctuations of normalization
	model_params params = * (struct model_params *) params_ptr;

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
	double half_length1, center1, half_length2, center2;

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

	double abscissa1, abscissa2, abscissa3, abscissa4;
	position my_pos;
	for (int i = 0; i <= order-1; i++)
	{
		abscissa1 = half_length1 * (*local_xpts_ptr1).at(i);
		my_pos.coord1 = center1 + abscissa1;
	for (int j = 0; j <= order-1; j++)
	{
		//if (dim < 2 && j > 0) break;
		abscissa2 = half_length2 * (*local_xpts_ptr2).at(j);
		my_pos.coord2 = center2 + abscissa2;

		(*ens_avg_ptr).at(local_index(i,j)) = integrate_1D(&ens_avg_integrand, &my_pos, fluct_lower, fluct_upper, 0);  //0 indicates finite integration interval
		//(*ens_avg_ptr).at(i) = integrate_1D(&ens_avg_integrand, &my_pos, fluct_lower, fluct_upper, 0);  //0 indicates finite integration interval
//cout << my_pos.coord2 << endl;
	}
	}
}

double ens_avg_integrand(double fluct, void * params_ptr)
{
	position my_pos = * (struct position *) params_ptr;
vector<double>* temp_ptr;
temp_ptr = new vector<double> (dim);
(*temp_ptr).at(0) = my_pos.coord1;
if (dim > 1) (*temp_ptr).at(1) = my_pos.coord2;
	model_params my_model_params_copy = my_model_params;
	switch(parameter_to_vary)
	{
		case 3:
		my_model_params_copy.R = fluct;
		break;
		case 4:
		my_model_params_copy.N = fluct;
		break;
		case 5:
		my_model_params_copy.r0 = fluct;
		break;
		case 6:
		my_model_params_copy.eps_n_bar = fluct;
		break;
		case 7:
		my_model_params_copy.harmonic_n = fluct;
		break;
		case 8:
		my_model_params_copy.psi_n_bar = fluct;
		break;
	}

	double S_val = S_function(temp_ptr, &my_model_params_copy);
	double prob = pdf(fluct, &my_pdf_params);
//cout << "fluct = " << fluct << endl;

	return (prob*S_val);
	//return (prob*exp(-my_pos.coord1*my_pos.coord1));
}

double check_normalization()
{
	vector<double>* temp_ptr;
	temp_ptr = new vector<double> (int(pow(order,dim)));

	for (int j = 0; j <= order-1; j++)
	{
		for (int k = 0; k <= order-1; k++)
		{
			(*temp_ptr).at(local_index(j,k)) = (*xi_0pinf_ptr)[j]*(*ens_avg_ptr).at(local_index(j,k));
		}
	}

	//TWO full-infinite integration intervals
	//double result = qng_2d_vec(temp_ptr, 0., 0., -PI, PI, 1, 0);  //this changes with coords!!
	//double result = qng_2d_vec(temp_ptr, 0., 0., 0., 0., 2, 2);  //this changes with coords!!
	double result = do_integrations(temp_ptr);

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

	return (result);
}

void fill_fourier_transform_ens_avg()
{
	for (int i = 0;  i <= (*q_vec_ptr).size() - 1; i++)
	{
		vector<double>* temp_re_ptr;
		temp_re_ptr = new vector<double> (int(pow(order,dim)));
		vector<double>* temp_im_ptr;
		temp_im_ptr = new vector<double> (int(pow(order,dim)));
		vector<double>* q_ptr;
		q_ptr = new vector<double> (dim);
		vector<double>* x_ptr;
		x_ptr = new vector<double> (dim);
		(*q_ptr).at(0) = ((*q_vec_ptr).at(i)).qcoord1;
		if (dim > 1) (*q_ptr).at(1) = ((*q_vec_ptr).at(i)).qcoord2;

		for (int j = 0; j <= order-1; j++)
		{
		for (int k = 0; k <= order-1; k++)
		{
		//if (dim < 2 && k > 0) continue;
			(*x_ptr).at(0) = (*grid_ptr).at(local_index(j,k)).coord1;
			if (dim > 1) (*x_ptr).at(1) = (*grid_ptr).at(local_index(j,k)).coord2;
			(*temp_re_ptr).at(local_index(j,k)) = fourier_kernel_real(q_ptr, x_ptr, dim)*(*ens_avg_ptr).at(local_index(j,k));
			(*temp_im_ptr).at(local_index(j,k)) = fourier_kernel_imag(q_ptr, x_ptr, dim)*(*ens_avg_ptr).at(local_index(j,k));
		}
		}

		//TWO full-infinite integration intervals
		(*FT_ens_avg_ptr).at(i) = do_integrations(temp_re_ptr) + complex_i*do_integrations(temp_im_ptr);

		delete temp_re_ptr;
		delete temp_im_ptr;
		delete q_ptr;
		delete x_ptr;
	}
}

void get_ens_avg_norm()
{
	vector<double>* temp_ptr;
	temp_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* q_ptr;
	q_ptr = new vector<double> (dim);
	vector<double>* x_ptr;
	x_ptr = new vector<double> (dim);
	(*q_ptr).at(0) = 0.;
	if (dim > 1) (*q_ptr).at(1) = 0.;

		for (int j = 0; j <= order-1; j++)
		{
		for (int k = 0; k <= order-1; k++)
		{
		//if (dim < 2 && k > 0) continue;
		(*x_ptr).at(0) = (*grid_ptr).at(local_index(j,k)).coord1;
		if (dim > 1) (*x_ptr).at(1) = (*grid_ptr).at(local_index(j,k)).coord2;
		(*temp_ptr).at(local_index(j,k)) = fourier_kernel_real(q_ptr, x_ptr, dim)*(*ens_avg_ptr).at(local_index(j,k));
//cout << j << "\t" << k << "\t" << (*ens_avg_ptr).at(local_index(j,k)) << "\t" << fourier_kernel_real(q_ptr, x_ptr, dim) << endl;
	}
	}
	//TWO full-infinite integration intervals
	double temp = do_integrations(temp_ptr);
	//double temp = qng_2d_vec(temp_ptr, 0., 0., -PI, PI, 1, 0);
	ens_avg_norm = temp;

	delete temp_ptr;
	delete q_ptr;
	delete x_ptr;
}

void fill_ens_avg_FTmf_vector()
{
	for ( int i = 0;  i <= (*q_vec_ptr).size()-1; i++)
	{
		my_pos.q_index = i;
		(*FT_model_function2_ens_avg_ptr).at(i) = integrate_1D(&ens_avg_FT2_integrand, &my_pos, fluct_lower, fluct_upper, 0);  //finite integration interval
	}
}

double ens_avg_FT2_integrand(double fluct, void * params_ptr)
{
	position my_pos = * (struct position *) params_ptr;
	int i = my_pos.q_index;
	model_params my_model_params_copy = my_model_params;
	switch(parameter_to_vary)
	{
		case 3:
		my_model_params_copy.R = fluct;
		break;
		case 4:
		my_model_params_copy.N = fluct;
		break;
		case 5:
		my_model_params_copy.r0 = fluct;
		break;
		case 6:
		my_model_params_copy.eps_n_bar = fluct;
		break;
		case 7:
		my_model_params_copy.harmonic_n = fluct;
		break;
		case 8:
		my_model_params_copy.psi_n_bar = fluct;
		break;
	}

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
		(*q_ptr).at(0) = ((*q_vec_ptr).at(i)).qcoord1;
		if (dim > 1) (*q_ptr).at(1) = ((*q_vec_ptr).at(i)).qcoord2;
		if (dim > 2) (*q_ptr).at(2) = ((*q_vec_ptr).at(i)).qcoord3;
		if (dim > 3) (*q_ptr).at(3) = ((*q_vec_ptr).at(i)).qcoord4;

		for (int j = 0; j <= order-1; j++)
		{
		for (int k = 0; k <= order-1; k++)
		{
			//cout << "j = " << j << " and k = " << k << endl;
			(*x_ptr).at(0) = (*grid_ptr).at(local_index(j,k)).coord1;
			//cout << (*x_ptr).at(0) << endl;
			if (dim > 1) (*x_ptr).at(1) = (*grid_ptr).at(local_index(j,k)).coord2;
			//cout << "made it?" << endl;
			//(*temp_re_ptr).at(local_index(j,k)) = fourier_kernel_real(q_ptr, x_ptr, dim)*S_function((*x_ptr).at(0), (*x_ptr).at(1), params_ptr);
			(*temp_re_ptr).at(local_index(j,k)) = fourier_kernel_real(q_ptr, x_ptr, dim)*S_function(x_ptr, params_ptr);
			//cout << "j = " << j << " and k = " << k << endl;
			//(*temp_im_ptr).at(local_index(j,k)) = fourier_kernel_imag(q_ptr, x_ptr, dim)*S_function((*x_ptr).at(0), (*x_ptr).at(1), params_ptr);
			(*temp_im_ptr).at(local_index(j,k)) = fourier_kernel_imag(q_ptr, x_ptr, dim)*S_function(x_ptr, params_ptr);
			//cout << "j = " << j << " and k = " << k << endl;
		}
		}
		//TWO full-infinite integration intervals
		(*FT_model_function_ptr).at(i) = do_integrations(temp_re_ptr) + complex_i*do_integrations(temp_im_ptr);
		//(*FT_model_function_ptr).at(i) = qng_2d_vec(temp_re_ptr, 0., 0., -PI, PI, 1, 0) + complex_i*qng_2d_vec(temp_im_ptr, 0., 0., -PI, PI, 1, 0);
		//(*FT_model_function_ptr).at(i) = qng_2d_vec(temp_re_ptr, 0., 0., 0., 0., 2, 2) + complex_i*qng_2d_vec(temp_im_ptr, 0., 0., 0., 0., 2, 2);

		delete temp_re_ptr;
		delete temp_im_ptr;
		delete q_ptr;
		delete x_ptr;
}

vector<double> cartesian_to_polar(double x, double y)
{
	vector<double> result (2);

	result.at(0) = sqrt(x*x+y*y);	//r
	result.at(1) = atan2(y,x);	//phi

	return (result);
}

vector<double> polar_to_cartesian(double r, double phi)
{
	vector<double> result (2);

	result.at(0) = r*cos(phi);	//x
	result.at(1) = r*sin(phi);	//y

	return (result);
}

/*vector<double> get_corr_fn_moments(int k1, int k2)
{
	double num1, num2, den1, den2;
	vector<double>* num1_ptr;
	num1_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* num2_ptr;
	num2_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* den1_ptr;
	den1_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* den2_ptr;
	den2_ptr = new vector<double> (int(pow(order,dim)));
	for (int i=0; i<int(pow(order,dim)); i++)
	{
		double temp_qr = (*q_vec_ptr).at(i).qcoord1;
		double temp_qphi = (*q_vec_ptr).at(i).qcoord2;
		double factor1 = (k1 == 1) ? temp_qr * cos(temp_qphi) : temp_qr * sin(temp_qphi);  //only works in dim = 2, will have to generalize eventually
		double factor2 = (k2 == 1) ? temp_qr * cos(temp_qphi) : temp_qr * sin(temp_qphi);
		(*num1_ptr).at(i) = factor1 * factor2 * ((*Cbar_ptr).at(i) - 1.);
		(*num2_ptr).at(i) = factor1 * factor2 * ((*C_ens_avg_ptr).at(i) - 1.);
		(*den1_ptr).at(i) = (*Cbar_ptr).at(i) - 1.;
		(*den2_ptr).at(i) = (*C_ens_avg_ptr).at(i) - 1.;
	}

	//num1 = qng_2d_vec(num1_ptr, 0., 0., -PI, PI, 1, 0);  //these will change as well with coords and dim!!!
	//num2 = qng_2d_vec(num2_ptr, 0., 0., -PI, PI, 1, 0);
	//den1 = qng_2d_vec(den1_ptr, 0., 0., -PI, PI, 1, 0);
	//den2 = qng_2d_vec(den2_ptr, 0., 0., -PI, PI, 1, 0);
	num1 = do_integrations(num1_ptr);
	num2 = do_integrations(num2_ptr);
	den1 = do_integrations(den1_ptr);
	den2 = do_integrations(den2_ptr);

	cout << "num1 = " << num1 << endl;
	cout << "den1 = " << den1 << endl;
	cout << "num2 = " << num2 << endl;
	cout << "den2 = " << den2 << endl;

	vector<double> result (2);
	result.at(0) = num1/den1;
	result.at(1) = num2/den2;
	return (result);
}*/

vector<double> get_corr_fn_moments(int k1, int k2)
{
	double num1, num2, den1, den2;
	vector<double>* num1_ptr;
	num1_ptr = new vector<double> (int(pow(nq_points,dim)));
	vector<double>* num2_ptr;
	num2_ptr = new vector<double> (int(pow(nq_points,dim)));
	vector<double>* den1_ptr;
	den1_ptr = new vector<double> (int(pow(nq_points,dim)));
	vector<double>* den2_ptr;
	den2_ptr = new vector<double> (int(pow(nq_points,dim)));
	for (int i=0; i<int(pow(nq_points,dim)); i++)
	{
		double temp_qr = (*q_vec_ptr).at(i).qcoord1;
		double temp_qphi = (*q_vec_ptr).at(i).qcoord2;
		double factor1 = (k1 == 1) ? temp_qr * cos(temp_qphi) : temp_qr * sin(temp_qphi);  //only works in dim = 2, will have to generalize eventually
		double factor2 = (k2 == 1) ? temp_qr * cos(temp_qphi) : temp_qr * sin(temp_qphi);
		(*num1_ptr).at(i) = factor1 * factor2 * ((*Cbar_ptr).at(i) - 1.);
		(*num2_ptr).at(i) = factor1 * factor2 * ((*C_ens_avg_ptr).at(i) - 1.);
		(*den1_ptr).at(i) = (*Cbar_ptr).at(i) - 1.;
		(*den2_ptr).at(i) = (*C_ens_avg_ptr).at(i) - 1.;
	}

	//num1 = qng_2d_vec(num1_ptr, 0., 0., -PI, PI, 1, 0);  //these will change as well with coords and dim!!!
	//num2 = qng_2d_vec(num2_ptr, 0., 0., -PI, PI, 1, 0);
	//den1 = qng_2d_vec(den1_ptr, 0., 0., -PI, PI, 1, 0);
	//den2 = qng_2d_vec(den2_ptr, 0., 0., -PI, PI, 1, 0);
	//num1 = do_integrations(num1_ptr);
	//num2 = do_integrations(num2_ptr);
	//den1 = do_integrations(den1_ptr);
	//den2 = do_integrations(den2_ptr);

num1 = q_int2d_alt(num1_ptr);
num2 = q_int2d_alt(num2_ptr);
den1 = q_int2d_alt(den1_ptr);
den2 = q_int2d_alt(den2_ptr);

	cout << "num1 = " << num1 << endl;
	cout << "den1 = " << den1 << endl;
	cout << "num2 = " << num2 << endl;
	cout << "den2 = " << den2 << endl;

	vector<double> result (2);
	result.at(0) = num1/den1;
	result.at(1) = num2/den2;
	return (result);
}

vector<double> get_corr_fn_moments()
{
	double num1 = 0., num2 = 0., den1 = 0., den2 = 0.;
	vector<double>* num1_ptr;
	num1_ptr = new vector<double> (int(pow(nq_points,dim)));
	vector<double>* num2_ptr;
	num2_ptr = new vector<double> (int(pow(nq_points,dim)));
	vector<double>* den1_ptr;
	den1_ptr = new vector<double> (int(pow(nq_points,dim)));
	vector<double>* den2_ptr;
	den2_ptr = new vector<double> (int(pow(nq_points,dim)));
	for (int i=0; i<nq_points; i++)
	{
		double temp_q = (*q_vec_ptr).at(i).qcoord1;
		double factor1 = temp_q;
		(*num1_ptr).at(i) = factor1 * factor1 * ((*Cbar_ptr).at(i) - 1.);
		(*num2_ptr).at(i) = factor1 * factor1 * ((*C_ens_avg_ptr).at(i) - 1.);
		(*den1_ptr).at(i) = (*Cbar_ptr).at(i) - 1.;
		(*den2_ptr).at(i) = (*C_ens_avg_ptr).at(i) - 1.;
	}

	for (int i=0; i<nq_points; i++)
	{
		num1 += (*num1_ptr).at(i);
		num2 += (*num2_ptr).at(i);
		den1 += (*den1_ptr).at(i);
		den2 += (*den2_ptr).at(i);
	}

	num1 -= 0.5*((*num1_ptr).at(0)+(*num1_ptr).at(nq_points-1));
	num2 -= 0.5*((*num2_ptr).at(0)+(*num2_ptr).at(nq_points-1));
	den1 -= 0.5*((*den1_ptr).at(0)+(*den1_ptr).at(nq_points-1));
	den2 -= 0.5*((*den2_ptr).at(0)+(*den2_ptr).at(nq_points-1));

	num1 *= q_interval;
	num2 *= q_interval;
	den1 *= q_interval;
	den2 *= q_interval;

	cout << "num1 = " << num1 << endl;
	cout << "den1 = " << den1 << endl;
	cout << "num2 = " << num2 << endl;
	cout << "den2 = " << den2 << endl;

	vector<double> result (2);
	result.at(0) = num1/den1;
	result.at(1) = num2/den2;
	return (result);
}

/*vector<double> get_corr_fn_moments()
{
double limit = 1.0;
	double num1 = 0., num2 = 0., den1 = 0., den2 = 0.;
	vector<double>* num1_ptr;
	num1_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* num2_ptr;
	num2_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* den1_ptr;
	den1_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* den2_ptr;
	den2_ptr = new vector<double> (int(pow(order,dim)));
	for (int i=0; i<order; i++)
	{
		double temp_q = (*q_vec_ptr).at(i).qcoord1;
		double factor1 = temp_q;
		(*num1_ptr).at(i) = factor1 * factor1 * ((*Cbar_ptr).at(i) - 1.);
		(*num2_ptr).at(i) = factor1 * factor1 * ((*C_ens_avg_ptr).at(i) - 1.);
		(*den1_ptr).at(i) = (*Cbar_ptr).at(i) - 1.;
		(*den2_ptr).at(i) = (*C_ens_avg_ptr).at(i) - 1.;
//		cout << "(1) = " << (*num1_ptr).at(i) << endl
//			<< "(2) = " << (*num2_ptr).at(i) << endl
//			<< "(3) = " << (*den1_ptr).at(i) << endl
//			<< "(4) = " << (*den2_ptr).at(i) << endl;
	}

	num1 = qng_1d_vec_alt(num1_ptr, q_lower, q_upper, 0);
	num2 = qng_1d_vec_alt(num2_ptr, q_lower, q_upper, 0);
	den1 = qng_1d_vec_alt(den1_ptr, q_lower, q_upper, 0);
	den2 = qng_1d_vec_alt(den2_ptr, q_lower, q_upper, 0);
	//num1 = do_integrations(num1_ptr);
	//num2 = do_integrations(num2_ptr);
	//den1 = do_integrations(den1_ptr);
	//den2 = do_integrations(den2_ptr);
	cout << "num1 = " << num1 << endl;
	cout << "den1 = " << den1 << endl;
	cout << "num2 = " << num2 << endl;
	cout << "den2 = " << den2 << endl;

	vector<double> result (2);
	result.at(0) = num1/den1;
	result.at(1) = num2/den2;
	return (result);
}*/

double do_integrations(vector<double>* params_ptr)
{
	double result = 0.;

	switch (coords)
	{
		case 0:
		switch (dim)
		{
			case 1:
			result = qng_1d_vec(params_ptr, 0., 0., 2);
			break;
			case 2:
			result = qng_2d_vec(params_ptr, 0., 0., 0., 0., 2, 2);
			break;
		}
		break;
		case 1:
		switch (dim)
		{
			case 1:
			result = qng_1d_vec(params_ptr, 0., 0., 2);
			break;
			case 2:
			result = qng_2d_vec(params_ptr, 0., 0., -PI, PI, 1, 0);
			break;
		}
		break;
		case 2:
		switch (dim)
		{
			case 1:
			result = qng_1d_vec(params_ptr, 0., 0., 2);
			break;
			case 2:
			result = qng_2d_vec(params_ptr, 0., 0., -PI, PI, 1, 0);
			break;
		}
		break;
	}

	return (result);
}

int local_index(int i, int j)
{
	int local_index;

	switch (dim)
	{
	case 1:
		local_index = i;
	break;
	case 2:
		local_index = index(i,j);
	break;
	}

	return local_index;
}

bool is_multiple(int x, int y)
{
	//x is number which might be a multiple of y
	return (double(x/y) == double(x)/double(y));
}

//End of file
