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

using namespace std;

#include "declarations_2.3.h"	//always include this one first
#include "corr_fn_fit_2D.h"
#include "emission_function.h"
#include "pdfs.h"
#include "outfilename_2.2.h"
#include "qng_1d_vec_2.2.h"
#include "qng_2d_vec_2.2.h"
#include "integration_routines_2.2.h"
#include "gauss.h"
#include "fourier_kernels.h"
#include "source_variances.h"

int
main (int argc, char *argv[])
{
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
/*	Cbar_ptr = new vector<double> (int(pow(nq_points,dim)));
	C_ens_avg_ptr = new vector<double> (int(pow(nq_points,dim)));
	q_vec_ptr = new vector<qposition> (int(pow(nq_points,dim)));
	FT_ens_avg_ptr = new vector<complex<double> > (int(pow(nq_points,dim)));
	FT_model_function_ptr = new vector<complex<double> > (int(pow(nq_points,dim)));
	FT_ens_avg2_ptr = new vector<double> (int(pow(nq_points,dim)));
	FT_model_function2_ens_avg_ptr = new vector<double> (int(pow(nq_points,dim)));*/
	Cbar_ptr = new vector<double> (nq_points1*nq_points2);
	C_ens_avg_ptr = new vector<double> (nq_points1*nq_points2);
	q1_pts_ptr = new vector<double> (nq_points1*nq_points2);
	q2_pts_ptr = new vector<double> (nq_points1*nq_points2);
	q_vec_ptr = new vector<qposition> (nq_points1*nq_points2);
	FT_ens_avg_ptr = new vector<complex<double> > (nq_points1*nq_points2);
	FT_model_function_ptr = new vector<complex<double> > (nq_points1*nq_points2);
	FT_ens_avg2_ptr = new vector<double> (nq_points1*nq_points2);
	FT_model_function2_ens_avg_ptr = new vector<double> (nq_points1*nq_points2);
	}
else if (program_function == 1 || program_function == 2) {
	Cbar_ptr = new vector<double> (int(pow(order,dim)));
	C_ens_avg_ptr = new vector<double> (int(pow(order,dim)));
	q_vec_ptr = new vector<qposition> (int(pow(order,dim)));
	FT_ens_avg_ptr = new vector<complex<double> > (int(pow(order,dim)));
	FT_ens_avg2_ptr = new vector<double> (int(pow(order,dim)));
	FT_model_function2_ens_avg_ptr = new vector<double> (int(pow(order,dim)));
	}

cout << "q_interval1 = " << q_interval1 << ", q_interval2 = " << q_interval2 << endl;

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
	string corrfnsqmomsfile = output_filename + "_q-moments.dat";
	string corrfnsradiifile = output_filename + "_radii.dat";
	time_t start, end, start2, end2;
	time(&start);

//if (1) return (0);

	//set parameters in model function and pdf function
	initialize_parameters();

	//if ( !fexists( paramsfile.data() ) ) print_params_to_output(paramsfile, &my_pdf_params, &my_model_params);
	print_params_to_output(paramsfile, &my_pdf_params, &my_model_params);
	ofstream corrfnsoutput (corrfnsfile.data());
	ofstream corrfnsqmomsoutput (corrfnsqmomsfile.data());
	ofstream corrfnsradiioutput (corrfnsradiifile.data());

cout << "Checkpoint 1" << endl;

if (coords >= 2) { cerr << "Only Cartesian and polar coordinates presently supported!  Ending run..." << endl;
				return (1);
			}

	//set the points needed to perform the integrations
	set_gaussian_points();
cout << "Checkpoint 2" << endl;

//cout << "The check integral is: " << setprecision(20) << check_qng() << endl;
//cout << "The exact answer is " << PI << endl;

	//set all values of q needed for desired computations
	fill_q_vector();

cout << "Checkpoint 2b" << endl;

	//define position vectors over coordinates corresponding to Gaussian points needed for integrations
	fill_pos_coords();
cout << "Checkpoint 3" << endl;

	//define a vector over spatial coordinates for doing remaining integrals
	fill_ens_avg_vector();
cout << "Checkpoint 4" << endl;

//if (1) return (0);

cout << "normalization is " << check_normalization() << endl;

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
for (int i = 0; i <= nq_points1-1; i++)
{
for (int j = 0; j <= nq_points2-1; j++)
{
	(*Cbar_ptr).at(nq_points2*i+j) = 1 + (*FT_ens_avg2_ptr).at(nq_points2*i+j) / (ens_avg_norm*ens_avg_norm);
	(*C_ens_avg_ptr).at(nq_points2*i+j) = 1 + (*FT_model_function2_ens_avg_ptr).at(nq_points2*i+j) / (ens_avg_norm*ens_avg_norm);
	cout << (*q_vec_ptr).at(nq_points2*i+j).qcoord1 << "\t" << (*q_vec_ptr).at(nq_points2*i+j).qcoord2 << "\t"
		<< (*Cbar_ptr).at(nq_points2*i+j) << "\t" << (*C_ens_avg_ptr).at(nq_points2*i+j) << endl;
	corrfnsoutput << (*q_vec_ptr).at(nq_points2*i+j).qcoord1 << "\t" << (*q_vec_ptr).at(nq_points2*i+j).qcoord2 << "\t"
		<< (*Cbar_ptr).at(nq_points2*i+j) << "\t" << (*C_ens_avg_ptr).at(nq_points2*i+j) << endl;
}
corrfnsoutput << endl;
}
cout << "Checkpoint 8" << endl;
double a, b, c;
if (program_function == 1 || program_function == 2 || program_function == 0)
{
	vector<double> result1 = get_corr_fn_moments(1, 1);
	vector<double> result2 = get_corr_fn_moments(1, 2);	// = get_corr_fn_moments(2, 1) by symmetry
	vector<double> result3 = get_corr_fn_moments(2, 2);

//this part only works for 2D ==> generalize!!!
double det0 = result1[0]*result3[0] - result2[0]*result2[0];
double det1 = result1[1]*result3[1] - result2[1]*result2[1];
	cout << "RESULTS (calculated q-moments):" << endl;
	cout << result1.at(0) << "     " << result1.at(1) << endl;
	cout << result2.at(0) << "     " << result2.at(1) << endl;
	cout << result3.at(0) << "     " << result3.at(1) << endl;
	cout << "RADII from q-moments:" << endl;
	cout << "Cbar:" << endl;
a=0.5*result3[0]/det0;
b=-0.5*result2[0]/det0;
c=0.5*result1[0]/det0;
	cout << a << "     " << b << endl;		//Cbar
	cout << b << "     " << c << endl;		//Cbar
cout << "Resulting determinant = " << (a)*(c) - (b*b) << endl;
cout << "Eigenvalues = (" << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ", " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ")" << endl;
cout << endl;
/***************************************************************************************/
//store radii information in output file
//"QM" = "q-moments"
corrfnsradiioutput << "_QM_CBAR_R2X_: " << a << endl;
corrfnsradiioutput << "_QM_CBAR_R2XY_: " << b << endl;
corrfnsradiioutput << "_QM_CBAR_R2Y_: " << c << endl;
corrfnsradiioutput << "_QM_CBAR_DET_: " << a*c - b*b << endl;
corrfnsradiioutput << "_QM_CBAR_EIG_1_: " << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << "_QM_CBAR_EIG_2_: " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << endl;
/***************************************************************************************/
cout << "Cavgd:" << endl;
a=0.5*result3[1]/det1;
b=-0.5*result2[1]/det1;
c=0.5*result1[1]/det1;
	cout << a << "     " << b << endl;		//Cavgd
	cout << b << "     " << c << endl;		//Cavgd
cout << "Resulting determinant = " << (0.5*result3[1]/det1)*(0.5*result1[1]/det1) - (-0.5*result2[1]/det1)*(-0.5*result2[1]/det1) << endl;
cout << "Eigenvalues = (" << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ", " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ")" << endl;
	corrfnsqmomsoutput << result1.at(0) << "\t" << result1.at(1) << endl;
	corrfnsqmomsoutput << result2.at(0) << "\t" << result2.at(1) << endl;
	corrfnsqmomsoutput << result3.at(0) << "\t" << result3.at(1) << endl;
corrfnsqmomsoutput << endl;
	corrfnsqmomsoutput << 0.5*result3[0]/det0 << "\t" << -0.5*result2[0]/det0 << endl;
	corrfnsqmomsoutput << -0.5*result2[0]/det0 << "\t" << 0.5*result1[0]/det0 << endl;
corrfnsqmomsoutput << endl;
	corrfnsqmomsoutput << 0.5*result3[1]/det1 << "\t" << -0.5*result2[1]/det1 << endl;
	corrfnsqmomsoutput << -0.5*result2[1]/det1 << "\t" << 0.5*result1[1]/det1 << endl;

/***************************************************************************************/
//store radii information in output file
//"QM" = "q-moments"
corrfnsradiioutput << "_QM_CAVGD_R2X_: " << a << endl;
corrfnsradiioutput << "_QM_CAVGD_R2XY_: " << b << endl;
corrfnsradiioutput << "_QM_CAVGD_R2Y_: " << c << endl;
corrfnsradiioutput << "_QM_CAVGD_DET_: " << a*c - b*b << endl;
corrfnsradiioutput << "_QM_CAVGD_EIG_1_: " << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << "_QM_CAVGD_EIG_2_: " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << endl;
/***************************************************************************************/
//vector<double> result = get_corr_fn_moments();
//	cout << "RESULTS:" << endl;
//	cout << result.at(0) << "\t" << result.at(1) << endl;
}

cout << endl << endl;



/********************************************************************/
/**************************SOURCE VARIANCES**************************/
/********************************************************************/
vector<double> exactresult_fromEA(3);
exactresult_fromEA = get_source_variances_fromEA();
cout << "Results from source variances:" << endl;
cout << "exactresult_fromEA[0] = " << exactresult_fromEA[0] << endl
     << "exactresult_fromEA[1] = " << exactresult_fromEA[1] << endl
     << "exactresult_fromEA[2] = " << exactresult_fromEA[2] << endl;
cout << "Resulting determinant = " << (exactresult_fromEA[0])*(exactresult_fromEA[2]) - (exactresult_fromEA[1])*(exactresult_fromEA[1]) << endl;
a=exactresult_fromEA[0];
b=exactresult_fromEA[1];
c=exactresult_fromEA[2];
cout << "Eigenvalues = (" << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ", " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ")" << endl;
/***************************************************************************************/
//store radii information in output file
//"SV" = "source variances"
corrfnsradiioutput << "_SVCBAR_R2X_: " << a << endl;
corrfnsradiioutput << "_SVCBAR_R2XY_: " << b << endl;
corrfnsradiioutput << "_SVCBAR_R2Y_: " << c << endl;
corrfnsradiioutput << "_SVCBAR_DET_: " << a*c - b*b << endl;
corrfnsradiioutput << "_SVCBAR_EIG_1_: " << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << "_SVCBAR_EIG_2_: " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << endl;
/***************************************************************************************/
double inverse_factor = 1./(ens_avg_norm*ens_avg_norm*((*C_ens_avg_ptr).at(0.) - 1.));
vector<double> exactresult_fromEBE(3);
exactresult_fromEBE = get_source_variances_fromEBE();
cout << "exactresult_fromEBE[0] = " << inverse_factor*exactresult_fromEBE[0] << endl
     << "exactresult_fromEBE[1] = " << inverse_factor*exactresult_fromEBE[1] << endl
     << "exactresult_fromEBE[2] = " << inverse_factor*exactresult_fromEBE[2] << endl;
a=inverse_factor*exactresult_fromEBE[0];
b=inverse_factor*exactresult_fromEBE[1];
c=inverse_factor*exactresult_fromEBE[2];
cout << "CHECK: inverse_factor = " << inverse_factor << endl;
cout << "Resulting determinant = " << a*c - b*b << endl;
cout << "Eigenvalues = (" << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ", " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ")" << endl;
/***************************************************************************************/
//store radii information in output file
//"SV" = "source variances"
corrfnsradiioutput << "_SVCAVGD_R2X_: " << a << endl;
corrfnsradiioutput << "_SVCAVGD_R2XY_: " << b << endl;
corrfnsradiioutput << "_SVCAVGD_R2Y_: " << c << endl;
corrfnsradiioutput << "_SVCAVGD_DET_: " << a*c - b*b << endl;
corrfnsradiioutput << "_SVCAVGD_EIG_1_: " << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << "_SVCAVGD_EIG_2_: " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << endl;
/***************************************************************************************/
cout << endl << endl;





corrfnsoutput.close();
corrfnsqmomsoutput.close();

/************************************************************************/
//TEST:
cout << "Gaussian fit to correlation functions (Cbar):" << endl;
vector<double>* GFCbar_ptr;
GFCbar_ptr = corr_fn_fit_2D(Cbar_ptr, q1_pts_ptr, q2_pts_ptr);
cout << "R2x = " << (*GFCbar_ptr)[0] << " +/- " << (*GFCbar_ptr)[1] << " fm^2" << endl;
cout << "R2y = " << (*GFCbar_ptr)[2] << " +/- " << (*GFCbar_ptr)[3] << " fm^2" << endl;
cout << "R2xy = " << (*GFCbar_ptr)[4] << " +/- " << (*GFCbar_ptr)[5] << " fm^2" << endl;
cout << "Resulting determinant = " << ((*GFCbar_ptr)[0])*((*GFCbar_ptr)[2]) - ((*GFCbar_ptr)[4])*((*GFCbar_ptr)[4]) << endl;
a=(*GFCbar_ptr)[0];
b=(*GFCbar_ptr)[4];
c=(*GFCbar_ptr)[2];
cout << "Eigenvalues = (" << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ", " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ")" << endl;
/***************************************************************************************/
//store radii information in output file
//"GF" = "Gaussian fit"
corrfnsradiioutput << "_GF_CBAR_R2X_: " << a << endl;
corrfnsradiioutput << "_GF_CBAR_R2XY_: " << b << endl;
corrfnsradiioutput << "_GF_CBAR_R2Y_: " << c << endl;
corrfnsradiioutput << "_GF_CBAR_DET_: " << a*c - b*b << endl;
corrfnsradiioutput << "_GF_CBAR_EIG_1_: " << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << "_GF_CBAR_EIG_2_: " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << endl;
/***************************************************************************************/
cout << endl << endl;

cout << "Gaussian fit to correlation functions (Cavgd):" << endl;
vector<double>* GFCavgd_ptr;
GFCavgd_ptr = corr_fn_fit_2D(C_ens_avg_ptr, q1_pts_ptr, q2_pts_ptr);
cout << "R2x = " << (*GFCavgd_ptr)[0] << " +/- " << (*GFCavgd_ptr)[1] << " fm^2" << endl;
cout << "R2y = " << (*GFCavgd_ptr)[2] << " +/- " << (*GFCavgd_ptr)[3] << " fm^2" << endl;
cout << "R2xy = " << (*GFCavgd_ptr)[4] << " +/- " << (*GFCavgd_ptr)[5] << " fm^2" << endl;
cout << "Resulting determinant = " << ((*GFCavgd_ptr)[0])*((*GFCavgd_ptr)[2]) - ((*GFCavgd_ptr)[4])*((*GFCavgd_ptr)[4]) << endl;
a=(*GFCavgd_ptr)[0];
b=(*GFCavgd_ptr)[4];
c=(*GFCavgd_ptr)[2];
cout << "Eigenvalues = (" << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ", " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << ")" << endl;
/***************************************************************************************/
//store radii information in output file
//"GF" = "Gaussian fit"
corrfnsradiioutput << "_GF_CAVGD_R2X_: " << a << endl;
corrfnsradiioutput << "_GF_CAVGD_R2XY_: " << b << endl;
corrfnsradiioutput << "_GF_CAVGD_R2Y_: " << c << endl;
corrfnsradiioutput << "_GF_CAVGD_DET_: " << a*c - b*b << endl;
corrfnsradiioutput << "_GF_CAVGD_EIG_1_: " << 0.5*(a + c - sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << "_GF_CAVGD_EIG_2_: " << 0.5*(a + c + sqrt(a*a + 4.*b*b - 2.*a*c + c*c)) << endl;
corrfnsradiioutput << endl;
/***************************************************************************************/
/************************************************************************/

corrfnsradiioutput.close();

	time(&end);
	return (0);
}

/************************************************************************/

void initialize_parameters()
{
	//Define default parameters for probability density function here
	//my_pdf_params.sigma = sigma_cmd;
	//my_pdf_params.mean = mean_cmd;
	my_pdf_params.mean_vec_ptr = &params_mean;
	my_pdf_params.sigma_vec_ptr = &params_sigma;
	my_pdf_params.fluct_lower_ptr = new vector<double> (params_mean.size());
	my_pdf_params.fluct_upper_ptr = new vector<double> (params_mean.size());

	//Define default parameters for EBE emission function here
	my_model_params.R = R_cmd;
	my_model_params.N = N_cmd;
	my_model_params.r0 = r0_cmd;
	my_model_params.eps_n_bar = eps_n_bar_cmd;
	my_model_params.harmonic_n = harmonic_n_cmd;
	my_model_params.psi_n_bar = psi_n_bar_cmd;

	//integration limits for ensemble-averaging (assumes uniform distribution)
	for (int i=0; i < params_mean.size(); i++)
	{
	(*(my_pdf_params.fluct_lower_ptr))[i] = (*(my_pdf_params.mean_vec_ptr))[i] - (*(my_pdf_params.sigma_vec_ptr))[i];
	(*(my_pdf_params.fluct_upper_ptr))[i] = (*(my_pdf_params.mean_vec_ptr))[i] + (*(my_pdf_params.sigma_vec_ptr))[i];
	//(*(my_pdf_params.fluct_lower_ptr))[i] = 0.;
	//(*(my_pdf_params.fluct_upper_ptr))[i] = 0.5;
	//cout << "i = " << i << ": " << (*(my_pdf_params.fluct_lower_ptr))[i] << "\t" << (*(my_pdf_params.fluct_upper_ptr))[i] << endl;
	}
}

/************************************************************************/

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
			center2 = 0.;
			q1_lower = q_lower;
			q2_lower = (nq_points2 == 1) ? 0. : -PI;
			q1_upper = q_upper;
			q2_upper = (nq_points2 == 1) ? 0. : PI;
			break;
		case 2:
			//temp_q_pts_ptr1 = xiq_0pinf_ptr;	//infinite limits don't work so well...
			temp_q_pts_ptr1 = xi_ptr;			//using finite limits instead...
			temp_q_pts_ptr2 = xi_ptr;
			half_length2 = PI;
			q1_lower = q_lower;
			q2_lower = (nq_points2 == 1) ? 0. : -PI;
			q1_upper = q_upper;
			q2_upper = (nq_points2 == 1) ? 0. : PI;
			break;
	}



	if (program_function == 0) {
		q1_interval = (nq_points1 == 1)? q1_upper - q1_lower : (q1_upper - q1_lower) / double(nq_points1-1);
		q2_interval = (nq_points2 == 1)? q2_upper - q2_lower : (q2_upper - q2_lower) / double(nq_points2-1);

	for (int i=0; i<nq_points1*nq_points2; i++)
		(*q1_pts_ptr).at(i) = q1_lower + double(i % nq_points1)*q1_interval;
/***********************************************************************************************************************/
//FIXING BUGS!!!!!!!!!!!
//USED TO BE THIS:
/*
if (dim >= 2) {
	for (int i=0; i<nq_points1*nq_points2; i++)
		(*q2_pts_ptr).at(i) = q2_lower + double(i / nq_points2)*q2_interval;
*/
/***********************************************************************************************************************/
if (dim >= 2) {
	for (int i=0; i<nq_points1*nq_points2; i++)
		(*q2_pts_ptr).at(i) = q2_lower + double(i / nq_points1)*q2_interval;
}
	for (int i=0; i<nq_points1*nq_points2; i++)
		{
			(*q_vec_ptr).at(i).qcoord1 = q1_lower + double(i % nq_points1)*q1_interval;
/***********************************************************************************************************************/
//FIXING BUGS!!!!!!!!!!!
//USED TO BE THIS:
/*
if (dim >= 2) {
	for (int i=0; i<nq_points1*nq_points2; i++)
		(*q2_pts_ptr).at(i) = q2_lower + double(i / nq_points2)*q2_interval;
*/
/***********************************************************************************************************************/
			if (dim >= 2) (*q_vec_ptr).at(i).qcoord2 = q2_lower + double(i / nq_points1)*q2_interval;
		}
	}
	else if (program_function == 1 || program_function == 2) {
	for (int i=0; i<int(pow(order,dim)); i++)
		{
			(*q_vec_ptr).at(i).qcoord1 = center1 + half_length1 * (*temp_q_pts_ptr1).at(i % order);
			if (dim >= 2) (*q_vec_ptr).at(i).qcoord2 = center2 + half_length2 * (*temp_q_pts_ptr2).at(i / order);
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
	//set_intervalvec(&intervalvec);
	intervalvec[0] = 0;
//cout << "CHECK: intervalvec[0] = " << intervalvec[0] << endl;
cout << "the fluctuation limits are " << (*my_pdf_params.fluct_lower_ptr)[0] << " and " << (*my_pdf_params.fluct_upper_ptr)[0] << endl;
cout << (*(my_pdf_params.fluct_lower_ptr))[0] << "\t" << (*(my_pdf_params.fluct_upper_ptr))[0] << endl;
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
		(*ens_avg_ptr).at(local_index(i,j)) = integrate(&ens_avg_integrand, &my_pos, my_pdf_params.fluct_lower_ptr, my_pdf_params.fluct_upper_ptr, &intervalvec);
		//(*ens_avg_ptr).at(i) = integrate_1D(&ens_avg_integrand, &my_pos, fluct_lower, fluct_upper, 0);  //0 indicates finite integration interval
//cout << my_pos.coord1 << "\t" << my_pos.coord2 << "\t\t";
//cout << i << "\t" << j << "\t" << (*ens_avg_ptr).at(local_index(i,j)) << endl;
	}
	}
}

double ens_avg_integrand(vector<double>* fluct_vec_ptr, void * params_ptr)
{
	position my_pos = * (struct position *) params_ptr;
vector<double>* temp_ptr;
temp_ptr = new vector<double> (dim);
(*temp_ptr).at(0) = my_pos.coord1;
if (dim > 1) (*temp_ptr).at(1) = my_pos.coord2;
	model_params my_model_params_copy = my_model_params;
for (int i=0; i<(*fluct_vec_ptr).size(); i++)
{
	switch(paramsvec2vary[i])
	{
		case 3:
		my_model_params_copy.R = (*fluct_vec_ptr)[i];
		break;
		case 4:
		my_model_params_copy.N = (*fluct_vec_ptr)[i];
		break;
		case 5:
		my_model_params_copy.r0 = (*fluct_vec_ptr)[i];
		break;
		case 6:
		my_model_params_copy.eps_n_bar = (*fluct_vec_ptr)[i];
		break;
		case 7:
		my_model_params_copy.harmonic_n = (*fluct_vec_ptr)[i];
		break;
		case 8:
		my_model_params_copy.psi_n_bar = (*fluct_vec_ptr)[i];
		break;
	}
}

	double S_val = S_function(temp_ptr, &my_model_params_copy);
	double prob = pdf(fluct_vec_ptr, &my_pdf_params);
	
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
//cout << j << "\t" << k << "\t" << (*xi_0pinf_ptr)[j] << "\t" << (*ens_avg_ptr).at(local_index(j,k)) << endl;
		}
	}

	//TWO full-infinite integration intervals
	//double result = qng_2d_vec(temp_ptr, 0., 0., -PI, PI, 1, 0);  //this changes with coords!!
	//double result = qng_2d_vec(temp_ptr, 0., 0., 0., 0., 2, 2);  //this changes with coords!!
	double result = do_integrations(temp_ptr);

	delete temp_ptr;

	return (result);
}

double test_integrand(vector<double>* x, void * params_ptr)
{
	double local_x = (*x)[0];
	double local_y = (*x)[1];
	return (exp(-local_x*local_x-local_y*local_y));
	//return (1.);
}

double check_qng()
{
	vector<double> test_fup (2);
	test_fup[0] = 0.;
	test_fup[1] = 0.;
	vector<double> test_fdown (2);
	test_fdown[0] = 0.;
	test_fdown[1] = 0.;
	vector<int> intervalveclocal (2);
	intervalveclocal[0] = 2;
	intervalveclocal[1] = 2;
	//double result = integrate_1D(&test_integrand, &my_pos, 0, 0, 2);  //full-infinite integration interval
	double result = integrate(&test_integrand, &my_pos, &test_fdown, &test_fup, &intervalveclocal);  //full-infinite integration interval

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
cout << "Entering fill_ens_avg_FTmf_vector() now..." << endl;
	//set_intervalvec(&intervalvec);
	intervalvec[0] = 0;
//(*my_pdf_params.fluct_lower_ptr)[0] = 0.;
//(*my_pdf_params.fluct_upper_ptr)[0] = 0.5;
cout << (*(my_pdf_params.fluct_lower_ptr))[0] << "\t" << (*(my_pdf_params.fluct_upper_ptr))[0] << endl;
	for ( int i = 0;  i <= (*q_vec_ptr).size()-1; i++)
	{
		my_pos.q_index = i;

		(*FT_model_function2_ens_avg_ptr).at(i) = integrate(&ens_avg_FT2_integrand, &my_pos, my_pdf_params.fluct_lower_ptr, my_pdf_params.fluct_upper_ptr, &intervalvec);
//cout << "(*FT_model_function2_ens_avg_ptr).at(" << i << ") = " << (*FT_model_function2_ens_avg_ptr).at(i) << endl;
		//(*FT_model_function2_ens_avg_ptr).at(i) = integrate_1D(&ens_avg_FT2_integrand, &my_pos, fluct_lower, fluct_upper, 0);  //finite integration interval
if (is_multiple(i,10))cout << "\tFinished loop " << i << "..." << endl;
	}
cout << "Leaving fill_ens_avg_FTmf_vector() now..." << endl;
}

double ens_avg_FT2_integrand(vector<double>* fluct_vec_ptr, void * params_ptr)
{
	position my_pos = * (struct position *) params_ptr;
	int i = my_pos.q_index;
	model_params my_model_params_copy = my_model_params;
for (int i=0; i<(*fluct_vec_ptr).size(); i++)
{
	switch(paramsvec2vary[i])
	{
		case 3:
		my_model_params_copy.R = (*fluct_vec_ptr)[i];
		break;
		case 4:
		my_model_params_copy.N = (*fluct_vec_ptr)[i];
		break;
		case 5:
		my_model_params_copy.r0 = (*fluct_vec_ptr)[i];
		break;
		case 6:
		my_model_params_copy.eps_n_bar = (*fluct_vec_ptr)[i];
		break;
		case 7:
		my_model_params_copy.harmonic_n = (*fluct_vec_ptr)[i];
		break;
		case 8:
		my_model_params_copy.psi_n_bar = (*fluct_vec_ptr)[i];
		break;
	}
}

	fill_fourier_transform_model_function(i, &my_model_params_copy);
	//fourier transform of given EBE emission function is stored as

	//function of q in FT_model_function_ptr

	double prob = pdf(fluct_vec_ptr, &my_pdf_params);
	//probability of given value of fluctuation
//cout << "prob = " << prob << "\t" << "(*FT_model_function_ptr).at(" << i << ") = " << (*FT_model_function_ptr).at(i) << endl;

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
			//cout << "j = " << j << " and k = " << k << " " << (*q_ptr)[0]<< " " << (*q_ptr)[1] << " " << (*x_ptr)[0] << " " << (*x_ptr)[1] << endl;
			//(*temp_im_ptr).at(local_index(j,k)) = fourier_kernel_imag(q_ptr, x_ptr, dim)*S_function((*x_ptr).at(0), (*x_ptr).at(1), params_ptr);
			(*temp_im_ptr).at(local_index(j,k)) = fourier_kernel_imag(q_ptr, x_ptr, dim)*S_function(x_ptr, params_ptr);
			//cout << "j = " << j << " and k = " << k << endl;
		}
		}
		//TWO full-infinite integration intervals
//cout << do_integrations(temp_re_ptr) << "\t" << do_integrations(temp_im_ptr) << endl;
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

vector<double> get_corr_fn_moments(int k1, int k2)
{
	vector<double>* num1_ptr;
	num1_ptr = new vector<double> (nq_points1*nq_points2);
	vector<double>* num2_ptr;
	num2_ptr = new vector<double> (nq_points1*nq_points2);
	vector<double>* den1_ptr;
	den1_ptr = new vector<double> (nq_points1*nq_points2);
	vector<double>* den2_ptr;
	den2_ptr = new vector<double> (nq_points1*nq_points2);

	for (int i=0; i<nq_points1; i++)
	{
	for (int j=0; j<nq_points2; j++)
	{
		double temp_qr = (*q_vec_ptr).at(nq_points2*i+j).qcoord1;
		double temp_qphi = (*q_vec_ptr).at(nq_points2*i+j).qcoord2;
		double factor1 = (k1 == 1) ? temp_qr * cos(temp_qphi) : temp_qr * sin(temp_qphi);  //only works in dim = 2, will have to generalize eventually
		double factor2 = (k2 == 1) ? temp_qr * cos(temp_qphi) : temp_qr * sin(temp_qphi);
		double jacobian = temp_qr;
		//double jacobian = 1.;
		(*num1_ptr).at(nq_points2*i+j) = jacobian * factor1 * factor2 * ((*Cbar_ptr).at(nq_points2*i+j) - 1.);
		(*num2_ptr).at(nq_points2*i+j) = jacobian * factor1 * factor2 * ((*C_ens_avg_ptr).at(nq_points2*i+j) - 1.);
		(*den1_ptr).at(nq_points2*i+j) = jacobian * ((*Cbar_ptr).at(nq_points2*i+j) - 1.);
		(*den2_ptr).at(nq_points2*i+j) = jacobian * ((*C_ens_avg_ptr).at(nq_points2*i+j) - 1.);
	}
	}

	double num1 = 0., num2 = 0., den1 = 0., den2 = 0.;

	for (int i = 0; i <= nq_points1-1; i++)
	{
	double num1x = 0., num2x = 0., den1x = 0., den2x = 0.;
		for (int j = 0; j <= nq_points2-1; j++)
		{
			double factor1 = 1.;
			double factor2 = 1.;
			if (i == 0 || i == nq_points1-1) factor1 = 0.5;
			if (j == 0 || j == nq_points2-1) factor2 = 0.5;
			num1x += factor1*factor2*(*num1_ptr).at(nq_points2*i+j);
			num2x += factor1*factor2*(*num2_ptr).at(nq_points2*i+j);
			den1x += factor1*factor2*(*den1_ptr).at(nq_points2*i+j);
			den2x += factor1*factor2*(*den2_ptr).at(nq_points2*i+j);
		}
num1 += num1x;
num2 += num2x;
den1 += den1x;
den2 += den2x;
	}
	num1 *= q_interval1 * q_interval2;
	num2 *= q_interval1 * q_interval2;
	den1 *= q_interval1 * q_interval2;
	den2 *= q_interval1 * q_interval2;

	cout << "num1 = " << num1 << endl;
	cout << "den1 = " << den1 << endl;
	cout << "num2 = " << num2 << endl;
	cout << "den2 = " << den2 << endl;
cout << endl;

	vector<double> result (2);
	result.at(0) = num1/den1;
	result.at(1) = num2/den2;
	return (result);
}

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


/*void set_intervalvec(vector<int>* intervalvecptr)
{
	for (int i = 0; i < (*intervalvecptr).size(); i++)
		{
		int pdfswitch = (pdfswitchvec.size() == 1) ? pdfswitchvec[0] : pdfswitchvec[i];
			switch (pdfswitch)
			{
				case 0:		//uniform distribution, integrate from fluct_upper to fluct_lower (finite)
				(*intervalvecptr)[i] = 0;
				break;
				case 1:		//Gaussian distribution, integrate from -inf to +inf (full-infinite)
				(*intervalvecptr)[i] = 2;
				break;
				case 2:		//Bessel-Gaussian distribution, integrate from 0 to +inf (semi-infinite)
				(*intervalvecptr)[i] = 1;
				break;
			}
		}
}*/

bool is_multiple(int x, int y)
{
	//x is number which might be a multiple of y
	return (double(x/y) == double(x)/double(y));
}

//End of file
