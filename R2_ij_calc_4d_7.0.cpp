//Program: R2_s_calc
//Author: Christopher Plumberg
//Date: November 26, 2012
//Modified from integration code written by Dick Furnstahl
//
//Comments: this algorithm only uses constant integration limits!
//May need to generalize to variable boundaries in future.
//
//Does not implement tau integration, since it is trivially (done analytically) in current model

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

const int order = 43;
const int half_order = (order - 1)/2;
const double limit = 5.;
const double PI = 3.14159265358979323846264338327950;
const double M1_4PI = 1./(4.*PI);
const double m_pion = 0.139;  //in GeV
const int angle_format = 2;	//0 == plot vs. Phi_K
				//1 == plot vs. Phi_K - psi_3_bar
				//2 == plot vs. Phi_K - Psi_3

int ecc_calc_flag=0;
double eps_2_bar_cmd;	//defined from 0..0.5
double eps_3_bar_cmd;	//defined from 0..0.5
double v_2_bar_cmd;	//defined from 0..0.5
double v_3_bar_cmd;	//defined from 0..0.5
double psi_2_bar_cmd;	//defined as fraction of PI
double psi_3_bar_cmd;	//defined as fraction of PI
double tau_f_cmd;
double Delta_tau_cmd;
double K_perp_cmd;
double eta_0_cmd;
double Delta_eta_cmd;
double eta_f_cmd;
double Rad_cmd;
double T0_cmd;
double flow_angle_cmd;
const double eps_2_bar_default = 0.;
const double eps_3_bar_default = 0.;
const double v_2_bar_default = 0.;
const double v_3_bar_default = 0.;
const double psi_2_bar_default = 0.;
const double psi_3_bar_default = 0.;
const double tau_f_default = 10.;  //in fm
const double Delta_tau_default = 1.;
const double K_perp_default = 0.5;  //in GeV
const double eta_0_default = 0.;
const double Delta_eta_default = 2.3;
const double eta_f_default = 0.6;
const double Rad_default = 4.5;
const double flow_angle_default = 0.;
const double T0_default = 0.12;
double r_lower=0., r_upper=4.*limit, phi_lower=-PI, phi_upper=PI;
double eta_lower=-limit, eta_upper=limit, tau_lower=0., tau_upper=4.*limit;
double xi[half_order];
double wi[half_order+1];

struct my_params
{
	double alpha;
	double beta;
};

typedef struct
{
	double r, phi, eta, tau;
	double S_val, norm;
	double xs, xo, tc, zc;
	double xs2, xo2, tc2, zc2;
	double x_os, x_ot, x_oz, x_st, x_sz, x_tz;
} Emissionfunction_data;

struct phys_params
{
	double M_perp, T_0, eta_0, Y_rapidity;
	double Rad, eta_f, tau_f, K_perp;
	double Delta_tau, Delta_eta, beta_perp, beta_long;

	double v_3_bar, eps_3_bar, psi_3_bar, Phi_K, v_2_bar, eps_2_bar, psi_2_bar;

	//the following parameters are useful for computational reasons
	double cos_psi, sin_psi, cos_3_psi, sin_3_psi, cos_2_psi, sin_2_psi;
	double M_1_2del_tau2, M_1_2del_eta2, M_1_2Rad2;
};

int index(int i, int j, int k);
int index(int i, int j, int k, int l);

//#include "declarations.h"
//#include "qng_mod.h"
#include "qng_3d_vec.h"
#include "generate_vector_2.4.h"
//#include "filename_new.h"
#include "outfilename_1.0.h"
#include "qng_1d_vec.h"
#include "gauss.h"

//function prototype(s)
phys_params get_parameters();
double S_function (double r, double phi, double eta, double tau, void *);
double eta_t(double r, double phi, void *);
double x2_s_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x_s_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x2_o_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x_o_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double t_c_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double t2_c_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double z2_c_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double z_c_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x_os_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x_ot_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x_oz_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x_st_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x_sz_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double x_tz_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val);
double test_function (double x, double y, double z, double t, void * params_ptr);
bool is_multiple(int x, int y);
double calc_flow_angle(void * params_ptr);
double S_no_tau(double r, double phi, double eta, double Phi_K, void * params_ptr);
double flowangle_integ1a(double r, double phi, double eta, double Phi_K, void * params_ptr);
double flowangle_integ1b(double r, double phi, double eta, double Phi_K, void * params_ptr);
double flowangle_integ2(double r, double phi, double eta, double Phi_K, void * params_ptr);

int
main (int argc, char *argv[])
{

//take arguments from command-line
if (argc == 1) { //if no arguments included in command-line, run with defaults
	eps_2_bar_cmd = eps_2_bar_default;
	eps_3_bar_cmd = eps_3_bar_default;
	v_2_bar_cmd = v_2_bar_default;
	v_3_bar_cmd = v_3_bar_default;
	psi_2_bar_cmd = psi_2_bar_default;
	psi_3_bar_cmd = psi_3_bar_default;
	tau_f_cmd = tau_f_default;
	Delta_tau_cmd = Delta_tau_default;
	K_perp_cmd = K_perp_default;
	eta_0_cmd = eta_0_default;
	Delta_eta_cmd = Delta_eta_default;
	eta_f_cmd = eta_f_default;
	Rad_cmd = Rad_default;
	ecc_calc_flag = 0;  //as default, don't bother calculating eccentricity stuff
	T0_cmd = T0_default;
	flow_angle_cmd = flow_angle_default;
}
else if (is_multiple(argc , 2)) {
	//if incorrect number of arguments, end program and return (1)
	cerr << "Incorrect number of arguments: expected (function_call) (int param_key) (double param_val)" << endl
	     << "param_key: 1 - eps_2_bar, 2 - eps_3_bar, 3 - v_2_bar, 4 - v_3_bar" << endl
	     << "5 - psi_2_bar, 6 - psi_3_bar, 7 - tau_f, 8 - Delta_tau, 9 - K_perp" << endl
	     << "10 - eta_0, 11 - Delta_eta, 12 - eta_f, 13 - Rad, 14 - ecc_calc_flag, 15 - T0" << endl;

	return (1);
}
else {
	eps_2_bar_cmd = eps_2_bar_default;
	eps_3_bar_cmd = eps_3_bar_default;
	v_2_bar_cmd = v_2_bar_default;
	v_3_bar_cmd = v_3_bar_default;
	psi_2_bar_cmd = psi_2_bar_default;
	psi_3_bar_cmd = psi_3_bar_default;
	tau_f_cmd = tau_f_default;
	Delta_tau_cmd = Delta_tau_default;
	K_perp_cmd = K_perp_default;
	eta_0_cmd = eta_0_default;
	Delta_eta_cmd = Delta_eta_default;
	eta_f_cmd = eta_f_default;
	Rad_cmd = Rad_default;
	ecc_calc_flag = 0;  //as default, don't bother calculating eccentricity stuff
	T0_cmd = T0_default;
	flow_angle_cmd = flow_angle_default;

int arg_index = 0;
do {
	switch (atoi(argv[arg_index+1]))
	{
		case 1:
			eps_2_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 2:
			eps_3_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 3:
			v_2_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 4:
			v_3_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 5:
			psi_2_bar_cmd = atof(argv[arg_index+2])*PI;
			break;
		case 6:
			psi_3_bar_cmd = atof(argv[arg_index+2])*PI;
			break;
		case 7:
			tau_f_cmd = atof(argv[arg_index+2]);
			break;
		case 8:
			Delta_tau_cmd = atof(argv[arg_index+2]);
			break;
		case 9:
			K_perp_cmd = atof(argv[arg_index+2]);
			break;
		case 10:
			eta_0_cmd = atof(argv[arg_index+2]);
			break;
		case 11:
			Delta_eta_cmd = atof(argv[arg_index+2]);
			break;
		case 12:
			eta_f_cmd = atof(argv[arg_index+2]);
			break;
		case 13:
			Rad_cmd = atof(argv[arg_index+2]);
			break;
		case 14:
			ecc_calc_flag = atof(argv[arg_index+2]);
			break;
		case 15:
			T0_cmd = atof(argv[arg_index+2]);
			break;
		case 16:
			flow_angle_cmd = atof(argv[arg_index+2]);
			break;
	}

//	cout << "Using eps_2_bar = " << eps_2_bar_cmd << ", eps_3_bar = " << eps_3_bar_cmd << "," << endl
//	     << "v_2_bar = " << v_2_bar_cmd << ", v_3_bar = " << v_3_bar_cmd << "," << endl
//	     << "psi_2_bar = " << psi_2_bar_cmd << ", psi_3_bar = " << psi_3_bar_cmd << endl << endl;
//cerr << "Taking K_perp_cmd = " << K_perp_cmd << endl;
arg_index += 2;  //advance to the next two input arguments
  } while (arg_index < argc-1);
}

  double norm;
  //double tau_interval = (tau_upper-tau_lower)/double(order);
  double tau_center = (tau_upper+tau_lower)*0.5;
  double tau_half_length = (tau_upper-tau_lower)*0.5;
  int n_points = 50;
  phys_params initial_params_copy = get_parameters();
  my_params test_params;
  test_params.alpha = 1.;
  test_params.beta = 0.;
  double beta_perp = initial_params_copy.beta_perp;
  double beta_long = initial_params_copy.beta_long;
  double K_perp = initial_params_copy.K_perp;
  string output_filename = outputfilename();  //trying a possible speed-up for code
  time_t start, end, start2, end2;
  string plumbergteststring;

	//calculate zeroes of legendre polynomials (abscissae) and corresponding weights
	double xpts[order];
	double wts[order];
	//double xi[half_order];
	//double wi[half_order+1];
	gauss(order,3,-1,1,xpts,wts);
	//cout << "i \t xi[i] \t wi[i]" << endl;
	for(int i = 0; i < half_order; i++) {
		xi[i] = xpts[i+half_order+1];
		wi[i] = wts[i+half_order+1];
		//cout << i << "\t" << xi[i] << "\t" << wi[i] << endl;
	}
	wi[half_order] = wts[half_order];

  double psi_3_calculated;
  double Phi_K_a = -PI;
  double Phi_K_b = PI;  //use these def'ns for angle_format==0 as default
  if (angle_format == 1) {
		Phi_K_a = -PI+initial_params_copy.psi_3_bar;
		Phi_K_b = PI+initial_params_copy.psi_3_bar;
	}
  else if (angle_format == 2) {
		//psi_3_calculated = calc_flow_angle(&initial_params_copy)/3.;
		psi_3_calculated = flow_angle_cmd;
		Phi_K_a = -PI+psi_3_calculated;
		Phi_K_b = PI+psi_3_calculated;
		cout << "The flow angle = " << psi_3_calculated << endl;
	}
  double Phi_K_interval = (Phi_K_b-Phi_K_a)/double(n_points);

cout<< "Computed Gaussian abscissae and corresponding weights..." << endl;

	time(&start2);

	double x_s_final, x_o_final, z_c_final, t_c_final;
	double x2_s_final, x2_o_final, t2_c_final, z2_c_final;
	double x_os_final, x_ot_final, x_oz_final, x_st_final, x_sz_final, x_tz_final;
	vector<double>* f1_ptr;
	f1_ptr = new vector<double> (order);
	vector<double>* f2_ptr;
	f2_ptr = new vector<double> (order);
	vector<double>* f3_ptr;
	f3_ptr = new vector<double> (order);
	vector<double>* f4_ptr;
	f4_ptr = new vector<double> (order);
	vector<double>* f5_ptr;
	f5_ptr = new vector<double> (order);
	vector<double>* f6_ptr;
	f6_ptr = new vector<double> (order);
	vector<double>* f7_ptr;
	f7_ptr = new vector<double> (order);
	vector<double>* f8_ptr;
	f8_ptr = new vector<double> (order);
	vector<double>* f9_ptr;
	f9_ptr = new vector<double> (order);
	vector<double>* f10_ptr;
	f10_ptr = new vector<double> (order);
	vector<double>* f11_ptr;
	f11_ptr = new vector<double> (order);
	vector<double>* f12_ptr;
	f12_ptr = new vector<double> (order);
	vector<double>* f13_ptr;
	f13_ptr = new vector<double> (order);
	vector<double>* f14_ptr;
	f14_ptr = new vector<double> (order);
	vector<double>* norm_ptr;
	norm_ptr = new vector<double> (order);

	cout << "The main output filename is " << output_filename + "_main.dat" << endl;
	cout << "The spacetime data is contained in " << output_filename + "_spacetime_data.dat" << endl;
	cout << "Initial parameters written to " << output_filename + "_init_params.dat" << endl;

	string outputfilename1 = output_filename + "_main.dat";
	string outputfilename2 = output_filename + "_spacetime_data.dat";
	string outputfilename3 = output_filename + "_init_params.dat";
	string outputfilename4 = output_filename + "_x2os_and_xos.dat";

	ofstream output (outputfilename1.data());  //save data
	ofstream output2 (outputfilename2.data());  //save data
	print_params_to_output (outputfilename3, &initial_params_copy);
	ofstream output3 ("normfile2.txt");
	//ofstream output4 (outputfilename4.data());

	output << "# Output stored in " << outputfilename1 << endl;
	output << "# For plotting exact R2_ij" << endl;
	output << "# npoints=" << n_points << ", eps3bar=" << initial_params_copy.eps_3_bar << ", v_3_bar="
	       << initial_params_copy.v_3_bar << ", psi_3_bar=" << initial_params_copy.psi_3_bar << endl;
	output << "# Integration limits are (r=0.." << limit << ", eta=" << -limit << ".." << limit << ", phi=-PI..PI, Phi_K=-PI..PI)" << endl;
	output << "# order of integration = " << order << endl;
	output << "# Phi_K-psi_3_bar \t R2_s_final \t R2_o_final \t R2_l_final \t R2_os_final \t R2_sl_final \t R2_ol_final" << endl;

	output2 << "# Phi_K-psi_3_bar \t x_s_final \t x_o_final \t z_c_final \t t_c_final \t x2_s_final \t x2_o_final \t t2_c_final \t z2_c_final \t ";
	output2 << "x_os_final \t x_ot_final \t x_oz_final \t x_st_final \t x_sz_final \t x_tz_final" << endl;

//added 09/23/2013 to study 0th harmonic difference between R^2_o and R^2_s
//double term1_avg = 0., term2_avg = 0., term3_avg = 0., term4_avg = 0.;
double term1_avga = 0., term2_avga = 0., term3_avga = 0., term4_avga = 0.;
double term1_avgb = 0., term2_avgb = 0., term3_avgb = 0., term4_avgb = 0.;

	for (int Phi_K_index = 0; Phi_K_index <= n_points; Phi_K_index++) {  //beginning of Phi_K loop

	double Phi_K_i = Phi_K_a + Phi_K_interval * Phi_K_index;
	initial_params_copy.Phi_K = Phi_K_i;

	time(&start);

	for (int l = 0; l <= order-1; l++) {  //beginning of tau loop

	//double tau_l = tau_lower + tau_interval * l;
	//if i'm using qng_1d_vec to do the tau integration below, i need to space the tau's correctly!!!
	double tau_abscissa;
	if (l < half_order) {tau_abscissa = -tau_half_length * xi[l];}
	else if (l == half_order) {tau_abscissa = 0.;}
	else {tau_abscissa = tau_half_length * xi[order-1-l];}
	double tau_l = tau_center + tau_abscissa;

	vector<Emissionfunction_data>* vec_ptr;
	vec_ptr = new vector<Emissionfunction_data> (order*order*order);

	//generate S function values over 3 spatial coordinates r \t phi, eta (in this order!!!) for given values of tau
	generate_vector(&S_function, &initial_params_copy, vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, tau_l);

	(*norm_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 1);
	//cout << "For Phi_K_l = " << Phi_K_l << " and norm = " << norm << endl;
	//output2 << Phi_K_l << "\t" << norm << endl;

	for (int i = 0; i <= order-1; i++){
	for (int j = 0; j <= order-1; j++){
	for (int k = 0; k <= order-1; k++){
			//note order of x2_s and x_s arguments!!!!!!!!!!!!!!!!  eta first, then phi!!!!!!!
			(*vec_ptr).at(index(i,j,k)).xs2 = x2_s_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).xs = x_s_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).xo2 = x2_o_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).xo = x_o_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).tc2 = t2_c_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).tc = t_c_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).zc2 = z2_c_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).zc = z_c_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).x_os = x_os_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).x_ot = x_ot_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).x_oz = x_oz_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).x_st = x_st_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).x_sz = x_sz_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
			(*vec_ptr).at(index(i,j,k)).x_tz = x_tz_integrand((*vec_ptr).at(index(i,j,k)).r, (*vec_ptr).at(index(i,j,k)).eta,
									(*vec_ptr).at(index(i,j,k)).phi, initial_params_copy.Phi_K,
									tau_l, (*vec_ptr).at(index(i,j,k)).S_val);
		}
		}
		}

	(*f1_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 2);  //x2s
	(*f2_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 3);  //xs
	(*f3_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 4);  //x2o
	(*f4_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 5);  //xo
	(*f5_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 6);  //t
	(*f6_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 7);  //t2
	(*f7_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 8);  //z
	(*f8_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 9);  //z2
	(*f9_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 10);  //xo xs
	(*f10_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 11);  //xo t
	(*f11_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 12);  //xo z
	(*f12_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 13);  //xs t
	(*f13_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 14);  //xs z
	(*f14_ptr).at(l) = qng_3d_vec (vec_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 15);  //t z

	delete vec_ptr;
	}  //end of tau loop

	time(&end);

	norm = qng_1d_vec(norm_ptr, tau_lower, tau_upper);
	output3 << norm << endl;

	x2_s_final = qng_1d_vec(f1_ptr, tau_lower, tau_upper)/norm;
	x_s_final = qng_1d_vec(f2_ptr, tau_lower, tau_upper)/norm;
	x2_o_final = qng_1d_vec(f3_ptr, tau_lower, tau_upper)/norm;
	x_o_final = qng_1d_vec(f4_ptr, tau_lower, tau_upper)/norm;
	t_c_final = qng_1d_vec(f5_ptr, tau_lower, tau_upper)/norm;
	t2_c_final = qng_1d_vec(f6_ptr, tau_lower, tau_upper)/norm;
	z_c_final = qng_1d_vec(f7_ptr, tau_lower, tau_upper)/norm;
	z2_c_final = qng_1d_vec(f8_ptr, tau_lower, tau_upper)/norm;
	x_os_final = qng_1d_vec(f9_ptr, tau_lower, tau_upper)/norm;
	x_ot_final = qng_1d_vec(f10_ptr, tau_lower, tau_upper)/norm;
	x_oz_final = qng_1d_vec(f11_ptr, tau_lower, tau_upper)/norm;
	x_st_final = qng_1d_vec(f12_ptr, tau_lower, tau_upper)/norm;
	x_sz_final = qng_1d_vec(f13_ptr, tau_lower, tau_upper)/norm;
	x_tz_final = qng_1d_vec(f14_ptr, tau_lower, tau_upper)/norm;

//added to get the 0th harmonics from R2o and R2s
if (Phi_K_index < n_points) {
//just name quadratic averages "term a" and
//linear terms squared "term b"  -->
//term1_avg += t2_c_final - t_c_final*t_c_final;
//term2_avg += x_ot_final - x_o_final*t_c_final;
//term3_avg += x2_o_final - x_o_final*x_o_final;
//term4_avg += x2_s_final - x_s_final*x_s_final;
term1_avga += t2_c_final;
term2_avga += x_ot_final;
term3_avga += x2_o_final;
term4_avga += x2_s_final;
term1_avgb += t_c_final*t_c_final;
term2_avgb += x_o_final*t_c_final;
term3_avgb += x_o_final*x_o_final;
term4_avgb += x_s_final*x_s_final;
}

	double R2_s_final = x2_s_final - x_s_final*x_s_final;

	double R2_o_final = x2_o_final - 2.*beta_perp*x_ot_final + beta_perp*beta_perp*t2_c_final
				- (x_o_final-beta_perp*t_c_final) * (x_o_final-beta_perp*t_c_final);

	double R2_l_final = z2_c_final - 2.*beta_long*x_tz_final + beta_long*beta_long*t2_c_final
				- (z_c_final-beta_long*t_c_final) * (z_c_final-beta_long*t_c_final);

	double R2_os_final = x_os_final - beta_perp*x_st_final - x_s_final*(x_o_final-beta_perp*t_c_final);

	double R2_sl_final = x_sz_final - beta_long*x_st_final - x_s_final*(z_c_final-beta_long*t_c_final);

	double R2_ol_final = x_oz_final - beta_long*x_ot_final - beta_perp*x_tz_final + beta_perp*beta_long*t2_c_final
				- (x_o_final-beta_perp*t_c_final)*(z_c_final-beta_long*t_c_final);

	//this if statement keeps the plot range of Phi_K - psi_3_bar within -PI..PI
	double angle;

	if (angle_format == 2) {
		if (Phi_K_i - psi_3_calculated < -PI) {angle = Phi_K_i - psi_3_calculated + 2.*PI;}
		else {angle = Phi_K_i - psi_3_calculated;}
		cout << "Phi_K - Psi_3 = " << angle << endl;
	}
	else if (angle_format == 1) {
		if (Phi_K_i - initial_params_copy.psi_3_bar < -PI) {angle = Phi_K_i - initial_params_copy.psi_3_bar + 2.*PI;}
		else {angle = Phi_K_i - initial_params_copy.psi_3_bar;}
		cout << "Phi_K - psi_3_bar = " << angle << endl;
	}
	else if (angle_format == 0) { angle = Phi_K_i; 
					cout << "Phi_K = " << angle << endl;}

	cout << "norm = " << norm << endl;
	cout << "R2_s_final = " << R2_s_final << endl;
	cout << "R2_o_final = " << R2_o_final << endl;
	cout << "R2_l_final = " << R2_l_final << endl;
	cout << "R2_os_final = " << R2_os_final << endl;
	cout << "R2_sl_final = " << R2_sl_final << endl;
	cout << "R2_ol_final = " << R2_ol_final << endl;
	cout << "The " << Phi_K_index << "th loop took " << difftime(end,start) << " seconds" << endl << endl;

	output << angle << "\t" << R2_s_final << "\t" << R2_o_final << "\t" << R2_l_final
	       << "\t" << R2_os_final << "\t" << R2_sl_final << "\t" << R2_ol_final << endl;
	output2 << angle << "\t" << x_s_final << "\t" << x_o_final << "\t" << z_c_final
		<< "\t" << t_c_final << "\t" << x2_s_final << "\t" << x2_o_final << "\t"
		<< t2_c_final << "\t" << z2_c_final << "\t" << x_os_final << "\t" << x_ot_final
		<< "\t" << x_oz_final << "\t" << x_st_final << "\t" << x_sz_final << "\t" << x_tz_final << endl;
//	output4 << angle << "\t" << x_os_final << "\t" << - x_s_final*x_o_final
//		<< "\t" << - beta_perp*x_st_final << "\t" << beta_perp*x_s_final*t_c_final << endl;

	}  //end of Phi_K loop

//added to get the average over Phi_K, i.e., the 0th harmonics of R2o and R2s
term1_avga /= double(n_points);
term1_avga *= beta_perp*beta_perp;
term2_avga /= double(n_points);
term2_avga *= -2.*beta_perp;
term3_avga /= double(n_points);
term4_avga /= double(n_points);
term1_avgb /= double(n_points);
term1_avgb *= beta_perp*beta_perp;
term2_avgb /= double(n_points);
term2_avgb *= -2.*beta_perp;
term3_avgb /= double(n_points);
term4_avgb /= double(n_points);
cout << setprecision(8) << "K_perp = " << K_perp << ": " << term1_avga << "\t" << term2_avga << "\t" << term3_avga << "\t" << term4_avga << "\t" << term1_avgb << "\t" << term2_avgb << "\t" << term3_avgb << "\t" << term4_avgb << "\t" << term4_avga-term4_avgb << "\t" << (term1_avga + term2_avga + term3_avga)-(term1_avgb + term2_avgb + term3_avgb) << endl;
cerr << "Finished run with K_perp = " << K_perp << endl;

	time(&end2);
	cout << "That took " << difftime(end2,start2) << " seconds" << endl;

	output.close();
	output2.close();
	output3.close();
	//output4.close();
	cout << "Data stored in " << outputfilename1 << ", " << outputfilename2 << " and " << outputfilename3 << endl;
	
  return (0);
}

/************************************************************************/

phys_params get_parameters()
{
	struct phys_params inits;

	//Define parameters for function to be integrated here
	inits.T_0 = T0_cmd;
	inits.Y_rapidity = 0.;
	inits.Rad = Rad_cmd;
	inits.eta_f = eta_f_cmd;
	inits.tau_f = tau_f_cmd;
	inits.K_perp = K_perp_cmd;
	inits.M_perp = sqrt(m_pion*m_pion + inits.K_perp*inits.K_perp);
	inits.eta_0 = eta_0_cmd;
	inits.Delta_tau = Delta_tau_cmd;
	inits.Delta_eta = Delta_eta_cmd;
	inits.beta_perp = inits.K_perp/inits.M_perp;
	//inits.beta_long = 0.;
	inits.beta_long = 0.;

	inits.eps_2_bar = eps_2_bar_cmd;
	inits.eps_3_bar = eps_3_bar_cmd;
	inits.v_2_bar = v_2_bar_cmd;
	inits.v_3_bar = v_3_bar_cmd;
	inits.psi_2_bar = psi_2_bar_cmd;
	inits.psi_3_bar = psi_3_bar_cmd;

	inits.cos_2_psi = cos(2.*inits.psi_2_bar);
	inits.sin_2_psi = sin(2.*inits.psi_2_bar);
	inits.cos_3_psi = cos(3.*inits.psi_3_bar);
	inits.sin_3_psi = sin(3.*inits.psi_3_bar);
	inits.M_1_2del_tau2 = 1./(2.*inits.Delta_tau*inits.Delta_tau);
	inits.M_1_2del_eta2 = 1./(2.*inits.Delta_eta*inits.Delta_eta);
	inits.M_1_2Rad2 = 1./(2.*inits.Rad*inits.Rad);

	return (inits);
}

/************************************************************************/

double S_function (double r, double phi, double eta, double tau, void * params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	double M_perp = params.M_perp;
	double T_0 = params.T_0;
	double eta_0 = params.eta_0;
	double Y_rapidity = params.Y_rapidity;
	double Rad = params.Rad;
	double tau_f = params.tau_f;
	double K_perp = params.K_perp;
	double Delta_tau = params.Delta_tau;
	double Delta_eta = params.Delta_eta;
	double Phi_K = params.Phi_K;

	double v_2_bar = params.v_2_bar;
	double psi_2_bar = params.psi_2_bar;
	double eps_2_bar = params.eps_2_bar;
	double v_3_bar = params.v_3_bar;
	double eps_3_bar = params.eps_3_bar;
	double psi_3_bar = params.psi_3_bar;

	//Expression for function is long and complicated,
	//so break it up into constituent terms
	double term0, term1, term2, term3, term4;
	term0 = (tau-tau_f)*(tau-tau_f) / (2.*Delta_tau*Delta_tau);
	//term0 = 0.;
	term1 = (eta-eta_0)*(eta-eta_0) / (2.*Delta_eta*Delta_eta);
	//term1 = 0.;
	term2 = (M_perp/T_0) * cosh(eta - Y_rapidity) * cosh(eta_t(r, phi, &params));
	//term2 = 0.;
	term3 = (K_perp/T_0) * (cos(phi) * cos(Phi_K) + sin(phi) * sin(Phi_K)) * sinh(eta_t(r, phi, &params));
	//term3 = 0.;
	term4 = (r*r)/(2.*Rad*Rad) * (1.+2.*eps_3_bar*( cos(3.*phi) * cos(3.*psi_3_bar) + sin(3.*phi) * sin(3.*psi_3_bar) )    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					+2.*eps_2_bar*( cos(2.*phi) * cos(2.*psi_2_bar) + sin(2.*phi) * sin(2.*psi_2_bar) ));  //changed minus signs to plus signs /*5/16/2013*/
	//term4 = 0.;

	return (exp(-term0 - term1 - term2 + term3 - term4));
	//return (1.);
}


double eta_t(double r, double phi, void * params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	double v_2_bar = params.v_2_bar;
	double psi_2_bar = params.psi_2_bar;
	double eta_f = params.eta_f;
	double Rad = params.Rad;
	double v_3_bar = params.v_3_bar;
	double psi_3_bar = params.psi_3_bar;

	return (eta_f * (r/Rad) * (1.+2.*v_3_bar*( cos(3.*phi) * cos(3.*psi_3_bar) + sin(3.*phi) * sin(3.*psi_3_bar) )
					+2.*v_2_bar*( cos(2.*phi) * cos(2.*psi_2_bar) + sin(2.*phi) * sin(2.*psi_2_bar) )));
}

double x2_s_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_s = r*sin(phi-Phi_K);

	return (r*S_val*x_s*x_s);
}

double x_s_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_s = r*sin(phi-Phi_K);

	return (r*S_val*x_s);
}

double x2_o_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_o = r*cos(phi-Phi_K);

	return (r*S_val*x_o*x_o);
}

double x_o_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_o = r*cos(phi-Phi_K);

	return (r*S_val*x_o);
}

double t2_c_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double t_c = tau*cosh(eta);

	return (r*S_val*t_c*t_c);
}

double t_c_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double t_c = tau*cosh(eta);

	return (r*S_val*t_c);
}

double z2_c_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double z_c = tau*sinh(eta);

	return (r*S_val*z_c*z_c);
}

double z_c_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double z_c = tau*sinh(eta);

	return (r*S_val*z_c);
}

double x_os_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_s = r*sin(phi-Phi_K);
	double x_o = r*cos(phi-Phi_K);

	return (r*S_val*x_o*x_s);
}

double x_ot_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_o = r*cos(phi-Phi_K);
	double t_c = tau*cosh(eta);

	return (r*S_val*x_o*t_c);
}

double x_oz_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_o = r*cos(phi-Phi_K);
	double z_c = tau*sinh(eta);

	return (r*S_val*x_o*z_c);
}

double x_st_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_s = r*sin(phi-Phi_K);
	double t_c = tau*cosh(eta);

	return (r*S_val*x_s*t_c);
}

double x_sz_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double x_s = r*sin(phi-Phi_K);
	double z_c = tau*sinh(eta);

	return (r*S_val*x_s*z_c);
}

double x_tz_integrand(double r, double eta, double phi, double Phi_K, double tau, double S_val)
{
	double z_c = tau*sinh(eta);
	double t_c = tau*cosh(eta);

	return (r*S_val*z_c*t_c);
}

double calc_flow_angle(void * params_ptr)
{
	//save original K_perp
	phys_params params = * (struct phys_params *) params_ptr;
	double K_perp_saved = params.K_perp;

	//changes from previous version:
	//	was K_T-dependent flow angle
	//	now needs to return K_T-independent flow angle
	vector<double>* K_T_vec1_ptr;
	K_T_vec1_ptr = new vector<double> (order);
	vector<double>* K_T_vec2_ptr;
	K_T_vec2_ptr = new vector<double> (order);

	double K_T_abscissa;
	//limits of numerical K_T integration to approximate
	//the actual semi-infinite interval of integration
	double K_T_lower = 0., K_T_upper = 3.5;
	double K_T_center = 0.5 * (K_T_lower + K_T_upper);
	double K_T_half_length = 0.5 * (K_T_upper - K_T_lower);
	for (int K_T_index = 0; K_T_index <= order-1; K_T_index++) {  //beginning of K_T loop
		if (K_T_index < half_order) {K_T_abscissa = -K_T_half_length * xi[K_T_index];}
		else if (K_T_index == half_order) {K_T_abscissa = 0.;}
		else {K_T_abscissa = K_T_half_length * xi[order-1-K_T_index];}
		double K_T_l = K_T_center + K_T_abscissa;
		//double K_T_l = xpts_hi[K_T_index];

	//change (*params_ptr).K_perp to K_T_l
	//cout << "K_T_l = " << K_T_l << endl;
	params.K_perp = K_T_l;

	vector<double>* vec1_ptr;
	vec1_ptr = new vector<double> (order);
	vector<double>* vec2_ptr;
	vec2_ptr = new vector<double> (order);
	//vector<double>* vec3_ptr;
	//vec3_ptr = new vector<double> (order);
	vector<Emissionfunction_data>* bigvec1_ptr;
	bigvec1_ptr = new vector<Emissionfunction_data> (order*order*order);
	vector<Emissionfunction_data>* bigvec2_ptr;
	bigvec2_ptr = new vector<Emissionfunction_data> (order*order*order);

	for (int l = 0; l <= order-1; l++) {  //beginning of Phi_K loop

	double Phi_K_abscissa;
	if (l < half_order) {Phi_K_abscissa = -PI * xi[l];}
	else if (l == half_order) {Phi_K_abscissa = 0.;}
	else {Phi_K_abscissa = PI * xi[order-1-l];}
	double Phi_K_l = 0. + Phi_K_abscissa;
	//cout << "l = " << l << ": Phi_K_abscissa = " << Phi_K_abscissa << endl;

	//generate S function values over 3 spatial coordinates r, phi, eta (in this order!!!) for given values of tau
	generate_vector(&flowangle_integ1a, &params, bigvec1_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, Phi_K_l);
	generate_vector(&flowangle_integ1b, &params, bigvec2_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, Phi_K_l);

	(*vec1_ptr).at(l) = qng_3d_vec(bigvec1_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 0);
	(*vec2_ptr).at(l) = qng_3d_vec(bigvec2_ptr, r_lower, r_upper, phi_lower, phi_upper, eta_lower, eta_upper, 0);
	//cout << "(*vec1_ptr).at(l) = " << (*vec1_ptr).at(l) << endl
	//	<< "(*vec2_ptr).at(l) = " << (*vec2_ptr).at(l) << endl;

	}  //end of Phi_K loop

	(*K_T_vec1_ptr).at(K_T_index) = K_T_l * qng_1d_vec(vec1_ptr, -PI, PI);
	(*K_T_vec2_ptr).at(K_T_index) = K_T_l * qng_1d_vec(vec2_ptr, -PI, PI);
	//cout << "(*K_T_vec1_ptr).at(K_T_index = " << K_T_index << ") = " << (*K_T_vec1_ptr).at(K_T_index) << endl;
	//cout << "(*K_T_vec2_ptr).at(K_T_index = " << K_T_index << ") = " << (*K_T_vec2_ptr).at(K_T_index) << endl << endl;

	delete vec1_ptr;
	delete vec2_ptr;
	//delete vec3_ptr;
	delete bigvec1_ptr;
	delete bigvec2_ptr;

//cerr << "Finished loop " << K_T_index+1 << " of " << order << endl;
	}  //end of K_T loop

	double integ1a = qng_1d_vec(K_T_vec1_ptr, K_T_lower, K_T_upper);
	double integ1b = qng_1d_vec(K_T_vec2_ptr, K_T_lower, K_T_upper);
	cout << setprecision(20) << "integ1a = " << integ1a << endl;
	cout << "integ1b = " << integ1b << endl;

complex<double> finalresult (integ1a, integ1b);

	params.K_perp = K_perp_saved;  //restore original value of K_perp

	return (arg(finalresult));
}

double S_no_tau(double r, double phi, double eta, double Phi_K, void * params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	//double M_perp = params.M_perp;
	//^^^this is different for every K_perp!  should be calculated directly
	//from new K_perp and m_pion!!!
	double T_0 = params.T_0;
	double eta_0 = params.eta_0;
	double Y_rapidity = params.Y_rapidity;
	double Rad = params.Rad;
	double tau_f = params.tau_f;
	double K_perp = params.K_perp;
	double M_perp = sqrt(K_perp*K_perp + m_pion*m_pion);
	double Delta_tau = params.Delta_tau;
	double Delta_eta = params.Delta_eta;
	//double Phi_K = params.Phi_K;

	double v_2_bar = params.v_2_bar;
	double psi_2_bar = params.psi_2_bar;
	double eps_2_bar = params.eps_2_bar;
	double v_3_bar = params.v_3_bar;
	double eps_3_bar = params.eps_3_bar;
	double psi_3_bar = params.psi_3_bar;

	//Expression for function is long and complicated,
	//so break it up into constituent terms
	double term0, term1, term2, term3, term4;
	//term0 = (tau-tau_f)*(tau-tau_f) / (2.*Delta_tau*Delta_tau);
	//term0 = 0.;
	term1 = (eta-eta_0)*(eta-eta_0) / (2.*Delta_eta*Delta_eta);
	//term1 = 0.;
	term2 = (M_perp/T_0) * cosh(eta - Y_rapidity) * cosh(eta_t(r, phi, &params));
	//term2 = 0.;
	term3 = (K_perp/T_0) * (cos(phi) * cos(Phi_K) + sin(phi) * sin(Phi_K)) * sinh(eta_t(r, phi, &params));
	//term3 = 0.;
	term4 = (r*r)/(2.*Rad*Rad) * (1.+2.*eps_3_bar*( cos(3.*phi) * cos(3.*psi_3_bar) + sin(3.*phi) * sin(3.*psi_3_bar) )    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					+2.*eps_2_bar*( cos(2.*phi) * cos(2.*psi_2_bar) + sin(2.*phi) * sin(2.*psi_2_bar) ));  //changed minus signs to plus signs /*5/16/2013*/
	//term4 = 0.;

	return (exp(- term1 - term2 + term3 - term4));
	//return (1.);
	//return (cos(3.*Phi_K));
}
double flowangle_integ1a(double r, double phi, double eta, double Phi_K, void * params_ptr)
{
	return (r*cos(3.*Phi_K)*S_no_tau(r, phi, eta, Phi_K, params_ptr));
}

double flowangle_integ1b(double r, double phi, double eta, double Phi_K, void * params_ptr)
{
	return (r*sin(3.*Phi_K)*S_no_tau(r, phi, eta, Phi_K, params_ptr));
}

double flowangle_integ2(double r, double phi, double eta, double Phi_K, void * params_ptr)
{
	return (r*S_no_tau(r, phi, eta, Phi_K, params_ptr));
}

double test_function (double x, double y, double z, double t, void * params_ptr)
{
	my_params params = * (struct my_params *) params_ptr;

	double alpha = params.alpha;
	double beta = params.beta;

	//return ( exp(-alpha*(0*x*x+0*y*y+0*z*z+t*t)+beta) );
	//return (x*x+y*y+z*z);
	//return (1.);
	return (exp(cos(x*x*x*x+y*y*y+z*z+t)));
}

int index(int i, int j, int k)
{
	return (order*order*i+order*j+k);
}

int index(int i, int j, int k, int l)
{
	return (order*order*order*i+order*order*j+order*k+l);
}

bool is_multiple(int x, int y)
{
	//x is number which might be a multiple of y
	return (double(x/y) == double(x)/double(y));
}

//End of file
