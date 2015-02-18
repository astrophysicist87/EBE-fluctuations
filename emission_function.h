#ifndef EMISSION_FUNCTION_H
#define EMISSION_FUNCTION_H

using namespace std;

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
	double factor1 = (abs(eps_n_bar) > 0.5) ? 0. : sqrt(1-4.*eps_n_bar*eps_n_bar);
	//double factor1 = 1.;
	double factor2 = (abs(eps_n_bar) > 0.5) ? 0. : exp(-(r*r)/(2.*R*R)*(1.+2.*eps_n_bar*cos(harmonic_n*(phi-psi_n_bar))));
//if (abs(eps_n_bar) > 0.5) cout << "HIT UPPER LIMIT!!!" << endl;

	const double result = factor1 * ( N / (2.*PI*R*R) ) * factor2;
	//const double result = sqrt(1-4.*eps_n_bar*eps_n_bar) * ( N / (2.*PI*R*R) ) * exp(-(r*r)/(2.*R*R)*(1.+2.*eps_n_bar*cos(harmonic_n*(phi-psi_n_bar))));
//if (my_flag) cout << exp(-(r*r)/(2.*R*R)*(1.+2.*eps_n_bar*cos(harmonic_n*(phi-psi_n_bar)))) << endl;

	return (result);
	//return (1.);
	//return (exp(-r*r/50.));
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
}
*/

#endif
