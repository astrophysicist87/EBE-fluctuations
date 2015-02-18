#ifndef PDFS_H
#define PDFS_H

using namespace std;

double pdf(vector<double>* fluct_vec_ptr, void * params_ptr);
double uniformPDF(int flucti, vector<double>* fluct_vec_ptr, void * params_ptr);
double gaussianPDF(int flucti, vector<double>* fluct_vec_ptr, void * params_ptr);
double besselgaussianPDF(int flucti, vector<double>* fluct_vec_ptr, void * params_ptr);
double quadraticPDF(int i, vector<double>* fluct_vec_ptr, void * params_ptr);
double gaussianaltPDF(int i, vector<double>* fluct_vec_ptr, void * params_ptr);
double theta_function(double x);

double pdf(vector<double>* fluct_vec_ptr, void * params_ptr)
{
//allow for different distributions for different fluctuating parameters
//if only one pdf_switch specified, use same pdf for all parameters
//otherwise, specify each individually
double result = 1.;
for (int i = 0; i < (*fluct_vec_ptr).size(); i++)
{
//i indexes the parameter currently being varied, unless all pdfs same
int pdf_switch = (pdfswitchvec.size() == 1) ? pdfswitchvec[0] : pdfswitchvec[i];
	switch(pdf_switch)
	{
		case 0:		//uniform
			result *= uniformPDF(i, fluct_vec_ptr, params_ptr);
			break;
		case 1:		//Gaussian
			result *= gaussianPDF(i, fluct_vec_ptr, params_ptr);
			break;
		case 2:		//Bessel-Gaussian
			result *= besselgaussianPDF(i, fluct_vec_ptr, params_ptr);
			break;
		case 3:		//quadratic
			result *= quadraticPDF(i, fluct_vec_ptr, params_ptr);
			break;
		case 4:		//Gaussian alternative which vanishes at center +/- sigma
			result *= gaussianaltPDF(i, fluct_vec_ptr, params_ptr);
			break;
	}
}
	return (result);
}

double uniformPDF(int i, vector<double>* fluct_vec_ptr, void * params_ptr)
{
	pdf_params params = * (struct pdf_params *) params_ptr;

	const double width = (*(params.sigma_vec_ptr))[i];
	const double center = (*(params.mean_vec_ptr))[i];
	const double result = ( abs((*fluct_vec_ptr)[i]-center) > width ) ? 0. : 1./(2.*width);

	return (result);
}

double gaussianPDF(int i, vector<double>* fluct_vec_ptr, void * params_ptr)
{
	pdf_params params = * (struct pdf_params *) params_ptr;

	const double width = (*(params.sigma_vec_ptr))[i];
	const double center = (*(params.mean_vec_ptr))[i];
	const double x = (*fluct_vec_ptr)[i];
	const double result = (1./(sqrt(2.*PI)*width)) * exp(-(x-center)*(x-center)/(2.*width*width));

	return (result);
}

double besselgaussianPDF(int i, vector<double>* fluct_vec_ptr, void * params_ptr)
{
	pdf_params params = * (struct pdf_params *) params_ptr;

	const double width = (*(params.sigma_vec_ptr))[i];
	const double center = (*(params.mean_vec_ptr))[i];
	const double x = (*fluct_vec_ptr)[i];
	double factor = (x*center/(width*width) > 700.) ? largest_double : gsl_sf_bessel_I0(x*center/(width*width));
	//const double result = (x/(width*width)) *  gsl_sf_bessel_I0(x*center/(width*width)) * exp(-(x*x + center*center)/(2.*width*width));
	const double result = (x/(width*width)) * factor * exp(-(x*x + center*center)/(2.*width*width));

	return (result);
}

double quadraticPDF(int i, vector<double>* fluct_vec_ptr, void * params_ptr)
{
	pdf_params params = * (struct pdf_params *) params_ptr;

	const double width = (*(params.sigma_vec_ptr))[i];
	const double center = (*(params.mean_vec_ptr))[i];
	const double x = (*fluct_vec_ptr)[i];
	const double result = theta_function(x-center+width) * theta_function(center+width-x) * (-3./(4.*width*width*width)) * (x-center+width) * (x-center-width);

	return (result);
}

double gaussianaltPDF(int i, vector<double>* fluct_vec_ptr, void * params_ptr)
{
	pdf_params params = * (struct pdf_params *) params_ptr;

	const double width = (*(params.sigma_vec_ptr))[i];
	const double center = (*(params.mean_vec_ptr))[i];
	const double x = (*fluct_vec_ptr)[i];
	const double result = theta_function(x-center+width) * theta_function(center+width-x) * (exp(0.5)/(2.*width*width*width))
				* (x-center+width) * (center+width-x) * exp(-(x-center)*(x-center)/(2.*width*width));

	return (result);
}

double theta_function(double x)
{
	return ((x>=0.) ? 1. : 0.);
}

//End of file

#endif
