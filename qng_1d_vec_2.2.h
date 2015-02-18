#ifndef QNG_1D_VEC_2_2_H
#define QNG_1D_VEC_2_2_H

using namespace std;

double qng_1d_vec (vector<double>* vector_ptr, double a1, double b1, int interval1);

double
qng_1d_vec (vector<double>* vector_ptr, double a1, double b1, int interval1)  	//lower and upper limits for x integration
{

	double abscissa1 = 0., half_length1 = 0., center1 = 0.;
	double result_int = 0.;
	vector<double>* local_xpts_ptr1;
	vector<double>* local_wts_ptr1;

	switch (interval1)
	{
		case 0:	//finite interval (a,b)
			half_length1 = 0.5 * (b1 - a1);
			center1 = 0.5 * (a1 + b1);
			local_xpts_ptr1 = xi_ptr;
			local_wts_ptr1 = wi_ptr;
			break;
		case 1:	//half-infinite interval (0,inf)
			half_length1 = 1.;
			center1 = 0.;
			local_xpts_ptr1 = xi_0pinf_ptr;
			local_wts_ptr1 = wi_0pinf_ptr;
			break;
		case 2:	//full-infinite interval (-inf,inf)
			half_length1 = 1.;
			center1 = 0.;
			local_xpts_ptr1 = xi_minfpinf_ptr;
			local_wts_ptr1 = wi_minfpinf_ptr;
			break;
	}

	for (int i = 0; i <= order-1; i++)
	{
		result_int += (*local_wts_ptr1)[i] * (*vector_ptr)[i];
	}
	result_int *= half_length1;

	return (result_int);
}

/*double qng_1d_vec_alt (vector<double>* vector_ptr, double a1, double b1, int interval1);

double
qng_1d_vec_alt (vector<double>* vector_ptr, double a1, double b1, int interval1)  	//lower and upper limits for x integration
{

	double abscissa1 = 0., half_length1 = 0., center1 = 0.;
	double result_int = 0.;
	vector<double>* local_xpts_ptr1;
	vector<double>* local_wts_ptr1;

	switch (interval1)
	{
		case 0:	//finite interval (a,b)
			half_length1 = 0.5 * (b1 - a1);
			center1 = 0.5 * (a1 + b1);
			local_xpts_ptr1 = xi_ptr;
			local_wts_ptr1 = wi_ptr;
			break;
		case 1:	//half-infinite interval (0,inf)
			half_length1 = 1.;
			center1 = 0.;
			local_xpts_ptr1 = xiq_0pinf_ptr;
			local_wts_ptr1 = wiq_0pinf_ptr;
			break;
		case 2:	//full-infinite interval (-inf,inf)
			half_length1 = 1.;
			center1 = 0.;
			local_xpts_ptr1 = xiq_minfpinf_ptr;
			local_wts_ptr1 = wiq_minfpinf_ptr;
			break;
	}

	for (int i = 0; i <= order-1; i++)
	{
		result_int += (*local_wts_ptr1)[i] * (*vector_ptr)[i];
	}
	result_int *= half_length1;

	return (result_int);
}*/

#endif
