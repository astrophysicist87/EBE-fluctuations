#ifndef INTEGRATION_ROUTINES_2_2_H
#define INTEGRATION_ROUTINES_2_2_H

using namespace std;

double integrate_1D(double (*function_original) (double, void *), void * params_ptr,
		double a1, double b1, int interval);
double integrate_2D(double (*function_original) (double, double, void *), void * params_ptr,
		double a1, double b1, double a2, double b2, int interval1, int interval2);
double integrate(double(*function_original) (vector<double> *, void *), void * params_ptr,
		vector<double> * lower_limits_ptr, vector<double> * upper_limits_ptr, vector<int> * intervals);

double integrate_1D(double (*function_original) (double, void *), void * params_ptr,
		double a1, double b1, int interval)
{
	double abscissa = 0., half_length = 0., center = 0.;
	double result_int = 0.;
	vector<double>* local_xpts_ptr;
	vector<double>* local_wts_ptr;

	switch (interval)
	{
		case 0:	//finite interval (a,b)
			half_length = 0.5 * (b1 - a1);
			center = 0.5 * (a1 + b1);
			local_xpts_ptr = xi_ptr;
			local_wts_ptr = wi_ptr;
			break;
		case 1:	//half-infinite interval (0,inf)
			half_length = 1.;
			center = 0.;
			local_xpts_ptr = xi_0pinf_ptr;
			local_wts_ptr = wi_0pinf_ptr;
			break;
		case 2:	//full-infinite interval (-inf,inf)
			half_length = 1.;
			center = 0.;
			local_xpts_ptr = xi_minfpinf_ptr;
			local_wts_ptr = wi_minfpinf_ptr;
			break;
	}

	for (int i = 0; i <= order-1; i++)
	{
		abscissa = half_length * (*local_xpts_ptr)[i];
		result_int += half_length * (*local_wts_ptr)[i] * function_original(center + abscissa, params_ptr);
	}

	return (result_int);
}

double integrate_2D(double (*function_original) (double, double, void *), void * params_ptr,
		double a1, double b1, double a2, double b2, int interval1, int interval2)
{
	double abscissa1 = 0., half_length1 = 0., center1 = 0.;
	double abscissa2 = 0., half_length2 = 0., center2 = 0.;
	double result_int = 0.;
	vector<double>* local_xpts_ptr1;
	vector<double>* local_wts_ptr1;
	vector<double>* local_xpts_ptr2;
	vector<double>* local_wts_ptr2;

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

	switch (interval2)
	{
		case 0:	//finite interval (a,b)
			half_length2 = 0.5 * (b2 - a2);
			center2 = 0.5 * (a2 + b2);
			local_xpts_ptr2 = xi_ptr;
			local_wts_ptr2 = wi_ptr;
			break;
		case 1:	//half-infinite interval (0,inf)
			half_length2 = 1.;
			center2 = 0.;
			local_xpts_ptr2 = xi_0pinf_ptr;
			local_wts_ptr2 = wi_0pinf_ptr;
			break;
		case 2:	//full-infinite interval (-inf,inf)
			half_length2 = 1.;
			center2 = 0.;
			local_xpts_ptr2 = xi_minfpinf_ptr;
			local_wts_ptr2 = wi_minfpinf_ptr;
			break;
	}

	for (int i = 0; i <= order-1; i++)
	{
		abscissa1 = half_length1 * (*local_xpts_ptr1)[i];
		for (int j = 0; j <= order-1; j++)
		{
			abscissa2 = half_length2 * (*local_xpts_ptr2)[j];
			result_int += (*local_wts_ptr1)[i] * (*local_wts_ptr2)[j] * function_original(center1 + abscissa1, center2 + abscissa2, params_ptr);
		}
	}

	return (result_int);
}

double integrate(double(*function_original) (vector<double> *, void *), void * params_ptr,
		vector<double> * lower_limits_ptr, vector<double> * upper_limits_ptr, vector<int> * intervals)
{
	double result = 0.;
	const int dimint = (*lower_limits_ptr).size();		//the dimension of the integral we need to do
	vector< vector<double> * > all_local_xpts_ptrs (dimint);
	vector< vector<double> * > all_local_wts_ptrs (dimint);

	vector<double> halflengths (dimint), centers (dimint);	//vectors to hold the halflengths and center of each respective integration interval

	for (int i = 0; i < dimint; i++)
	{
//cout << "(*intervals)[" << i << "] = " << (*intervals)[i] << endl;
		switch ((*intervals)[i])
		{
		case 0:	//finite interval (a,b)
			halflengths[i] = 0.5 * ((*upper_limits_ptr)[i] - (*lower_limits_ptr)[i]);
			centers[i] = 0.5 * ((*upper_limits_ptr)[i] + (*lower_limits_ptr)[i]);
//cout << "went to case 0!" << endl;
			all_local_xpts_ptrs[i] = xi_ptr;
			all_local_wts_ptrs[i] = wi_ptr;
			break;
		case 1:	//half-infinite interval (0,inf)
			halflengths[i] = 1.;
			centers[i] = 0.;
//cout << "went to case 1!" << endl;
			all_local_xpts_ptrs[i] = xi_0pinf_ptr;
			all_local_wts_ptrs[i] = wi_0pinf_ptr;
			break;
		case 2:	//full-infinite interval (-inf,inf)
			halflengths[i] = 1.;
			centers[i] = 0.;
//cout << "went to case 2!" << endl;
			all_local_xpts_ptrs[i] = xi_minfpinf_ptr;
			all_local_wts_ptrs[i] = wi_minfpinf_ptr;
			break;
		}
	}
	for (int i = 0; i < pow(order, dimint); i++)
	{
		double wt_factor = 1.;
		vector<int> indices = convert_base(i, order, dimint);	//treats all points for multi-dimensional integration as one giant vector
									//indices vector now holds index value of loop variable i in order**dimint grid
									//total number of points is order**dimint
		//if (i == 10*order) cout << indices[0] << "\t" << indices[1] << endl;
		vector<double> xlocs (dimint);
		for (int j = 0; j < dimint; j++)
		{
			//cout << "INSIDE integrate(): centers["<< j << "] = " << centers[j] << " and halflengths["<< j << "] = " << halflengths[j] << endl;
			xlocs[j] = centers[j] + halflengths[j] * (*(all_local_xpts_ptrs[j]))[indices[j]];
			wt_factor *= halflengths[j] * (*(all_local_wts_ptrs[j]))[indices[j]];
		}
		result += wt_factor * function_original(&xlocs, params_ptr);
	}

	return (result);
}

//End of file

#endif
