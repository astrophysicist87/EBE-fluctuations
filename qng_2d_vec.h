double qng_2d_vec (vector<double>* vector_ptr, double a1, double b1, double a2, double b2);
double qng_2d_vec (vector<double>* vector_ptr,
			double a1, double b1, double a2, double b2, int interval1, int interval2);
double q_int2d_vec_alt(vector<double>* vector_ptr,
			double a1, double b1, double a2, double b2, int interval1, int interval2);

double qng_2d_vec (vector<double>* vector_ptr,
        double a1, double b1,  	//lower and upper limits for x integration
	double a2, double b2  	//lower and upper limits for y integration
	)
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	const double half_length2 = 0.5 * (b2 - a2);
	const double center2 = 0.5 * (a2 + b2);
	double result_int = 0.;

	for (int i = 0; i <= half_order; i++) {
		double result_int_x = 0.;
	for (int j = 0; j <= half_order; j++) {
		int i_c = order-1-i;
		int j_c = order-1-j;
		const double fval
			= (*vector_ptr)[index(i_c,j_c)] + (*vector_ptr)[index(i,j_c)] + (*vector_ptr)[index(i_c,j)] + (*vector_ptr)[index(i,j)];
		
		double factor = 1.;
		if (i == half_order) {factor /= 2.;}
		if (j == half_order) {factor /= 2.;}

		result_int_x += fval * (*wi_ptr)[i] * (*wi_ptr)[j] * factor;
	}
		result_int += result_int_x;
	}
	result_int *= half_length1 * half_length2;

	return (result_int);
}

double qng_2d_vec(vector<double>* vector_ptr,
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
			//result_int += (*local_wts_ptr1)[i] * (*local_wts_ptr2)[j] * function_original(center1 + abscissa1, center2 + abscissa2, params_ptr);
			result_int += (*local_wts_ptr1)[i] * (*local_wts_ptr2)[j] * (*vector_ptr)[index(i,j)];
		}
	}
	result_int *= half_length1 * half_length2;

	return (result_int);
}

double q_int2d_alt(vector<double>* vector_ptr)
{
	double result_int = 0., result_int_x = 0.;

	for (int i = 0; i <= nq_points-1; i++)
	{
		for (int j = 0; j <= nq_points-1; j++)
		{
			double factor1 = 1.;
			double factor2 = 1.;
			if (i == 0 || i == nq_points-1) factor1 = 0.5;
			if (j == 0 || j == nq_points-1) factor2 = 0.5;
			result_int_x += factor1*factor2*(*vector_ptr).at(nq_points*i+j);
		}
result_int += result_int_x;
	}
	result_int *= q_interval * q_interval;

	return (result_int);
}

//End of file
