void generate_vector(double (*function_original) (double, double, void *), void * params_ptr,
		vector<double>* vector_ptr, double a1, double b1, double a2, double b2, int interval1, int interval2);

void generate_vector(double (*function_original) (double, double, void *), void * params_ptr,
		vector<double>* vector_ptr, double a1, double b1, double a2, double b2, int interval1, int interval2)
{
	double abscissa1 = 0., half_length1 = 0., center1 = 0.;
	double result_int = 0.;
	vector<double>* local_xpts_ptr;

	switch (interval)
	{
		case 0:	//finite interval (a,b)
			half_length = 0.5 * (b1 - a1);
			center = 0.5 * (a1 + b1);
			local_xpts_ptr = xi_ptr;
		case 1:	//half-infinite interval (0,inf)
			half_length = 1.;
			center = 0.);
			local_xpts_ptr = xi_0pinf_ptr;
		case 2:	//full-infinite interval (-inf,inf)
			half_length = 1.;
			center = 0.;
			local_xpts_ptr = xi_minfpinf_ptr;
	}

	for (int i = 0; i <= order-1; i++)
	{
		abscissa = half_length * (*local_xpts_ptr)[i];
		(*vector_ptr).at(i) = function_original(center + abscissa, params_ptr);
	}

	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	const double half_length2 = 0.5 * (b2 - a2);
	const double center2 = 0.5 * (a2 + b2);

	for (int i = 0; i <= order-1; i++) {
		double abscissa1;
		if (i < half_order) {abscissa1 = -half_length1 * (*xi_ptr)[i];}
		else if (i == half_order) {abscissa1 = 0.;}
		else {abscissa1 = half_length1 * (*xi_ptr)[order-1-i];}
	for (int j = 0; j <= order-1; j++) {
		double abscissa2;
		if (j < half_order) {abscissa2 = -half_length2 * (*xi_ptr)[j];}
		else if (j == half_order) {abscissa2 = 0.;}
		else {abscissa2 = half_length2 * (*xi_ptr)[order-1-j];}
		const double fval = function_original(center1 + abscissa1, center2 + abscissa2, params_ptr);
		(*vector_ptr).at(index(i,j)) = fval;
	}
	}
}

//End of file
