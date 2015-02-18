double qng_2d_vec (vector<double>* vector_ptr, double a1, double b1, double a2, double b2);

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

		result_int_x += fval * wi[i] * wi[j] * factor;
	}
		result_int += result_int_x;
	}
	result_int *= half_length1 * half_length2;

	return (result_int);
}

//End of file
