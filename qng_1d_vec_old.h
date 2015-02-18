double qng_1d_vec (vector<double>* vector_ptr, double a1, double b1);

double
qng_1d_vec (vector<double>* vector_ptr, double a1, double b1)  	//lower and upper limits for x integration
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	double result_int = 0.;

	for (int i = 0; i <= half_order; i++) {
		const double fval = (*vector_ptr).at(order-1-i) + (*vector_ptr).at(i);
			//this part sums up all function values at different abscissae which receive the same weights
		double factor = 1.;
		if (i == half_order) {factor /= 2.;}
		result_int += fval * wi[i] * factor;
	}
	result_int *= half_length1;

	return (result_int);
}
