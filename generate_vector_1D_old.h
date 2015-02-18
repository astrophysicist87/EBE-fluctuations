void generate_vector(double (*function_original) (double, void *), void * params_ptr,
		vector<double>* vector_ptr, double a1, double b1);
void generate_vector(double (*function_original) (double, position *),
		position * position_ptr, vector<double>* vector_ptr, double a1, double b1);

void generate_vector(double (*function_original) (double, void *), void * params_ptr,
		vector<double>* vector_ptr, double a1, double b1)
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	double result_int = 0.;

	for (int i = 0; i <= order-1; i++) {
		double abscissa1;
		if (i < half_order) {abscissa1 = -half_length1 * xi[i];}
		else if (i == half_order) {abscissa1 = 0.;}
		else {abscissa1 = half_length1 * xi[order-1-i];}
		const double fval = function_original(center1 + abscissa1, params_ptr);
		(*vector_ptr).at(i) = fval;
	}
}

void generate_vector(double (*function_original) (double, position *),
		position * position_ptr, vector<double>* vector_ptr, double a1, double b1)
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	double result_int = 0.;

	for (int i = 0; i <= order-1; i++) {
		double abscissa1;
		if (i < half_order) {abscissa1 = -half_length1 * xi[i];}
		else if (i == half_order) {abscissa1 = 0.;}
		else {abscissa1 = half_length1 * xi[order-1-i];}
		(*vector_ptr).at(i) = function_original(center1 + abscissa1, position_ptr);
	}
}

//End of file
