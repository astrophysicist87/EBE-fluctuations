void generate_vector(double (*function_original) (double, double, void *), void * params_ptr,
		vector<double>* vector_ptr, double a1, double b1, double a2, double b2);
void generate_vector(double (*function_original) (double, double, position *),
		position * position_ptr, vector<double>* vector_ptr, double a1, double b1);

void generate_vector(double (*function_original) (double, double, void *), void * params_ptr,
		vector<double>* vector_ptr, double a1, double b1, double a2, double b2)
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	const double half_length2 = 0.5 * (b2 - a2);
	const double center2 = 0.5 * (a2 + b2);

	for (int i = 0; i <= order-1; i++) {
		double abscissa1;
		if (i < half_order) {abscissa1 = -half_length1 * xi[i];}
		else if (i == half_order) {abscissa1 = 0.;}
		else {abscissa1 = half_length1 * xi[order-1-i];}
	for (int j = 0; j <= order-1; j++) {
		double abscissa2;
		if (j < half_order) {abscissa2 = -half_length2 * xi[j];}
		else if (j == half_order) {abscissa2 = 0.;}
		else {abscissa2 = half_length2 * xi[order-1-j];}
		const double fval = function_original(center1 + abscissa1, center2 + abscissa2, params_ptr);
		(*vector_ptr).at(index(i,j)) = fval;
	}
	}
}

void generate_vector(double (*function_original) (double, double, position *),
		position * position_ptr, vector<double>* vector_ptr, double a1, double b1, double a2, double b2)
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	const double half_length2 = 0.5 * (b2 - a2);
	const double center2 = 0.5 * (a2 + b2);

	for (int i = 0; i <= order-1; i++) {
		double abscissa1;
		if (i < half_order) {abscissa1 = -half_length1 * xi[i];}
		else if (i == half_order) {abscissa1 = 0.;}
		else {abscissa1 = half_length1 * xi[order-1-i];}
	for (int j = 0; j <= order-1; j++) {
		double abscissa2;
		if (j < half_order) {abscissa2 = -half_length2 * xi[j];}
		else if (j == half_order) {abscissa2 = 0.;}
		else {abscissa2 = half_length2 * xi[order-1-j];}
		const double fval = function_original(center1 + abscissa1, center2 + abscissa2, position_ptr);
		(*vector_ptr).at(index(i,j)) = fval;
	}
	}
}

//End of file
