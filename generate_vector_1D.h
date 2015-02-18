void generate_vector(double (*function_original) (double, void *), void * params_ptr,
		vector<double>* vector_ptr, double a1, double b1, int interval);

void generate_vector(double (*function_original) (double, void *), void * params_ptr,
		vector<double>* vector_ptr, double a1, double b1, int interval)
{
	double abscissa = 0., half_length = 0., center = 0.;
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
}

//End of file
