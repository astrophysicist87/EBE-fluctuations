double qng_3d_vec (vector<Emissionfunction_data>* vector_ptr,
		double a1, double b1, double a2, double b2, double a3, double b3, int int_option);

double
qng_3d_vec (vector<Emissionfunction_data>* vector_ptr,
        double a1, double b1,  	//lower and upper limits for x integration
	double a2, double b2,  	//lower and upper limits for y integration
	double a3, double b3,  	//lower and upper limits for z integration
	int int_option)
{
	const double half_length1 = 0.5 * (b1 - a1);
	const double center1 = 0.5 * (a1 + b1);
	const double half_length2 = 0.5 * (b2 - a2);
	const double center2 = 0.5 * (a2 + b2);
	const double half_length3 = 0.5 * (b3 - a3);
	const double center3 = 0.5 * (a3 + b3);
	double result_int = 0.;

	for (int i = 0; i <= half_order; i++) {
		double result_int_x = 0.;
	for (int j = 0; j <= half_order; j++) {
		double result_int_y = 0.;
	for (int k = 0; k <= half_order; k++) {
		int i_c = order-1-i;
		int j_c = order-1-j;
		int k_c = order-1-k;
		double fval;

		switch(int_option)
		{
			case 0:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].S_val + (*vector_ptr)[index(i,j_c,k_c)].S_val
					+ (*vector_ptr)[index(i_c,j,k_c)].S_val + (*vector_ptr)[index(i,j,k_c)].S_val
					+ (*vector_ptr)[index(i_c,j_c,k)].S_val + (*vector_ptr)[index(i,j_c,k)].S_val
					+ (*vector_ptr)[index(i_c,j,k)].S_val + (*vector_ptr)[index(i,j,k)].S_val;
				break;
			case 1:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].norm + (*vector_ptr)[index(i,j_c,k_c)].norm
					+ (*vector_ptr)[index(i_c,j,k_c)].norm + (*vector_ptr)[index(i,j,k_c)].norm
					+ (*vector_ptr)[index(i_c,j_c,k)].norm + (*vector_ptr)[index(i,j_c,k)].norm
					+ (*vector_ptr)[index(i_c,j,k)].norm + (*vector_ptr)[index(i,j,k)].norm;
				break;
			case 2:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].xs2 + (*vector_ptr)[index(i,j_c,k_c)].xs2
					+ (*vector_ptr)[index(i_c,j,k_c)].xs2 + (*vector_ptr)[index(i,j,k_c)].xs2
					+ (*vector_ptr)[index(i_c,j_c,k)].xs2 + (*vector_ptr)[index(i,j_c,k)].xs2
					+ (*vector_ptr)[index(i_c,j,k)].xs2 + (*vector_ptr)[index(i,j,k)].xs2;
				break;
			case 3:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].xs + (*vector_ptr)[index(i,j_c,k_c)].xs
					+ (*vector_ptr)[index(i_c,j,k_c)].xs + (*vector_ptr)[index(i,j,k_c)].xs
					+ (*vector_ptr)[index(i_c,j_c,k)].xs + (*vector_ptr)[index(i,j_c,k)].xs
					+ (*vector_ptr)[index(i_c,j,k)].xs + (*vector_ptr)[index(i,j,k)].xs;
				break;
			case 4:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].xo2 + (*vector_ptr)[index(i,j_c,k_c)].xo2
					+ (*vector_ptr)[index(i_c,j,k_c)].xo2 + (*vector_ptr)[index(i,j,k_c)].xo2
					+ (*vector_ptr)[index(i_c,j_c,k)].xo2 + (*vector_ptr)[index(i,j_c,k)].xo2
					+ (*vector_ptr)[index(i_c,j,k)].xo2 + (*vector_ptr)[index(i,j,k)].xo2;
				break;
			case 5:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].xo + (*vector_ptr)[index(i,j_c,k_c)].xo
					+ (*vector_ptr)[index(i_c,j,k_c)].xo + (*vector_ptr)[index(i,j,k_c)].xo
					+ (*vector_ptr)[index(i_c,j_c,k)].xo + (*vector_ptr)[index(i,j_c,k)].xo
					+ (*vector_ptr)[index(i_c,j,k)].xo + (*vector_ptr)[index(i,j,k)].xo;
				break;
			case 6:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].tc + (*vector_ptr)[index(i,j_c,k_c)].tc
					+ (*vector_ptr)[index(i_c,j,k_c)].tc + (*vector_ptr)[index(i,j,k_c)].tc
					+ (*vector_ptr)[index(i_c,j_c,k)].tc + (*vector_ptr)[index(i,j_c,k)].tc
					+ (*vector_ptr)[index(i_c,j,k)].tc + (*vector_ptr)[index(i,j,k)].tc;
				break;
			case 7:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].tc2 + (*vector_ptr)[index(i,j_c,k_c)].tc2
					+ (*vector_ptr)[index(i_c,j,k_c)].tc2 + (*vector_ptr)[index(i,j,k_c)].tc2
					+ (*vector_ptr)[index(i_c,j_c,k)].tc2 + (*vector_ptr)[index(i,j_c,k)].tc2
					+ (*vector_ptr)[index(i_c,j,k)].tc2 + (*vector_ptr)[index(i,j,k)].tc2;
				break;
			case 8:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].zc + (*vector_ptr)[index(i,j_c,k_c)].zc
					+ (*vector_ptr)[index(i_c,j,k_c)].zc + (*vector_ptr)[index(i,j,k_c)].zc
					+ (*vector_ptr)[index(i_c,j_c,k)].zc + (*vector_ptr)[index(i,j_c,k)].zc
					+ (*vector_ptr)[index(i_c,j,k)].zc + (*vector_ptr)[index(i,j,k)].zc;
				break;
			case 9:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].zc2 + (*vector_ptr)[index(i,j_c,k_c)].zc2
					+ (*vector_ptr)[index(i_c,j,k_c)].zc2 + (*vector_ptr)[index(i,j,k_c)].zc2
					+ (*vector_ptr)[index(i_c,j_c,k)].zc2 + (*vector_ptr)[index(i,j_c,k)].zc2
					+ (*vector_ptr)[index(i_c,j,k)].zc2 + (*vector_ptr)[index(i,j,k)].zc2;
				break;
			case 10:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].x_os + (*vector_ptr)[index(i,j_c,k_c)].x_os
					+ (*vector_ptr)[index(i_c,j,k_c)].x_os + (*vector_ptr)[index(i,j,k_c)].x_os
					+ (*vector_ptr)[index(i_c,j_c,k)].x_os + (*vector_ptr)[index(i,j_c,k)].x_os
					+ (*vector_ptr)[index(i_c,j,k)].x_os + (*vector_ptr)[index(i,j,k)].x_os;
				break;
			case 11:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].x_ot + (*vector_ptr)[index(i,j_c,k_c)].x_ot
					+ (*vector_ptr)[index(i_c,j,k_c)].x_ot + (*vector_ptr)[index(i,j,k_c)].x_ot
					+ (*vector_ptr)[index(i_c,j_c,k)].x_ot + (*vector_ptr)[index(i,j_c,k)].x_ot
					+ (*vector_ptr)[index(i_c,j,k)].x_ot + (*vector_ptr)[index(i,j,k)].x_ot;
				break;
			case 12:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].x_oz + (*vector_ptr)[index(i,j_c,k_c)].x_oz
					+ (*vector_ptr)[index(i_c,j,k_c)].x_oz + (*vector_ptr)[index(i,j,k_c)].x_oz
					+ (*vector_ptr)[index(i_c,j_c,k)].x_oz + (*vector_ptr)[index(i,j_c,k)].x_oz
					+ (*vector_ptr)[index(i_c,j,k)].x_oz + (*vector_ptr)[index(i,j,k)].x_oz;
				break;
			case 13:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].x_st + (*vector_ptr)[index(i,j_c,k_c)].x_st
					+ (*vector_ptr)[index(i_c,j,k_c)].x_st + (*vector_ptr)[index(i,j,k_c)].x_st
					+ (*vector_ptr)[index(i_c,j_c,k)].x_st + (*vector_ptr)[index(i,j_c,k)].x_st
					+ (*vector_ptr)[index(i_c,j,k)].x_st + (*vector_ptr)[index(i,j,k)].x_st;
				break;
			case 14:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].x_sz + (*vector_ptr)[index(i,j_c,k_c)].x_sz
					+ (*vector_ptr)[index(i_c,j,k_c)].x_sz + (*vector_ptr)[index(i,j,k_c)].x_sz
					+ (*vector_ptr)[index(i_c,j_c,k)].x_sz + (*vector_ptr)[index(i,j_c,k)].x_sz
					+ (*vector_ptr)[index(i_c,j,k)].x_sz + (*vector_ptr)[index(i,j,k)].x_sz;
				break;
			case 15:
				fval = 	  (*vector_ptr)[index(i_c,j_c,k_c)].x_tz + (*vector_ptr)[index(i,j_c,k_c)].x_tz
					+ (*vector_ptr)[index(i_c,j,k_c)].x_tz + (*vector_ptr)[index(i,j,k_c)].x_tz
					+ (*vector_ptr)[index(i_c,j_c,k)].x_tz + (*vector_ptr)[index(i,j_c,k)].x_tz
					+ (*vector_ptr)[index(i_c,j,k)].x_tz + (*vector_ptr)[index(i,j,k)].x_tz;
				break;
		}

		double factor = 1.;
		if (i == half_order) {factor /= 2.;}
		if (j == half_order) {factor /= 2.;}
		if (k == half_order) {factor /= 2.;}

		result_int_y += fval * wi[i] * wi[j] * wi[k] * factor;
	}
		result_int_x += result_int_y;
	}
		result_int += result_int_x;
	}
	result_int *= half_length1 * half_length2 * half_length3;

	return (result_int);
}

/*
int index(int i, int j, int k)
{
	return (order*order*i+order*j+k);
}
*/
