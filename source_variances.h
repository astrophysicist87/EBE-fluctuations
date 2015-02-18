#ifndef SOURCE_VARIANCES_H
#define SOURCE_VARIANCES_H

using namespace std;

double normS(vector<double>*, void *);
double xS(vector<double>*, void *);
double yS(vector<double>*, void *);
double x2S(vector<double>*, void *);
double y2S(vector<double>*, void *);
double xyS(vector<double>*, void *);

double normS(vector<double>* x_ptr, void * params_ptr)
{
	double jac;
	switch (coords)
	{
		case 0:
			jac = 1.;
			break;
		case 1:
			jac = (*x_ptr).at(0);
			break;
	}

	return (jac*S_function(x_ptr,params_ptr));
}

double xS(vector<double>* x_ptr, void * params_ptr)
{
	double x;
	double jac;
	double Sval;
	switch (coords)
	{
		case 0:
			x=(*x_ptr).at(0);
			jac = 1.;
			break;
		case 1:
			x=(*x_ptr).at(0) * cos((*x_ptr).at(1));
			jac = (*x_ptr).at(0);
			break;
	}

	return (jac*x*S_function(x_ptr,params_ptr));
}

double yS(vector<double>* x_ptr, void * params_ptr)
{
	double y;
	double jac;
	switch (coords)
	{
		case 0:
			y=(*x_ptr).at(1);
			jac = 1.;
			break;
		case 1:
			y=(*x_ptr).at(0) * sin((*x_ptr).at(1));
			jac = (*x_ptr).at(0);
			break;
	}

	return (jac*y*S_function(x_ptr,params_ptr));
}

double x2S(vector<double>* x_ptr, void * params_ptr)
{
	double x;
	double jac;
	switch (coords)
	{
		case 0:
			x=(*x_ptr).at(0);
			jac = 1.;
			break;
		case 1:
			x=(*x_ptr).at(0) * cos((*x_ptr).at(1));
			jac = (*x_ptr).at(0);
			break;
	}

	return (jac*x*x*S_function(x_ptr,params_ptr));
}

double y2S(vector<double>* x_ptr, void * params_ptr)
{
	double y;
	double jac;
	switch (coords)
	{
		case 0:
			y=(*x_ptr).at(1);
			jac = 1.;
			break;
		case 1:
			y=(*x_ptr).at(0) * sin((*x_ptr).at(1));
			jac = (*x_ptr).at(0);
			break;
	}

	return (jac*y*y*S_function(x_ptr,params_ptr));
}

double xyS(vector<double>* x_ptr, void * params_ptr)
{
	double x, y;
	double jac;
	switch (coords)
	{
		case 0:
			x=(*x_ptr).at(0);
			y=(*x_ptr).at(1);
			jac = 1.;
			break;
		case 1:
			x=(*x_ptr).at(0) * cos((*x_ptr).at(1));
			y=(*x_ptr).at(0) * sin((*x_ptr).at(1));
			jac = (*x_ptr).at(0);
			break;
	}

	return (jac*x*y*S_function(x_ptr,params_ptr));
}

vector<double> get_source_variances(void * params_ptr)
{
	vector<double> result (3);
	vector<int> local_intervalsvec(2);

	local_intervalsvec[0] = 1;	//r-integral
	local_intervalsvec[1] = 0;	//phi-integral

	//assume polar coordinates for timebeing
	vector<double> local_lower_limits_vec (2);
	vector<double> local_upper_limits_vec (2);

	local_lower_limits_vec[0] = 0.;
	local_lower_limits_vec[1] = -PI;
	local_upper_limits_vec[0] = 0.;
	local_upper_limits_vec[1] = PI;

	//EA stands for "ensemble average"
	//EBE stands for "event-by-event"
	double norm = integrate(&normS, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double xEA = integrate(&xS, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double yEA = integrate(&yS, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double x2EA = integrate(&x2S, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double y2EA = integrate(&y2S, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double xyEA = integrate(&xyS, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);

	result.at(0) = x2EA/norm - (xEA/norm)*(xEA/norm);
	result.at(1) = xyEA/norm - (xEA/norm)*(yEA/norm);
	result.at(2) = y2EA/norm - (yEA/norm)*(yEA/norm);

	return (result);
}

vector<double> get_source_variances_fromEA()
{
	vector<double> result (3);

	vector<double>* temp1_ptr;
	temp1_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* temp2_ptr;
	temp2_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* temp3_ptr;
	temp3_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* temp4_ptr;
	temp4_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* temp5_ptr;
	temp5_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* temp6_ptr;
	temp6_ptr = new vector<double> (int(pow(order,dim)));
	vector<double>* x_ptr;
	x_ptr = new vector<double> (dim);
	double x, y;
	double jac;

		for (int j = 0; j <= order-1; j++)
		{
		for (int k = 0; k <= order-1; k++)
		{
		//if (dim < 2 && k > 0) continue;
		(*x_ptr).at(0) = (*grid_ptr).at(local_index(j,k)).coord1;
		if (dim > 1) (*x_ptr).at(1) = (*grid_ptr).at(local_index(j,k)).coord2;
	switch (coords)
	{
		case 0:
			x=(*x_ptr).at(0);
			y=(*x_ptr).at(1);
			jac = 1.;
			break;
		case 1:
			x=(*x_ptr).at(0) * cos((*x_ptr).at(1));
			y=(*x_ptr).at(0) * sin((*x_ptr).at(1));
			jac = (*x_ptr).at(0);
			break;
	}
		(*temp1_ptr).at(local_index(j,k)) = jac*(*ens_avg_ptr).at(local_index(j,k));		//norm
		(*temp2_ptr).at(local_index(j,k)) = jac*x*(*ens_avg_ptr).at(local_index(j,k));		//<x>
		(*temp3_ptr).at(local_index(j,k)) = jac*y*(*ens_avg_ptr).at(local_index(j,k));		//<y>
		(*temp4_ptr).at(local_index(j,k)) = jac*x*x*(*ens_avg_ptr).at(local_index(j,k));		//<x2>
		(*temp5_ptr).at(local_index(j,k)) = jac*x*y*(*ens_avg_ptr).at(local_index(j,k));		//<xy>
		(*temp6_ptr).at(local_index(j,k)) = jac*y*y*(*ens_avg_ptr).at(local_index(j,k));		//<y2>
	}
	}
	double norm = do_integrations(temp1_ptr);
//cout << "CHECK IN SVEA: norm = " << norm << endl;
	double xfromEA = do_integrations(temp2_ptr);
//cout << "CHECK IN SVEA: xfromEA = " << norm << endl;
	double yfromEA = do_integrations(temp3_ptr);
	double x2fromEA = do_integrations(temp4_ptr);
//cout << "CHECK IN SVEA: x2fromEA = " << norm << endl;
	double xyfromEA = do_integrations(temp5_ptr);
	double y2fromEA = do_integrations(temp6_ptr);

	delete temp1_ptr;
	delete temp2_ptr;
	delete temp3_ptr;
	delete temp4_ptr;
	delete temp5_ptr;
	delete temp6_ptr;
	delete x_ptr;

	result.at(0) = x2fromEA/norm - (xfromEA/norm)*(xfromEA/norm);
	result.at(1) = xyfromEA/norm - (xfromEA/norm)*(yfromEA/norm);
	result.at(2) = y2fromEA/norm - (yfromEA/norm)*(yfromEA/norm);

	return (result);
}

/*
double xSnormed(vector<double>* x_ptr, void * params_ptr)
{
	double Sval = (abs(S_function(x_ptr, params_ptr)) <= 1e-10) ? 1. : S_function(x_ptr, params_ptr);
	return (xS(x_ptr, params_ptr)/Sval);
}
double ySnormed(vector<double>* x_ptr, void * params_ptr)
{
	double Sval = (abs(S_function(x_ptr, params_ptr)) <= 1e-10) ? 1. : S_function(x_ptr, params_ptr);
	return (yS(x_ptr, params_ptr)/Sval);
}
double x2Snormed(vector<double>* x_ptr, void * params_ptr)
{
	double Sval = (abs(S_function(x_ptr, params_ptr)) <= 1e-10) ? 1. : S_function(x_ptr, params_ptr);
	return (x2S(x_ptr, params_ptr)/Sval);
}
double xySnormed(vector<double>* x_ptr, void * params_ptr)
{
	double Sval = (abs(S_function(x_ptr, params_ptr)) <= 1e-10) ? 1. : S_function(x_ptr, params_ptr);
	return (xyS(x_ptr, params_ptr)/Sval);
}
double y2Snormed(vector<double>* x_ptr, void * params_ptr)
{
	double Sval = (abs(S_function(x_ptr, params_ptr)) <= 1e-10) ? 1. : S_function(x_ptr, params_ptr);
	return (y2S(x_ptr, params_ptr)/Sval);
}

vector<double> get_EAsource_variances(void * params_ptr)
{
	vector<double> result (3);
	vector<int> local_intervalsvec(2);

	local_intervalsvec[0] = 1;	//r-integral
	local_intervalsvec[1] = 0;	//phi-integral

//assume polar coordinates for timebeing
	vector<double> local_lower_limits_vec (2);
	vector<double> local_upper_limits_vec (2);

	local_lower_limits_vec[0] = 0.;
	local_lower_limits_vec[1] = -PI;
	local_upper_limits_vec[0] = 0.;
	local_upper_limits_vec[1] = PI;

//EA stands for "ensemble average"
	double xEAnorm = integrate(&xSnormed, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double yEAnorm = integrate(&ySnormed, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double x2EAnorm = integrate(&x2Snormed, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double y2EAnorm = integrate(&y2Snormed, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double xyEAnorm = integrate(&xySnormed, params_ptr, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);

	result.at(0) = x2EAnorm - (xEAnorm)*(xEAnorm);
	result.at(1) = xyEAnorm - (xEAnorm)*(yEAnorm);
	result.at(2) = y2EAnorm - (yEAnorm)*(yEAnorm);

	return (result);
}
*/

double x2varEBE_integrand(vector<double>* fluct_vec_ptr, void * params_ptr)
{
	vector<int> local_intervalsvec(2);

	local_intervalsvec[0] = 1;	//r-integral
	local_intervalsvec[1] = 0;	//phi-integral

	//assume polar coordinates for timebeing
	vector<double> local_lower_limits_vec (2);
	vector<double> local_upper_limits_vec (2);

	local_lower_limits_vec[0] = 0.;
	local_lower_limits_vec[1] = -PI;
	local_upper_limits_vec[0] = 0.;
	local_upper_limits_vec[1] = PI;

	model_params my_model_params_copy = my_model_params;
	for (int i=0; i<(*fluct_vec_ptr).size(); i++)
	{
		switch(paramsvec2vary[i])
		{
			case 3:
			my_model_params_copy.R = (*fluct_vec_ptr)[i];
			break;
			case 4:
			my_model_params_copy.N = (*fluct_vec_ptr)[i];
			break;
			case 5:
			my_model_params_copy.r0 = (*fluct_vec_ptr)[i];
			break;
			case 6:
			my_model_params_copy.eps_n_bar = (*fluct_vec_ptr)[i];
			break;
			case 7:
			my_model_params_copy.harmonic_n = (*fluct_vec_ptr)[i];
			break;
			case 8:
			my_model_params_copy.psi_n_bar = (*fluct_vec_ptr)[i];
			break;
		}
	}

	double prob = pdf(fluct_vec_ptr, &my_pdf_params);

	//EBE stands for "event-by-event"
	double normEBE = integrate(&normS, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
//cout << "CHECK IN SVEBE: normEBE = " << normEBE << endl;
	double xEBE = integrate(&xS, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
//cout << "CHECK IN SVEBE: xEBE = " << normEBE << endl;
	double x2EBE = integrate(&x2S, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
//cout << "CHECK IN SVEBE: x2EBE = " << normEBE << endl;

	double result = x2EBE*normEBE - xEBE*xEBE;	//from extra factor of N^2_k in definition of EBE R^2_{ij}

	return (prob*result);
}

double xyvarEBE_integrand(vector<double>* fluct_vec_ptr, void * params_ptr)
{
	vector<int> local_intervalsvec(2);

	local_intervalsvec[0] = 1;	//r-integral
	local_intervalsvec[1] = 0;	//phi-integral

	//assume polar coordinates for timebeing
	vector<double> local_lower_limits_vec (2);
	vector<double> local_upper_limits_vec (2);

	local_lower_limits_vec[0] = 0.;
	local_lower_limits_vec[1] = -PI;
	local_upper_limits_vec[0] = 0.;
	local_upper_limits_vec[1] = PI;

	model_params my_model_params_copy = my_model_params;
	for (int i=0; i<(*fluct_vec_ptr).size(); i++)
	{
		switch(paramsvec2vary[i])
		{
			case 3:
			my_model_params_copy.R = (*fluct_vec_ptr)[i];
			break;
			case 4:
			my_model_params_copy.N = (*fluct_vec_ptr)[i];
			break;
			case 5:
			my_model_params_copy.r0 = (*fluct_vec_ptr)[i];
			break;
			case 6:
			my_model_params_copy.eps_n_bar = (*fluct_vec_ptr)[i];
			break;
			case 7:
			my_model_params_copy.harmonic_n = (*fluct_vec_ptr)[i];
			break;
			case 8:
			my_model_params_copy.psi_n_bar = (*fluct_vec_ptr)[i];
			break;
		}
	}

	double prob = pdf(fluct_vec_ptr, &my_pdf_params);

	//EBE stands for "event-by-event"
	double normEBE = integrate(&normS, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double xEBE = integrate(&xS, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double yEBE = integrate(&yS, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double xyEBE = integrate(&xyS, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);

	double result = xyEBE*normEBE - xEBE*yEBE;	//from extra factor of N^2_k in definition of EBE R^2_{ij}

	return (prob*result);
}

double y2varEBE_integrand(vector<double>* fluct_vec_ptr, void * params_ptr)
{
	vector<int> local_intervalsvec(2);

	local_intervalsvec[0] = 1;	//r-integral
	local_intervalsvec[1] = 0;	//phi-integral

	//assume polar coordinates for timebeing
	vector<double> local_lower_limits_vec (2);
	vector<double> local_upper_limits_vec (2);

	local_lower_limits_vec[0] = 0.;
	local_lower_limits_vec[1] = -PI;
	local_upper_limits_vec[0] = 0.;
	local_upper_limits_vec[1] = PI;

	model_params my_model_params_copy = my_model_params;
	for (int i=0; i<(*fluct_vec_ptr).size(); i++)
	{
		switch(paramsvec2vary[i])
		{
			case 3:
			my_model_params_copy.R = (*fluct_vec_ptr)[i];
			break;
			case 4:
			my_model_params_copy.N = (*fluct_vec_ptr)[i];
			break;
			case 5:
			my_model_params_copy.r0 = (*fluct_vec_ptr)[i];
			break;
			case 6:
			my_model_params_copy.eps_n_bar = (*fluct_vec_ptr)[i];
			break;
			case 7:
			my_model_params_copy.harmonic_n = (*fluct_vec_ptr)[i];
			break;
			case 8:
			my_model_params_copy.psi_n_bar = (*fluct_vec_ptr)[i];
			break;
		}
	}

	double prob = pdf(fluct_vec_ptr, &my_pdf_params);

	//EBE stands for "event-by-event"
	double normEBE = integrate(&normS, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double yEBE = integrate(&yS, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);
	double y2EBE = integrate(&y2S, &my_model_params_copy, &local_lower_limits_vec, &local_upper_limits_vec, &local_intervalsvec);

	double result = y2EBE*normEBE - yEBE*yEBE;	//from extra factor of N^2_k in definition of EBE R^2_{ij}

	return (prob*result);
}

vector<double> get_source_variances_fromEBE()
{
	intervalvec[0] = 0;	//finite integration range
	vector<double> result (3);

//note that first call of &intervalvec is currently just a dummy argument - not used for anything here
	result.at(0) = integrate(&x2varEBE_integrand, &intervalvec, my_pdf_params.fluct_lower_ptr, my_pdf_params.fluct_upper_ptr, &intervalvec);
	result.at(1) = integrate(&xyvarEBE_integrand, &intervalvec, my_pdf_params.fluct_lower_ptr, my_pdf_params.fluct_upper_ptr, &intervalvec);
	result.at(2) = integrate(&y2varEBE_integrand, &intervalvec, my_pdf_params.fluct_lower_ptr, my_pdf_params.fluct_upper_ptr, &intervalvec);

	return (result);
}

//End of file

#endif
