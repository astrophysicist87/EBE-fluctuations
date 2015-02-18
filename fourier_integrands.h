/*************************************************************************/
double fourier_kernel_real(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (dim)
	{
		case 1:
			return FI_real_1D(q_ptr, x_ptr);
			break;
		case 2:
			return FI_real_2D(q_ptr, x_ptr);
			break;
		case 3:
			return FI_real_3D(q_ptr, x_ptr);
			break;
		case 4:
			return FI_real_4D(q_ptr, x_ptr);
			break;
	}
}
/*************************************************************************/

double FI_real_1D(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			return (cos(q*x));
			break;
		case 2:	//cylindrical
			return (cos(q*x));
			break;
		case 3:	//tau, eta, r, phi
			return (cos(q*x));
			break;
	}
}

double FI_real_2D(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			//return (cos(q*r*cos(q_phi - phi)));
			break;
		case 2:	//cylindrical
			break;
		case 3:	//tau, eta, r, phi
			break;
	}
}

double FI_real_3D(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			break;
		case 2:	//cylindrical
			break;
		case 3:	//tau, eta, r, phi
			break;
	}
}

double FI_real_4D(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			break;
		case 2:	//cylindrical
			break;
		case 3:	//tau, eta, r, phi
			break;
	}
}

/*************************************************************************/
double fourier_kernel_imag(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (dim)
	{
		case 1:
			return FI_imag_1D(q_ptr, x_ptr);
			break;
		case 2:
			return FI_imag_2D(q_ptr, x_ptr);
			break;
		case 3:
			return FI_imag_3D(q_ptr, x_ptr);
			break;
		case 4:
			return FI_imag_4D(q_ptr, x_ptr);
			break;
	}
}
/*************************************************************************/

double FI_imag_1D(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (coords)
	{
		case 0:	//cartesian
			return (sin(q*x));
			break;
		case 1:	//spherical
			return (sin(q*x));
			break;
		case 2:	//cylindrical
			return (sin(q*x));
			break;
		case 3:	//tau, eta, r, phi
			return (sin(q*x));
			break;
	}
}

double FI_imag_2D(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			break;
		case 2:	//cylindrical
			break;
		case 3:	//tau, eta, r, phi
			break;
	}
}

double FI_imag_3D(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			break;
		case 2:	//cylindrical
			break;
		case 3:	//tau, eta, r, phi
			break;
	}
}

double FI_imag_4D(vector<double>* q_ptr, vector<double>* x_ptr, int dim)
{
	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			break;
		case 2:	//cylindrical
			break;
		case 3:	//tau, eta, r, phi
			break;
	}
}

//End of file
