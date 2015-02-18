#ifndef FOURIER_KERNELS_H
#define FOURIER_KERNELS_H

using namespace std;

double FI_real_1D(vector<double>* q_ptr, vector<double>* x_ptr);
double FI_real_2D(vector<double>* q_ptr, vector<double>* x_ptr);
double FI_real_3D(vector<double>* q_ptr, vector<double>* x_ptr);
double FI_imag_1D(vector<double>* q_ptr, vector<double>* x_ptr);
double FI_imag_2D(vector<double>* q_ptr, vector<double>* x_ptr);
double FI_imag_3D(vector<double>* q_ptr, vector<double>* x_ptr);

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
	}
}
/*************************************************************************/

double FI_real_1D(vector<double>* q_ptr, vector<double>* x_ptr)
{
	double q, x;

	switch (coords)
	{
		case 0:	//cartesian
			q = (*q_ptr).at(0);
			x = (*x_ptr).at(0);
			return (cos(q*x));
			break;
		case 1:	//spherical
			q = (*q_ptr).at(0);
			x = (*x_ptr).at(0);
			return (cos(q*x));
			break;
		case 2:	//cylindrical
			q = (*q_ptr).at(0);
			x = (*x_ptr).at(0);
			return (cos(q*x));
			break;
	}
}

double FI_real_2D(vector<double>* q_ptr, vector<double>* x_ptr)
{
	double q, qx, qy, q_phi;
	double r, x, y, phi;

	switch (coords)
	{
		case 0:	//cartesian
			qx = (*q_ptr).at(0);
			qy = (*q_ptr).at(1);
			x = (*x_ptr).at(0);
			y = (*x_ptr).at(1);
			return (cos(qx*x+qy*y));
			break;
		case 1:	//spherical
			q = (*q_ptr).at(0);
			q_phi = (*q_ptr).at(1);
			r = (*x_ptr).at(0);
			phi = (*x_ptr).at(1);
			return (r*cos(q*r*cos(q_phi-phi)));
			break;
		case 2:	//cylindrical
			q = (*q_ptr).at(0);
			q_phi = (*q_ptr).at(1);
			r = (*x_ptr).at(0);
			phi = (*x_ptr).at(1);
			return (r*cos(q*r*cos(q_phi-phi)));
			break;
	}
}

double FI_real_3D(vector<double>* q_ptr, vector<double>* x_ptr)
{
	double q, x;

	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			break;
		case 2:	//cylindrical
			break;
	}
}

double FI_real_4D(vector<double>* q_ptr, vector<double>* x_ptr)
{
	double q, x;

	switch (coords)
	{
		case 0:	//cartesian
			return (cos(q*x));
			break;
		case 1:	//spherical
			break;
		case 2:	//cylindrical
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
	}
}
/*************************************************************************/

double FI_imag_1D(vector<double>* q_ptr, vector<double>* x_ptr)
{
	double q, x;

	switch (coords)
	{
		case 0:	//cartesian
			q = (*q_ptr).at(0);
			x = (*x_ptr).at(0);
			return (sin(q*x));
			break;
		case 1:	//spherical
			q = (*q_ptr).at(0);
			x = (*x_ptr).at(0);
			return (sin(q*x));
			break;
		case 2:	//cylindrical
			q = (*q_ptr).at(0);
			x = (*x_ptr).at(0);
			return (sin(q*x));
			break;
	}
}

double FI_imag_2D(vector<double>* q_ptr, vector<double>* x_ptr)
{
	double q, qx, qy, q_phi;
	double r, x, y, phi;

	switch (coords)
	{
		case 0:	//cartesian
			qx = (*q_ptr).at(0);
			qy = (*q_ptr).at(1);
			x = (*x_ptr).at(0);
			y = (*x_ptr).at(1);
			return (sin(qx*x+qy*y));
			break;
		case 1:	//spherical
			q = (*q_ptr).at(0);
			q_phi = (*q_ptr).at(1);
			r = (*x_ptr).at(0);
			phi = (*x_ptr).at(1);
			return (r*sin(q*r*cos(q_phi-phi)));
			break;
		case 2:	//cylindrical
			q = (*q_ptr).at(0);
			q_phi = (*q_ptr).at(1);
			r = (*x_ptr).at(0);
			phi = (*x_ptr).at(1);
			return (r*sin(q*r*cos(q_phi-phi)));
			break;
	}
}

double FI_imag_3D(vector<double>* q_ptr, vector<double>* x_ptr)
{
	double q, qx, qy, qz, q_phi;
	double r, x, y, z, phi;

	switch (coords)
	{
		case 0:	//cartesian
			qx = (*q_ptr).at(0);
			qy = (*q_ptr).at(1);
			x = (*x_ptr).at(0);
			y = (*x_ptr).at(1);
			return (sin(qx*x+qy*y+qz*z));
		case 1:	//spherical
			break;
		case 2:	//cylindrical
			break;
	}
}

//End of file

#endif
