const int order = 83;
const int half_order = (order - 1)/2;
const double limit = 5.;
const double PI = 3.14159265358979323846264338327950;
const double M1_4PI = 1./(4.*PI);
const int coords = 0;
const int dim = 1;
const complex<double> complex_i (0.,1.);
//double r_lower=0., r_upper=3.*limit, phi_lower=-PI, phi_upper=PI;
//double eta_lower=-limit, eta_upper=limit, tau_lower=0., tau_upper=limit;

double r_lower=0., r_upper=4.*limit, phi_lower=-PI, phi_upper=PI;
double eta_lower=-limit, eta_upper=limit, tau_lower=0., tau_upper=4.*limit;
double t_lower=-limit, t_upper=limit, x_lower=-4.*limit, x_upper=4.*limit;
double y_lower=-limit, y_upper=limit, z_lower=-4.*limit, z_upper=4.*limit;
double theta_lower=0., theta_upper=PI;
double fluct_lower=1., fluct_upper=9.;
double q_lower = 0., q_upper = 2.;
const int nq_points = 50;
double q_interval = (q_upper - q_lower)/double(nq_points);
//double xi[half_order];
//double wi[half_order+1];
double Phi_K_a = -PI, Phi_K_b = PI;
double Phi_K_interval = (Phi_K_b-Phi_K_a)/double(order);
double phi_interval = 0.;

static double xi[half_order];
static double wi[half_order+1];

struct pdf_params
{
	double mean;
	double sigma;
	double norm;
};

struct model_params
{
	double mean;
	double sigma;
	double norm;
};

struct position
{
	double t;
	double x;
	double y;
	double z;
	double r;
	double phi;
	double theta;
	double tau;
	double eta;
	double q_index;
};

pdf_params my_pdf_params;
model_params my_model_params;
position my_pos;
position lower_limit_vec;
position upper_limit_vec;

vector<position>* grid_ptr;
vector<double>* ens_avg_ptr;
vector<complex<double> >* FT_ens_avg_ptr;
vector<complex<double> >* FT_model_function_ptr;
vector<double>* FT_ens_avg2_ptr;
vector<double>* FT_model_function2_ens_avg_ptr;
vector<double>* Cbar_ptr;
vector<double>* C_ens_avg_ptr;

int index(int i, int j)
{
	return (order*i+j);
}

int index(int i, int j, int k)
{
	return (order*order*i+order*j+k);
}

int index(int i, int j, int k, int l)
{
	return (order*order*order*i+order*order*j+order*k+l);
}

//End of file
