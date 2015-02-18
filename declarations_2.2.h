const int order = 63;
const int half_order = (order - 1)/2;
//const double limit = 4.5;
const double PI = 3.14159265358979323846264338327950;
const double M1_4PI = 1./(4.*PI);
const double largest_double = std::numeric_limits<double>::max();
const int coords = 1;
const int dim = 2;
const complex<double> complex_i (0.,1.);
int pdf_switches[] = {0};	//0 - uniform distribution
				//1 - Gaussian distribution
				//2 - Bessel-Gaussian distribution
				//3 - quadratic distribution
				//4 - alternative Gaussian distribution
int params2vary[] = {6};
double pmean[] = {0.25};
double psigma[] = {0.2};
//int params2vary[] = {3,4};
//double pmean[] = {5., 1.};
//double psigma[] = {1., 0.5};
int vec_of_intervals[] = {0};

double sigma_cmd, mean_cmd, R_cmd, N_cmd, r0_cmd;
double eps_n_bar_cmd, harmonic_n_cmd, psi_n_bar_cmd;
double sigma_default = 0.5;
double mean_default = 0.0;
double R_default = 5.;
double N_default = 1.;
double r0_default = 0.;
double eps_n_bar_default = 0.;
double harmonic_n_default = 2.;
double psi_n_bar_default = 0.;

bool my_flag = false;

double fluct_lower, fluct_upper;
double q_lower = 0., q_upper = 0.75;
//const int nq_points = 51;
const int nq_points1 = 51;
const int nq_points2 = 51;
double q_interval1 = (nq_points1 == 1)? q_upper - q_lower : (q_upper - q_lower)/double(nq_points1-1);
double q_interval2 = (nq_points2 == 1)? 2.*PI : (PI - (-PI))/double(nq_points2-1);
double Phi_K_a = -PI, Phi_K_b = PI;
double Phi_K_interval = (Phi_K_b-Phi_K_a)/double(order);
double phi_lower=-PI, phi_upper=PI;
double theta_lower=0., theta_upper=PI;
double phi_interval = (phi_upper-phi_lower)/double(order);

const int program_function = 0;	//0 - compute correlation functions at nq_points and output (useful for plotting)
				//1 - extract normalized curvature from correlation functions at origin
				//2 - compute q-moments of correlation functions and use to extract real HBT radii

double ens_avg_norm = 0.;

vector<double>* xi_ptr;
vector<double>* wi_ptr;
vector<double>* xi_0pinf_ptr;
vector<double>* wi_0pinf_ptr;
vector<double>* xi_minfpinf_ptr;
vector<double>* wi_minfpinf_ptr;
vector<double>* xiq_ptr;
vector<double>* wiq_ptr;
vector<double>* xiq_0pinf_ptr;
vector<double>* wiq_0pinf_ptr;
vector<double>* xiq_minfpinf_ptr;
vector<double>* wiq_minfpinf_ptr;

struct pdf_params
{
	vector<double>* mean_vec_ptr;
	vector<double>* sigma_vec_ptr;
	vector<double>* fluct_lower_ptr;
	vector<double>* fluct_upper_ptr;
};

struct model_params
{
	double mean;
	double sigma;
	double norm;
	double r0;
	double R;
	double N;
	double eps_n_bar;
	double harmonic_n;
	double psi_n_bar;
};

struct position
{
	double coord1, coord2, coord3, coord4;
	double q_index;
};

struct qposition
{
	double qcoord1, qcoord2, qcoord3, qcoord4;
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
vector<qposition>* q_vec_ptr;

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

//need these to define vector of parameters to fluctuate
template <typename T, size_t N>
T* begin(T(&arr)[N]) { return &arr[0]; }
template <typename T, size_t N>
T* end(T(&arr)[N]) { return &arr[0]+N; }

std::vector<int> paramsvec2vary(begin(params2vary), end(params2vary));
std::vector<int> pdfswitchvec(begin(pdf_switches), end(pdf_switches));
//std::vector<int> intervalvec (paramsvec2vary.size());
std::vector<int> intervalvec (begin(vec_of_intervals), end(vec_of_intervals));
std::vector<double> params_mean(begin(pmean), end(pmean));
std::vector<double> params_sigma(begin(psigma), end(psigma));

vector<int> convert_base(int number, int base)
{
	const int length = 1 + int(log(number)/log(base));
	vector<int> result(length);
	for (int dim = 1; dim <= length; dim++) result[dim-1] = (number%int(pow(base,dim)))/int(pow(base,dim-1));
	//for (int dim = 1; dim <= length; dim++) {result[length-dim] = (number%int(pow(base,dim)))/int(pow(base,dim-1));  cout << result[length-dim] << " ";}

	return (result);
}

vector<int> convert_base(int number, int base, int length)
{
	//const int length = 1 + int(log(number)/log(base));
	vector<int> result(length);
	for (int dim = 1; dim <= length; dim++) result[dim-1] = (number%int(pow(base,dim)))/int(pow(base,dim-1));
	//for (int dim = 1; dim <= length; dim++) {result[length-dim] = (number%int(pow(base,dim)))/int(pow(base,dim-1));  cout << result[length-dim] << " ";}

	return (result);
}

//End of file
