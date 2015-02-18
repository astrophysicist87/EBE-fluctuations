#ifndef OUTFILENAME_2_2_H
#define OUTFILENAME_2_2_H

using namespace std;

void print_params_to_output(string outfilename, void * params_ptr1, void * params_ptr2);
bool fexists(const char *filename);

void print_params_to_output(string outfilename, void * params_ptr1, void * params_ptr2)
{
	pdf_params params1 = * (struct pdf_params *) params_ptr1;
	model_params params2 = * (struct model_params *) params_ptr2;

	ofstream outputparams (outfilename.data());

	outputparams << "order = " << order << endl;
	outputparams << "q_lower = " << q_lower << endl;
	outputparams << "q_upper = " << q_upper << endl;
	outputparams << "nq_points1 = " << nq_points1 << endl;
	outputparams << "nq_points2 = " << nq_points2 << endl;

	for (int i=0; i<pdfswitchvec.size(); i++)
	{
		outputparams << "pdfswitchvec(" << i << ") = " << pdfswitchvec[i] << endl;
	}

	outputparams << "PDF parameters:" << endl;
	//outputparams << "parameter to vary = " << parameter_to_vary << endl;
	for (int i=0; i<(*(params1.mean_vec_ptr)).size(); i++)
	{
	outputparams << "parameter to vary(" << i << ") = " << paramsvec2vary[i] << endl;
	outputparams << "sigma(" << i << ") = " << (*(params1.sigma_vec_ptr))[i] << endl;
	outputparams << "mean(" << i << ") = " << (*(params1.mean_vec_ptr))[i] << endl;
	}
	cout << endl;

	outputparams << "Model parameters:" << endl;
	outputparams << "R = " << params2.R << endl;
	outputparams << "N = " << params2.N << endl;
	outputparams << "r0 = " << params2.r0 << endl;
	outputparams << "eps_n_bar = " << params2.eps_n_bar << endl;
	outputparams << "harmonic_n = " << params2.harmonic_n << endl;
	outputparams << "psi_n_bar = " << params2.psi_n_bar << endl << endl;

	outputparams.close();
}

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

//End of file

#endif
