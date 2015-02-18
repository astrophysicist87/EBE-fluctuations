void print_params_to_output(string outfilename, void * params_ptr1, void * params_ptr2);
bool fexists(const char *filename);

void print_params_to_output(string outfilename, void * params_ptr1, void * params_ptr2)
{
	pdf_params params1 = * (struct pdf_params *) params_ptr1;
	model_params params2 = * (struct model_params *) params_ptr2;

	ofstream outputparams (outfilename.data());

	outputparams << "PDF parameters:" << outfilename << endl;
	outputparams << "parameter to vary = " << parameter_to_vary << endl;
	outputparams << "sigma = " << params1.sigma << endl;
	outputparams << "mean = " << params1.mean << endl << endl;

	outputparams << "Model parameters:" << outfilename << endl;
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
