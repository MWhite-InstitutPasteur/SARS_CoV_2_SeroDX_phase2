#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>
#include "randlib.h"
#include <omp.h>
#include <vector>
#include <algorithm>

using namespace std;


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//          //                                        //
//   ####   //  Setting up objects,                   //
//  ##  ##  //  declaring array sizes                 //
//  ##  ##  //                                        //
//  ##  ##  //                                        //
//   ####   //                                        //
//          //                                        //
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
// 0.1 Declare global variables of sizes of arrays

#define N_mcmc 10000000        // Number of MCMC iterations: indexed by mc

#define N_adapt 600000        // Number of MCMC iterations in adaptive phase
#define N_tune_start 10000    // Number of MCMC iterations in adaptive phase
#define N_tune_end 500000

#define N_data_cols 33        // Number of columns in data frame
#define N_part 194             // Number of participants to read in data from
#define N_negcon 270           // Number of negative control participants to read in data from
#define N_t 15                // Maximum number of data points per participant
#define N_loc_par 8           // Number of individual-level parameters
#define N_glob_par 8          // Number of population-level parameters
#define LARGE 1e12            // large number needed for priors

#define log2 0.69314718055995
#define sqrt2 1.414214

/////////////////////////////////////////////////////////////////	
// 0.2 Create structure to hold data for participant n
//     and local parameter estimates

struct part_n
{
	//////////////////////////////////////
	// Covariate information

	int site;                // 1 = Bichat; 2 = Strasbourg; 3 = Cochin; 4 = Thailand; 5 = Peru; 6 = EFS
	int status;              // negative = 1; positive = 2
	int symptoms;            // mild = 1; severe = 2


	//////////////////////////////////////
	// Antibody data

	int N_sam;                   // Number of samples of antibody data

	vector<double> AB;
	vector<double> tt;

	vector<double> lAB;

	double AB_min;   // Minimum antibody level for the individual
	double AB_max;   // Maximum antibody level for the individual

	double lAB_min;  // Log minimum antibody level for the individual
	double lAB_max;  // Log maximum antibody level for the individual

	double lNC_max;  // Log of maximum background antibody level (determined by negative controls)

	//////////////////////////////////////
	// Individual-level parameters

	double Ab_bg;         // background antibody level
	double beta;          // boost in ASCs - vaccine dose 1
	double tau;           // delay in boosting of antibody responses
	double t_delta;       // half-life of memory B cells
	double t_short;       // half-life of short-lived ASCs
	double t_long;        // half-life of long-lived ASCs
	double t_IgG;         // half-life of IgG molecules
	double rho;           // proportion of short-lived ASCs

	double lAb_bg;        // log(boost in antibody levels)
	double lbeta;         // log(boost in ASCs) 
	double ltau;          // delay in boosting of antibody responses
	double lt_delta;      // half-life of memory B cells
	double lt_short;      // half-life of short-lived ASCs
	double lt_long;       // half-life of long-lived ASCs
	double lt_IgG;        // half-life of IgG molecules
	double logitrho;      // logit proportion of short-lived ASCs - vaccine dose 1

	double r_delta;       // drug decay rate
	double r_short;       // drug decay rate
	double r_long;        // drug decay rate
	double r_IgG;         // drug decay rate


	//////////////////////////////////////
	// Likelihood

	double data_like;     // data likelihood
	double mix_like;      // mixed-effects likelihood

	double lhood;         // individual-level likelihood
};


/////////////////////////////////////////////////////////////////
// 0.3 Create structure for global parameters to be estimated

struct params
{
	/////////////////////////////////////////////////////////////
	// Population-level parameters describing mixed effects

	double mu_par[N_glob_par];
	double tau_par[N_glob_par];


	/////////////////////////////////////////////////////////////
	// Maximum and minimum antibody levels

	double AB_min;   // Global minimum antibody level
	double AB_max;   // Global maximum antibody level

	double lAB_min;  // Log global minimum antibody level
	double lAB_max;  // Log global maximum antibody level


	/////////////////////////////////////////////////////////////
	// Parameter for observational error

	double sig_obs;        // precision of observational error
	double log_sig_obs;

	double sig_obs_scale;


	/////////////////////////////////////////////////////////////
	// Log likelihood and prior

	double loglike;
	double prior;


	/////////////////////////////////////////////////////////////
	// Prior distributions

	double prior_MM[N_glob_par];
	double prior_MM_CV[N_glob_par];
	double prior_SIG[N_glob_par];
	double prior_SIG_CV[N_glob_par];

	double prior_LN_MM[N_glob_par];
	double prior_LN_MM_CV[N_glob_par];
	double prior_LN_SIG[N_glob_par];
	double prior_LN_SIG_CV[N_glob_par];

	double prior_mu[N_glob_par];
	double prior_tau[N_glob_par];
	double prior_k[N_glob_par];
	double prior_theta[N_glob_par];


	/////////////////////////////////////////////////////////////
	// Individual-level parameter book-keeping

	double Y_par[N_glob_par];
	double Ymu2_par[N_glob_par];
};


/////////////////////////////////////////////////////////////////
// 0.4 Individual-level structure for MCMC tuning

struct part_n_MCMC
{
	float par_vec[N_loc_par];                             // Parameter vector (in float format for setgmn) (lAB_0, rr)

	float par_vec_test[N_loc_par];                        // Test parameter vector for MCMC update (in float format for setgmn)
	float work[N_loc_par];                                // Scratch vector for setgmn

	double par_S1[N_loc_par];                             // Sum of parameters
	double par_S2[N_loc_par][N_loc_par];                  // Sum of product of pairs

	float COV_MAT[N_loc_par][N_loc_par];                  // covariance matrix (in float format for setgmn)
	float COV_MAT_dummy[N_loc_par][N_loc_par];            // dummy covariance matrix: setgmn gives back sqrt(COV_MAT) or similar so we feed it a dummy

	float GMN_parm[(N_loc_par)*(N_loc_par + 3) / 2 + 1];  // array for setgmn output

	int denom;                                            // denominator for tracking SD calculations

	double step_scale;                                    // scalar for tuning acceptance rate

	int accept;                                           // number of accepted steps
};


////////////////////////////////////////////////////
// 0.5 Initialise functions

double data_like_n(part_n* p, params* theta);
double mix_like_n(part_n* p, params* theta);
double global_prior(params* priors);
double local_prior(part_n* p);
double rm_scale(double step_scale, int step, int N_step_adapt, double log_prob);
double gammln(const double xx);


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//        //                                            //
//   ##   //  Initialise main object, read in data and  //
//  ###   //  fill out objects                          //
//   ##   //                                            //
//   ##   //                                            // 
//  ####  //                                            //
//        //                                            //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// 1.1 Initialise main object - need to choose whether
//     to run in console or as a .exe

int main(int argc, char** argv)
{

	// do we have the correct command line?
	if (argc != 5)
	{
		std::cout << "Incorrect command line.\n";
		return 0;
	}

	char* AB_input_File        = argv[1];
	char* global_output_File   = argv[2];
	char* local_output_File    = argv[3];
	double long_half           = atof(argv[4]);

	//////////////////////////////////////////////////////
	// 1.2 Declare seed, buffer for writing to and clock

	setall(time(NULL), 7);

	int cl = clock();


	///////////////////////////////////////////////////////
	// 1.3 Read in antibody data (infected individuals)

	std::ifstream AB_Stream(AB_input_File);

	if (AB_Stream.fail())
	{
		std::cout << "Failure reading in data." << endl;
	}

	vector<vector<double>> AB_data_read;

	AB_data_read.resize(N_part+ N_negcon);
	for (int i = 0; i < (N_part + N_negcon); i++)
	{
		AB_data_read[i].resize(N_data_cols);
	}

	for (int i = 0; i< (N_part + N_negcon); i++)
	{
		for (int j = 0; j<N_data_cols; j++)
		{
			AB_Stream >> AB_data_read[i][j];
		}
	}

	AB_Stream.close();


	//////////////////////////////////////////////////////
	// 1.4 Create global parameter objects

	params theta, theta_p1;


	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	// Priors on global parameter

	////////////////////////////
	// Ab_bg; background Ab level

	theta.prior_MM[0]    = 0.001;    // Mean value of the parameter in the population
	theta.prior_MM_CV[0] = 0.5;      // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[0]    = 0.005;    // Between Person standard deviation
	theta.prior_SIG_CV[0] = 0.33;     // Coefficient of variation in Between Person sd


	////////////////////////////
	// beta_mild; B cell boost

	theta.prior_MM[1]     = 0.001;     // Mean value of the parameter in the population
	theta.prior_MM_CV[1]  = 2;        // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[1]    = 0.01;     // Between Person standard deviation
	theta.prior_SIG_CV[1] = 0.5;      // Coefficient of variation in Between Person sd


	///////////////////////////
	// tau;  half-life of memory B cells

	theta.prior_MM[2]    = 9.6;           // Mean value of the parameter in the population
	theta.prior_MM_CV[2] = 0.25;          // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[2]    = 4.2;           // Between Person standard deviation
	theta.prior_SIG_CV[2] = 0.33;          // Coefficient of variation in Between Person sd


	///////////////////////////
	// t_delta;  half-life of memory B cells

	theta.prior_MM[3]    = 2.0;          // Mean value of the parameter in the population
	theta.prior_MM_CV[3] = 0.25;          // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[3]    = 1.8;           // Between Person standard deviation
	theta.prior_SIG_CV[3] = 0.33;          // Coefficient of variation in Between Person sd


    ///////////////////////////
	// t_short;  half-life of short-lived component

	theta.prior_MM[4]     = 2.5;          // Mean value of the parameter in the population
	theta.prior_MM_CV[4]  = 0.25;          // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[4]    = 1.3;           // Between Person standard deviation
	theta.prior_SIG_CV[4] = 0.33;          // Coefficient of variation in Between Person sd


	////////////////////////////
	// t_long;  half-life of long-lived component

	theta.prior_MM[5]     = long_half;        // Mean value of the parameter in the population
	theta.prior_MM_CV[5]  = 0.05;           // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[5]    = 100.0;        // Between Person standard deviation
	theta.prior_SIG_CV[5] = 0.33;          // Coefficient of variation in Between Person sd


	////////////////////////////
	// t_IgG;  half-life of IgG

	theta.prior_MM[6]     = 21.0;          // Mean value of the parameter in the population
	theta.prior_MM_CV[6]  = 0.02;          // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[6]    = 3.0;           // Between Person standard deviation
	theta.prior_SIG_CV[6] = 0.33;          // Coefficient of variation in Between Person sd


	////////////////////////////
	// rho;  proportion of short-lived component at start

	theta.prior_MM[7]     = 0.9;           // Mean value of the parameter in the population
	theta.prior_MM_CV[7]  = 0.075;           // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[7]    = 0.05;           // Between Person standard deviation
	theta.prior_SIG_CV[7] = 0.25;          // Coefficient of variation in Between Person sd
	

	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	// Transformation of priors 

	for (int p = 0; p < 7; p++)
	{
		theta.prior_LN_MM[p]     = log( theta.prior_MM[p] / sqrt(1.0 + pow(theta.prior_SIG[p] / theta.prior_MM[p], 2.0)) );   
		theta.prior_LN_MM_CV[p]  = theta.prior_MM_CV[p];

		theta.prior_LN_SIG[p]    = sqrt(log( 1.0 + pow(theta.prior_SIG[p] / theta.prior_MM[p], 2.0) ));          
		theta.prior_LN_SIG_CV[p] = theta.prior_SIG_CV[p];                                                       
	}


	theta.prior_LN_MM[7]    = 2.3046092;
	theta.prior_LN_MM_CV[7] = 2.5*theta.prior_MM_CV[7];

	theta.prior_LN_SIG[7]    = 0.5386774;
	theta.prior_LN_SIG_CV[7] = 2.5*theta.prior_SIG_CV[7];


	for (int p = 0; p < N_glob_par; p++)
	{
		theta.prior_mu[p]  = theta.prior_LN_MM[p];
		theta.prior_tau[p] = 1.0 / pow(theta.prior_LN_MM[p] * theta.prior_LN_MM_CV[p], 2.0);

		theta.prior_k[p]     = 1.0 / pow(2.0*theta.prior_LN_SIG_CV[p], 2.0);
		theta.prior_theta[p] = pow(2.0*theta.prior_LN_SIG_CV[p] / theta.prior_LN_SIG[p], 2.0);
	}

	for (int p = 0; p < N_glob_par; p++)
	{
		theta.mu_par[p]  = genunf(0.9, 1.1)*theta.prior_mu[p];
		theta.tau_par[p] = 0.1*genunf(0.9, 1.1)*theta.prior_k[p] * theta.prior_theta[p];
	}

	for (int p = 0; p < N_glob_par; p++)
	{
		cout << theta.mu_par[p] << "\t" << theta.tau_par[p] << endl;
	}
	cout << endl;


	theta.sig_obs = genunf(0.9, 1.1);      // Precision of observational error

	theta.log_sig_obs = log(theta.sig_obs);

	theta.sig_obs_scale = 1.0;


	theta.AB_min = 1.95e-6;
	theta.AB_max = 0.02;

	theta.lAB_min = log(theta.AB_min);
	theta.lAB_max = log(theta.AB_max);


	//////////////////////////////////////////////////////////
	// 1.5 Create individual-level objects for participant n

	part_n* part;
	part = new part_n[N_part + N_negcon];


	//////////////////////////////////////////////////////////
	// 1.6.1 Infected individuals

	for (int n = 0; n<(N_part+N_negcon); n++)
	{
		part[n].site     = AB_data_read[n][0];
		part[n].status   = AB_data_read[n][1];
		part[n].symptoms = AB_data_read[n][2];

		/////////////////////////////////////////////////
		// Fill antibody data

		part[n].N_sam = 0;

		for (int j = 0; j<N_t; j++)
		{
			if (AB_data_read[n][3 + j] > -0.5)
			{
				part[n].tt.push_back(AB_data_read[n][3 + j]);
				part[n].AB.push_back(AB_data_read[n][3 + N_t + j]);

				part[n].lAB.push_back(log(AB_data_read[n][3 + N_t + j]));

				part[n].N_sam = part[n].N_sam + 1;
			}
		}


		part[n].AB_min = 1e10;
		part[n].AB_max = 1e-10;


		for (int j = 0; j < part[n].N_sam; j++)
		{
			if (part[n].AB[j] < part[n].AB_min)
			{
				part[n].AB_min = part[n].AB[j];
			}
		}


		for (int j = 0; j < part[n].N_sam; j++)
		{
			if (part[n].AB[j] > part[n].AB_max)
			{
				part[n].AB_max = part[n].AB[j];
			}
		}

		//cout << n << "\t" << part[n].AB_min << "\t" << part[n].AB_max << endl;


		part[n].lAB_min = log(part[n].AB_min);
		part[n].lAB_max = log(part[n].AB_max);


		/////////////////////////////////////////////////
		// Randomly assign individual-level parameters

		if (part[n].status == 2)
		{
			part[n].Ab_bg   = genunf(2e-5, 3e-5);
			part[n].beta    = genunf(0.00001, 0.00005);
			part[n].tau     = genunf(0.1, 5.0);
			part[n].t_delta = genunf(5.0, 10.0);
			part[n].t_short = genunf(5.0, 10.0);
			part[n].t_long  = genunf(500.0, 1500.0);
			part[n].t_IgG   = genunf(10.0, 30.0);
			part[n].rho     = genunf(0.8, 0.9);

			part[n].lAb_bg   = log(part[n].Ab_bg);
			part[n].ltau     = log(part[n].tau);
			part[n].lbeta    = log(part[n].beta);
			part[n].lt_delta = log(part[n].t_delta);
			part[n].lt_short = log(part[n].t_short);
			part[n].lt_long  = log(part[n].t_long);
			part[n].lt_IgG   = log(part[n].t_IgG);

			part[n].r_delta = log2 / part[n].t_delta;
			part[n].r_short = log2 / part[n].t_short;
			part[n].r_long  = log2 / part[n].t_long;
			part[n].r_IgG   = log2 / part[n].t_IgG;

			part[n].logitrho = log(part[n].rho / (1.0 - part[n].rho));
		}


		if (part[n].status == 1)
		{
			part[n].Ab_bg   = genunf(2e-5, 3e-5);
			part[n].beta    = 0.01;
			part[n].tau     = 2.0;
			part[n].t_delta = 2.0;
			part[n].t_short = 10.0;
			part[n].t_long  = 100.0;
			part[n].t_IgG   = 10.0;
			part[n].rho     = 0.2;

			part[n].lAb_bg   = log(part[n].Ab_bg);
			part[n].ltau     = log(part[n].tau);
			part[n].lbeta    = log(part[n].beta);
			part[n].lt_delta = log(part[n].t_delta);
			part[n].lt_short = log(part[n].t_short);
			part[n].lt_long  = log(part[n].t_long);
			part[n].lt_IgG   = log(part[n].t_IgG);

			part[n].r_delta = log2 / part[n].t_delta;
			part[n].r_short = log2 / part[n].t_short;
			part[n].r_long  = log2 / part[n].t_long;
			part[n].r_IgG   = log2 / part[n].t_IgG;

			part[n].logitrho = log(part[n].rho / (1.0 - part[n].rho));
		}

		/////////////////////////////////////////////////
		// Calculate individual-level likelihood

		part[n].data_like = data_like_n(&part[n], &theta);

		part[n].mix_like  = mix_like_n(&part[n], &theta);

		part[n].lhood     = part[n].data_like + part[n].mix_like;
	
	}

	AB_data_read.clear();



	double NC_max = 0.0;


	for (int n = 0; n < (N_part + N_negcon); n++)
	{
		if (part[n].status == 1)
		{
			if (part[n].AB[0] > NC_max)
			{
				NC_max = part[n].AB[0];
			}
		}
	}


	for (int n = 0; n < (N_part + N_negcon); n++)
	{
		part[n].lNC_max = log(1.1*NC_max);
	}
	   


	//////////////////////////////////////////////////////
	// 1.7 Initialise adaptive MCMC object for individual-level parameters
	//     One object for each participant.

	part_n_MCMC* part_MCMC;
	part_MCMC = new part_n_MCMC[N_part + N_negcon];

	for (int n = 0; n<(N_part + N_negcon); n++)
	{
		///////////////////////////////////
		// Parameter vector for MVN update

		part_MCMC[n].par_vec[0] = part[n].lAb_bg;
		part_MCMC[n].par_vec[1] = part[n].lbeta;
		part_MCMC[n].par_vec[2] = part[n].ltau;
		part_MCMC[n].par_vec[3] = part[n].lt_delta;
		part_MCMC[n].par_vec[4] = part[n].lt_short;
		part_MCMC[n].par_vec[5] = part[n].lt_long;
		part_MCMC[n].par_vec[6] = part[n].lt_IgG;
		part_MCMC[n].par_vec[7] = part[n].logitrho;


		/////////////////////////////
		// Initialise diagonal covariance matrix

		for (int p = 0; p<N_loc_par; p++)
		{
			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].COV_MAT[p][q] = 0.0;
			}
		}

		part_MCMC[n].COV_MAT[0][0] = 0.2*0.2;
		part_MCMC[n].COV_MAT[1][1] = 0.2*0.2;
		part_MCMC[n].COV_MAT[2][2] = 0.2*0.2;
		part_MCMC[n].COV_MAT[3][3] = 0.2*0.2;
		part_MCMC[n].COV_MAT[4][4] = 0.2*0.2;
		part_MCMC[n].COV_MAT[5][5] = 0.2*0.2;
		part_MCMC[n].COV_MAT[6][6] = 0.2*0.2;
		part_MCMC[n].COV_MAT[7][7] = 0.2*0.2;

		/////////////////////////////
		// Counting moments

		for (int p = 0; p<N_loc_par; p++)
		{
			part_MCMC[n].par_S1[p] = part_MCMC[n].par_vec[p];

			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].par_S2[p][q] = part_MCMC[n].par_vec[p] * part_MCMC[n].par_vec[q];
			}
		}

		part_MCMC[n].denom = 1;


		/////////////////////////////
		// Set up dummy covariance matrix including
		// step-size scaling

		part_MCMC[n].step_scale = 0.1;

		for (int p = 0; p<N_loc_par; p++)
		{
			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].COV_MAT_dummy[p][q] = part_MCMC[n].step_scale*part_MCMC[n].COV_MAT[p][q];
			}
		}

		part_MCMC[n].accept = 0.0;
	}


	////////////////////////////////////////////////////////
	// 1.8 Book-keeping

	for (int p = 0; p < N_glob_par; p++)
	{
		theta.Y_par[p]    = 0.0;
		theta.Ymu2_par[p] = 0.0;

		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			theta.Y_par[p]    = theta.Y_par[p] + part_MCMC[n].par_vec[p];
			theta.Ymu2_par[p] = theta.Ymu2_par[p] + (part_MCMC[n].par_vec[p] - theta.mu_par[p])*(part_MCMC[n].par_vec[p] - theta.mu_par[p]);
		}
	}

	theta_p1 = theta;


	////////////////////////////////////////////////////////
	// 1.9 Create objects for updating local parameters

	part_n* part_p1;
	part_p1 = new part_n[N_part + N_negcon];

	for (int n = 0; n<(N_part + N_negcon); n++)
	{
		part_p1[n] = part[n];
	}



	///////////////////////////////////////////////////////////////////////////
	// 1.10 Test output of likelihood

	for (int n = 0; n < (N_part + N_negcon); n++)
	{
		cout << n << "\t" << "Ab_bg: " << part[n].lAb_bg << "\t" << "beta: " << part[n].lbeta << "\t"
			              << "tau: " << part[n].ltau << "\t" << "t_delta: " << part[n].lt_delta << "\t" << "t_short: " << part[n].lt_short << "\t" << "t_long: " << part[n].lt_long << "\t" << "t_IgG: " << part[n].lt_IgG << "\t"
			              << "rho: " << part[n].logitrho << "\t" << "logL: " << part[n].lhood << endl;

		//cout << "N_sam = " << part[n].N_sam << "\t" << part[n].AB[0] << "\t" << part[n].lAB[0] << endl;

		//system("PAUSE");
	}


	
	///////////////////////////////////////////////////////////////////////////
	// 1.11 Initialise parameters for MCMC likelihood, Robbins-Munro 
	//      acceptance and output

	double loglike = global_prior(&theta);
	for (int n = 0; n<(N_part + N_negcon); n++)
	{
		loglike = loglike + part[n].lhood;
	}

	double log_prob, loglike_p1;

	double log_loc_prob[N_part + N_negcon];

	int glob_out = max(2, (int)((int)N_mcmc) / 10000);
	int loc_out = max(2, (int)((int)N_mcmc) / 1000);


	vector<double> randomU(N_part + N_negcon);

	vector<double> loglike_vec_p1(N_part + N_negcon);


	///////////////////////////////////////////////////////////////////////////
	// 1.12 Open file for output and write first line

	std::cout << "START MCMC" << endl;

	cout << 0 << "\t";
	for (int p = 0; p < N_glob_par; p++)
	{
		cout << theta.mu_par[p] << "\t";
	}
	for (int p = 0; p < N_glob_par; p++)
	{
		cout << theta.tau_par[p] << "\t";
	}
	cout << theta.sig_obs << "\t" << loglike << "\t" << global_prior(&theta) << endl;


	std::ofstream global_MCMC_Stream(global_output_File);
	
	for (int p = 0; p < N_glob_par; p++)
	{
		global_MCMC_Stream << theta.mu_par[p] << "\t";
	}
	for (int p = 0; p < N_glob_par; p++)
	{
		global_MCMC_Stream << theta.tau_par[p] << "\t";
	}
	global_MCMC_Stream << theta.sig_obs << "\t" << loglike << "\t" << global_prior(&theta) << endl;


	std::ofstream local_MCMC_Stream(local_output_File);

	for (int n = 0; n<(N_part + N_negcon); n++)
	{
		local_MCMC_Stream << part[n].Ab_bg << "\t" << part[n].beta << "\t" <<
			                 part[n].tau << "\t" << part[n].t_delta << "\t" << part[n].t_short << "\t" << part[n].t_long << "\t" << part[n].t_IgG << "\t" <<
			                 part[n].rho << "\t" << part[n].lhood << "\t";
	}
	local_MCMC_Stream << endl;


	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	//          //                                         //
	//   ####   //  Begin MCMC fitting procedure           //
	//  ##  ##  //                                         //
	//     ##   //                                         //
	//    ##    //                                         //
	//   #####  //                                         //
	//          //                                         //
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////

	for (int mc = 1; mc<N_mcmc; mc++)
	{
		/////////////////////////////////////////////////
		/////////////////////////////////////////////////
		//       //                                    //
		//  2.1  //  UPDATE STAGE 1: INDIVIDUAL-LEVEL  //
		//       //  Metropolis-Hastings sampler       //
		//       //                                    //
		/////////////////////////////////////////////////
		/////////////////////////////////////////////////

		/////////////////////////////////////////////
		// 2.1.1. Proposal step

		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			////////////////////////////////////////////////
			// Update COV_MAT_dummay

			for (int p = 0; p < N_loc_par; p++)
			{
				for (int q = 0; q < N_loc_par; q++)
				{
					part_MCMC[n].COV_MAT_dummy[p][q] = part_MCMC[n].step_scale*part_MCMC[n].COV_MAT[p][q];
				}
			}


			///////////////////////////////////////////////
			// Multi-variate Normal proposal step

			setgmn(part_MCMC[n].par_vec, *part_MCMC[n].COV_MAT_dummy, N_loc_par, part_MCMC[n].GMN_parm);

			genmn(part_MCMC[n].GMN_parm, part_MCMC[n].par_vec_test, part_MCMC[n].work);

			part_p1[n].lAb_bg    = part_MCMC[n].par_vec_test[0];

			if( part[n].status == 2)
			{
				part_p1[n].lbeta    = part_MCMC[n].par_vec_test[1];
				part_p1[n].ltau     = part_MCMC[n].par_vec_test[2];
				part_p1[n].lt_delta = part_MCMC[n].par_vec_test[3];
				part_p1[n].lt_short = part_MCMC[n].par_vec_test[4];
				part_p1[n].lt_long  = part_MCMC[n].par_vec_test[5];
				part_p1[n].lt_IgG   = part_MCMC[n].par_vec_test[6];
				part_p1[n].logitrho = part_MCMC[n].par_vec_test[7];
			}

			randomU[n] = genunf(0.0, 1.0);
		}


		/////////////////////////////////////////////
		// 2.1.2. Update step
		
		//#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			////////////////////////////////////////////////////////
			// 2.1.2.1. Only proceed if allowable parameters proposed

			if (local_prior(&part_p1[n]) > -0.5*LARGE)
			{

				part_p1[n].Ab_bg = exp(part_p1[n].lAb_bg);

				if (part[n].status == 2)
				{
					part_p1[n].beta = exp(part_p1[n].lbeta);

					part_p1[n].tau = exp(part_p1[n].ltau);

					part_p1[n].t_delta = exp(part_p1[n].lt_delta);
					part_p1[n].r_delta = log2 / part_p1[n].t_delta;

					part_p1[n].t_short = exp(part_p1[n].lt_short);
					part_p1[n].r_short = log2 / part_p1[n].t_short;

					part_p1[n].t_long = exp(part_p1[n].lt_long);
					part_p1[n].r_long = log2 / part_p1[n].t_long;

					part_p1[n].t_IgG = exp(part_p1[n].lt_IgG);
					part_p1[n].r_IgG = log2 / part_p1[n].t_IgG;

					part_p1[n].rho = exp(part_p1[n].logitrho) / (1.0 + exp(part_p1[n].logitrho));
				}

				part_p1[n].data_like = data_like_n(&part_p1[n], &theta);

				part_p1[n].mix_like = mix_like_n(&part_p1[n], &theta);

				part_p1[n].lhood = part_p1[n].data_like + part_p1[n].mix_like;


				double log_prob_n = part_p1[n].lhood - part[n].lhood;

				log_loc_prob[n] = _finite(log_prob_n) ? std::min(log_prob_n, 0.0) : -LARGE;


				////////////////////////////////////////
				// 2.1.2.2. Update if necessary

				if (log(randomU[n]) < log_loc_prob[n])
				{
					part[n] = part_p1[n];

					for (int p = 0; p < N_loc_par; p++)
					{
						part_MCMC[n].par_vec[p] = part_MCMC[n].par_vec_test[p];
					}

					part_MCMC[n].accept = part_MCMC[n].accept + 1;
				}


				////////////////////////////////////////
				// 2.1.3. Adjust step-size with Robbins-Monro
				//        Only do this for a local step within allowed range

				if (mc < N_adapt)
				{
					part_MCMC[n].step_scale = rm_scale(part_MCMC[n].step_scale, mc, N_adapt, log_loc_prob[n]);
				}
			}


			////////////////////////////////////////////////////////////
			// Running account of sums and sums of squares

			if (mc < N_tune_end)
			{
				for (int p = 0; p < N_loc_par; p++)
				{
					part_MCMC[n].par_S1[p] = part_MCMC[n].par_S1[p] + part_MCMC[n].par_vec[p];
					
					for (int q = 0; q < N_loc_par; q++)
					{
						part_MCMC[n].par_S2[p][q] = part_MCMC[n].par_S2[p][q] + part_MCMC[n].par_vec[p] * part_MCMC[n].par_vec[q];
					}
				}

				part_MCMC[n].denom = part_MCMC[n].denom + 1;
			}


			////////////////////////////////////////////////////////////
			// 2.1.4. TUNING STAGE 1

			/////////////////////////////////
			// Update covariance matrix

			if ((mc >= N_tune_start) && (mc < N_tune_end))
			{
				for (int p = 0; p < N_loc_par; p++)
				{
					for (int q = 0; q < N_loc_par; q++)
					{
						if (part_MCMC[n].accept / part_MCMC[n].denom > 0.01)
						{
							part_MCMC[n].COV_MAT[p][q] = part_MCMC[n].par_S2[p][q] / (part_MCMC[n].denom) - part_MCMC[n].par_S1[p] * part_MCMC[n].par_S1[q] / (part_MCMC[n].denom*part_MCMC[n].denom);
						}
					}
				}
			}
		}

		//////////////////////////////////////////////////////////////
		// 2.1.5. Update the total likelihood given the local updates

		loglike = global_prior(&theta);
		
		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			loglike = loglike + part[n].lhood;
		}


		///////////////////////////////////////////////////
		///////////////////////////////////////////////////
		//       //                                      //
		//  2.2  //  UPDATE STAGE 2 - POPULATION-LEVEL   //
		//       //  Gibbs sampler                       // 
		//       //                                      //
		///////////////////////////////////////////////////
		///////////////////////////////////////////////////


		///////////////////////////////////////////////////
		// Ab_bg 

		theta.Y_par[0] = 0.0;

		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			theta.Y_par[0] = theta.Y_par[0] + part_MCMC[n].par_vec[0];
		}

		theta.mu_par[0] = gennor( (theta.prior_tau[0] * theta.prior_mu[0] + theta.tau_par[0] * theta.Y_par[0]) / (theta.prior_tau[0] + (N_part + N_negcon)* theta.tau_par[0]),
			                      1.0 / sqrt(theta.prior_tau[0] + (N_part + N_negcon)* theta.tau_par[0]));

		theta.Ymu2_par[0] = 0.0;

		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			theta.Ymu2_par[0] = theta.Ymu2_par[0] + (part_MCMC[n].par_vec[0] - theta.mu_par[0])*(part_MCMC[n].par_vec[0] - theta.mu_par[0]);
		}

		theta.tau_par[0] = gengam( 1.0 / theta.prior_theta[0] + 0.5*theta.Ymu2_par[0],
			                       0.5*(N_part + N_negcon) + theta.prior_k[0]);


		///////////////////////////////////////////////////
		// beta

		theta.Y_par[1] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Y_par[1] = theta.Y_par[1] + part_MCMC[n].par_vec[1];
		}

		theta.mu_par[1] = gennor( (theta.prior_tau[1] * theta.prior_mu[1] + theta.tau_par[1] * theta.Y_par[1]) / (theta.prior_tau[1] + N_part*theta.tau_par[1]),
				                  1.0 / sqrt(theta.prior_tau[1] + N_part*theta.tau_par[1]));

		theta.Ymu2_par[1] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Ymu2_par[1] = theta.Ymu2_par[1] + (part_MCMC[n].par_vec[1] - theta.mu_par[1])*(part_MCMC[n].par_vec[1] - theta.mu_par[1]);
		}

		theta.tau_par[1] = gengam( 1.0 / theta.prior_theta[1] + 0.5*theta.Ymu2_par[1],
				                   0.5*N_part + theta.prior_k[1] );


		///////////////////////////////////////////////////
		// tau 

		theta.Y_par[2] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Y_par[2] = theta.Y_par[2] + part_MCMC[n].par_vec[2];
		}

		theta.mu_par[2] = gennor( (theta.prior_tau[2] * theta.prior_mu[2] + theta.tau_par[2] * theta.Y_par[2]) / (theta.prior_tau[2] + N_part * theta.tau_par[2]),
			                      1.0 / sqrt(theta.prior_tau[2] + N_part * theta.tau_par[2]));

		theta.Ymu2_par[2] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Ymu2_par[2] = theta.Ymu2_par[2] + (part_MCMC[n].par_vec[2] - theta.mu_par[2])*(part_MCMC[n].par_vec[2] - theta.mu_par[2]);
		}

		theta.tau_par[2] = gengam( 1.0 / theta.prior_theta[2] + 0.5*theta.Ymu2_par[2],
			                       0.5*N_part + theta.prior_k[2]);


		///////////////////////////////////////////////////
		// delta 

		theta.Y_par[3] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Y_par[3] = theta.Y_par[3] + part_MCMC[n].par_vec[3];
		}

		theta.mu_par[3] = gennor( (theta.prior_tau[3] * theta.prior_mu[3] + theta.tau_par[3] * theta.Y_par[3]) / (theta.prior_tau[3] + N_part * theta.tau_par[3]),
			                      1.0 / sqrt(theta.prior_tau[3] + N_part * theta.tau_par[3]));

		theta.Ymu2_par[3] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Ymu2_par[3] = theta.Ymu2_par[3] + (part_MCMC[n].par_vec[3] - theta.mu_par[3])*(part_MCMC[n].par_vec[3] - theta.mu_par[3]);
		}

		theta.tau_par[3] = gengam( 1.0 / theta.prior_theta[3] + 0.5*theta.Ymu2_par[3],
			                       0.5*N_part + theta.prior_k[3]);


		///////////////////////////////////////////////////
		// t_short

		theta.Y_par[4] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Y_par[4] = theta.Y_par[4] + part_MCMC[n].par_vec[4];
		}

		theta.mu_par[4] = gennor( (theta.prior_tau[4] * theta.prior_mu[4] + theta.tau_par[4] * theta.Y_par[4]) / (theta.prior_tau[4] + N_part*theta.tau_par[4]),
				                   1.0 / sqrt(theta.prior_tau[4] + N_part *theta.tau_par[4]) );

		theta.Ymu2_par[4] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Ymu2_par[4] = theta.Ymu2_par[4] + (part_MCMC[n].par_vec[4] - theta.mu_par[4])*(part_MCMC[n].par_vec[4] - theta.mu_par[4]);
		}

		theta.tau_par[4] = gengam( 1.0/theta.prior_theta[4] + 0.5*theta.Ymu2_par[4],
			                        0.5*N_part + theta.prior_k[4]);
		

		///////////////////////////////////////////////////
		// t_long

		theta.Y_par[5] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Y_par[5] = theta.Y_par[5] + part_MCMC[n].par_vec[5];
		}

		theta.mu_par[5] = gennor( (theta.prior_tau[5] * theta.prior_mu[5] + theta.tau_par[5] * theta.Y_par[5]) / (theta.prior_tau[5] + N_part*theta.tau_par[5]),
			                       1.0 / sqrt(theta.prior_tau[5] + N_part*theta.tau_par[5]) );

		theta.Ymu2_par[5] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Ymu2_par[5] = theta.Ymu2_par[5] + (part_MCMC[n].par_vec[5] - theta.mu_par[5])*(part_MCMC[n].par_vec[5] - theta.mu_par[5]);
		}

		theta.tau_par[5] = gengam( 1.0 / theta.prior_theta[5] + 0.5*theta.Ymu2_par[5],
			                        0.5*N_part + theta.prior_k[5] );


		///////////////////////////////////////////////////
		// t_IgG

		theta.Y_par[6] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Y_par[6] = theta.Y_par[6] + part_MCMC[n].par_vec[6];
		}

		theta.mu_par[6] = gennor( (theta.prior_tau[6] * theta.prior_mu[6] + theta.tau_par[6] * theta.Y_par[6]) / (theta.prior_tau[6] + N_part*theta.tau_par[6]),
			                       1.0 / sqrt(theta.prior_tau[6] + N_part *theta.tau_par[6]));

		theta.Ymu2_par[6] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Ymu2_par[6] = theta.Ymu2_par[6] + (part_MCMC[n].par_vec[6] - theta.mu_par[6])*(part_MCMC[n].par_vec[6] - theta.mu_par[6]);
		}


		theta.tau_par[6] = gengam( 1.0 / theta.prior_theta[6] + 0.5*theta.Ymu2_par[6],
			                        0.5*N_part + theta.prior_k[6]);


		///////////////////////////////////////////////////
		// rho

		theta.Y_par[7] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Y_par[7] = theta.Y_par[7] + part_MCMC[n].par_vec[7];
		}

		theta.mu_par[7] = gennor( (theta.prior_tau[7] * theta.prior_mu[7] + theta.tau_par[7] * theta.Y_par[7]) / (theta.prior_tau[7] + N_part*theta.tau_par[7]),
			                       1.0 / sqrt(theta.prior_tau[7] + N_part*theta.tau_par[7]) );

		theta.Ymu2_par[7] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Ymu2_par[7] = theta.Ymu2_par[7] + (part_MCMC[n].par_vec[7] - theta.mu_par[7])*(part_MCMC[n].par_vec[7] - theta.mu_par[7]);
		}

		theta.tau_par[7] = gengam( 1.0 / theta.prior_theta[7] + 0.5*theta.Ymu2_par[7],
			                       0.5*N_part + theta.prior_k[7] );


		theta_p1 = theta;


		//////////////////////////////////////////////////////////////
		// 2.1.5. Update the total likelihood given the local updates

		loglike = global_prior(&theta);

		#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			//part[n].data_like = data_like_n(&part[n], &theta);

			part[n].mix_like = mix_like_n(&part[n], &theta);

			part[n].lhood = part[n].data_like + part[n].mix_like;

			loglike = loglike + part[n].lhood;
		}

		

		///////////////////////////////////////////////////
		///////////////////////////////////////////////////
		//       //                                      //
		//  2.3  //  UPDATE STAGE 3 - OBS-VARIANCE       //
		//       //  Metropolis-Hastings update          // 
		//       //                                      //
		///////////////////////////////////////////////////
		///////////////////////////////////////////////////

		///////////////////////////////////////////////
		// 2.3.1 Propose new parameters

		#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			part_p1[n] = part[n];
		}

		///////////////////////////////////////////////
		// Normal proposal step

		theta_p1.sig_obs = gennor(theta.sig_obs, 0.1*theta.sig_obs_scale);
	
		//////////////////////////////////////////////
		// 2.3.2 If parameters are acceptable, attempt update

		if (global_prior(&theta_p1) > -0.5*LARGE)
		{
			theta_p1.log_sig_obs = log(theta_p1.sig_obs);

			#pragma omp parallel for schedule(dynamic,4)
			for (int n = 0; n<(N_part + N_negcon); n++)
			{
				part_p1[n].data_like = data_like_n(&part_p1[n], &theta_p1);
				
				part_p1[n].mix_like = mix_like_n(&part_p1[n], &theta_p1);

				part_p1[n].lhood = part_p1[n].data_like + part_p1[n].mix_like;

				loglike_vec_p1[n] = part_p1[n].lhood;
			}


			loglike_p1 = global_prior(&theta_p1);
			for (int n = 0; n<(N_part + N_negcon); n++)
			{
				loglike_p1 = loglike_p1 + loglike_vec_p1[n];
			}

			const double log_prob0 = loglike_p1 - loglike;
			log_prob = _finite(log_prob0) ? std::min(log_prob0, 0.0) : -LARGE;

			if (log(genunf(0, 1)) < log_prob)
			{
				loglike = loglike_p1;
				theta = theta_p1;

				for (int n = 0; n<(N_part + N_negcon); n++)
				{
					//part[n].lhood = loglike_vec_p1[n];

					part[n] = part_p1[n];
				}
			}


			////////////////////////////////////////
			// 2.3.4 Adjust step-size with Robbins-Monro

			if (mc < N_adapt)
			{
				theta.sig_obs_scale = rm_scale(theta.sig_obs_scale, mc + 1, N_adapt, log_prob);
			}

		}


		//////////////////////////////////////////////////////////////
		// 2.1.5. Update the total likelihood given the local updates

		loglike = global_prior(&theta);

		for (int n = 0; n < (N_part + N_negcon); n++)
		{
			loglike = loglike + part[n].lhood;
		}


		/////////////////////////////////////////
		/////////////////////////////////////////
		//       //                            //
		//  2.3  //  Output results to file    //
		//       //                            //
		/////////////////////////////////////////
		/////////////////////////////////////////

		if (mc%(glob_out*10) == 0)
		{
			cout << 100 * ((double)mc) / ((double)N_mcmc) << "% complete." << "\t";
			for (int p = 0; p < N_glob_par; p++)
			{
				cout << theta.mu_par[p] << "\t";
			}
			for (int p = 0; p < N_glob_par; p++)
			{
				cout << theta.tau_par[p] << "\t";
			}
			cout << theta.sig_obs << "\t" << loglike << "\t" << global_prior(&theta) << endl;
		}

		if (mc%glob_out == 0)
		{
			for (int p = 0; p < N_glob_par; p++)
			{
				global_MCMC_Stream << theta.mu_par[p] << "\t";
			}
			for (int p = 0; p < N_glob_par; p++)
			{
				global_MCMC_Stream << theta.tau_par[p] << "\t";
			}
			global_MCMC_Stream << theta.sig_obs << "\t" << loglike << "\t" << global_prior(&theta) << endl;
		}



		if (mc%loc_out == 0)
		{
			for (int n = 0; n<(N_part + N_negcon); n++)
			{
				local_MCMC_Stream << part[n].Ab_bg << "\t" << part[n].beta << "\t" <<
					                 part[n].tau << "\t" << part[n].t_delta << "\t" << part[n].t_short << "\t" << part[n].t_long << "\t" << part[n].t_IgG << "\t" <<
					                 part[n].rho << "\t" << part[n].lhood << "\t";
			}
			local_MCMC_Stream << endl;
		}

	}


	//////////////////////////////////////
	// 2.4. Calculate and output acceptance rates

	cout << "local acceptance:  " << part_MCMC[0].accept << "\t" << (double(part_MCMC[0].accept)) / (double(N_mcmc)) << endl;

	cout << "Time taken: " << cl << endl;

	return 0;
}



///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//          //                                       //
//   ####   //  Functions                            //
//  ##  ##  //                                       //
//     ##   //                                       //
//  ##  ##  //                                       //
//   ####   //                                       //
//          //                                       // 
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// 3.1 Individual-level likelihood


double data_like_n(part_n* p, params* theta)
{
	///////////////////////////////////////
	// Model-predicted antibody titre

	vector<double> model_AB(p->N_sam), model_lAB(p->N_sam);

	if (p->status == 1)
	{
		for (int j = 0; j < p->N_sam; j++)
		{
			model_AB[j] = p->Ab_bg;
		}
	}

	if (p->status == 2)
	{
		double part_1 = p->rho / ((p->r_short - p->r_IgG)*(p->r_short - p->r_delta));
		double part_2 = (1.0 - p->rho) / ((p->r_long - p->r_IgG)*(p->r_long - p->r_delta));
		double part_3 = (p->r_IgG + p->r_short*(p->rho - 1.0) - p->r_long*p->rho) / ((p->r_short - p->r_IgG)*(p->r_long - p->r_IgG)*(p->r_IgG - p->r_delta));
		double part_4 = -(p->r_delta + p->r_short*(p->rho - 1.0) - p->r_long*p->rho) / ((p->r_short - p->r_delta)*(p->r_long - p->r_delta)*(p->r_IgG - p->r_delta));


		for (int j = 0; j < p->N_sam; j++)
		{
			model_AB[j] = p->Ab_bg;

			if (p->tt[j] > p->tau)
			{
				model_AB[j] = model_AB[j] +
					          p->beta*( part_1*exp(-p->r_short * (p->tt[j] - p->tau)) + part_2*exp(-p->r_long * (p->tt[j] - p->tau)) +
						                part_3*exp(-p->r_IgG * (p->tt[j] - p->tau))   + part_4*exp(-p->r_delta * (p->tt[j] - p->tau)) );
			}
		}
	}


	for (int j = 0; j < p->N_sam; j++)
	{
		model_lAB[j] = log(model_AB[j]);
	}


	double loglike = 0.0;

	for (int j = 0; j < p->N_sam; j++)
	{
		/////////////////////////////////////////////////////////////
		// CASE 1: Data lies between assay minimum and maximum

		if ((p->AB[j] > theta->AB_min) && (p->AB[j] < theta->AB_max))
		{
			loglike = loglike - 0.9189385 - theta->log_sig_obs - p->lAB[j] - 0.5*(p->lAB[j] - model_lAB[j])*(p->lAB[j] - model_lAB[j]) / (theta->sig_obs*theta->sig_obs)
				- log(0.5*erf((theta->lAB_max - model_lAB[j]) / (sqrt2*theta->sig_obs)) - 0.5*erf((theta->lAB_min - model_lAB[j]) / (sqrt2*theta->sig_obs)));
			
			//cout << "AB = " << p->AB[j] << "    model = " << model_AB[j] << endl;

			//cout << -0.9189385 - theta->log_sig_obs - p->lAB[j] - 0.5*(p->lAB[j] - model_lAB[j])*(p->lAB[j] - model_lAB[j]) / (theta->sig_obs*theta->sig_obs) << endl;

			//cout << - log(0.5*erf((theta->lAB_max - model_lAB[j]) / (sqrt2*theta->sig_obs)) - 0.5*erf((theta->lAB_min - model_lAB[j]) / (sqrt2*theta->sig_obs))) << endl;

			//cout << 0.5*erf( (theta->lAB_max - model_lAB[j]) / (sqrt2*theta->sig_obs) )  << endl;

			//cout << 0.5*erf((theta->lAB_min - model_lAB[j]) / (sqrt2*theta->sig_obs)) << endl;

			//cout << "CASE 1: " << loglike << endl;
			//cout << endl;
			
		}

		/////////////////////////////////////////////////////////////
		// CASE 2: Data equals assay maximum

		if (p->AB[j] > 0.999999*theta->AB_max)
		{
			loglike = loglike + log(0.5 - 0.5*erf((p->lAB[j] - model_lAB[j]) / (sqrt2*theta->sig_obs)));

			//cout << "AB = " << p->AB[j] << "    model = " << model_AB[j] << endl;

			//cout << "CASE 2: " << loglike << endl;
			//cout << endl;
		}

		/////////////////////////////////////////////////////////////
		// CASE 3: Data equals assay minimum

		if (p->AB[j] < 1.000001*theta->AB_min)
		{
			loglike = loglike + log(0.5 + 0.5*erf((p->lAB[j] - model_lAB[j]) / (sqrt2*theta->sig_obs)));

			//cout << "AB = " << p->AB[j] << "    model = " << model_AB[j] << endl;

			//cout << "CASE 3: " << loglike << endl;
			//cout << endl;
		}
	}


	return loglike;
}



double mix_like_n(part_n* p, params* theta)
{
	double mixlike = 0.0;

	//////////////////////////////////////////
	// Mixed-effects component of likelihood

	if (p->status == 2)
	{
		mixlike = -6.43257 + 0.5*log(theta->tau_par[0] * theta->tau_par[1] * theta->tau_par[2] * theta->tau_par[3] * theta->tau_par[4] * theta->tau_par[5] * theta->tau_par[6] * theta->tau_par[7])
			                - 0.5*theta->tau_par[0] * (p->lAb_bg   - theta->mu_par[0]) * (p->lAb_bg   - theta->mu_par[0])
			                - 0.5*theta->tau_par[1] * (p->lbeta    - theta->mu_par[1]) * (p->lbeta    - theta->mu_par[1])
			                - 0.5*theta->tau_par[2] * (p->ltau     - theta->mu_par[2]) * (p->ltau     - theta->mu_par[2])
			                - 0.5*theta->tau_par[3] * (p->lt_delta - theta->mu_par[3]) * (p->lt_delta - theta->mu_par[3])
			                - 0.5*theta->tau_par[4] * (p->lt_short - theta->mu_par[4]) * (p->lt_short - theta->mu_par[4])
			                - 0.5*theta->tau_par[5] * (p->lt_long  - theta->mu_par[5]) * (p->lt_long  - theta->mu_par[5])
		 	                - 0.5*theta->tau_par[6] * (p->lt_IgG   - theta->mu_par[6]) * (p->lt_IgG   - theta->mu_par[6])
			                - 0.5*theta->tau_par[7] * (p->logitrho - theta->mu_par[7]) * (p->logitrho - theta->mu_par[7]);
	}
	
	if (p->status == 1)
	{
		mixlike = - 0.9189385 + 0.5*log(theta->tau_par[0])
			                  - 0.5*theta->tau_par[0] * (p->lAb_bg - theta->mu_par[0]) * (p->lAb_bg - theta->mu_par[0]);
	}
	
	return mixlike;
}



///////////////////////////////////////////////////////
// 3.2 Individual-level prior likelihood: excludes ineligible steps

double local_prior(part_n* p)
{
	if (p->lAb_bg < -13.14768) { return -LARGE; } // log(0.1*1.95e-5) = -13.14768 
	if (p->lAb_bg > p->lNC_max) { return -LARGE; } // Greater than the maximum of the neg cons

	if( p->status == 2)
	{
		//if (p->ltau < 0.0) { return -LARGE; }
		if (p->ltau > 2.995732) { return -LARGE; }  // log(20) = 2.995732

		if (p->lt_short < 0.0) { return -LARGE; }

		if (p->lt_short > p->lt_long) { return -LARGE; }

		if (p->lt_IgG > 4.094345) { return -LARGE; }  // log(60) = 4.094345

		//if (p->rho < 0.0) { return -LARGE; }
		//if (p->rho > 1.0) { return -LARGE; }

		//if (p->logitrho < -1.386294) { return -LARGE; }   // log( 0.2/(1-0.2) ) = - 1.386294
	}

	return 0.0;
}


///////////////////////////////////////////////////////
// 3.3 Population-level prior likelihood: excludes ineligible steps

double global_prior(params* priors)
{
	double logprior = 0.0;

	for (int p = 0; p < N_glob_par; p++)
	{
		///////////////////////////
		// Normal prior on mean_par[p] 

		logprior = logprior - 0.9189385 + 0.5*log(priors->prior_tau[p]) - 0.5*priors->prior_tau[p] * (priors->mu_par[p] - priors->prior_mu[p])*(priors->mu_par[p] - priors->prior_mu[p]);


		///////////////////////////
		// Gamma prior on tau_par[p]

		logprior = logprior + (priors->prior_k[p] - 1.0)*log(priors->tau_par[p]) - priors->tau_par[p] / priors->prior_theta[p] - priors->prior_k[p] *log(priors->prior_theta[p]) - gammln(priors->prior_k[p]);
	}


	///////////////////////////
	// Gamma prior on tau_obs

	if (priors->sig_obs < 0.0) { return -LARGE; }
	if (priors->sig_obs > 2.0) { return -LARGE; }

	double mean_sig_obs = 1.0;// 0.01;
	double sd_sig_obs = 1.0;// 0.02;
	double k_sig_obs = (mean_sig_obs / sd_sig_obs)*(mean_sig_obs / sd_sig_obs);
	double theta_sig_obs = sd_sig_obs*sd_sig_obs / mean_sig_obs;

	double log_prior_sig_obs = (k_sig_obs - 1)*log(priors->sig_obs) - priors->sig_obs / theta_sig_obs - k_sig_obs*log(theta_sig_obs) - gammln(k_sig_obs);

//	double log_prior_sig_obs = (priors->prior_k_tau_obs - 1.0)*log(priors->tau_obs) - priors->tau_obs / priors->prior_theta_tau_obs - priors->prior_k_tau_obs*log(priors->prior_theta_tau_obs) - gammln(priors->prior_k_tau_obs);

	return logprior + log_prior_sig_obs;
}


/////////////////////////////////////////////////////////
// 3.4 Robbins-Monro algorithm for setting acceptance rate.
//
// Robbins-Munroe stochastic approximation adjustment
//	return adjustment
//	i iteration number
//	a scaling factor for adjustment size
//	m number of iterations when adjustment halves from value at i=0
//	d outcome at iteration i (e.g. calculated acceptance probability (min of MH ratio and 1),
//								or acceptance or no acceptance)
//	p desired mean for d (e.g. desired acceptance probability)

double rm_scale(double step_scale, int step, int N_step_adapt, double log_prob)
{
	double dd = exp(log_prob);
	if (dd < -30) { dd = 0.0; }
	dd = std::min(dd, 1.0);

	double rm_temp = (dd - 0.23) / ((double(step) + 1) / (0.01*(double(N_step_adapt)) + 1));

	double out = step_scale*exp(rm_temp);

	out = std::max(out, 0.02);
	out = std::min(out, 5.0);

	return out;
}


/////////////////////////////////////////////////////////
// 3.6 Log gamma function, based on gamma.h from NRC3

double gammln(const double xx) {
	int j;
	double x, tmp, y, ser;
	static const double cof[14] = { 57.1562356658629235, -59.5979603554754912,
		14.1360979747417471, -0.491913816097620199, 0.339946499848118887e-4,
		0.465236289270485756e-4, -0.983744753048795646e-4, 0.158088703224912494e-3,
		-0.210264441724104883e-3, 0.217439618115212643e-3, -0.164318106536763890e-3,
		0.844182239838527433e-4, -0.261908384015814087e-4, 0.368991826595316234e-5 };
	if (xx <= 0) throw("bad arg in gammln");
	y = x = xx;
	tmp = x + 5.242187500000000;
	tmp = (x + 0.5)*log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j<14; j++) ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005*ser / x);
}

