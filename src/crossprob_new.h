//
//  crossprob.h
//  
//
//  Created by Raymond on 2017/8/9.
//
//

#ifndef crossprob_h
#define crossprob_h

#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <numeric>

#include "stdio.h"
#include <iterator>
#include <fstream>
#include <string>
#include <math.h>
#include <set>
#include <cmath>
#include <random>


#include "one_sided_noncrossing_probability.h"
#include "one_sided_noncrossing_probability_n2logn.h"
#include "one_sided_noncrossing_probability_n2.h"
#include "two_sided_noncrossing_probability.h"
#include "read_boundaries_file.h"
#include "string_utils.h"

using namespace std;

//static void print_usage();

vector <double> TheoreticalDF (vector <double> pmf);

vector <double> EmpiricalDF (vector <double> obs);

vector <double> BinomialPF(int trial, double prob);

vector <double> BinomialDF(int trial, double prob);

vector <double> MixDF (vector <double> obs);

vector <double> LeftMixDF (vector <double> points, double eps);

vector <double> LeftDF (vector <double> obs, double eps);

double Mixed_KS_stat (vector <double> obs);

int upperlowerbounds_discontinuous(vector <double> points, int sample_size, double critical, double eps);

int upperlowerbounds(double n, double x);

double ecdf_noncrossing_probability_loop(int n, const vector<double>& g_steps, const vector<double>& h_steps, bool use_fft, int loop);

//static int cont_ks_distribution(long n);
double cont_ks_distribution(long n);

void upperlowerbounds_alternative(double n, double x);







//int integer_crossing_times();
//
//static int handle_command_line_arguments(int argc, char* argv[]);
//
//int main(int argc, char* argv[]);


#endif /* crossprob_h */
