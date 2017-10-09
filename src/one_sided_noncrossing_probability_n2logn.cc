#include <vector>
#include <iostream>
#include <stdexcept>

#include "one_sided_noncrossing_probability_n2logn.h"
#include "common.h"
#include "poisson_pmf.h"
#include "fftwconvolver.h"
#include "aligned_mem.h"
#include "string_utils.h"

using namespace std;


vector<double> poisson_lower_noncrossing_probability_n2logn(int n, double intensity, const vector<double>& lower_bound_steps)
{
    DoubleBuffer<double> buffers(n+1, 0.0);
    buffers.get_src()[0] = 1.0;

    FFTWConvolver fftconvolver(n+1);

    double* pmf = allocate_aligned_doubles(n+1);
    PoissonPMFGenerator pmfgen(n+1);

    int n_steps = lower_bound_steps.size();
    double prev_step_location = 0.0;

    for (int i = 0; i < n_steps; ++i) {
        pmfgen.compute_pmf(n-i+1, intensity*(lower_bound_steps[i]-prev_step_location), pmf);
        fftconvolver.convolve_same_size(n-i+1, pmf, &buffers.get_src()[i], &buffers.get_dest()[i]);
        //cout << "==== i: " << i << " ==================================\n";
        //cout << "src: " << vector_to_string(buffers.get_src());
        //vector<double> tmp(pmf, pmf+n+1);
        //cout << "pmf: " << vector_to_string(tmp);
        //cout << "dest: " << vector_to_string(buffers.get_dest());

        buffers.get_dest()[i] = 0.0;
        buffers.get_src()[i] = 0.0; // Not strictly necessary. This just keeps the arrays cleaner when printing.
        prev_step_location = lower_bound_steps[i];
        buffers.flip();
    }
    pmfgen.compute_pmf(n-n_steps+1, intensity*(1.0-prev_step_location), pmf);
    fftconvolver.convolve_same_size(n+1-n_steps, pmf, &buffers.get_src()[n_steps], &buffers.get_dest()[n_steps]);

    free_aligned_mem(pmf);
    return buffers.get_dest();
}

double ecdf_lower_noncrossing_probability_n2logn(int n, const vector<double>& lower_bound_steps)
{
    if ((int)lower_bound_steps.size() > n) {
        stringstream ss;
        ss << "Empirical CDF must cross lower boundary g(t) since g(1)==" << lower_bound_steps.size() << " > n and the number of samples is n. h_steps:\n";
        throw runtime_error(ss.str() + vector_to_string(lower_bound_steps));
    }

    vector<double> poisson_nocross_probabilities = poisson_lower_noncrossing_probability_n2logn(n, n, lower_bound_steps);
    return poisson_nocross_probabilities[n] / poisson_pmf(n, n);
}

double ecdf_upper_noncrossing_probability_n2logn(int n, const vector<double>& upper_bound_steps)
{
    if ((int)upper_bound_steps.size() < n) {
        stringstream ss;
        ss << "empirical CDF must cross upper boundary h(t) since h(1)==" << upper_bound_steps.size() << " < n and the number of samples is n.";
        throw runtime_error(ss.str());
    }

    vector<double> symmetric_steps(n, 0.0);
    for (int i = n-upper_bound_steps.size(); i < n; ++i) {
        symmetric_steps[i] = 1.0 - upper_bound_steps[upper_bound_steps.size() - 1 - i];
    }

    return ecdf_lower_noncrossing_probability_n2logn(n, symmetric_steps);
}
