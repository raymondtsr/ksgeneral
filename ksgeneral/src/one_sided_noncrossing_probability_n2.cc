#include <vector>
#include <iostream>
#include <stdexcept>

#include "one_sided_noncrossing_probability_n2.h"
#include "common.h"
#include "poisson_pmf.h"
#include "fftwconvolver.h"
#include "aligned_mem.h"
#include "string_utils.h"

using namespace std;

vector<double> poisson_lower_noncrossing_probability_n2(int n, double intensity, const vector<double>& lower_bound_steps, int jump_size)
{
    assert(jump_size <= n);
    DoubleBuffer<double> buffers(n+1, 0.0);
    DoubleBuffer<double> minibuffers(jump_size, 0.0);
    buffers.get_src()[0] = 1.0;

    FFTWConvolver fftconvolver(n+1);

    double* pmf = allocate_aligned_doubles(n+1);
    PoissonPMFGenerator pmfgen(n+1);

    double* tmp = allocate_aligned_doubles(n+1);

    int n_steps = lower_bound_steps.size();
    double I_prev_location = 0.0;
    int I_prev = -1;
    int I = I_prev+jump_size;
    while (true) {
        // cout << "I: " << I << " I_prev: " << I_prev << endl;
        pmfgen.compute_pmf(n-I_prev, intensity*(lower_bound_steps[I]-I_prev_location), pmf);
        fftconvolver.convolve_same_size(n-I_prev, pmf, &buffers.get_src()[I_prev+1], tmp);
        fill(&buffers.get_dest()[0], &buffers.get_dest()[I+1], 0.0);
        copy(&tmp[I-I_prev], &tmp[n-I_prev], &buffers.get_dest()[I+1]);

        copy(&buffers.get_src()[I_prev+1], &buffers.get_src()[I+1], &minibuffers.get_src()[0]);
        double i_prev_location = I_prev_location;
        for (int i = I_prev+1; i < I; i++) {
            pmfgen.compute_pmf(I-i+1, intensity*(lower_bound_steps[i]-i_prev_location), pmf);
            fftconvolver.convolve_same_size(I-i+1, pmf, &minibuffers.get_src()[i-I_prev-1], &minibuffers.get_dest()[i-I_prev-1]);

            double prob_exit_now = minibuffers.get_dest()[i-I_prev-1];
            double lambda = intensity*(lower_bound_steps[I]- lower_bound_steps[i]);
            for (int j = I+1; j < n+1; ++j) {
                buffers.get_dest()[j] -= prob_exit_now * pmfgen.evaluate_pmf(lambda, j-i);
            }

            minibuffers.get_dest()[i-I_prev-1] = 0.0;
            minibuffers.get_src()[i-I_prev-1] = 0.0;
            minibuffers.flip();
            i_prev_location = lower_bound_steps[i];
        }
        // cout << "I: " << I << " I_prev: " << I_prev << endl;
        I_prev = I;
        // cout << "I: " << I << " I_prev: " << I_prev << endl;
        I_prev_location = lower_bound_steps[I];
        buffers.flip();
        if (I == (n_steps-1)) {
            break;
        }

        // cout << "I: " << I << " I_prev: " << I_prev << endl;
        I = min(I+jump_size, n_steps-1);
        // cout << "I: " << I << " I_prev: " << I_prev << endl;
    }
    pmfgen.compute_pmf(n-n_steps+1, intensity*(1.0-I_prev_location), pmf);
    fftconvolver.convolve_same_size(n-n_steps+1, pmf, &buffers.get_src()[n_steps], &buffers.get_dest()[n_steps]);
    fill(&buffers.get_dest()[0], &buffers.get_dest()[n_steps], 0.0);

    free_aligned_mem(tmp);
    free_aligned_mem(pmf);
    return buffers.get_dest();
}

double ecdf_lower_noncrossing_probability_n2(int n, const vector<double>& lower_bound_steps)
{
    if ((int)lower_bound_steps.size() > n) {
        stringstream ss;
        ss << "Empirical CDF must cross lower boundary g(t) since g(1)==" << lower_bound_steps.size() << " > n and the number of samples is n. h_steps:\n";
        throw runtime_error(ss.str() + vector_to_string(lower_bound_steps));
    }
    // Asymptotically any k in the range [logn, n/logn] should give optimal results as n goes to infinity.
    // Setting k=c*sqrt(n) and minimizing the asymptotic runtime, we obtain k=sqrt(2*n),
    // however, empirically 5*sqrt(n) gives better results.
    int k = min(5*sqrt(n),n/log2(n));
    vector<double> poisson_nocross_probabilities = poisson_lower_noncrossing_probability_n2(n, n, lower_bound_steps, k);
    return poisson_nocross_probabilities[n] / poisson_pmf(n, n);
}
// For n=10000, best results k=400...600

double ecdf_upper_noncrossing_probability_n2(int n, const vector<double>& upper_bound_steps)
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

    return ecdf_lower_noncrossing_probability_n2(n, symmetric_steps);
}
