
#ifndef __one_sided_noncrossing_probabiity_n2_h__
#define __one_sided_noncrossing_probabiity_n2_h__

#include <vector>

std::vector<double> poisson_lower_noncrossing_probability_n2(int n, double intensity, const std::vector<double>& lower_bound_steps, int jump_size);
double ecdf_lower_noncrossing_probability_n2(int n, const std::vector<double>& lower_bound_steps);
double ecdf_upper_noncrossing_probability_n2(int n, const std::vector<double>& upper_bound_steps);

#endif
