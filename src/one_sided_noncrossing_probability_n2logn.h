
#ifndef __one_sided_noncrossing_probabiity_n2logn_h__
#define __one_sided_noncrossing_probabiity_n2logn_h__

#include <vector>

double ecdf_lower_noncrossing_probability_n2logn(int n, const std::vector<double>& lower_bound_steps);
double ecdf_upper_noncrossing_probability_n2logn(int n, const std::vector<double>& upper_bound_steps);

#endif
