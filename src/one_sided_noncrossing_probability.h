#ifndef __one_sided_noncrossing_probabiity_h__
#define __one_sided_noncrossing_probabiity_h__

#include <vector>

double ecdf_lower_noncrossing_probability(int n, const std::vector<double>& lower_bound_steps);
double ecdf_upper_noncrossing_probability(int n, const std::vector<double>& upper_bound_steps);

#endif
