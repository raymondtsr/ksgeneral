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

//static void print_usage()
//{
//    cout << "SYNOPSIS\n";
//    cout << "    crossprob ecdf <n> 'default.txt' [--no-fft]\n";
//    cout << endl;
//    cout << "DESCRIPTION\n";
//    cout << "    Let g(t), h(t):[0,1] -> R be two functions such that g(t) <= h(t). This program computes\n";
//    cout << "    the probability that a Poisson process or an empirical CDF will cross these boundaries.\n";
//    cout << "    For more details see https://github.com/mosco/crossing-probability\n";
//    cout << endl;
//    cout << "    crossprob ecdf <n> <boundary-functions-file> [--no-fft]\n";
//    cout << "        Computes the probability that an empirical CDF F^(t) will cross either the lower\n";
//    cout << "        boundary g(t) or the upper boundary h(t) at some point t, where F^(t) is the\n";
//    cout << "        empirical CDF of n samples drawn uniformly from the interval [0,1].\n";
//    cout << "        i.e. Letting X_1, ..., X_n ~ U[0,1], we have F^(t) = (number of X_i <= t) / n\n";
//    cout << "        One application is to compute the complementary(cdf) of the two-sided one-sample Kolmogorov-Smirnov (K-S) test statistic Dn when the underlying cdf is (dis)continuous\n";
//    cout << endl;
//    cout << "OPTIONS\n";
//    cout << "    <n>\n";
//    cout << "        We input the size of the samples in the K-S test here\n";
//    cout << endl;
//    cout << "    'default.txt'\n";
//    cout << "        This file describes the boundary functions g(t) and h(t).\n";
//    cout << "        It must contain 2 lines of monotone non-decreasing comma-separated numbers between 0 and 1\n";
//    cout << endl;
//    cout << "        In the two-sided one-sample K-S test when the underlying cdf is (dis)continuous:\n";
//    cout << "        Line 1: the i-th number in this list is inf{t in [0,1] : g(t) > i-1} = F(bi) as shown in (6) on page 5 of the manuscript by Dimitrova et al. (2016)\n";
//    cout << "        Line 2: the i-th number in this list is sup{t in [0,1] : h(t) < i} = F(ai-) as shown in (6) on page 5 of the manuscript by Dimitrova et al. (2016)\n";
//    cout << endl;
//    cout << "        Example:\n";
//    cout << "            0.3, 0.7, 0.9, 1\n";
//    cout << "            0, 0, 0.15, 0.5, 0.8\n";
//    cout << endl;
//    cout << "    --no-fft\n";
//    cout << "        Revert to algorithm [2] instead of [1].\n";
//    cout << endl;
//    cout << "NOTES\n";
//    cout << "    The two-sided crossing commands ecdf implement the O(n^2 log n) procedure described in [1] \n";
//    cout << "    unless the '--no-fft' flag is used, in which case the O(n^3) procedure of [2] is used instead.\n";
//    cout << endl;
//    cout << "REFERENCES\n";
//    cout << "    [1] A. Moscovich, B. Nadler, Fast calculation of boundary crossing probabilities for\n";
//    cout << "        Poisson processes (2016), http://arxiv.org/abs/1503.04363\n";
//    cout << "    [2] E. Khmaladze, E. Shinjikashvili, Calculation of noncrossing probabilities for Poisson\n";
//    cout << "        processes and its corollaries, Adv. Appl. Probab. 33 (2001) 702-716, http://doi.org/10.1239/aap/\n";
//}

/**********************************************************************************************************/
/*
 This function, TheoreticalDF(), calculates the cumulative sum of a given vector.
 The input, pmf, is a vector.
 The output is the vector containing the cumulative sum of the vector input.
 */
vector <double> TheoreticalDF (vector <double> pmf)
{

    double cumulative_sum = 0.0;
    vector <double> TDF(pmf.size());

    for(int i = 0; i < pmf.size(); ++i){

        cumulative_sum += pmf[i];
        TDF[i] = cumulative_sum;

    }
    return TDF;
}

/*========================================================================*/

/*
 This function, Empirical(), calculates the empirical cumulative distribution function (EDF) of a given random sample.
 Note that the ECDF starts with 0
 The input, obs, is the vector containing a random sample from a given distribution
 The output is the vector containing the EDF of the random sample.
 */
vector <double> EmpiricalDF (vector <double> obs)
{

    vector <double> observed = obs;

    set<double> s;
    for (int i = 0; i < obs.size(); ++i){
        s.insert(obs[i]);
    }
    obs.assign(s.begin(), s.end());

    vector <double> counting(obs.size());

    for (int j = 0; j < obs.size(); ++j){

        for (int i = 0; i < observed.size(); ++i){
            if (observed[i] == obs[j]){
                counting[j] += 1.0;
            }
            else
            {
                counting[j] += 0.0;
            }
        }
    }


    vector <double> masses(counting.size());


    for (int i = 0; i < counting.size(); ++i){
        masses[i] = counting[i]/observed.size();
    }


    vector <double>::iterator it;
    it = masses.begin();
    masses.insert(it, 0.0);

    vector <double> EDF = TheoreticalDF(masses);
    return EDF;
}

/*========================================================================*/

/*
 This function, BinomialPF(), defines binomial probability mass function (PMF).
 This is used in the numerical illustration of the complementary cdf of Dn for binomial(n, p) distribution.
 Note that the PMF starts with 0.
 The first input, trial, is the number of Bernoulli trials, that is, the parameter n in the binomial(n, p) distribution
 The second input, p, is the success probability of the Bernoulli trials, that is, the parameter p in the binomial(n, p) distribution
 */
vector <double> BinomialPF(int trial, double prob)
{

    int num = trial + 1;
    vector <double> PF(num);
    for (int i = 0; i < num; ++i){
        PF[i] = exp(lgamma(num) - lgamma(num - i) - lgamma(i + 1)) * pow(prob, i) * pow(1 - prob, trial - i);
    }

    vector <double>::iterator it;
    it = PF.begin();
    PF.insert(it, 0.0);
    return PF;
}

/*========================================================================*/

/*
 This function, BinomialDF(), defines the binomial cumulative distribution function (CDF).
 This is used in the numerical illustration of the complementary cdf of Dn for binomial(n, p) distribution.
 The first input, trial, is the number of Bernoulli trials, that is, the parameter n in the binomial(n, p) distribution
 The second input, p, is the sucess probability of the Bernoulli trials, that is, the parameter p in the binomial (n, p) distribution
 */

vector <double> BinomialDF(int trial, double prob)
{


    vector <double> Bin_CDF = TheoreticalDF(BinomialPF(trial, prob));
    Bin_CDF[Bin_CDF.size()-1]= 1.0;
    return Bin_CDF;

}

/*========================================================================*/

/*
 This function, MixDF(), defines the mixed CDF, discrete CDF or continuous CDF, F(x).
 The input, obs, represents x
 The output is the function defined by the users, F(x).

 Examples in Dimitrova et al. (2016) are shown:
 e.g.,
 the mixed cdf $F_{Y}(y)$ in the reinsurance example in (23),
 the Binomial (3, 0.5) distribution in Table 3,
 the Binomial (7, 0.5) distribution in Table 4,
 the Binomial (15, 0.5) distribution in Table 5,
 the Discrete Uniform (1, 10) distribution in Table 2.
 */

vector <double> MixDF (vector <double> obs)
{

    vector <double> observed = obs;

    set<double> s;
    for (int i = 0; i < obs.size(); ++i){
        s.insert(obs[i]);
    }
    obs.assign(s.begin(), s.end());


    vector <double> DF(obs.size());

    /*========================================================================*/

    /* The distribution in the reinsurance example in (23) */

//    for (int i = 0; i < obs.size(); ++i){
//        if (obs[i] < 0.0){
//            DF[i] = 0.0;
//        }
//        else if (obs[i] < log(2.5)){
//            DF[i] = 1 - 0.5 * exp(-1.0 * obs[i]);
//        }
//        else
//        {
//            DF[i] = 1.0;
//        }
//    }

    /*========================================================================*/

    /* The Binomial (3, 0.5) distribution in Table 3 */

//    vector <double> Binom_pmf = BinomialDF(3, 0.5);
//    for (int i = 0; i < obs.size(); ++i){
//        if (obs[i] < 0.0){
//            DF[i] = 0.0;
//        }
//        else if (obs[i] < 1.0){
//            DF[i] = Binom_pmf[1];
//        }
//        else if (obs[i] < 2.0){
//            DF[i] = Binom_pmf[2];
//        }
//        else if (obs[i] < 3.0){
//            DF[i] = Binom_pmf[3];
//        }
//        else
//        {
//            DF[i] = 1.0;
//        }
//    }

    /*========================================================================*/

    /* The Binomial (7, 0.5) distribution in Table 4 */

    //vector <double> Binom_pmf = BinomialDF(7, 0.5);
    //for (int i = 0; i < obs.size(); ++i){
    //    if (obs[i] < 0.0){
    //        DF[i] = 0.0;
    //    }
    //    else if (obs[i] < 1.0){
    //        DF[i] = Binom_pmf[1];
    //    }
    //    else if (obs[i] < 2.0){
    //        DF[i] = Binom_pmf[2];
    //    }
    //    else if (obs[i] < 3.0){
    //        DF[i] = Binom_pmf[3];
    //    }
    //    else if (obs[i] < 4.0){
    //        DF[i] = Binom_pmf[4];
    //    }
    //    else if (obs[i] < 5.0){
    //        DF[i] = Binom_pmf[5];
    //    }
    //    else if (obs[i] < 6.0){
    //        DF[i] = Binom_pmf[6];
    //    }
    //    else if (obs[i] < 7.0){
    //        DF[i] = Binom_pmf[7];
    //    }
    //    else
    //    {
    //        DF[i] = 1.0;
    //    }
    //}

    /*========================================================================*/

    /* The Binomial (15, 0.5) distribution in Table 5 */

    vector <double> Binom_pmf = BinomialDF(15, 0.5);
    for (int i = 0; i < obs.size(); ++i){
        if (obs[i] < 0.0){
            DF[i] = 0.0;
        }
        else if (obs[i] < 1.0){
            DF[i] = Binom_pmf[1];
        }
        else if (obs[i] < 2.0){
            DF[i] = Binom_pmf[2];
        }
        else if (obs[i] < 3.0){
            DF[i] = Binom_pmf[3];
        }
        else if (obs[i] < 4.0){
            DF[i] = Binom_pmf[4];
        }
        else if (obs[i] < 5.0){
            DF[i] = Binom_pmf[5];
        }
        else if (obs[i] < 6.0){
            DF[i] = Binom_pmf[6];
        }
        else if (obs[i] < 7.0){
            DF[i] = Binom_pmf[7];
       }
        else if (obs[i] < 8.0){
            DF[i] = Binom_pmf[8];
        }
        else if (obs[i] < 9.0){
            DF[i] = Binom_pmf[9];
        }
        else if (obs[i] < 10.0){
            DF[i] = Binom_pmf[10];
        }
        else if (obs[i] < 11.0){
            DF[i] = Binom_pmf[11];
        }
        else if (obs[i] < 12.0){
            DF[i] = Binom_pmf[12];
        }
        else if (obs[i] < 13.0){
            DF[i] = Binom_pmf[13];
        }
        else if (obs[i] < 14.0){
            DF[i] = Binom_pmf[14];
        }
        else if (obs[i] < 15.0){
            DF[i] = Binom_pmf[15];
        }
        else
        {
            DF[i] = 1.0;
        }
    }

    /*========================================================================*/

    /* An example of continuous distribution: the Exponential(1) distribution */

    //for (int i = 0; i < obs.size(); ++i){
    //    if (obs[i] < 0.0){
    //        DF[i] = 0.0;
    //    }
    //    else
    //    {
    //        DF[i] = 1 - exp(-1.0 * obs[i]);
    //    }
    //}

    /*========================================================================*/

    /* The Discrete Uniform (1, 10) distribution in Table 2 */

    //    for (int i = 0; i < obs.size(); ++i){
    //        if (obs[i] < 1.0){
    //            DF[i] = 0.0;
    //        }
    //        else if (obs[i] < 2.0){
    //            DF[i] = 0.1;
    //        }
    //        else if (obs[i] < 3.0){
    //            DF[i] = 0.2;
    //        }
    //        else if (obs[i] < 4.0){
    //            DF[i] = 0.3;
    //        }
    //        else if (obs[i] < 5.0){
    //            DF[i] = 0.4;
    //        }
    //        else if (obs[i] < 6.0){
    //            DF[i] = 0.5;
    //        }
    //        else if (obs[i] < 7.0){
    //            DF[i] = 0.6;
    //        }
    //        else if (obs[i] < 8.0){
    //            DF[i] = 0.7;
    //        }
    //        else if (obs[i] < 9.0){
    //            DF[i] = 0.8;
    //        }
    //        else if (obs[i] < 10.0){
    //            DF[i] = 0.9;
    //        }
    //        else
    //        {
    //            DF[i] = 1.0;
    //        }
    //    }


    return DF;
}

/*========================================================================*/

/*
 This function, LeftMixDF(), calculates the left limit of the mixed CDF at given points.
 The first input, points, is a vector containing points for which the left limit of the mixed CDF is calculated.
 The second input, eps, is a positive number close to zero.
 The output is the vector containing the left limit of the mixed CDF at given points.
 */
vector <double> LeftMixDF (vector <double> points, double eps)
{

    vector <double> LLimit(points.size());
    vector <double> LPoints(points.size());

    for (int i = 0; i < points.size(); ++i){

        LPoints[i] = points[i] - eps;


    }

    LLimit = MixDF(LPoints);
    return LLimit;
}

/*========================================================================*/

/*
 This function, LeftDF(), calculates the left limits of the mixed CDF at members in the random sample.
 The first input, obs, is a vector containing the random sample.
 The second input, eps, is a positive number close to zero.
 The output is the vector containing the left limit of the mixed CDF at members in the random sample.
 */
vector <double> LeftDF (vector <double> obs, double eps)
{

    vector <double> observed = obs;

    set<double> s;
    for (int i = 0; i < obs.size(); ++i){
        s.insert(obs[i]);
    }
    obs.assign(s.begin(), s.end());

    vector <double> LLimit(obs.size());
    vector <double> LPoints(obs.size());

    for (int i = 0; i < obs.size(); ++i){

        LPoints[i] = obs[i] - eps;

    }

    LLimit = MixDF(LPoints);
    return LLimit;
}

/*========================================================================*/

/*
 This function, Mixed_KS_stat(), calculates the value of the two-sided one-sample Kolmogorov-Smirnov (KS)
 test statistic Dn, possibly when the theoretical distribution is mixed, given a random sample and the
 prespecified underlying distribution F(x).

 The input, obs, is a vector containing the random sample assumed to follow the distribution F(x).
 The output is the value of the two-sided one-sample KS test statistic Dn.
 */
double Mixed_KS_stat (vector <double> obs)
{
    vector <double> observed = obs;

    set<double> s;
    for (int i = 0; i < obs.size(); ++i){
        s.insert(obs[i]);
    }
    obs.assign(s.begin(), s.end());



    vector <double> EDF = EmpiricalDF(observed);
    vector <double> TDF = MixDF(observed);
    vector <double> LTDF = LeftDF(observed, 0.0000000001);


    vector <double> Diff_left(obs.size());
    vector <double> Diff_right(obs.size());

    double Max_left;
    double Max_right;
    double d_n = 0.0;


    for (int i = 0; i < obs.size(); ++i){
        Diff_right[i] = fabs(EDF[i+1] - TDF[i]);
    }


    for (int i = 0; i < obs.size(); ++i){
        Diff_left[i] = fabs(EDF[i] - LTDF[i]);
    }

    for (int i = 0; i < obs.size(); ++i){
        sort(Diff_left.begin(), Diff_left.end());
        sort(Diff_right.begin(), Diff_right.end());

        Max_left = Diff_left[obs.size()-1];
        Max_right = Diff_right[obs.size()-1];
        d_n = max(Max_left, Max_right);
    }

    return d_n;
}

/*========================================================================*/

/*
 Given a discontinuous theoretical distribution F(x), this function calculates the F(a_{i}-) and F(b_{i}) in Step 2
 of Dimitrova et al. (2016) (cf., (6)).

 Since the KS test for mixed distribution is not distribution-free, the number of "else if" statments
 depends on the number of jumps in the mixed distribution F(x).

 For example, if there are 3 jumps in the mixed distribution F(x), one needs
 "else if ((a_bounds[i] > Left_DF[k]) && (a_bounds[i] <= Right_DF[k])){
 F_a[i] = Left_DF[k];
 }"

 and

 "else if ((b_bounds[i] >= Left_DF[k]) && (b_bounds[i] < Right_DF[k])){
 F_b[i] = Right_DF[k];
 }"

 up to k = 3-1 = 2.
 That is, the number of "else if" statements equals the number of jumps in the mixed distribution F(x).

 The first input, points, is a vector containing all points where F(x) has a jump.
 The second input, sample_size, is a number denoting the sample size, n, as shown in Step 1 of Dimitrova et al. (2016) (cf.,(3)).
 The third input, critical, is a number denoting the quantile, q, as shown in Step 1 of Dimitrova et al. (2016) (cf.,(3)).
 The fourth input, eps, is a positive number close to zero, \epsilon, as shown in Step 1 of Dimitrova et al. (2016) (cf.,(3)).

 The output is a text file containing values of $F(a_{i}-)$ and $F(b_{i})$ in Step 2 of Dimitrova et al. (2016) (cf.,(6)).
 */


/*========================================================================*/

int upperlowerbounds_discontinuous(vector <double> points, int sample_size, double critical, double eps)
{

    vector <double> Left_DF = LeftMixDF(points, 0.00000000001);
    vector <double> Right_DF = MixDF(points);


    double b_bounds = 0.0;
    double F_b = 0.0;

    b_bounds = 0.0 + critical * sample_size - eps * sample_size;
    b_bounds = b_bounds/sample_size;

    if (b_bounds < 0.0){
        F_b = 0.0;
    }
    else if ((b_bounds >= Left_DF[0]) && (b_bounds < Right_DF[0])){
        F_b = Right_DF[0];
    }
    else if ((b_bounds >= Left_DF[1]) && (b_bounds < Right_DF[1])){
        F_b = Right_DF[1];
    }
        else if ((b_bounds >= Left_DF[2]) && (b_bounds < Right_DF[2])){
            F_b = Right_DF[2];
        }
        else if ((b_bounds >= Left_DF[3]) && (b_bounds < Right_DF[3])){
            F_b = Right_DF[3];
        }
        else if ((b_bounds >= Left_DF[4]) && (b_bounds < Right_DF[4])){
            F_b = Right_DF[4];
        }
        else if ((b_bounds >= Left_DF[5]) && (b_bounds < Right_DF[5])){
            F_b = Right_DF[5];
        }
        else if ((b_bounds >= Left_DF[6]) && (b_bounds < Right_DF[6])){
            F_b = Right_DF[6];
        }
        else if ((b_bounds >= Left_DF[7]) && (b_bounds < Right_DF[7])){
            F_b = Right_DF[7];
        }
        else if ((b_bounds >= Left_DF[8]) && (b_bounds < Right_DF[8])){
            F_b = Right_DF[8];
        }
        else if ((b_bounds >= Left_DF[9]) && (b_bounds < Right_DF[9])){
            F_b = Right_DF[9];
        }
        else if ((b_bounds >= Left_DF[10]) && (b_bounds < Right_DF[10])){
            F_b = Right_DF[10];
        }
        else if ((b_bounds >= Left_DF[11]) && (b_bounds < Right_DF[11])){
            F_b = Right_DF[11];
        }
        else if ((b_bounds >= Left_DF[12]) && (b_bounds < Right_DF[12])){
            F_b = Right_DF[12];
        }
        else if ((b_bounds >= Left_DF[13]) && (b_bounds < Right_DF[13])){
            F_b = Right_DF[13];
        }
        else if ((b_bounds >= Left_DF[14]) && (b_bounds < Right_DF[14])){
            F_b = Right_DF[14];
        }
        else if ((b_bounds >= Left_DF[15]) && (b_bounds < Right_DF[15])){
            F_b = Right_DF[15];
        }
    //        else if ((b_bounds >= Left_DF[16]) && (b_bounds < Right_DF[16])){
    //            F_b = Right_DF[16];
    //        }
    //        else if ((b_bounds >= Left_DF[17]) && (b_bounds < Right_DF[17])){
    //            F_b = Right_DF[17];
    //        }
    //        else if ((b_bounds >= Left_DF[18]) && (b_bounds < Right_DF[18])){
    //            F_b = Right_DF[18];
    //        }
    //        else if ((b_bounds >= Left_DF[19]) && (b_bounds < Right_DF[19])){
    //            F_b = Right_DF[19];
    //        }
    //else if ((b_bounds >= Left_DF[20]) && (b_bounds < Right_DF[20])){
    //    F_b = Right_DF[20];
    //}
    //else if ((b_bounds >= Left_DF[21]) && (b_bounds < Right_DF[21])){
    //    F_b = Right_DF[21];
    //}
    else
    {
        F_b = (0 + 0.0)/sample_size + critical;
    }

    if (F_b >= 1.0){
        F_b = 1.0;
    }

    ofstream myfile3;

    myfile3.open("Boundary_Crossing_Time.txt");

    myfile3 << fixed << setprecision(20) << F_b;



    for (int i = 1; i < sample_size; ++i){

        b_bounds = (i + 0.0) + critical * sample_size - eps * sample_size;
        b_bounds = b_bounds/sample_size;



        if (b_bounds < 0.0){
            F_b = 0.0;
        }
        else if ((b_bounds >= Left_DF[0]) && (b_bounds < Right_DF[0])){
            F_b = Right_DF[0];
        }
        else if ((b_bounds >= Left_DF[1]) && (b_bounds < Right_DF[1])){
            F_b = Right_DF[1];
        }
                else if ((b_bounds >= Left_DF[2]) && (b_bounds < Right_DF[2])){
                    F_b = Right_DF[2];
                }
                else if ((b_bounds >= Left_DF[3]) && (b_bounds < Right_DF[3])){
                    F_b = Right_DF[3];
                }
                else if ((b_bounds >= Left_DF[4]) && (b_bounds < Right_DF[4])){
                    F_b = Right_DF[4];
                }
                else if ((b_bounds >= Left_DF[5]) && (b_bounds < Right_DF[5])){
                    F_b = Right_DF[5];
                }
                else if ((b_bounds >= Left_DF[6]) && (b_bounds < Right_DF[6])){
                    F_b = Right_DF[6];
                }
                else if ((b_bounds >= Left_DF[7]) && (b_bounds < Right_DF[7])){
                    F_b = Right_DF[7];
                }
                else if ((b_bounds >= Left_DF[8]) && (b_bounds < Right_DF[8])){
                    F_b = Right_DF[8];
                }
                else if ((b_bounds >= Left_DF[9]) && (b_bounds < Right_DF[9])){
                    F_b = Right_DF[9];
                }
                else if ((b_bounds >= Left_DF[10]) && (b_bounds < Right_DF[10])){
                    F_b = Right_DF[10];
                }
                else if ((b_bounds >= Left_DF[11]) && (b_bounds < Right_DF[11])){
                    F_b = Right_DF[11];
                }
                else if ((b_bounds >= Left_DF[12]) && (b_bounds < Right_DF[12])){
                    F_b = Right_DF[12];
                }
                else if ((b_bounds >= Left_DF[13]) && (b_bounds < Right_DF[13])){
                    F_b = Right_DF[13];
                }
                else if ((b_bounds >= Left_DF[14]) && (b_bounds < Right_DF[14])){
                    F_b = Right_DF[14];
                }
                else if ((b_bounds >= Left_DF[15]) && (b_bounds < Right_DF[15])){
                    F_b = Right_DF[15];
                }
        //        else if ((b_bounds >= Left_DF[16]) && (b_bounds < Right_DF[16])){
        //            F_b = Right_DF[16];
        //        }
        //        else if ((b_bounds >= Left_DF[17]) && (b_bounds < Right_DF[17])){
        //            F_b = Right_DF[17];
        //        }
        //        else if ((b_bounds >= Left_DF[18]) && (b_bounds < Right_DF[18])){
        //            F_b = Right_DF[18];
        //        }
        //        else if ((b_bounds >= Left_DF[19]) && (b_bounds < Right_DF[19])){
        //            F_b = Right_DF[19];
        //        }
        //else if ((b_bounds >= Left_DF[20]) && (b_bounds < Right_DF[20])){
        //    F_b = Right_DF[20];
        //}
        //else if ((b_bounds >= Left_DF[21]) && (b_bounds < Right_DF[21])){
        //    F_b = Right_DF[21];
        //}
        else
        {
            F_b = (i + 0.0)/sample_size + critical;
        }

        if (F_b >= 1.0){
            F_b = 1.0;
        }

        myfile3 << fixed << setprecision(20) << ", " << F_b;


    }

    myfile3.close();


    /**********************************************************************************/
    double a_bounds = 0.0;
    double F_a = 0.0;

    a_bounds = (0 + 1.0)/sample_size - critical + eps;

    if (a_bounds <= 0.0){
        F_a = 0.0;
    }
    else if ((a_bounds > Left_DF[0]) && (a_bounds <= Right_DF[0])){
        F_a = Left_DF[0];
    }
    else if ((a_bounds > Left_DF[1]) && (a_bounds <= Right_DF[1])){
        F_a = Left_DF[1];
    }
        else if ((a_bounds > Left_DF[2]) && (a_bounds <= Right_DF[2])){
            F_a = Left_DF[2];
        }
        else if ((a_bounds > Left_DF[3]) && (a_bounds <= Right_DF[3])){
            F_a = Left_DF[3];
        }
        else if ((a_bounds > Left_DF[4]) && (a_bounds <= Right_DF[4])){
            F_a = Left_DF[4];
        }
        else if ((a_bounds > Left_DF[5]) && (a_bounds <= Right_DF[5])){
            F_a = Left_DF[5];
        }
        else if ((a_bounds > Left_DF[6]) && (a_bounds <= Right_DF[6])){
            F_a = Left_DF[6];
        }
        else if ((a_bounds > Left_DF[7]) && (a_bounds <= Right_DF[7])){
            F_a = Left_DF[7];
        }
        else if ((a_bounds > Left_DF[8]) && (a_bounds <= Right_DF[8])){
            F_a = Left_DF[8];
        }
        else if ((a_bounds > Left_DF[9]) && (a_bounds <= Right_DF[9])){
            F_a = Left_DF[9];
        }
        else if ((a_bounds > Left_DF[10]) && (a_bounds <= Right_DF[10])){
            F_a = Left_DF[10];
        }
        else if ((a_bounds > Left_DF[11]) && (a_bounds <= Right_DF[11])){
            F_a = Left_DF[11];
        }
        else if ((a_bounds > Left_DF[12]) && (a_bounds <= Right_DF[12])){
            F_a = Left_DF[12];
        }
        else if ((a_bounds > Left_DF[13]) && (a_bounds <= Right_DF[13])){
            F_a = Left_DF[13];
        }
        else if ((a_bounds > Left_DF[14]) && (a_bounds <= Right_DF[14])){
            F_a = Left_DF[14];
        }
        else if ((a_bounds > Left_DF[15]) && (a_bounds <= Right_DF[15])){
            F_a = Left_DF[15];
        }
    //        else if ((a_bounds > Left_DF[16]) && (a_bounds <= Right_DF[16])){
    //            F_a = Left_DF[16];
    //        }
    //        else if ((a_bounds > Left_DF[17]) && (a_bounds <= Right_DF[17])){
    //            F_a = Left_DF[17];
    //        }
    //        else if ((a_bounds > Left_DF[18]) && (a_bounds <= Right_DF[18])){
    //            F_a = Left_DF[18];
    //        }
    //        else if ((a_bounds > Left_DF[19]) && (a_bounds <= Right_DF[19])){
    //            F_a = Left_DF[19];
    //        }
    //else if ((a_bounds > Left_DF[20]) && (a_bounds <= Right_DF[20])){
    //    F_a = Left_DF[20];
    //}
    //else if ((a_bounds > Left_DF[21]) && (a_bounds <= Right_DF[21])){
    //    F_a = Left_DF[21];
    //}
    else
    {
        F_a = (0 + 1.0)/sample_size - critical;
    }

    if (F_a < 0.0){
        F_a = 0.0;
    }

    myfile3.open("Boundary_Crossing_Time.txt", ios::app);
    myfile3 << "\n";


    myfile3 << fixed << setprecision(20) << F_a;

    for (int i = 1; i < sample_size; ++i){

        a_bounds = (i + 1.0)/sample_size - critical + eps;

        if (a_bounds <= 0.0){
            F_a = 0.0;
        }
        else if ((a_bounds > Left_DF[0]) && (a_bounds <= Right_DF[0])){
            F_a = Left_DF[0];
        }
        else if ((a_bounds > Left_DF[1]) && (a_bounds <= Right_DF[1])){
            F_a = Left_DF[1];
        }
                else if ((a_bounds > Left_DF[2]) && (a_bounds <= Right_DF[2])){
                    F_a = Left_DF[2];
                }
                else if ((a_bounds > Left_DF[3]) && (a_bounds <= Right_DF[3])){
                    F_a = Left_DF[3];
                }
                else if ((a_bounds > Left_DF[4]) && (a_bounds <= Right_DF[4])){
                    F_a = Left_DF[4];
                }
                else if ((a_bounds > Left_DF[5]) && (a_bounds <= Right_DF[5])){
                    F_a = Left_DF[5];
                }
                else if ((a_bounds > Left_DF[6]) && (a_bounds <= Right_DF[6])){
                    F_a = Left_DF[6];
                }
                else if ((a_bounds > Left_DF[7]) && (a_bounds <= Right_DF[7])){
                    F_a = Left_DF[7];
                }
                else if ((a_bounds > Left_DF[8]) && (a_bounds <= Right_DF[8])){
                    F_a = Left_DF[8];
                }
                else if ((a_bounds > Left_DF[9]) && (a_bounds <= Right_DF[9])){
                    F_a = Left_DF[9];
                }
                else if ((a_bounds > Left_DF[10]) && (a_bounds <= Right_DF[10])){
                    F_a = Left_DF[10];
                }
                else if ((a_bounds > Left_DF[11]) && (a_bounds <= Right_DF[11])){
                    F_a = Left_DF[11];
                }
                else if ((a_bounds > Left_DF[12]) && (a_bounds <= Right_DF[12])){
                    F_a = Left_DF[12];
                }
                else if ((a_bounds > Left_DF[13]) && (a_bounds <= Right_DF[13])){
                    F_a = Left_DF[13];
                }
                else if ((a_bounds > Left_DF[14]) && (a_bounds <= Right_DF[14])){
                    F_a = Left_DF[14];
                }
                else if ((a_bounds > Left_DF[15]) && (a_bounds <= Right_DF[15])){
                    F_a = Left_DF[15];
                }
        //        else if ((a_bounds > Left_DF[16]) && (a_bounds <= Right_DF[16])){
        //            F_a = Left_DF[16];
        //        }
        //        else if ((a_bounds > Left_DF[17]) && (a_bounds <= Right_DF[17])){
        //            F_a = Left_DF[17];
        //        }
        //        else if ((a_bounds > Left_DF[18]) && (a_bounds <= Right_DF[18])){
        //            F_a = Left_DF[18];
        //        }
        //        else if ((a_bounds > Left_DF[19]) && (a_bounds <= Right_DF[19])){
        //            F_a = Left_DF[19];
        //        }
        //else if ((a_bounds > Left_DF[20]) && (a_bounds <= Right_DF[20])){
        //    F_a = Left_DF[20];
        //}
        //else if ((a_bounds > Left_DF[21]) && (a_bounds <= Right_DF[21])){
        //    F_a = Left_DF[21];
        //}
        else
        {
            F_a = (i + 1.0)/sample_size - critical;
        }

        if (F_a < 0.0){
            F_a = 0.0;
        }

        myfile3 << fixed << setprecision(20) << ", " << F_a;

    }

    myfile3.close();

    return 0;
}

/*========================================================================*/

/*
 This function generates the upper bound [(i - 1)/n + q] and the lower bound [i/n - q] for the uniform order statistics, U(i),
 when the distribution function F(x) is continuous (cf., (25) of Dimitrova et al. (2016)).

 The first input, n, is the sample size, n in (25) of Dimitrova et al. (2016).
 The second input, x, is the quantile, q in (25) of Dimitrova et al. (2016).

 The output is a text file containing the upper bound and the lower bound for the uniform order statistics, U(i),
 when the distribution function F(x) is continuous.
 */
/*========================================================================*/


int upperlowerbounds(double n, double x)
{

    double line_1 = 0.0;///////////////////////////
    double line_2 = 0.0;///////////////////////////


    line_1 = (0 + 0.0)/n + x;
    if (line_1 >= 1.0){
        line_1 = 1.0;
    }



    ofstream myfile2;

    myfile2.open("Boundary_Crossing_Time.txt");

    myfile2 << fixed << setprecision(16) << line_1;

    for (int i = 1; i < n; ++i){

        line_1 = (i + 0.0)/n + x;
        if (line_1 >= 1.0){
            line_1 = 1.0;
        }

        myfile2 << fixed << setprecision(16) << ", " << line_1;


    }

    myfile2.close();

    line_2 = (0 + 1.0)/n - x;
    if (line_2 < 0.0){
        line_2 = 0.0;
    }


    myfile2.open("Boundary_Crossing_Time.txt", ios::app);
    myfile2 << "\n";


    myfile2 << fixed << setprecision(16) << line_2;

    for (int i = 1; i < n; ++i){

        line_2 = (i + 1.0)/n - x;
        if (line_2 < 0.0){
            line_2 = 0.0;
        }

        myfile2 << fixed << setprecision(16) << ", " << line_2;


    }

    myfile2.close();

    return 0;
}

/*========================================================================*/


void upperlowerbounds_alternative(double n, double x)
{
    
    double line_1 = 0.0;///////////////////////////
    double line_2 = 0.0;///////////////////////////
    
    
    line_1 = (0 + 0.0)/n + x;
    if (line_1 >= 1.0){
        line_1 = 1.0;
    }
    
    
    
    ofstream myfile2;
    
    myfile2.open("Boundary_Crossing_Time.txt");
    
    myfile2 << fixed << setprecision(20) << line_1;
    
    for (int i = 1; i < n; ++i){
        
        line_1 = (i + 0.0)/n + x;
        if (line_1 >= 1.0){
            line_1 = 1.0;
        }
        
        myfile2 << fixed << setprecision(20) << ", " << line_1;
        
        
    }
    
    myfile2.close();
    
    line_2 = (0 + 1.0)/n - x;
    if (line_2 < 0.0){
        line_2 = 0.0;
    }
    
    
    myfile2.open("Boundary_Crossing_Time.txt", ios::app);
    myfile2 << "\n";
    
    
    myfile2 << fixed << setprecision(20) << line_2;
    
    for (int i = 1; i < n; ++i){
        
        line_2 = (i + 1.0)/n - x;
        if (line_2 < 0.0){
            line_2 = 0.0;
        }
        
        myfile2 << fixed << setprecision(20) << ", " << line_2;
        
        
    }
    
    myfile2.close();
    
}



/**********************************************************************************/
/**********************************************************************************/
double ecdf_noncrossing_probability_loop(int n, const vector<double>& g_steps, const vector<double>& h_steps, bool use_fft, int loop)
{
    double result_;
    for (int i = 0; i < loop; ++i){

        result_ = ecdf_noncrossing_probability(n, g_steps, h_steps, use_fft);

    }

    return result_;

}

/**********************************************************************************/
/**********************************************************************************/

//int integer_crossing_times()
//{
//    int types;
//    double inputs;
//    int objectives;
//    vector <double> vector_input1;
//    vector <double> vector_input2;
//    vector <double> vector_input3;
//
//    cout << "Please enter the distribution type: 1 for Continuous Distribution, 2 for Discontinuous Distributions: " << endl;
//    cin >> types;
//    if (types == 1){
//        double sample_sizes;
//        double d_n = 0.0;
//
//        cout << "Please enter the sample sizes: " << endl;
//        cin >> sample_sizes;
//
//        cout << "Please enter the objective: 1 for K-S Complementary Distribution, 2 for P-Values: " << endl;
//        cin >> objectives;
//
//        if (objectives == 1){
//            cout << "Please enter the quantile: " << endl;
//            cin >> d_n;
//
//            //d_n = sqrt(4.0/sample_sizes);
//        }
//        else if (objectives == 2){
//            cout << "Please enter the observed values and finish by inputting -1: " << endl;
//            while((cin >> inputs) && inputs != -1)
//                vector_input2.push_back(inputs);
//
//            d_n = Mixed_KS_stat(vector_input2);
//            cout << "The K_S test statistic is: " << fixed << setprecision(8) << d_n << endl;
//        }
//
//
//        upperlowerbounds(sample_sizes, d_n);
//
//    }
//    else if (types == 2){
//
//        double sample_sizes = 0.0;
//        double critical = 0.0;
//
//
//
//        cout << "Please enter the objective: 1 for K-S Complementary Distribution, 2 for P-Values: " << endl;
//        cin >> objectives;
//
//        if (objectives == 1){
//            cout << "Please enter the sample sizes and quantile: " << endl;
//            cin >> sample_sizes >> critical;
//        }
//        else if (objectives == 2){
//            cout << "Please enter the sample sizes: " << endl;
//            cin >> sample_sizes;
//
//            cout << "Please enter the observed values and finish by inputting -1: " << endl;
//            while((cin >> inputs) && inputs != -1)
//                vector_input2.push_back(inputs);
//
//            critical = Mixed_KS_stat(vector_input2);
//            cout << "The K_S test statistic is: " << fixed << setprecision(8) << critical << endl;
//        }
//
//        //        cout << "Please enter the points where there are jumps in F(x) and finish by inputting -1: " << endl;
//        //        while((cin >> inputs) && inputs != -1)
//        //            vector_input3.push_back(inputs);
//
//        vector_input3 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
//        //vector_input3 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//
//        //vector_input3 = {0.0, log(2.5)};
//
//        vector <double> result = MixDF(vector_input3);
//
//        upperlowerbounds_discontinuous(vector_input3, sample_sizes, critical, 0.0000000001);
//
//        //
//    }
//
//
//    char * dir = getcwd(NULL, 0); // Platform-dependent, outputting the current directory
//    printf("Current dir: %s\n", dir);
//
//    return 0;
//
//}
/********************************  NEW  **********************************************/
//static int cont_ks_distribution(long n){
double cont_ks_distribution(long n){


    pair<vector<double>, vector<double> > bounds = read_boundaries_file("Boundary_Crossing_Time.txt");
    const vector<double>& lower_bound_steps = bounds.first;
    const vector<double>& upper_bound_steps = bounds.second;

    bool use_fft = true;


//    cout << fixed << setprecision(16) << "Probability: " << 1.0 - ecdf_noncrossing_probability(n, lower_bound_steps, upper_bound_steps, use_fft) << endl;

    return 1.0 - ecdf_noncrossing_probability(n, lower_bound_steps, upper_bound_steps, use_fft);
//    return 0;

}




/**********************************************************************************/
//static int handle_command_line_arguments(int argc, char* argv[])
//{
//    if (!((argc == 4) || (argc == 5))) {
//        print_usage();
//        throw runtime_error("Expecting 4 or 5 command line arguments!");
//    }
//
//
//    string command = string(argv[1]);
//    long n = string_to_long(argv[2]);
//    if (n < 0) {
//        print_usage();
//        throw runtime_error("n must be non-negative!");
//    }
//    string filename = string(argv[3]);
//    bool use_fft = true;
//    if (argc == 5) {
//        if (string(argv[4]) == "--no-fft") {
//            use_fft = false;
//        } else {
//            print_usage();
//            throw runtime_error("If a 5th command line argument is provided, it must be '--no-fft'");
//        }
//    }
//
///**********************************************************************************/
//
//
//    pair<vector<double>, vector<double> > bounds = read_boundaries_file(filename);
//    const vector<double>& lower_bound_steps = bounds.first;
//    const vector<double>& upper_bound_steps = bounds.second;
//
//    if (command == "ecdf") {
//        verify_boundary_is_valid(lower_bound_steps);
//        verify_boundary_is_valid(upper_bound_steps);
//
//        clock_t tStart = clock();
//
//
//        cout << fixed << setprecision(16) << "Probability: " << 1.0 - ecdf_noncrossing_probability(n, lower_bound_steps, upper_bound_steps, use_fft) << "\nTime taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
//    }
//    else {
//        print_usage();
//        throw runtime_error("Second command line argument must be 'ecdf'");
//    }
//
//    return 0;
//}
//
//int main(int argc, char* argv[])
//{
//    try {
//        integer_crossing_times();
//        handle_command_line_arguments(argc, argv);
//        return 0;
//    } catch (runtime_error& e) {
//        cout << "runtime_error exception caught:" << endl;
//        cout << e.what() << endl;
//        return 1;
//    } catch (ifstream::failure& e) {
//        cout << "ifstream::failure exception caught:" << endl;
//        cout << e.what() << endl;
//        return 2;
//    }
//}
