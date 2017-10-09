#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cassert>
#include "mm_malloc.h"

#include "one_sided_noncrossing_probability.h"

using namespace std;

#define FLOAT_PRECISION_BITS 80

#if FLOAT_PRECISION_BITS == 128
#include "string_utils.hh"
    extern "C" {
        #include <quadmath.h>
    }
    typedef __float128 FLOAT;
    #define EXP expq
    #define LOG logq
    #define GAMMA tgammaq
    #define LOG_GAMMA lgammaq
    #define POW powq
    #define FREXP frexpq
    #define LDEXP ldexpq
    #define MAX_EXP FLT128_MAX_EXP
    #define TO_STRING float128_to_string
    string float128_to_string(__float128 x)
    {
        static char s[1000];
        quadmath_snprintf(s, sizeof(s), "%.30Qg", x);
        return string(s);
    }
#elif FLOAT_PRECISION_BITS == 80
    typedef long double FLOAT;
    #include <cmath>
    #include <cfloat>
    #define EXP expl
    #define LOG logl
    #define GAMMA tgammal
    #define LOG_GAMMA lgammal
    #define POW powl
    #define FREXP frexpl
    #define LDEXP ldexpl
    #define MAX_EXP LDBL_MAX_EXP
    #define TO_STRING long_double_to_string
#elif FLOAT_PRECISION_BITS == 64
    typedef double FLOAT;
    #include <cmath>
    #include <cfloat>
    #define EXP exp
    #define LOG log
    #define GAMMA tgamma
    #define LOG_GAMMA lgamma
    #define POW pow
    #define FREXP frexp
    #define LDEXP ldexp
    #define MAX_EXP DBL_MAX_EXP
    #define TO_STRING double_to_string
#else
    #error FLOAT_PRECISION_BITS must be 64, 80 or 128.
#endif

// static string double_to_string(double x)
// {
//     stringstream ss;
//     ss.precision(17); // http://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
//     ss << x;
//     return ss.str();
// }

static string long_double_to_string(long double x)
{
    stringstream ss;
    ss.precision(22); // http://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
    ss << x;
    return ss.str();
}

class PolynomialTranslatedMonomials {
public:
    PolynomialTranslatedMonomials(int max_degree);
    ~PolynomialTranslatedMonomials();
    FLOAT get_multiplicative_coefficient(int degree) const;
    void set_multiplicative_coefficient(int degree, FLOAT multiplicative_coefficient);
    void set_additive_coefficient(int degree, FLOAT additive_coefficient);
    FLOAT evaluate(FLOAT x) const;
    void integrate();
    void ldexp_all_multiplicative_coefficients(int exp);
    friend ostream& operator<<(ostream& stream, const PolynomialTranslatedMonomials& poly);
    int degree;
private:
    FLOAT* __restrict__ multiplicative_coefficients;
    FLOAT* __restrict__ additive_coefficients;
};

PolynomialTranslatedMonomials::PolynomialTranslatedMonomials(int max_degree) :
    degree(0)
{
    assert(max_degree >= 0);
    multiplicative_coefficients = (FLOAT*)_mm_malloc(sizeof(FLOAT)*(max_degree+1), 32);
    additive_coefficients = (FLOAT*)_mm_malloc(sizeof(FLOAT)*(max_degree+1), 32);
    memset(multiplicative_coefficients, 0, sizeof(FLOAT)*(max_degree+1));
    memset(additive_coefficients, 0, sizeof(FLOAT)*(max_degree+1));
}

PolynomialTranslatedMonomials::~PolynomialTranslatedMonomials()
{
    _mm_free(multiplicative_coefficients);
    _mm_free(additive_coefficients);
}

FLOAT PolynomialTranslatedMonomials::get_multiplicative_coefficient(int degree) const
{
    assert(degree >= 0);
    assert(degree <= this->degree);
    return multiplicative_coefficients[degree];
}

void PolynomialTranslatedMonomials::set_multiplicative_coefficient(int degree, FLOAT multiplicative_coefficient)
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    multiplicative_coefficients[degree] = multiplicative_coefficient;
}

void PolynomialTranslatedMonomials::set_additive_coefficient(int degree, FLOAT additive_coefficient)
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    additive_coefficients[degree] = additive_coefficient;
}

FLOAT PolynomialTranslatedMonomials::evaluate(FLOAT x) const
{
    FLOAT result = 0.0;
    for (int i = 0; i < degree+1; ++i) {
        result += multiplicative_coefficients[i] * POW(x+additive_coefficients[i], i);
    }
    return result;
}

void PolynomialTranslatedMonomials::integrate()
{
    for (int i = degree+1; i >= 1; --i) {
        multiplicative_coefficients[i] = multiplicative_coefficients[i-1] / FLOAT(i);
        additive_coefficients[i] = additive_coefficients[i-1];
    }
    additive_coefficients[1] = 0.0;
    additive_coefficients[0] = 0.0;
    multiplicative_coefficients[0] = 0.0;
    ++degree;
}

void PolynomialTranslatedMonomials::ldexp_all_multiplicative_coefficients(int exp)
{
    for (int i = 0; i < degree+1; ++i) {
        multiplicative_coefficients[i] = LDEXP(multiplicative_coefficients[i], exp);
    }
}

// static void print_FLOAT_array(const FLOAT* arr, int n)
// {
//     for (int i = 0; i < n; ++i) {
//         cout << TO_STRING(arr[i]) << ", ";
//     }
//     cout << endl;
// }

ostream& operator<<(ostream& stream, const PolynomialTranslatedMonomials& poly)
{
    for (int i = poly.degree; i >= 0; --i) {
        FLOAT coef = poly.multiplicative_coefficients[i];
        if (coef >= 0) {
            stream << (i != poly.degree ? " + " : "") << TO_STRING(coef);
            
        } else {
            stream << (i != poly.degree ? " - " : "") << TO_STRING(-coef);
        }
        if (i >= 1) {
            stream << " (x + " << TO_STRING(poly.additive_coefficients[i]) << ")^" << i;
        }
    }
    return stream;
}

// static FLOAT minimum_multiplicative_coefficient(const PolynomialTranslatedMonomials& p)
// {
//     assert(p.degree >= 1);
//     FLOAT m = p.get_multiplicative_coefficient(0);
//     for (int i = 0; i < p.degree+1; ++i) {
//         FLOAT x = p.get_multiplicative_coefficient(i);
//         if (x < m) {
//             m = x;
//         }
//     }
//     return m;
// }

static FLOAT maximum_multiplicative_coefficient(const PolynomialTranslatedMonomials& p)
{
    assert(p.degree >= 1);
    FLOAT m = p.get_multiplicative_coefficient(0);
    for (int i = 0; i < p.degree+1; ++i) {
        FLOAT x = p.get_multiplicative_coefficient(i);
        if (x > m) {
            m = x;
        }
    }
    return m;
}

static int extract_exponent(FLOAT x)
{
    static int exponent;
    FREXP(x, &exponent);
    return exponent;
}

double ecdf_upper_noncrossing_probability(int n, const vector<double>& upper_bound_steps)
{
    const int NUM_ITERATIONS_BETWEEN_EXP_FIXES = 16;

    // cout << "Using translated polynomials! precision: " << FLOAT_PRECISION_BITS << endl;
    // cout << "sizeof(FLOAT) = " << sizeof(FLOAT) << endl;
    if ((int)upper_bound_steps.size() < n) {
        stringstream ss;
        ss << "empirical CDF must cross upper boundary h(t) since h(1)==" << upper_bound_steps.size() << " < n and the number of samples is n.";
        throw runtime_error(ss.str());
    }
    PolynomialTranslatedMonomials p(n+1);
    p.set_multiplicative_coefficient(0, 1.0);

    int total_exponent_delta = 0;
    for (int i = 0; i < (int)upper_bound_steps.size(); ++i) {
        //cout << "============= i == " << i << ": " << p << endl;
        p.integrate();
        p.set_additive_coefficient(1, -upper_bound_steps[i]);
        p.set_multiplicative_coefficient(0, -p.evaluate(upper_bound_steps[i]));

        //cout << "Before expfix: " << p << endl;

        if ((i % NUM_ITERATIONS_BETWEEN_EXP_FIXES) == 0) {
            FLOAT max_coef = maximum_multiplicative_coefficient(p);
            int highest_exponent = extract_exponent(max_coef);
            if (highest_exponent < 0) {
                //cout << "highest exponent: " << highest_exponent << endl;
                p.ldexp_all_multiplicative_coefficients(0 - highest_exponent);
                total_exponent_delta += 0 - highest_exponent;
            }
        }
    }
    //FLOAT max_coef = maximum_multiplicative_coefficient(p);
    //FLOAT min_coef = minimum_multiplicative_coefficient(p);
    FLOAT integral_result = p.evaluate(1);
    FLOAT log_integral_result = LOG(integral_result);
    FLOAT noncrossing_probability = EXP(LOG_GAMMA(n+1) + log_integral_result - total_exponent_delta*LOG(2));
    // cout << p << endl;
    // cout << "max_coef: " << TO_STRING(max_coef) << endl;
    // cout << "min_coef: " << TO_STRING(min_coef) << endl;
    // cout << "Integral result: (no exponent fix) " << TO_STRING(integral_result) << endl;
    // cout << "Exponent delta: " << total_exponent_delta << endl;
    // cout << "LOG_GAMMA(n+1) == " << TO_STRING(LOG_GAMMA(n+1)) << endl;
    // cout << "log_integral_result (no exponent fix) == " << TO_STRING(log_integral_result) << endl;
    // cout << "Final result: " << TO_STRING(noncrossing_probability) << endl;
    return noncrossing_probability;
}

double ecdf_lower_noncrossing_probability(int n, const vector<double>& lower_bound_steps)
{
    if ((int)lower_bound_steps.size() > n) {
        stringstream ss;
        ss << "Empirical CDF must cross lower boundary g(t) since g(1)==" << lower_bound_steps.size() << " > n and the number of samples is n.";
        throw runtime_error(ss.str());
        return 0;
    }

    vector<double> symmetric_steps(n, 0.0);
    for (int i = n-lower_bound_steps.size(); i < n; ++i) {
        symmetric_steps[i] = 1.0 - lower_bound_steps[lower_bound_steps.size() - 1 - i];
    }

    return ecdf_upper_noncrossing_probability(n, symmetric_steps);
}
