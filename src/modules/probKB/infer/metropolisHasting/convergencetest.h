#ifndef CONVERGENCETEST_H_SEP_30_2005
#define CONVERGENCETEST_H_SEP_30_2005

#include <cmath>
#include <cstring>
#include <iostream>
#include "../common/meanvariance.h"

// Grabbed from http://home.online.no/~pjacklam/notes/invnorm/impl/misra/
// normsinv.html:
//
// Lower tail quantile for standard normal distribution function.
//
// This function returns an approximation of the inverse cumulative
// standard normal distribution function.  I.e., given P, it returns
// an approximation to the X satisfying P = Pr{Z <= X} where Z is a
// random variable from the standard normal distribution.
//
// The algorithm uses a minimax approximation by rational functions
// and the result has a relative error whose absolute value is less
// than 1.15e-9.
//
// Author:    Peter J. Acklam
// (Javascript version by Alankar Misra @ Digital Sutras
// (alankar@digitalsutras.com))
// Time-stamp:  2003-05-05 05:15:14
// E-mail:    pjacklam@online.no
// WWW URL:    http://home.online.no/~pjacklam

// An algorithm with a relative error less than 1.15�10-9 in the
// entire region.

inline double NORMSINV(double p)
{
    double a[] = { -3.969683028665376e+01,  2.209460984245205e+02,
                   -2.759285104469687e+02,  1.383577518672690e+02,
                   -3.066479806614716e+01,  2.506628277459239e+00
                 };
    double b[] = { -5.447609879822406e+01,  1.615858368580409e+02,
                   -1.556989798598866e+02,  6.680131188771972e+01,
                   -1.328068155288572e+01
                 };

    double c[] = { -7.784894002430293e-03, -3.223964580411365e-01,
                   -2.400758277161838e+00, -2.549732539343734e+00,
                   4.374664141464968e+00,  2.938163982698783e+00
                 };

    double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                  2.445134137142996e+00, 3.754408661907416e+00
                 };

    // Define break-points.
    double plow  = 0.02425;
    double phigh = 1 - plow;

    // Rational approximation for lower region:
    if (p < plow) {
        double q = sqrt(-2 * log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }

    // Rational approximation for upper region:
    if (phigh < p) {
        double q  = sqrt(-2 * log(1 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
    // Rational approximation for central region:
    double q = p - 0.5;
    double r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}

/////////////////////////////////////////////////////////////////////////////


class ConvergenceTest {
public:
    // epsilonFrac_ is allowable error, as a fraction of the mean of the samples
    // gamma is the probability that you want to be within that error
    ConvergenceTest(const int &numChains, const double &gamma,
                    const double &epsilonFrac)
        : numChains_(numChains), chainTotals_(new double[numChains]),
          totalOverAllChains_(0), numSamplesPerChain_(0),
          epsilonFrac_(epsilonFrac), gamma_(gamma) {
        memset(chainTotals_, 0, numChains_ * sizeof(double));
    }


    ~ConvergenceTest() {
        delete [] chainTotals_;
    }


    // values is array of size numChains_
    void appendNewValues(const double *const &values) {
        for (int i = 0; i < numChains_; i++) {
            chainTotals_[i] += values[i];
            totalOverAllChains_ += values[i];
        }
        numSamplesPerChain_++;
    }


    void appendNewValues(const vector<bool> &values) {
        for (int i = 0; i < numChains_; i++) {
            chainTotals_[i] += ((double)values[i]);
            totalOverAllChains_ += ((double)values[i]);
        }
        numSamplesPerChain_++;
    }


    // An estimate of how many samples are neccessary, per chain
    double getConvergenceScore() {
        //double m0 = numSamplesPerChain_;
        //double S = getS();
        //double sigma = sqrt(m0) * S;
        //double v = getV(sigma);
        //return v / numChains_;
        return getV(sqrt((double)numSamplesPerChain_) * getS()) / numChains_;
    }


    int getNumSamplesAdded() const {
        return numSamplesPerChain_;
    }


    static bool checkConvergenceOfAtLeast(ConvergenceTest *tests[],
                                          const int &numTests,
                                          const double &threshold,
                                          const double &fracConverged,
                                          const int& qid) {
        if (threshold <= 0) {
            return false;
        }
        double score = tests[qid-1]->getConvergenceScore();
        if (score > threshold) 
           return true;
        else
           return false;
    }


private:
    double getS() {
        meanVar_.reset();
        for (int i = 0; i < numChains_; i++) {
            meanVar_.appendValue(chainTotals_[i] / numSamplesPerChain_);
        }
        return sqrt(meanVar_.getVariance());
    }


    double getV(const double &sigma) {
        double invNorm = NORMSINV((gamma_ + 1) / 2);
        double val = invNorm * sigma / epsilonFrac_;
        return val * val;
    }


private:
    int numChains_;
    double *chainTotals_;
    double totalOverAllChains_;
    int numSamplesPerChain_;

    double epsilonFrac_;
    double gamma_;

    MeanVariance meanVar_;
};






#endif
