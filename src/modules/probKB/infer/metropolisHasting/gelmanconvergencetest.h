#ifndef GELMANCONVERGENCETEST_H_SEP_30_2005
#define GELMANCONVERGENCETEST_H_SEP_30_2005

#include "../common/meanvariance.h"

/**
 *
 */
class GelmanConvergenceTest {
public:
    GelmanConvergenceTest(const int &numChains) :
        numChains_(numChains),
        withinChainMeanVars_(new MeanVariance[numChains]),
        numSamples_(0) {}

    ~GelmanConvergenceTest() {
        delete [] withinChainMeanVars_;
    }


    // values is array of size numChains_
    void appendNewValues(const double *const &values) {
        for (int i = 0; i < numChains_; i++) {
            withinChainMeanVars_[i].appendValue(values[i]);
        }
        numSamples_++;
    }

    // values is array of size numChains_
    void appendNewValues(const vector<bool>& values) {
        for (int i = 0; i < numChains_; i++) {
            withinChainMeanVars_[i].appendValue((double)values[i]);
        }
        numSamples_++;
    }


    double getConvergenceScore() {
        betweenChainsMeanVar_.reset();
        double totalW = 0;
        for (int i = 0; i < numChains_; i++) {
            betweenChainsMeanVar_.appendValue(withinChainMeanVars_[i].getMean());
            totalW += withinChainMeanVars_[i].getVariance();
        }
        int numValues = withinChainMeanVars_[0].getNumValues();

        double B = betweenChainsMeanVar_.getVariance() * numValues;
        double W = totalW / numChains_;

        // score as stated in "Probability and Statistics", DeGroot and Schervish
        double score = B / W;
        return score;
    }


    int getNumSamplesAdded() {
        return numSamples_;
    }


    static bool checkConvergenceOfAll(GelmanConvergenceTest *tests[],
                                      const int &numTests,const double &threshold_, const int& qid) {
        // threshold as stated in "Probability and Statistics", DeGroot & Schervish
        double threshold = 1 + 0.44 * threshold_;
        double score = tests[qid-1]->getConvergenceScore(); 
        if (score > threshold)
           return false;
        else return true;
    }


private:
    int numChains_;
    MeanVariance *withinChainMeanVars_;
    int numSamples_;
    MeanVariance betweenChainsMeanVar_;

};


#endif
