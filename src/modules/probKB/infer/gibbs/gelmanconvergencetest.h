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
    void appendNewValues(const vector<double>& values, int samplesPerTest) {
        for (int i = 0; i < numChains_; i++) {
            withinChainMeanVars_[i].appendValue(values[i]);
        }
        numSamples_ += samplesPerTest;
    }

    // values is array of size numChains_
    void appendNewValues(const vector<bool>& values, int samplesPerTest) {
        for (int i = 0; i < numChains_; i++) {
            withinChainMeanVars_[i].appendValue((double)values[i]);
        }
        numSamples_ += samplesPerTest;
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
                                      const int &numTests) {
        // threshold as stated in "Probability and Statistics", DeGroot & Schervish
        double threshold = 1 + 0.44 * tests[0]->getNumSamplesAdded();
        double maxScore = -1;
        int numbad      = 0;

        for (int f = 0; f < numTests; f++) {
            double score = tests[f]->getConvergenceScore();

            if (!finite(score)) {
                numbad++;
                continue;
            }

            if (score > threshold) {
                return false;
            }

            if (score > maxScore) {
                maxScore = score;
            }
        }

        if (numbad == numTests) {
            return false;
        }

        // at this point all scores are less than the threshold
        return true;
    }


private:
    int numChains_;
    MeanVariance *withinChainMeanVars_;
    int numSamples_;
    MeanVariance betweenChainsMeanVar_;

};


#endif
