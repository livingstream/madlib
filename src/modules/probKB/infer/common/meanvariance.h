#ifndef MEANVARIANCE_H_SEP_30_2005
#define MEANVARIANCE_H_SEP_30_2005


class MeanVariance {
public:
    MeanVariance() : totalX_(0), totalXSquared_(0), numValues_(0) {}

    void reset() {
        totalX_ = 0;
        totalXSquared_ = 0;
        numValues_ = 0;
    }

    void appendValue(const double &x) {
        totalX_ += x;
        totalXSquared_ += x * x;
        numValues_++;
    }

    int getNumValues() const        {
        return numValues_;
    }
    void setNumValues(const int &n) {
        numValues_ = n;
    }
    double getMean() const          {
        return totalX_ / numValues_;
    }

    double getVariance() {
        // mean(xsquared) - (mean(X))^2
        double meanXSquared = totalXSquared_ / numValues_;
        double meanX = totalX_ / numValues_;
        return (((double)numValues_) / (numValues_ - 1)) * (meanXSquared - meanX * meanX);
    }

private:
    double totalX_;
    double totalXSquared_;
    int numValues_;
};

#endif
