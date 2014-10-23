#ifndef VOSE_H_
#define VOSE_H_
#include <cstdlib>
#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <numeric>
#include <queue>
using namespace std;
/**
  * This is the Vose Alias algorithm. It creates two huge arrays
  * to help in bias sampling.
  *
  */
class Vose {
  private:
    vector<double> probability;
    vector<int> alias;

  public:
    Vose() {
      this->probability = vector<double>();
      this->alias = vector<int>();
    }
    Vose(vector<double> prob, vector<int> _alias):
      probability(prob), alias(_alias) {}
    Vose& operator=(const Vose& rhs) {
      this->probability = rhs.probability;
      this->alias = rhs.alias;
      return *this;
    }
    Vose(bool normalize, vector<double> _scores)/*: scores(_scores)*/ {
      build(normalize, _scores);
    }
    /**
      * Take a vector of scores that are positive real values.
      */
    void build(bool normalize, vector<double> scores) {
      probability.resize(scores.size());
      alias.resize(scores.size());

      double denom = 0.0L;
      double average = 1.0 / (double)scores.size();
      if (normalize) {
        denom = std::accumulate(scores.begin(), scores.end(), 0.0);
        for(size_t i=0; i<scores.size(); i++) {
           scores[i] *= (1/denom);
        }
      }
      //log_trace("build::denom == %lf", denom);

      deque<int> small;
      deque<int> large;

      if (scores.size() == 1) {
        scores[0] = 1.0L;
      }
      else {
        for (size_t i = 0; i != scores.size(); ++i) {
          if (scores[i] >= average)
            large.push_back((int)i);
          else
            small.push_back((int)i);
        }
        while (!large.empty() && !small.empty()) {
          //if(counter++ % 1000 == 0) printf("::%d::", counter); // debug
          
          int less = small.back(); small.pop_back(); 
          int more = large.back(); large.pop_back(); 

          probability[less] = (double)scores[less] * (double)scores.size();
          alias[less] = more;

          scores[more] += scores[less] - average;
          
          if (scores[more] >= average) 
            large.push_back(more);
          else
            small.push_back(more);
        }
        //log_trace("small.size: %d, large.size: %d", small.size(), large.size());

        while (!small.empty()) { 
          probability[small.back()] = 1.0l;
          small.pop_back();
        }
        while (!large.empty()) { 
          probability[large.back()] = 1.0l;
          large.pop_back();
        }
      }

    }

  float RandDouble() {
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  }
  /** This needs a randomInt and a randomDouble */
  int next(/* Pass in random numbers? */ ) {
    int randint = rand();
    int column = randint % (int)probability.size();   
    bool coinToss = RandDouble() < probability[column];
    int result = (coinToss)?column: alias[column];
    return result;
  }

};
#endif
