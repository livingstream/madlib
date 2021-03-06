#ifndef GROUNDCLAUSE_H_JUN_26_2005
#define GROUNDCLAUSE_H_JUN_26_2005
#include<vector>
#include<iostream>
using namespace std;

class GroundClause {
public:
    vector<int> gndPreds; // 4 + 4*n bytes (n is no. of preds)
    size_t component;
    double wt_; // 8 bytes
    size_t getNumGroundPredicates() const {
        return gndPreds.size();
    }
    int getGroundPredicateIndex(const size_t &i) const {
        if(i>=gndPreds.size()) cout << "out of bound exception" << endl;
        return gndPreds[i];
    }
    bool getGroundPredicateSense(const size_t &i) const {
        if(i>=gndPreds.size()) cout << "out of bound exception" << endl;
        return (gndPreds[i] > 0);
    }

};

#endif
