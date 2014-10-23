#ifndef VARIABLESTATE_H_
#define VARIABLESTATE_H_
#include "groundclause.h"
#include <vector>
#include <stdlib.h>     /* abs */

class VariableState {
public:
    VariableState(size_t numAtoms_) {
        numAtoms = numAtoms_;
        gndClauses_ = new vector<GroundClause *>;
    }
    ~VariableState() { }

    void init() {
        occurence_.resize(2 * numAtoms + 1);
        for (size_t i = 0; i < (*gndClauses_).size(); i++) {
            GroundClause *gndClause = getGndClause(i);
            size_t numGndPreds = gndClause->gndPreds.size();
            for (size_t j = 0; j < numGndPreds; j++) {
                int lit = getAtomInClause(j, i);
                int litIdx = 2 * abs(lit) - (lit > 0);
                occurence_[litIdx].push_back(i);
            }
        }
    }

    vector<size_t> &getNegOccurenceVector(const size_t &atomIdx) {
        size_t litIdx = 2 * atomIdx;
        return getOccurenceVector(litIdx);
    }

    vector<size_t> &getPosOccurenceVector(const size_t &atomIdx) {
        size_t litIdx = 2 * atomIdx - 1;
        return getOccurenceVector(litIdx);
    }

    vector<size_t> &getOccurenceVector(const size_t &idx) {
        return occurence_[idx];
    }

    size_t getNumTrueLits(const size_t &clauseIdx) {
        return numTrueLits_[clauseIdx];
    }

    int getAtomInClause(const size_t &atomIdxInClause, const size_t &clauseIdx) {
        return getGndClause(clauseIdx)->gndPreds.at(atomIdxInClause);
    }

    GroundClause *getGndClause(const size_t &index) {
        return (*gndClauses_)[index];
    }

public:
    size_t numAtoms;
    //vector<GroundClause*>* gndClauses_;
    vector<GroundClause *> *gndClauses_;
    // Number of true literals in each clause
    vector<size_t> numTrueLits_;
    // Indexed as occurence_[2*abs(lit) - (lit > 0)][occurence_num]
    vector<vector<size_t> > occurence_;
};

#endif /*VARIABLESTATE_H_*/
