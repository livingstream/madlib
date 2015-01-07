#include <vector>
#include <stdlib.h>
#include <math.h>
#include "../../common/timer.h"
#include "../../state/variablestate.h"
#include "../../state/groundclause.h"
#include "../convergencetest.h"
#include "../gelmanconvergencetest.h"
#include "../gibbsparams.h"
#include <iostream>
using namespace std;
class Gibbs {
public:
    size_t numAtoms;
    size_t numClauses;
    bool warmStart;
    VariableState *state;
    vector<bool> affectedGndPredFlag;
    vector<vector<bool> > truthValues;
    vector<vector<double> > wtsWhenFalse;
    vector<vector<double> > wtsWhenTrue;
    vector<double> numTrue;
    vector<vector<int> > numTrueLits;
    // No. of chains which MCMC will use
    size_t numChains;
    // Min. no. of burn-in steps MCMC will take per chain
    int burnMinSteps;
    // Max. no. of burn-in steps MCMC will take per chain
    int burnMaxSteps;
    // Min. no. of sampling steps MCMC will take per chain
    int minSteps;
    // Max. no. of sampling steps MCMC will take per chain
    int maxSteps;
    // Max. no. of seconds MCMC should run
    int maxSeconds;

    // Gamma used by convergence test
    double gamma;
    // Epsilon used by convergence test
    double epsilonError;
    // Fraction of samples needed to converge
    double fracConverged;
    // Number of samples between checking for convergence
    int samplesPerTest;
    // Convergence test for burning in
    GelmanConvergenceTest **burnConvergenceTests;
    // Convergence test for sampling
    ConvergenceTest **gibbsConvergenceTests;

    Gibbs(size_t inNumAtoms, size_t inNumClauses, VariableState *inState, bool inWarmStart) {
        GibbsParams *params = new GibbsParams();
        numAtoms = inNumAtoms;
        numClauses = inNumClauses;
        state = inState;
        warmStart = inWarmStart;
        numChains = params->numChains;
        burnMinSteps = params->burnMinSteps;
        burnMaxSteps = params->burnMaxSteps;
        minSteps = params->minSteps;
        maxSteps = params->maxSteps;
        maxSeconds = params->maxSeconds;
        gamma = params->gamma;
        epsilonError = params->epsilonError;
        fracConverged = params->fracConverged;
        samplesPerTest = params->samplesPerTest;
        delete params;
        init();
    }

    ~Gibbs() {
        deleteConvergenceTests();
    }

    static void  *callMemberFunction(void *arg) {
        return ((Gibbs *)arg)->infer();
    }

    void init() {
        randomInitGndPredsTruthValues();
        initTruthValuesAndWts();
        initNumTrueLits();
        initNumTrue();
        initConvergenceTests();
    }
    void warmInit() {
        initTruthValuesAndWts();
        initNumTrueLits();
    }

    void initTruthValuesAndWts() {
        wtsWhenFalse.resize(numAtoms);
        wtsWhenTrue.resize(numAtoms);
        affectedGndPredFlag.resize(numAtoms, false);
        for (size_t i = 0; i < numAtoms; i++) {
            wtsWhenFalse[i].resize(numChains, 0);
            wtsWhenTrue[i].resize(numChains, 0);
        }

        numTrueLits.resize(numClauses);
        for (size_t i = 0; i < numClauses; i++) {
            numTrueLits[i].resize(numChains, 0);
        }
    }

    void randomInitGndPredsTruthValues() {
        truthValues.resize(numAtoms);
        for (size_t i = 0; i < numAtoms; i++) {
            truthValues[i].resize(numChains, false);
        }
        for (size_t c = 0; c < numChains; c++) {
            // Random tv for all not in blocks
            for (size_t i = 0; i < truthValues.size(); i++) {
                bool tv = genTruthValueForProb(0.5);
                // Truth values are stored differently for multi-chain
                truthValues[i][c] = tv;
            }
        }
    }


    bool genTruthValueForProb(const double &p) {
        if (p == 1.0) {
            return true;
        }
        if (p == 0.0) {
            return false;
        }
        bool r = random() <= p * RAND_MAX;
        return r;
    }

    /**
     * Computes the probability of a ground predicate in a chain.
     *
     * @param predIdx Index of predicate.
     * @param chainIdx Index of chain.
     * @param invTemp InvTemp used in simulated tempering.
     *
     * @return Probability of predicate.
     */
    double getProbabilityOfPred(const size_t &predIdx, const size_t &chainIdx,
                                const double &invTemp) {
        // Different for multi-chain
        return 1.0 /
               (1.0 + exp((wtsWhenFalse[predIdx][chainIdx] -
                           wtsWhenTrue[predIdx][chainIdx]) *
                          invTemp));
    }

    void performGibbsStep(const size_t &chainIdx, const bool &burningIn,
                          vector<size_t> &affectedGndPredIndices) {
        // Now go through all preds not in blocks
        for (size_t i = 0; i < numAtoms; i++) {
            bool newAssignment
                = genTruthValueForProb(getProbabilityOfPred(i, chainIdx, 1));

            // Truth values are stored differently for multi-chain
            bool truthValue = truthValues[i][chainIdx];
            // If gndPred is flipped, do updates & find all affected gndPreds
            if (newAssignment != truthValue) {
                truthValues[i][chainIdx] = newAssignment;
                affectedGndPredIndices.clear();
                std::fill(affectedGndPredFlag.begin(), affectedGndPredFlag.end(), false);
                gndPredFlippedUpdates(i, chainIdx, affectedGndPredIndices);
                updateWtsForGndPreds(affectedGndPredIndices, chainIdx);
            }

            if (!burningIn && newAssignment) {
                numTrue[i]++;
            }
        }

    }

    void updateWtsForGndPreds(vector<size_t> &gndPredIndices, const size_t &chainIdx) {
        // for each ground predicate whose MB has changed
        for (size_t g = 0; g < gndPredIndices.size(); g++) {
            double wtIfNoChange = 0, wtIfInverted = 0, wt;
            // Ground clauses in which this pred occurs
            vector<size_t> &negGndClauses =
                state->getNegOccurenceVector(gndPredIndices[g] + 1);
            vector<size_t> &posGndClauses =
                state->getPosOccurenceVector(gndPredIndices[g] + 1);
            size_t gndClauseIdx;
            bool sense;

            for (size_t i = 0; i < negGndClauses.size() + posGndClauses.size(); i++) {
                if (i < negGndClauses.size()) {
                    gndClauseIdx = negGndClauses[i];
                    sense = false;
                } else {
                    gndClauseIdx = posGndClauses[i - negGndClauses.size()];
                    sense = true;
                }

                GroundClause *gndClause = state->getGndClause(gndClauseIdx);
                wt = gndClause->wt_;
                int numSatLiterals = numTrueLits[gndClauseIdx][chainIdx];
                if (numSatLiterals > 1) {
                    if (wt > 0) {
                        wtIfNoChange += wt;
                        wtIfInverted += wt;
                    }
                } else if (numSatLiterals == 1) {
                    if (wt > 0) {
                        wtIfNoChange += wt;
                    }
                    bool truthValue = truthValues[gndPredIndices[g]][chainIdx];
                    if (truthValue == sense) {
                        if (wt < 0) {
                            wtIfInverted += fabs(wt);
                        }
                    } else {
                        if (wt > 0) {
                            wtIfInverted += wt;
                        }
                    }
                } else if (numSatLiterals == 0) {
                    // None satisfy, so when gndPred switch to its negative, it'll satisfy
                    if (wt > 0) {
                        wtIfInverted += wt;
                    } else if (wt < 0) {
                        wtIfNoChange += fabs(wt);
                    }
                }
            } // for each ground clause that gndPred appears in

            // Clause info is stored differently for multi-chain
            if (truthValues[gndPredIndices[g]][chainIdx]) {
                wtsWhenTrue[gndPredIndices[g]][chainIdx] = wtIfNoChange;
                wtsWhenFalse[gndPredIndices[g]][chainIdx] = wtIfInverted;
            } else {
                wtsWhenFalse[gndPredIndices[g]][chainIdx] = wtIfNoChange;
                wtsWhenTrue[gndPredIndices[g]][chainIdx] = wtIfInverted;
            }
        } // for each ground predicate whose MB has changed
    }


    void gndPredFlippedUpdates(const size_t &gndPredIdx, const size_t &chainIdx,
                               vector<size_t> &affectedGndPredIndices) {
        affectedGndPredIndices.push_back(gndPredIdx);

        vector<size_t> &negGndClauses =
            state->getNegOccurenceVector(gndPredIdx + 1);
        vector<size_t> &posGndClauses =
            state->getPosOccurenceVector(gndPredIdx + 1);
        size_t gndClauseIdx;
        GroundClause *gndClause;
        bool sense;

        // Find the Markov blanket of this ground predicate
        for (size_t i = 0; i < negGndClauses.size() + posGndClauses.size(); i++) {
            if (i < negGndClauses.size()) {
                gndClauseIdx = negGndClauses[i];
                sense = false;
            } else {
                gndClauseIdx = posGndClauses[i - negGndClauses.size()];
                sense = true;
            }
            gndClause = state->getGndClause(gndClauseIdx);

            if (truthValues[gndPredIdx][chainIdx] == sense) {
                numTrueLits[gndClauseIdx][chainIdx]++;
            } else {
                numTrueLits[gndClauseIdx][chainIdx]--;
            }

            for (size_t j = 0; j < gndClause->getNumGroundPredicates(); j++) {
                size_t predIndex = abs(gndClause->getGroundPredicateIndex(j)) - 1;
                if (!affectedGndPredFlag[ predIndex ]) {
                    affectedGndPredIndices.push_back(predIndex);
                    affectedGndPredFlag[ predIndex ] = true;
                }

            }
        }
    }

    void initNumTrueLits() {
        for (size_t i = 0; i < numClauses; i++) {
            GroundClause *gndClause = state->getGndClause(i);
            for (size_t j = 0; j < gndClause->getNumGroundPredicates(); j++) {
                const size_t atomIdx = abs(state->getAtomInClause(j, i)) - 1;
                const bool sense = gndClause->getGroundPredicateSense(j);
                for (size_t c = 0; c < numChains; c++) {
                    if (truthValues[atomIdx][c] == sense) {
                        numTrueLits[i][c]++;
                    }
                }
            }
        }
    }

    void initNumTrue() {
        numTrue.resize(numAtoms);
        for (size_t i = 0; i < numTrue.size(); i++) {
            numTrue[i] = 0;
        }
    }

    void setProbTrue(const size_t &predIdx, const double &p) {
        numTrue[predIdx] = p;
        //std::cout << "probability of predicate " << predIdx + 1 << " = " << p << endl;
    }

    void *infer() {
        if(warmStart) {
           warmInit();
        }
        Timer timer;
        bool burningIn = (burnMaxSteps > 0) ? true : false;
        double secondsElapsed = 0;
        double startTimeSec = timer.time();
        double currentTimeSec;

        vector<size_t> affectedGndPredIndices;

        for (size_t i = 0; i < numAtoms; i++) {
            affectedGndPredIndices.push_back(i);
        }
        for (size_t c = 0; c < numChains; c++) {
            updateWtsForGndPreds(affectedGndPredIndices, c);
        }
        affectedGndPredIndices.clear();
        // Sampling loop
        int sample = 0;
        size_t numSamplesPerPred = 0;
        bool done = false;
        while (!done) {
            ++sample;
            cout << "number of samples = " << sample << endl;
            //cout << "second elapsed = " << timer.time() - startTimeSec << endl;
            if (sample % samplesPerTest == 0) {
                currentTimeSec = timer.time();
                secondsElapsed = currentTimeSec - startTimeSec;
            }
            for (size_t c = 0; c < numChains; c++) {
                performGibbsStep(c, burningIn, affectedGndPredIndices);
                if (!burningIn) {
                    numSamplesPerPred++;
                }
            }

	    // Add current truth values to the convergence testers
	    for (size_t i = 0; i < numAtoms; i++) {
		    //WARNING: implicit cast from bool* to double*
		    if (burningIn) burnConvergenceTests[i]->appendNewValues(truthValues[i], 1);
		    else           gibbsConvergenceTests[i]->appendNewValues(truthValues[i], 1);
	    }

            if (sample % samplesPerTest != 0) {
                continue;
            }
            if (burningIn) {
                bool burnConverged = GelmanConvergenceTest::checkConvergenceOfAll(
                                         burnConvergenceTests, (int)numAtoms);
                if ((sample >= burnMinSteps && burnConverged)
                        || (burnMaxSteps >= 0 && sample >= burnMaxSteps)
                        || (maxSeconds > 0 && secondsElapsed >= maxSeconds)) {
                    burningIn = false;
                    sample = 0;
                }
            } else {
                bool gibbsConverged = ConvergenceTest::checkConvergenceOfAtLeast(
                                          gibbsConvergenceTests, (int)numAtoms, sample, fracConverged);
                if ((sample >= minSteps && gibbsConverged)
                        || (maxSteps >= 0 && sample >= maxSteps)
                        || (maxSeconds > 0 && secondsElapsed >= maxSeconds)) {
                    done = true;
                }
            }
        }
        // update gndPreds probability that it is true
        for (size_t i = 0; i < numAtoms; i++) {
            setProbTrue(i, numTrue[i] / static_cast<double>(numSamplesPerPred));
        }
        return NULL;
    }


    void initConvergenceTests() {
        burnConvergenceTests = new GelmanConvergenceTest*[numAtoms];
        gibbsConvergenceTests = new ConvergenceTest*[numAtoms];
        for (size_t i = 0; i < numAtoms; i++) {
            burnConvergenceTests[i]  = new GelmanConvergenceTest((int)numChains);
            gibbsConvergenceTests[i] = new ConvergenceTest((int)numChains, gamma, epsilonError);
        }
    }


    void deleteConvergenceTests() {
        for (size_t i = 0; i < numAtoms; i++) {
            delete burnConvergenceTests[i];
            delete gibbsConvergenceTests[i];
        }
        delete [] burnConvergenceTests;
        delete [] gibbsConvergenceTests;
    }

};
