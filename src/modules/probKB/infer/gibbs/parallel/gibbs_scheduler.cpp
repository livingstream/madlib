#include "gibbs_gist.h"
#include "gibbs_scheduler.h"
#include <stdint.h>

GibbsScheduler::GibbsScheduler(size_t chainId_, size_t inNumAtoms, size_t inNumClauses, GibbsGist *inState)
 : chainId(chainId_) {
   numAtoms = inNumAtoms;
   numClauses = inNumClauses;
   st = inState;
   init(false);
   /*stringstream ss;
   ss << "log" << chainId_;
   string str = ss.str();
   f.open(str.c_str(), ios::out | ios::trunc);
   
   f << "initialized truth value = " << loc_truthValues[0] << endl;
   f << "numTrueLits = " << loc_numTrueLits[0] << endl;*/     

}

GibbsScheduler::~GibbsScheduler()
{
}


void GibbsScheduler::init(bool warm)
{
   initTruthValuesAndWts(warm);
   randomInitGndPredsTruthValues(warm);
   for (size_t i = 0; i < numAtoms; i++) {
      loc_affectedGndPredIndices.push_back(i);
   }
   updateWtsForGndPreds(loc_affectedGndPredIndices);
   loc_affectedGndPredIndices.clear();
}

void GibbsScheduler::initTruthValuesAndWts(bool warm)
{
   loc_wtsWhenFalse.resize(numAtoms,0);
   loc_wtsWhenTrue.resize(numAtoms,0);
   loc_numTrue.resize(numAtoms,0);
   loc_numTrueTemp.resize(numAtoms,0);
   loc_affectedGndPredFlag.resize(numAtoms, false);
   if(!warm)
   loc_truthValues.resize(numAtoms, false);
   loc_numTrueLits.resize(numClauses, 0);
}

void GibbsScheduler::randomInitGndPredsTruthValues(bool warm)
{
   if(!warm)
   for (size_t i = 0; i < numAtoms; i++) {
     bool tv = genTruthValueForProb(0.5);
     loc_truthValues[i] = tv;
   }

   for (size_t i = 0; i < numClauses; i++) {
      GroundClause *gndClause = st->varst->getGndClause(i);
      for (size_t j = 0; j < gndClause->getNumGroundPredicates(); j++) {
         const size_t atomIdx = abs(st->varst->getAtomInClause(j, i)) - 1;
         const bool sense = gndClause->getGroundPredicateSense(j);
         if (loc_truthValues[atomIdx]== sense) {
             loc_numTrueLits[i]++;
         }
      }
   }
}


void *GibbsScheduler::callMemberFunction(void *arg)
{  
   ((GibbsScheduler *)arg)->performGibbsStep();
   pthread_exit(NULL);
}

bool GibbsScheduler::genTruthValueForProb(const double &p)
{
   if (p == 1.0) {
      return true;
   }
   if (p == 0.0) {
      return false;
   }
   bool r= random() <= p * RAND_MAX;
   return r;
}
 

/**
 * Computes the probability of a ground predicate in a chain.
 */
double GibbsScheduler::getProbabilityOfPred(const size_t &predIdx)
{
   return 1.0 / (1.0 +
                 exp(loc_wtsWhenFalse[predIdx] - loc_wtsWhenTrue[predIdx]));
}

void *GibbsScheduler::performGibbsStep()
{
   if(st->warmStart) { 
      init(true);
   }
   /*stringstream ss;
   ss << "log" << chainId_;
   string str = ss.str();
   f.open(str.c_str(), ios::out | ios::trunc);
   
   f << "initialized truth value = " << loc_truthValues[0] << endl;
   f << "numTrueLits = " << loc_numTrueLits[0] << endl;*/     

   for (size_t i = 0; i < numAtoms; i++) {
      loc_affectedGndPredIndices.push_back(i);
   }
   updateWtsForGndPreds(loc_affectedGndPredIndices);
   loc_affectedGndPredIndices.clear();

   int sample = 0;
   while(sample < st->gGla->para->samplesPerTest) {
     sample++; 
     for (size_t i = 0; i < numAtoms; i++) {
        bool newAssign = genTruthValueForProb(getProbabilityOfPred(i));
        bool truthValue = loc_truthValues[i];
        if (newAssign != truthValue) {
           loc_truthValues[i] = newAssign;
           loc_affectedGndPredIndices.clear();
           std::fill(loc_affectedGndPredFlag.begin(), loc_affectedGndPredFlag.end(), false);
           gndPredFlippedUpdates(i);
           updateWtsForGndPreds(loc_affectedGndPredIndices);
        }
        if (newAssign) {
            loc_numTrueTemp[i]+=1;
            if(!st->gGla->burnIn) {
              loc_numTrue[i]++;
            }
        }
     }
  }

  return NULL;
}

void GibbsScheduler::updateWtsForGndPreds(vector<size_t> &gndPredIndices)
{
   // for each ground predicate whose MB has changed
   for (size_t g = 0; g < gndPredIndices.size(); g++) {
      double wtIfNoChange = 0, wtIfInverted = 0, wt;
      // Ground clauses in which this pred occurs
      vector<size_t> &negGndClauses =
         st->varst->getNegOccurenceVector(gndPredIndices[g] + 1);
      vector<size_t> &posGndClauses =
         st->varst->getPosOccurenceVector(gndPredIndices[g] + 1);
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

         GroundClause *gndClause = st->varst->getGndClause(gndClauseIdx);
         wt = gndClause->wt_;
         int numSatLiterals = loc_numTrueLits[gndClauseIdx];
         if (numSatLiterals > 1) {
            if (wt > 0) {
               wtIfNoChange += wt;
               wtIfInverted += wt;
            }
         } else if (numSatLiterals == 1) {
            if (wt > 0) {
               wtIfNoChange += wt;
            }
            bool truthValue = loc_truthValues[gndPredIndices[g]];
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
         } else {
           cout << "unexpected error !!!!!!!!!!!!!!! with chainId = " << chainId << endl;
           exit(1);
         }
      } // for each ground clause that gndPred appears in

      // Clause info is stored differently for multi-chain
      if (loc_truthValues[gndPredIndices[g]]) {
         loc_wtsWhenTrue[gndPredIndices[g]] = wtIfNoChange;
         loc_wtsWhenFalse[gndPredIndices[g]] = wtIfInverted;
      } else {
         loc_wtsWhenFalse[gndPredIndices[g]] = wtIfNoChange;
         loc_wtsWhenTrue[gndPredIndices[g]] = wtIfInverted;
      }
   } // for each ground predicate whose MB has changed
}


void GibbsScheduler::gndPredFlippedUpdates(const size_t &gndPredIdx)
{
   loc_affectedGndPredIndices.push_back(gndPredIdx);

   vector<size_t> &negGndClauses =
      st->varst->getNegOccurenceVector(gndPredIdx + 1);
   vector<size_t> &posGndClauses =
      st->varst->getPosOccurenceVector(gndPredIdx + 1);
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
      gndClause = st->varst->getGndClause(gndClauseIdx);
      if (loc_truthValues[gndPredIdx] == sense) {
         loc_numTrueLits[gndClauseIdx]++;
      } else {
         loc_numTrueLits[gndClauseIdx]--;
      }

      for (size_t j = 0; j < gndClause->getNumGroundPredicates(); j++) {
         size_t predIndex = abs(gndClause->getGroundPredicateIndex(j)) - 1;
         if (!loc_affectedGndPredFlag[ predIndex ]) {
            loc_affectedGndPredIndices.push_back(predIndex);
            loc_affectedGndPredFlag[predIndex] = true;
         }

      }
   }
}
