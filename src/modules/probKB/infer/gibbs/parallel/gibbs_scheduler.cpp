#include "gibbs_gist.h"
#include "gibbs_scheduler.h"

GibbsScheduler::GibbsScheduler(size_t chainId_, size_t inNumAtoms, GibbsGist *inState)
{
   chainId = chainId_;
   numAtoms = inNumAtoms;
   st = inState;

   for (size_t i = 0; i < st->truthValues.size(); i++) {
     bool tv = genTruthValueForProb(0.5);
     st->truthValues[i][chainId] = tv;
   }

   for (size_t i = 0; i < st->numClauses; i++) {
      GroundClause *gndClause = st->varst->getGndClause(i);
      for (size_t j = 0; j < gndClause->getNumGroundPredicates(); j++) {
         const size_t atomIdx = abs(st->varst->getAtomInClause(j, i)) - 1;
         const bool sense = gndClause->getGroundPredicateSense(j);
            if (st->truthValues[atomIdx][chainId] == sense) {
               st->numTrueLits[i][chainId]++;
            }
      }
   }

   for (size_t i = 0; i < numAtoms; i++) {
      st->affectedGndPredIndices[chainId].push_back(i);
   }
   updateWtsForGndPreds(st->affectedGndPredIndices[chainId]);
   st->affectedGndPredIndices[chainId].clear();
}

GibbsScheduler::~GibbsScheduler()
{
}

void *GibbsScheduler::callMemberFunction(void *arg)
{
   return ((GibbsScheduler *)arg)->performGibbsStep();
}

bool GibbsScheduler::genTruthValueForProb(const double &p)
{
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
 */
double GibbsScheduler::getProbabilityOfPred(const size_t &predIdx)
{
   return 1.0 / (1.0 +
                 exp(st->wtsWhenFalse[predIdx][chainId] - st->wtsWhenTrue[predIdx][chainId]));
}

void *GibbsScheduler::performGibbsStep()
{
   int sample = 0;
   while(sample < st->gGla->para->samplesPerTest) {
     sample++; 
     for (size_t i = 0; i < numAtoms; i++) {
        bool newAssign = genTruthValueForProb(getProbabilityOfPred(i));

        // Truth values are stored differently for multi-chain
        bool truthValue = st->truthValues[i][chainId];
        // If gndPred is flipped, do updates & find all affected gndPreds
        if (newAssign != truthValue) {
           st->truthValues[i][chainId] = newAssign;
           st->affectedGndPredIndices[chainId].clear();
           std::fill(st->affectedGndPredFlag[chainId].begin(), st->affectedGndPredFlag[chainId].end(), false);
           gndPredFlippedUpdates(i);
           updateWtsForGndPreds(st->affectedGndPredIndices[chainId]);
        }
        if (newAssign) {
            st->numTrueTemp[i][chainId]++;
            if(!st->gGla->burnIn) {
              st->numTrue[i][chainId]++;
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
         int numSatLiterals = st->numTrueLits[gndClauseIdx][chainId];
         if (numSatLiterals > 1) {
            if (wt > 0) {
               wtIfNoChange += wt;
               wtIfInverted += wt;
            }
         } else if (numSatLiterals == 1) {
            if (wt > 0) {
               wtIfNoChange += wt;
            }
            bool truthValue = st->truthValues[gndPredIndices[g]][chainId];
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
           cout << "unexpected error !!!!!!!!!!!!!!!" << endl;
           exit(1);
         }
      } // for each ground clause that gndPred appears in

      // Clause info is stored differently for multi-chain
      if (st->truthValues[gndPredIndices[g]][chainId]) {
         st->wtsWhenTrue[gndPredIndices[g]][chainId] = wtIfNoChange;
         st->wtsWhenFalse[gndPredIndices[g]][chainId] = wtIfInverted;
      } else {
         st->wtsWhenFalse[gndPredIndices[g]][chainId] = wtIfNoChange;
         st->wtsWhenTrue[gndPredIndices[g]][chainId] = wtIfInverted;
      }
   } // for each ground predicate whose MB has changed
}


void GibbsScheduler::gndPredFlippedUpdates(const size_t &gndPredIdx)
{
   st->affectedGndPredIndices[chainId].push_back(gndPredIdx);

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

      if (st->truthValues[gndPredIdx][chainId] == sense) {
         st->numTrueLits[gndClauseIdx][chainId]++;
      } else {
         st->numTrueLits[gndClauseIdx][chainId]--;
      }

      for (size_t j = 0; j < gndClause->getNumGroundPredicates(); j++) {
         size_t predIndex = abs(gndClause->getGroundPredicateIndex(j)) - 1;
         if (!st->affectedGndPredFlag[chainId][ predIndex ]) {
            st->affectedGndPredIndices[chainId].push_back(predIndex);
            st->affectedGndPredFlag[chainId][predIndex] = true;
         }

      }
   }
}
