#include "gibbs_gist.h"

GibbsGist::GibbsGist(size_t numAtoms_, size_t numClauses_, VariableState *varst_)
{
   numAtoms = numAtoms_;
   numClauses = numClauses_;
   numChains = NUM_THREADS;
   varst = varst_;
   gGla = new GibbsCGLA(numAtoms_, numChains, this);
   init();
   for (int i = 0; i < NUM_THREADS; i++) {
      GibbsScheduler *ls = new GibbsScheduler(i, numAtoms, this);
      gibbsVec.push_back(ls);
   }
}

void GibbsGist::init()
{
   randomInitGndPredsTruthValues();
   initTruthValuesAndWts();
}

void GibbsGist::initTruthValuesAndWts()
{
   wtsWhenFalse.resize(numAtoms);
   wtsWhenTrue.resize(numAtoms);
   numTrue.resize(numAtoms);
   numTrueTemp.resize(numAtoms);
   probs.resize(numAtoms, 0);
   for (size_t i = 0; i < numAtoms; i++) {
      wtsWhenFalse[i].resize(numChains, 0);
      wtsWhenTrue[i].resize(numChains, 0);
      numTrue[i].resize(numChains, 0);
      numTrueTemp[i].resize(numChains, 0);
   }
   affectedGndPredIndices.resize(numChains);
   affectedGndPredFlag.resize(numChains);
   for (size_t i = 0; i < numChains; i++) {
      affectedGndPredFlag[i].resize(numAtoms, false);
   }
   numTrueLits.resize(numClauses);
   for (size_t i = 0; i < numClauses; i++) {
      numTrueLits[i].resize(numChains, 0);
   }
}

void GibbsGist::randomInitGndPredsTruthValues()
{
   truthValues.resize(numAtoms);
   for (size_t i = 0; i < numAtoms; i++) {
      truthValues[i].resize(numChains, false);
   }
}

void GibbsGist::infer()
{
   do {
      for (int i = 0; i < NUM_THREADS; i++) {
         pthread_create(&threads[i], 0, &GibbsScheduler::callMemberFunction, gibbsVec[i]);
      }

      for (int i = 0; i < NUM_THREADS; i++) {
         pthread_join(threads[i], NULL);
      }

   } while (gGla->shouldIterate());

   finalize();
}

void GibbsGist::finalize()
{
   for (size_t i = 0; i < numAtoms; i++)
      for (size_t j = 0; j < numChains; j++) {
         probs[i] += numTrue[i][j];
      }
   for (size_t i = 0; i < numAtoms; i++) {
      probs[i] = probs[i] / ((double)numChains * gGla->sample);
   }
}
