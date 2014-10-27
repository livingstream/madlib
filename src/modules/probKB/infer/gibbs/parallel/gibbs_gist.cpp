#include "gibbs_gist.h"

GibbsGist::GibbsGist(size_t numAtoms_, size_t numClauses_, VariableState *varst_)
{
   numAtoms = numAtoms_;
   numClauses = numClauses_;
   numChains = NUM_THREADS;
   varst = varst_;
   gGla = new GibbsCGLA(numAtoms_, numChains, this);
   initTruthValues();
   for (int i = 0; i < NUM_THREADS; i++) {
      GibbsScheduler *ls = new GibbsScheduler(i, numAtoms, numClauses, this);
      gibbsVec.push_back(ls);
   }
}


void GibbsGist::initTruthValues()
{
   numTrue.resize(numAtoms);
   numTrueTemp.resize(numAtoms);
   probs.resize(numAtoms, 0);
   for (size_t i = 0; i < numAtoms; i++) {
      numTrue[i].resize(numChains, 0);
      numTrueTemp[i].resize(numChains, 0);
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
         probs[i] += gibbsVec[j]->loc_numTrue[i];
      }
   for (size_t i = 0; i < numAtoms; i++) {
      probs[i] = probs[i] / ((double)numChains * gGla->sample);
      cout << "probability " << i << " = " << probs[i] << endl;
   }
}
