#include "gibbs_gist.h"
#include "gibbs_cgla.h"

GibbsCGLA::GibbsCGLA(size_t numAtoms_, size_t numChains_, GibbsGist *state_)
{
   numAtoms = numAtoms_;
   numChains = numChains_;
   done = false;
   sample = 0;
   secondsElapsed = 0;
   startTimeSec = timer.time();
   currentTimeSec = 0;
   para = new GibbsParams();
   state = state_;
   burnIn = (para->burnMaxSteps > 0) ? true : false;
   initConvergenceTests();
}

GibbsCGLA::~GibbsCGLA()
{
   deleteConvergenceTests();
}

/*void GibbsCGLA::merge(GibbsCGLA *cgla)
{
}*/

bool GibbsCGLA::shouldIterate()
{
   //cout << "sample = " << sample << endl;
   sample += para->samplesPerTest;
   currentTimeSec = timer.time();
   secondsElapsed = currentTimeSec - startTimeSec;

   for(size_t i = 0; i < numAtoms; i++) {
       for (size_t j = 0; j < numChains; j++) {
	    state->numTrueTemp[i][j] = state->gibbsVec[j]->loc_numTrueTemp[i];
            state->gibbsVec[j]->loc_numTrueTemp[i] = 0;
       }
   }

   // Add current truth values to the convergence testers
   for (size_t i = 0; i < numAtoms; i++) {
      //WARNING: implicit cast from bool* to double*
      if (burnIn) {
         burnConvergenceTests[i]->appendNewValues(state->numTrueTemp[i], para->samplesPerTest);
      } else {
         gibbsConvergenceTests[i]->appendNewValues(state->numTrueTemp[i], para->samplesPerTest);
      }
   }

   if (burnIn) {
      bool burnConverged = GelmanConvergenceTest::checkConvergenceOfAll(
                              burnConvergenceTests, (int)numAtoms);
      if ((sample >= para->burnMinSteps && burnConverged)
            || (para->burnMaxSteps >= 0 && sample >= para->burnMaxSteps)
            || (para->maxSeconds > 0 && secondsElapsed >= para->maxSeconds)) {
         burnIn = false;
         sample = 0;
      }
   } else {
      bool gibbsConverged = ConvergenceTest::checkConvergenceOfAtLeast(
                               gibbsConvergenceTests, (int)numAtoms, sample, para->fracConverged);
      if ((sample >= para->minSteps && gibbsConverged)
            || (para->maxSteps >= 0 && sample >= para->maxSteps)
            || (para->maxSeconds > 0 && secondsElapsed >= para->maxSeconds)) {
         return false;
      }
   }

   return true;
}

void GibbsCGLA::initConvergenceTests()
{
   burnConvergenceTests = new GelmanConvergenceTest*[numAtoms];
   gibbsConvergenceTests = new ConvergenceTest*[numAtoms];
   for (size_t i = 0; i < numAtoms; i++) {
      burnConvergenceTests[i]  = new GelmanConvergenceTest((int)numChains);
      gibbsConvergenceTests[i] = new ConvergenceTest((int)numChains, para->gamma, para->epsilonError);
   }
}

void GibbsCGLA::deleteConvergenceTests()
{
   for (size_t i = 0; i < numAtoms; i++) {
      delete burnConvergenceTests[i];
      delete gibbsConvergenceTests[i];
   }
   delete [] burnConvergenceTests;
   delete [] gibbsConvergenceTests;
}
