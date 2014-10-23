#ifndef METROPOLISPARAMS_H
#define METROPOLISPARAMS_H

/**
 * This struct holds parameters needed to run Gibbs sampling.
 *
 * @see MCMCParams
 */
struct MetropolisParams {
    MetropolisParams() {
        numChains = 10;
        burnMinSteps = 100;
        burnMaxSteps = 100;
        minSteps = -1;
        maxSteps = 1000;
        maxSeconds = -1;
        gamma = 1 - 0.05;
        epsilonError = 0.01;
        fracConverged = 0.95;
        samplesPerTest = 100;
    }
    // No. of chains which MCMC will use
    int numChains;
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

    double gamma;
    double epsilonError;
    double fracConverged;
    int    samplesPerTest;
};

#endif /*METROPOLISPARAMS_H*/
