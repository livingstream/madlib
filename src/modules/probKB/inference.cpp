/* ----------------------------------------------------------------------- *//**
 *
 * @file inference.cpp
 *
 *//* ----------------------------------------------------------------------- */
#define gibbs true /* Gibbs or Metropolis-Hasting */
#define parallel true /* parallel or sequential */

#include <iostream>
#include <dbconnector/dbconnector.hpp>
#include <modules/shared/HandleTraits.hpp>
#include <modules/probKB/infer/state/variablestate.h>

#if gibbs
    #if parallel
        #include <modules/probKB/infer/gibbs/parallel/gibbs_gist.h>
    #else 
        #include <modules/probKB/infer/gibbs/sequential/gibbs.cpp>
    #endif
#else 
    #if parallel
        #include <modules/probKB/infer/metropolisHasting/parallel/metropolis-hasting_gist.h>
    #else 
        #include <modules/probKB/infer/metropolisHasting/sequential/metropolis-hasting.cpp>
    #endif
#endif

#include "inference.hpp"
#include <unistd.h>
#include <map>
#include <stdlib.h>

namespace madlib {

namespace modules {

// Import names from other MADlib modules
using dbal::NoSolutionFoundException;

namespace inference {

// Use Eigen
using namespace dbal::eigen_integration;

// Internal functions
AnyType stateToResult(const Allocator &inAllocator,
                      const HandleMap < const ColumnVector,
                      TransparentHandle<double> > &intruth);

template <class Handle>
class inferenceTransitionState {
    template <class OtherHandle>
    friend class inferenceTransitionState;

public:
    inferenceTransitionState(const AnyType &inArray)
        : mStorage(inArray.getAs<Handle>()) {
        rebind(static_cast<uint32_t>(mStorage[2]), static_cast<uint32_t>(mStorage[3]));
    }

    /**
     * @brief Convert to backend representation
     *
     * We define this function so that we can use State in the
     * argument list and as a return type.
     */
    inline operator AnyType() const {
        return mStorage;
    }

    inline void initialize(const Allocator &inAllocator, uint32_t inNumAtoms, uint32_t inSize, int32_t inQid) {
        mStorage = inAllocator.allocateArray < double, dbal::AggregateContext,
        dbal::DoZero, dbal::ThrowBadAlloc > (arraySize(inNumAtoms, inSize));
        rebind(inNumAtoms, inSize);
        numAtoms = inNumAtoms;
        size = inSize;
        qid = inQid;
    }

    /**
     * @brief We need to support assigning the previous state
     */
    template <class OtherHandle>
    inferenceTransitionState &operator=(
        const inferenceTransitionState<OtherHandle> &inOtherState) {
        for (size_t i = 0; i < mStorage.size(); i++) {
            mStorage[i] = inOtherState.mStorage[i];
        }
        return *this;
    }

    /**
     * @brief Merge with another State object by copying the intra-iteration
     * fields
     */
    template <class OtherHandle>
    inferenceTransitionState &operator+=(
        const inferenceTransitionState<OtherHandle> &inOtherState) {
        if (mStorage.size() != inOtherState.mStorage.size()) {
            throw std::logic_error("Internal error: Incompatible transition states");
        }
        numRows += inOtherState.numRows;
        for (int i = 0; i < inOtherState.curSize; i++) {
            clauses[curSize + i] = inOtherState.clauses[i];
        }
        curSize += inOtherState.curSize;
        return *this;
    }

private:
    static inline uint32_t arraySize(const uint32_t numAtoms, const uint32_t size) {
        return 5 + numAtoms + size;
    }

    void rebind(uint32_t inNumAtoms, uint32_t inSize) {
        numRows.rebind(&mStorage[0]);
        curSize.rebind(&mStorage[1]);
        numAtoms.rebind(&mStorage[2]);
        size.rebind(&mStorage[3]);
        qid.rebind(&mStorage[4]);
        truth.rebind(&mStorage[5], inNumAtoms);
        clauses.rebind(&mStorage[5 + inNumAtoms], inSize);
    }
    Handle mStorage;

public:
    typename HandleTraits<Handle>::ReferenceToUInt64 numRows;
    typename HandleTraits<Handle>::ReferenceToInt64 curSize;
    typename HandleTraits<Handle>::ReferenceToUInt64 numAtoms;
    typename HandleTraits<Handle>::ReferenceToUInt64 size;
    typename HandleTraits<Handle>::ReferenceToInt64 qid;
    typename HandleTraits<Handle>::ColumnVectorTransparentHandleMap truth;
    typename HandleTraits<Handle>::ColumnVectorTransparentHandleMap clauses;
};


/**
 * @brief Compute the log likelihood and gradient vector for each tuple
 */
AnyType
gibbs_step_transition::run(AnyType &args)
{
    inferenceTransitionState<MutableArrayHandle<double> > state = args[0];
    MappedColumnVector clause = args[1].getAs<MappedColumnVector>();
    int clauseSize = static_cast<int>(clause.size());
    double weight = args[2].getAs<double>();
    int component = static_cast<uint32_t>(args[3].getAs<double>());
    if (state.numRows == 0) {
        state.initialize(*this, static_cast<uint32_t>(args[4].getAs<double>()),
                         static_cast<uint32_t>(args[5].getAs<double>()), static_cast<int32_t>(args[6].getAs<double>()));
    }
    /*if(state.numRows == 0) {
       s << state.curSize << "|" << state.clauses.size() <<"|"<< state.numAtoms;
       throw std::logic_error("error: transition states" + s.str());
    }*/
    state.numRows++;
    state.clauses[state.curSize] = clauseSize;
    state.clauses[state.curSize + 1] = weight;
    state.clauses[state.curSize + 2] = component;
    state.curSize += 3;
    for (int i = 0; i < clauseSize; i++) {
        state.clauses[state.curSize + i] = clause[i];
    }
    state.curSize += clauseSize;
    return state;
}

/**
 * @brief Perform the perliminary aggregation function: Merge transition states
 */
AnyType
gibbs_step_merge_states::run(AnyType &args)
{
    inferenceTransitionState<MutableArrayHandle<double> > stateLeft = args[0];
    inferenceTransitionState<ArrayHandle<double> > stateRight = args[1];
    // We first handle the trivial case where this function is called with one
    // of the states being the initial state
    if (stateLeft.numRows == 0) {
        return stateRight;
    } else if (stateRight.numRows == 0) {
        return stateLeft;
    }

    // Merge states together and return
    stateLeft += stateRight;
    return stateLeft;
}


AnyType
gibbs_step_final::run(AnyType &args)
{
    inferenceTransitionState<MutableArrayHandle<double> > state = args[0];
    if (state.numRows == 0) {
        return Null();
    }
    size_t i = 0;
    size_t numClauses = 0;
    std::map<int, int> indexMap;

    VariableState *varState = new VariableState(state.numAtoms);
    while (i < static_cast<size_t>(state.clauses.size())) {
        GroundClause *gc = new GroundClause();
        size_t clauseSize =  static_cast<size_t>(state.clauses[i]);
        gc->wt_ = state.clauses[i + 1];
        gc->component = static_cast<size_t>(state.clauses[i + 2]);
        i += 3;
        if (state.qid < 0) {
            for (size_t j = 0; j < clauseSize; j++) {
                gc->gndPreds.push_back(static_cast<int>(state.clauses[i + j]));
            }
        } else {
            for (size_t j = 0; j < clauseSize; j++) {
                int oriId = static_cast<int>(state.clauses[i + j]);
                if (indexMap.count(abs(oriId)) == 0) {
                    int newId = (int)indexMap.size() + 1;
                    indexMap[abs(oriId)] = newId;
                }
                gc->gndPreds.push_back(indexMap[abs(oriId)] * (oriId / abs(oriId)));
            }
        }
        i += clauseSize;
        numClauses += 1;
        varState->gndClauses_->push_back(gc);
    }

    varState->init();
    #if gibbs 
       #if parallel 
           GibbsGist instance(state.numAtoms, numClauses, varState);
       #else
           gibbs instance(state.numAtoms, numClauses, varState);
       #endif
    #else
       #if parallel 
           MetropolisGist instance(state.qid, state.numAtoms, numClauses, varState);
       #else
           Metropolis instance(state.qid, state.numAtoms, numClauses, varState);
       #endif
    #endif
    instance.infer();
    if (state.qid < 0) {
        for (i = 0; i < state.numAtoms; i++) {
            #if parallel
               state.truth[i] = instance.probs[i];
            #else
               state.truth[i] = instance.numTrue[i];
            #endif
        }
    } else {
        #if parallel
           state.truth[0] = instance.probs[indexMap[(int)state.qid] - 1];
        #else
           state.truth[0] = instance.numTrue[indexMap[(int)state.qid] - 1];
        #endif
    }
    return state;
}

AnyType
internal_gibbs_result::run(AnyType &args)
{
    inferenceTransitionState<ArrayHandle<double> > state = args[0];
    return stateToResult(*this, state.truth);
}

AnyType stateToResult(
    const Allocator &inAllocator,
    const HandleMap<const ColumnVector, TransparentHandle<double> > &intruth)
{
    MutableNativeColumnVector truth(
        inAllocator.allocateArray<double>(intruth.size()));
    truth = intruth;
    AnyType tuple;
    tuple << truth;
    return tuple;
}

} // namespace inference

} // namespace modules

} // namespace madlib
