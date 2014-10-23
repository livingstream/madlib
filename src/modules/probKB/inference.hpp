/* ----------------------------------------------------------------------- *//**
 *
 * @file gibbs.hpp
 *
 *//* ----------------------------------------------------------------------- */
/**
 * @brief gibbs: Transition function
 */
DECLARE_UDF(inference, gibbs_step_transition)

/**
 * @brief gibbs: State merge function
 */
DECLARE_UDF(inference, gibbs_step_merge_states)

/**
 * @brief gibbs: Final function
 */
DECLARE_UDF(inference, gibbs_step_final)

/**
 * @brief gibbs: Convert transition state to result tuple
 */
DECLARE_UDF(inference, internal_gibbs_result)
