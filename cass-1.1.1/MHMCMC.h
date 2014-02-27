/* Class to run the "vanilla" Metropolis-Hastings algorithm
 * in an MCMC framework. Step size is constant (temperature
 * does not e.g. decrease throughout), acceptance function is 
 * the commonly used
 *
 * P(Pd(s), Pd(s'), T) = |---- 1, if Pd(s') > Pd(s)
 *                       |
 *		         |---- exp(-(Pd(s) - Pd(s'))/T)
 *
 * which maximizes the function Pd(s).
 */

/* Idempotency.
 */
#ifndef MHMCMC_H
#define MHMCMC_H

/* Forward declarations.
 */
class CRandomMersenne;

/* Includes.
 */
#include "MCMC.h"

class MHMCMC : public MCMC{
	public:
		MHMCMC();
		MHMCMC(shared_ptr<Proposal> startState, float startTemp, int numSteps);
		virtual ~MHMCMC();
	
		/* Function to take a single MCMC step with tunable
		 * variance in step sizes.
		 */
		void takeStep(float stepSize, float stepVar);

		static void setPRNG(CRandomMersenne &rangeGenerator);

	protected:
		/* Function to evaluate moving to a proposed new state
		 * (acceptance function).
		 */
		virtual bool moveToState(shared_ptr<Proposal> newState);

		static shared_ptr<CRandomMersenne> pRandGen;

};

/* Idempotency end.
 */
#endif
