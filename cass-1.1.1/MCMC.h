/* Object to hold information about and perform a Markov Chain
 * Monte Carlo process for parameter optimization.
 *
 * IS A LITTLE BACKWARDS IN THAT IT'S TRYING TO MINIMIZE
 * PROBABILITY DENSITIES INSTEAD OF MAXIMIZING: THEY'RE
 * ENERGIES IN MY CASE.
 */

/* Idempotency.
 */
#ifndef MCMC_H
#define MCMC_H

/* Forward declarations.
 */
class Proposal;

/* Includes.
 */
#include <tr1/memory>

using namespace std;

/* Interface.
 */
class MCMC {
	public:
		MCMC();
		MCMC(shared_ptr<Proposal> startState, float startTemp, int numSteps);
		virtual ~MCMC();

		/* Making a single step.
		 *
		 * SHOULD ALSO DECREASE TEMPERATURE WHEN USING SIMULATED
		 * ANNEALING (HENCE THE VIRTUALIZATION).
		 */
		virtual void takeStep();

		/* Convenience function to do a full simulation at once.
		 */
		void run();

		/* Accessor methods.
		 */
		shared_ptr<Proposal> getCurrentState() const;
		shared_ptr<Proposal> getBestState() const;

	protected:
		
		/* Function to evaluate moving to a proposed new state
		 * (acceptance function).
		 */
		virtual bool moveToState(shared_ptr<Proposal> newState);

		/* Current position in parameter space.
		 */
		shared_ptr<Proposal> currentState;

		/* Best previous position in parameter space.
		 */
		shared_ptr<Proposal> bestState;

		/* Simulation parameters.
		 */
		float temp;
		int currentStep;
		int maxSteps;
		
	private:
};

/* Idempotency end.
 */
#endif
