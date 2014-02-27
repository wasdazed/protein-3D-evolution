/* Object to hold information about a proposed state
 * (ANY STATE?) in a Markov Chain Monte Carlo (MCMC)
 * process.
 */

/* Make source file idempotent (avoid multiple includes).
 */
#ifndef PROPOSAL_H
#define PROPOSAL_H

/* Forward declarations.
 */
class Parameter;

/* Includes
 */
#include <string>
#include <vector>
#include <tr1/memory>

using namespace std;

/* Interface.
 */
class Proposal{
	public:
		Proposal();
		Proposal(vector< shared_ptr<Parameter> > startParams);
		virtual ~Proposal();

		/* Propose a new state from the current one
		 * (with or without a specified variance).
		 */
		virtual shared_ptr<Proposal> nextStep(float stepSize);
		virtual shared_ptr<Proposal> nextStep(float stepSize, float stepVar);

		/* Propose a random state.
		 */
		virtual shared_ptr<Proposal> randomStep(float scaleFactor);

		/* Relevant properties for display.
		 */
		string toString();

		/* Accessor method for the probability
		 * density.
		 */
		float getProbDensity() const;

		/* Acessor method for floating-point representation
		 * of parameter values.
		 */
		vector<float> getFloatParams() const;

	protected:
		/* Parameters and probability density of that
		 * parameter set.
		 */
		vector< shared_ptr<Parameter> > params;
		float probDensity;

		/* Function for calculating the probability density of
		 * the current parameter set.
		 */
		virtual float calcProbDensity(); 

	private:
};

/* End idempotency statement.
 */
#endif
