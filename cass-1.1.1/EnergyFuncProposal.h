/* Class to make proposals about parameterizations for folding/binding
 * energy functions under my simulation framework.
 */

/* Idempotency.
 */
#ifndef ENERGYFUNCPROPOSAL_H
#define ENERGYFUNCPROPOSAL_H

/* Forward declarations.
 */

/* Includes.
 */
#include "Proposal.h"

class EnergyFuncProposal : public Proposal{
	public:
		EnergyFuncProposal();
		EnergyFuncProposal(vector< shared_ptr<Parameter> > startParams, vector<float> nativeSeqTermValues, vector< vector<float> > decoyTermValues, vector< vector<float> > randomSeqTermValues);
		virtual ~EnergyFuncProposal();

		/* Propose a new state from the current one
		 * with varying-size steps.
	 	*/
		shared_ptr<Proposal> nextStep(float stepSize);

		/* Propose a new state from the current one with
		 * step sizes varying with some variance.
		 */
		shared_ptr<Proposal> nextStep(float stepSize, float stepVar);

		/* Propose a random state based on some
		 * scaling of the parameter space.
		 */
		shared_ptr<Proposal> randomStep(float scaleFactor);

	private:
		/* Function for calculating the probability density, i.e.
		 * SOME_METRIC, for the current parameter set.
		 *
		 * METRIC WILL POSSIBLY CHANGE.
		 */
		float calcProbDensity(); 

		/* Unscaled values for energy terms for the native
		 * sequence in this conformation.
		 */
		vector<float> natSeqTermVals;

		/* Unscaled values for energy terms for the native
		 * sequence in one or more decoy conformations.
		 */
		vector< vector<float> > decoyTermVals;

		/* Unscaled values of energy terms for some number of 
		 * random sequences
		 */
		vector< vector<float> > randSeqTermVals;
};

/* Idempotency end.
 */
#endif
