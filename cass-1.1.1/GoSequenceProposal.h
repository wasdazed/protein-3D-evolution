/* Object to hold a particular proposed state (sequence + Go-model
 * energy) during an MCMC search of sequence space.
 */

/* Idempotency.
 */
#ifndef GOSEQUENCEPROPOSAL_H
#define GOSEQUENCEPROPOSAL_H

/* Forward declarations.
 */

/* Includes.
 */
#include "Proposal.h"

/* Interface.
 */
class GoSequenceProposal : public Proposal{
	public:
		GoSequenceProposal();
		GoSequenceProposal(vector< shared_ptr<Parameter> > startSequences, string nativeSequence, vector<float> nativeResidueContribs);
		virtual ~GoSequenceProposal();

		/* Propose a new state from the current one
		 * (with or without a specified variance).
		 */
		shared_ptr<Proposal> nextStep();
		shared_ptr<Proposal> nextStep(float stepSize);
		shared_ptr<Proposal> nextStep(float stepSize, float stepVar);

		/* Propose a random state.
		 */
		shared_ptr<Proposal> randomStep(float scaleFactor);

	protected:
	
		/* Parameters and probability density of that
		 * parameter set.
		 */
		vector<float> nativeResidueEngs;
		string nativeSequence;

		/* Function for calculating the probability density of
		 * the current parameter set.
		 */
		float calcProbDensity(); 
};

/* Idempotency end.
 */
#endif
