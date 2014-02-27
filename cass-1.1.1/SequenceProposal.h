/* Object to hold a particular proposed state (sequence + folding/binding
 * energy) during an MCMC search of sequence space.
 */

/* Idempotency.
 */
#ifndef SEQUENCEPROPOSAL_H
#define SEQUENCEPROPOSAL_H

/* Forward declarations.
 */
class EnergyModel;
class Structure;

/* Includes.
 */
#include "Proposal.h"

/* Interface.
 */
class SequenceProposal : public Proposal{
	public:
		SequenceProposal();
		SequenceProposal(vector< shared_ptr<Parameter> > startSequences, shared_ptr<Structure> proteinStructure, shared_ptr<EnergyModel> energyModel);
		SequenceProposal(vector< shared_ptr<Parameter> > startSequences, shared_ptr<Structure> proteinStructure, shared_ptr<Structure> ligandStructure, shared_ptr<EnergyModel> energyModel);
		virtual ~SequenceProposal();

		SequenceProposal(vector< shared_ptr<Parameter> > startSequences, shared_ptr<Structure> proteinStructure, shared_ptr<EnergyModel> energyModel, shared_ptr<Structure> invariantGeometry);

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
		shared_ptr<Structure> protStruct;
		shared_ptr<Structure> ligStruct;
		shared_ptr<Structure> geometrySource;
		shared_ptr<EnergyModel> engModel;

		/* Function for calculating the probability density of
		 * the current parameter set.
		 */
		float calcProbDensity(); 
};

/* Idempotency end.
 */
#endif
