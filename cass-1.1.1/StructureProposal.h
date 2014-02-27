/* Class to make proposals about angular (???) parametrizations
 * for the state of residues in a protein structure.
 */

/* Idempotency.
 */
#ifndef STRUCTUREPROPOSAL_H
#define STRUCTUREPROPOSAL_H

/* Forward declarations.
 */
class TwoBeadStructure;

/* Includes.
 */
#include "Proposal.h"

/* Definitions.
 */
#define USE_SIMPLE_ENERGY_FUNCTION false

/* Interface.
 */
class StructureProposal : public Proposal{

	public:
		/* Constructors and destructors.
		 */
		StructureProposal();
		StructureProposal(shared_ptr<TwoBeadStructure> inputStruct, string newSequence, float neighborDistance);
		StructureProposal(shared_ptr<TwoBeadStructure> inputStruct, vector<bool> movingResidues);

		virtual ~StructureProposal();

		/* Propose a new state from the current one
		 * with varying-size steps.
	 	*/
		shared_ptr<Proposal> nextStep(float stepSize);

		/* Propose a new state from the current one with
		 * step sizes varying with some variance.
		 */
		shared_ptr<Proposal> nextStep(float stepSize, float stepVar);

		/* Propose a random state within parameter space.
		 */
		shared_ptr<Proposal> randomStep();

		/* Returns the current full Structure.
		 */
		shared_ptr<TwoBeadStructure> getStructure() const;

	private:

		/* Constructor from list of parameters. 	
		 */
		StructureProposal(vector< shared_ptr<Parameter> > startParams, shared_ptr<TwoBeadStructure> scaffold, vector<bool> activeResidues, string novelSequence);

		/* Function that updates the old structure with current
		 * angle values.
		 */
		void updateProteinStructure();

		/* Function for calculating the probability density, i.e.
		 * SOME_METRIC, for the current parameter set.
		 *
		 * METRIC WILL POSSIBLY CHANGE. L-J VS LINEAR ENERGY TERMS?
		 * CAPPED OR NOT?
		 */
		float calcProbDensity();
		float calcProbDensity(bool simpleMethod);

		/* Protein structures.
		 */
		shared_ptr<TwoBeadStructure> pStartStruct;
		shared_ptr<TwoBeadStructure> pCurrentStruct; 

		/* Sequence to aspire to.
		 */
		string targetSeq;

		/* List of residues that can move around.
		 */
		vector<bool> active;
};

/* Idempotency end.
 */
#endif
