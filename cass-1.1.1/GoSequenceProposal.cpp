#include "GoSequenceProposal.h"
#include "Parameter.h"
#include "myUtils.h"
#include <iostream>

/* Default constructor, does nothing much.
 */
GoSequenceProposal::GoSequenceProposal():
	Proposal()
{
}

/* Constructor from just a protein structure. Implies use of folding
 * energy, rather than binding energy, as probability density.
 */
GoSequenceProposal::GoSequenceProposal(vector< shared_ptr<Parameter> > startSequences, string nativeSeq, vector<float> nativeResidueContribs):
	Proposal(startSequences)
{
		nativeSequence = nativeSeq;		
		nativeResidueEngs = nativeResidueContribs; 	
		probDensity = calcProbDensity();
}


/* Virtual destructor, required for inheritance.
 */
GoSequenceProposal::~GoSequenceProposal(){
}

/* Propose a new state from the current one, without a number of
 * mutations (or variation in the same) specified.
 */
shared_ptr<Proposal> GoSequenceProposal::nextStep(){
	shared_ptr<Proposal> newProp;
	vector< shared_ptr<Parameter> > newParams;
	unsigned int i;
	float stepSize = 1;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->nextValue(stepSize));
	}

	/* Propagate the new parameters.
	 */
	newProp = shared_ptr<Proposal>(new GoSequenceProposal(newParams, nativeSequence, nativeResidueEngs));

	return newProp;
}

/* Propose a new position in sequence space, a fixed number of mutations
 * away from the current one.
 *
 * Currently just takes a step of 1 anyway.
 */
shared_ptr<Proposal> GoSequenceProposal::nextStep(float stepSize){
	return nextStep();
}

/* Propose a new position in sequence space, with some variable number
 * of differences from the current one.
 *
 * Currently just takes a step of 1 anyway.
 */
shared_ptr<Proposal> GoSequenceProposal::nextStep(float stepSize, float stepVar){
	return nextStep();
}

/* Propose a random state, somewhere in sequence space.
 *
 * Scale factor meaningless, a really random sequence is returned.
*/
shared_ptr<Proposal> GoSequenceProposal::randomStep(float scaleFactor){
	shared_ptr<Proposal> newProp;
	vector< shared_ptr<Parameter> > newParams;
	unsigned int i;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->randomValue(scaleFactor));
	}

	/* Propagate the ligand if it was loaded, ignore it otherwise.
	 */
	newProp = shared_ptr<Proposal>(new GoSequenceProposal(newParams, nativeSequence, nativeResidueEngs));

	return newProp;

}

/* Function for calculating the probability density of
 * the current parameter set, under a frustrated Go model of sequence
 * space.
 *
 * Current state: experimental.
 */
float GoSequenceProposal::calcProbDensity(){
	float eng = 0.0;
	unsigned int i;
	float nearNatWeight = 0.5; // Scales contribution from near-native substitutions

	string mutSeq = params[0]->toString();

	for(i = 0; i < nativeSequence.length(); i++){
		if(mutSeq.substr(i,1).compare(nativeSequence.substr(i,1)) == 0){
			eng += nativeResidueEngs[i];
		}
		else if(sameBiochemClass(mutSeq[i], nativeSequence[i])){
			/* DEBUG
			 *
			cerr << "Near native state (" << mutSeq[i] << " vs " << nativeSequence[i] << "), getting half-weight." << endl;
			*/

			eng += nearNatWeight*nativeResidueEngs[i];
		}
	}
	
	eng *= -1.0; // We're maximizing the probability density, not a free energy measure

	return eng;
}
