#include "SequenceProposal.h"
#include "Parameter.h"
#include "EnergyModel.h"
#include "Structure.h"
#include "BastollaAugmentedModel.h"
#include <iostream>
#include <stdexcept>
#include <limits>

/* Default constructor, does nothing much.
 */
SequenceProposal::SequenceProposal():
	Proposal()
{
}

/* Constructor from just a protein structure. Implies use of folding
 * energy, rather than binding energy, as probability density.
 */
SequenceProposal::SequenceProposal(vector< shared_ptr<Parameter> > startSequences, shared_ptr<Structure> proteinStructure, shared_ptr<EnergyModel> energyModel):
	Proposal(startSequences)
{
	protStruct = proteinStructure->threadSequence(startSequences[0]->toString());
	engModel = energyModel;

	shared_ptr<Structure> pNothing;
	ligStruct = pNothing;

	try{
		probDensity = calcProbDensity();
	}
	/* On threading failure, assign "negative inifinity" as
	 * a score.
	 */
	catch(runtime_error& e){
		probDensity = -numeric_limits<float>::max();
	}
}

SequenceProposal::SequenceProposal(vector< shared_ptr<Parameter> > startSequences, shared_ptr<Structure> proteinStructure, shared_ptr<EnergyModel> energyModel, shared_ptr<Structure> invariantGeometry):
	Proposal(startSequences)
{
	protStruct = proteinStructure->threadSequence(startSequences[0]->toString());
	engModel = energyModel;

	shared_ptr<Structure> pNothing;
	ligStruct = pNothing;

	geometrySource = invariantGeometry;

	try{
		probDensity = calcProbDensity();
	}
	/* On threading failure, assign "negative inifinity" as
	 * a score.
	 */
	catch(runtime_error& e){
		probDensity = -numeric_limits<float>::max();
	}
}

/* Constructor from protein stucture and ligand structure. Implies use of
 * binding energy as probability density.
 */
SequenceProposal::SequenceProposal(vector< shared_ptr<Parameter> > startSequences, shared_ptr<Structure> proteinStructure, shared_ptr<Structure> ligandStructure, shared_ptr<EnergyModel> energyModel):
	Proposal(startSequences)
{
	protStruct = proteinStructure->threadSequence(startSequences[0]->toString());
	ligStruct = ligandStructure;
	engModel = energyModel;

	try{
		probDensity = calcProbDensity();
	}
	/* On threading failure, assign "negative inifinity" as
	 * a score.
	 */
	catch(runtime_error& e){
		probDensity = -numeric_limits<float>::max();
	}
}

/* Virtual destructor, required for inheritance.
 */
SequenceProposal::~SequenceProposal(){
}

/* Propose a new state from the current one, without a number of
 * mutations (or variation in the same) specified.
 */
shared_ptr<Proposal> SequenceProposal::nextStep(){
	shared_ptr<Proposal> newProp;
	vector< shared_ptr<Parameter> > newParams;
	unsigned int i;
	float stepSize = 1;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->nextValue(stepSize));
	}

	/* Propagate the ligand if it was loaded, ignore it otherwise.
	 */
	if(ligStruct.get() != 0){
		newProp = shared_ptr<Proposal>(new SequenceProposal(newParams, protStruct, ligStruct, engModel));
	}
	else{
		if(geometrySource.get() != 0){
			newProp = shared_ptr<Proposal>(new SequenceProposal(newParams, protStruct, engModel, geometrySource));
		}
		else{
			newProp = shared_ptr<Proposal>(new SequenceProposal(newParams, protStruct, engModel));
		}
	}

	return newProp;
}

/* Propose a new position in sequence space, a fixed number of mutations
 * away from the current one.
 *
 * Currently just takes a step of 1 mutation anyway.
 */
shared_ptr<Proposal> SequenceProposal::nextStep(float stepSize){
	return nextStep();
}

/* Propose a new position in sequence space, with some variable number
 * of differences from the current one.
 *
 * Currently just takes a step of 1 mutation anyway.
 */
shared_ptr<Proposal> SequenceProposal::nextStep(float stepSize, float stepVar){
	return nextStep();
}

/* Propose a random state, somewhere in sequence space.
 *
 * Scale factor meaningless, a really random sequence is returned.
*/
shared_ptr<Proposal> SequenceProposal::randomStep(float scaleFactor){
	shared_ptr<Proposal> newProp;
	vector< shared_ptr<Parameter> > newParams;
	unsigned int i;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->randomValue(scaleFactor));
	}

	/* Propagate the ligand if it was loaded, ignore it otherwise.
	 */
	if(ligStruct.get() != 0){
		newProp = shared_ptr<Proposal>(new SequenceProposal(newParams, protStruct, ligStruct, engModel));
	}
	else{
		if(geometrySource.get() != 0){
			newProp = shared_ptr<Proposal>(new SequenceProposal(newParams, protStruct, engModel, geometrySource));
		}
		else{
			newProp = shared_ptr<Proposal>(new SequenceProposal(newParams, protStruct, engModel));
		}
	}

	return newProp;

}

/* Function for calculating the probability density of
 * the current parameter set. If a ligand was provided at construction
 * time, it's the binding energy; if no ligand was provided, it's the
 * folding energy.
 *
 * Throws runtime_error if there isn't a protein structure to
 * perform calculations on (e.g. threading failed).
 */
float SequenceProposal::calcProbDensity(){
	float eng;

	/* Check for existance of structures before doing anything.
	 */
	if(protStruct.get() == 0){
		string message = "Protein structure doesn't exist! Threading failure?";
		throw runtime_error(message);
	}

	/* Calculate binding energy if the ligand was provided,
	 * folding gap score otherwise. 
	 */
	if(ligStruct.get() != 0){
		eng = abs(engModel->calcInteractionScore(protStruct, ligStruct));
	}
	else{
		if(geometrySource.get() == 0){
			eng = -engModel->calcFoldEnergyGap(protStruct);
		}
		else{
			eng = -engModel->calcFoldEnergyGap(protStruct, geometrySource);
		}
	}

	return eng;
}
