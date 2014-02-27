/* Includes.
 */
#include "MultipleStructureOrganism.h"
#include "Structure.h"
#include "EnergyModel.h"
#include <iostream>
#include <limits>

/* Default constructor, doesn't do much.
 */
MultipleStructureOrganism::MultipleStructureOrganism():
	Organism()
{
}

/* Constructor from values, calls upward a lot.
 */
MultipleStructureOrganism::MultipleStructureOrganism(float mutationRate, vector<string> codingSequences, vector< shared_ptr<Structure> > nativeStructures, vector< shared_ptr<Structure> > nativeLigands, vector< shared_ptr<Structure> > decoyLigands, shared_ptr<EnergyModel> engModel, float foldEnergyThreshold, float bindEnergyThreshold, shared_ptr<Structure> geometrySample):
	Organism(mutationRate, codingSequences, nativeStructures, engModel, foldEnergyThreshold)
{
	bindThreshold = bindEnergyThreshold;
	nativeLigs = nativeLigands;
	decoyLigs = decoyLigands;
	invariant = geometrySample;
}

/* Virtual destructor, required by gcc.
 */
MultipleStructureOrganism::~MultipleStructureOrganism(){
}

/* Copy constructor with possibilty of mutation.
 */
MultipleStructureOrganism::MultipleStructureOrganism(const MultipleStructureOrganism &org, bool evolve):
	Organism(org, evolve)
{
	/* Copy data members not taken care of in the superclass.
	 */
	bindThreshold = org.bindThreshold;
	nativeLigs = org.nativeLigs;
	decoyLigs = org.decoyLigs;
	invariant = org.invariant;
}

/* Copy constructor with possibility of mutation, if using
 * polymorphism.
 */
shared_ptr<Organism> MultipleStructureOrganism::Clone(bool evolve){
	return shared_ptr<Organism>(new MultipleStructureOrganism(*this, evolve));
}

/* Calculate all the folding energies, but avoid re-calculating
 * for known structures.
 */
vector<float> MultipleStructureOrganism::calcFoldingEnergies(){
	vector<float> engs; 
	unsigned int i;

	/* Go over all the native structures.
	 */
	for(i = 0; i < structs.size(); i++){
		/* Did we have a threading failure? 
		 */
		if(structs[i].get() == 0){
			/* DEBUG
			 *
			cerr << "Threading failure! No energy calculated." << endl;
			*/
			
			engs.push_back(numeric_limits<float>::max());
		}
		else{
			/* DEBUG
			 */


			/* This isn't a previously calculated energy.
			 */
			if(Organism::knownFoldEnergies[i].find(structs[i]->getSequence()) == Organism::knownFoldEnergies[i].end()){
				/* Use the specified energy model to estimate the folding energy.
				 */
				engs.push_back( pEngModel->calcFoldEnergyGap(structs[i], invariant) );

				/* Store seq-eng combo in appropriate part of folding energy
				 * container.
				 */
				Organism::knownFoldEnergies[i].insert(stringAndFloat(structs[i]->getSequence(),engs[i]));

				/* DEBUG
				 *
				 cerr << Organism::knownFoldEnergies[i].size() << " sequence energies now calculated." << endl;
				 */
			}
			/* We've done this already, retrieve the known energy.
			 */
			else{
				engs.push_back( Organism::knownFoldEnergies[i].find(structs[i]->getSequence())->second );
				
				/* DEBUG
				 *
				cerr << "Known folding energy for #" << i+1 << ": " << engs[i] << endl;
				*/
			}
		}
	}

	return engs;
}

/* Calculate all the native and decoy binding energies for
 * each native protein conformation.
 *
 * For N proteins with K native ligands and M decoy ligands, the 
 * output vector is organized as follows:
 *
 * 1..K native engs for protein 1, 1..M decoy engs for protein 1,
 * 1..K native engs for protein 2, 1..M decoy engs for protein 2,
 * ...,
 * 1..K native engs for protein N, 1..M decoy engs for protein N,
 * in that order.
 */
vector<float> MultipleStructureOrganism::calcBindingEnergies(){
	vector<float> engs;
	unsigned int i, j;
	unsigned int numValues;

	for(i = 0; i < structs.size(); i++){

		/* Make sure there wasn't a threading
		 * failure somewhere.
		 */
		if(structs[i].get() != 0){

			/* First the native ligands.
			*/
			for(j = 0; j < nativeLigs.size(); j++){
				engs.push_back( pEngModel->calcInteractionScore(structs[i], nativeLigs[j]) );
			}

			/* Then the decoy ones.
			*/
			for(j = 0; j < decoyLigs.size(); j++){
				engs.push_back( pEngModel->calcInteractionScore(structs[i], decoyLigs[j]) );
			}
		}
		else{
			numValues = nativeLigs.size()+decoyLigs.size();
			vector<float> maxes(numValues, numeric_limits<float>::max());
			engs.insert(engs.end(), maxes.begin(), maxes.end());
		}
	}

	return engs;
}

/* Calculate overall fitness. See header file for fitness definitions.
 */
void MultipleStructureOrganism::calcFitnessAndPhenotype(){

	unsigned int i,j, offset;
	float energy;

	vector<float> foldEngs = calcFoldingEnergies();
	vector<float> bindEngs = calcBindingEnergies();
	fitness = 0.5;
	phenotype = "Normal";

	/* DEBUG
	 *
	cerr << "Folding energies:";
	for(j = 0; j < foldEngs.size(); j++){
		cerr << "\t" << foldEngs[j];
	}
	cerr << endl;
	cerr << "Binding energies:";
	for(j = 0; j < bindEngs.size(); j++){
		cerr << "\t" << bindEngs[j];
	}
	cerr << endl;
	*/

	/* Check if native interactions are still maintained,
	 * and if there are any additional ones from decoy
	 * ligands.
	 */
	for(i = 0; i < structs.size(); i++){
		/* Folding failure, all interactions lost, fitness
		 * wiped out.
		 */
		if(foldEngs[i] > foldThreshold){

			/* DEBUG
			 *
			cerr << "Folding failure for protein " << (i+1) <<  " (" << foldEngs[i] <<  " > " << foldThreshold << ")." << endl;
			*/

			fitness = 0.0;
			phenotype = "Dead";
			return;
		}
		else{
			/* How are the native binding interactions?
			 */
			offset = i*(nativeLigs.size()+decoyLigs.size());
			for(j = 0; j < nativeLigs.size(); j++){
				energy = bindEngs[offset+j]; //?

				/* Binding failure, no fitness.
				 */
				if(energy > bindThreshold){

					/* DEBUG
					 *
					cerr << "Binding failure for ligand " << (j+1) << " (score " << energy << " > " << bindThreshold << ")." << endl;
					*/

					fitness = 0.0;
					phenotype = "Dead";
					return;
				}
			}

			/* Any decoys bind, too?
			 */
			offset = i*(nativeLigs.size()+decoyLigs.size())+nativeLigs.size();
			for(j = 0; j < decoyLigs.size(); j++){
				energy = bindEngs[offset+j]; //??

				// CURRENTLY 10% FITNESS PENALTY PER DECOY
				if(energy < bindThreshold){
					/* DEBUG
					 *
					cerr << "Binding decoy " << (j+1) << " (score " << energy << ")." << endl;
					*/

					fitness *= 0.9;
					phenotype = "Decoy";
				}
			}
		}
	}
}
