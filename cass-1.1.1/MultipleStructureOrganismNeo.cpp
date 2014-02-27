/* Includes.
 */
#include "MultipleStructureOrganismNeo.h"
#include "Structure.h"
#include "EnergyModel.h"
#include <iostream>
#include <limits>

/* Default constructor, doesn't do much.
 */
MultipleStructureOrganismNeo::MultipleStructureOrganismNeo():
	MultipleStructureOrganism()
{
}

/* Constructor from values, calls upward a lot.
 */
MultipleStructureOrganismNeo::MultipleStructureOrganismNeo(float mutationRate, vector<string> codingSequences, vector< shared_ptr<Structure> > nativeStructures, vector< shared_ptr<Structure> > nativeLigands, vector< shared_ptr<Structure> > decoyLigands, vector< shared_ptr<Structure> > novelLigands, shared_ptr<EnergyModel> engModel, float foldEnergyThreshold, float bindEnergyThreshold, shared_ptr<Structure> geometrySample):
	MultipleStructureOrganism(mutationRate, codingSequences, nativeStructures, nativeLigands, decoyLigands, engModel, foldEnergyThreshold, bindEnergyThreshold, geometrySample)
{
	neoLigs = novelLigands;
}

/* Virtual destructor, required by gcc.
 */
MultipleStructureOrganismNeo::~MultipleStructureOrganismNeo(){
}

/* Copy constructor with possibilty of mutation.
 */
MultipleStructureOrganismNeo::MultipleStructureOrganismNeo(const MultipleStructureOrganismNeo &org, bool evolve):
	MultipleStructureOrganism(org, evolve)
{
	/* Copy data members not taken care of in the superclass.
	 */
	neoLigs = org.neoLigs;
}

/* Copy constructor with possibility of mutation, if using
 * polymorphism.
 */
shared_ptr<Organism> MultipleStructureOrganismNeo::Clone(bool evolve){
	return shared_ptr<Organism>(new MultipleStructureOrganismNeo(*this, evolve));
}

/* Calculate all the native, decoy and novel binding energies for
 * each native protein conformation.
 *
 * For N proteins with K native ligands, M decoy ligands, and L novel ligands, 
 * the output vector is organized as follows:
 *
 * 1..K native engs for protein 1, 1..M decoy engs for protein 1, 1..L novel engs for protein 1,
 * 1..K native engs for protein 2, 1..M decoy engs for protein 2, 1..L novel engs for protein 2,
 * ...,
 * 1..K native engs for protein N, 1..M decoy engs for protein N, 1..L novel engs for protein N
 *
 * in that order. 
 *
 * On threading failure for a structure, all binding energies are FLOAT_MAX.
 */
vector<float> MultipleStructureOrganismNeo::calcBindingEnergies(){
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

			/* And finally the novel ones.
			 */
			for(j = 0; j < neoLigs.size(); j++){
				engs.push_back( pEngModel->calcInteractionScore(structs[i], neoLigs[j]) );
			}
		}
		else{
			numValues = nativeLigs.size()+decoyLigs.size()+neoLigs.size();
			vector<float> maxes(numValues, numeric_limits<float>::max());
			engs.insert(engs.end(), maxes.begin(), maxes.end());
		}
	}

	return engs;
}

/* Calculate overall fitness. See header file for fitness definitions.
 *
 * MAKE SURE THIS DOESN'T BLOW UP IN YOUR FACE...
 */
void MultipleStructureOrganismNeo::calcFitnessAndPhenotype(){

	unsigned int i,j, offset;
	float energy;

	vector<float> foldEngs = calcFoldingEnergies();
	vector<float> bindEngs = calcBindingEnergies();
	fitness = 0.5;
	phenotype = "Original";

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
	 * and if there are any additional ones from decoy/novel
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
			offset = i*(nativeLigs.size()+decoyLigs.size()+neoLigs.size());
			for(j = 0; j < nativeLigs.size(); j++){
				energy = bindEngs[offset+j];

				/* Binding failure, no fitness.
				 */
				if(energy > bindThreshold){

//					/* DEBUG
//					 *
					cout << "Binding failure for ligand " << (j+1) << " (score " << energy << " > " << bindThreshold << ")." << endl;
//					*/

					fitness = 0.0;
					phenotype = "Dead";
					return;
				}
			}

			/* Any decoys bind, too?
			 */
			offset = i*(nativeLigs.size()+decoyLigs.size()+neoLigs.size())+nativeLigs.size();
			for(j = 0; j < decoyLigs.size(); j++){
				energy = bindEngs[offset+j];

				// CURRENTLY 10% FITNESS PENALTY PER DECOY
				if(energy < bindThreshold){
//					/* DEBUG
//					 *
					cout << "Binding decoy " << (j+1) << " (score " << energy << ")." << endl;
//					*/

					fitness *= 0.9;
					phenotype = "Decoy";
				}
			}

			/* How about the novel ligands?
			 */
			offset = i*(nativeLigs.size()+decoyLigs.size()+neoLigs.size())+(nativeLigs.size()+decoyLigs.size());
			for(j = 0; j < decoyLigs.size(); j++){
				energy = bindEngs[offset+j];

				// CURRENTLY 5% FITNESS GAIN PER NOVEL LIGAND
				if(energy < bindThreshold){
//					/* DEBUG
//					 *
					cout << "Binding novel " << (j+1) << " (score " << energy << ")." << endl;
//					*/

					fitness *= 1.05;
					phenotype += "+Neo"+itos(j+1);
				}
			}
		}
	}
}
