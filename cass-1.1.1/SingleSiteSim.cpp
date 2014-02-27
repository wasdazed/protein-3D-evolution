#include "SingleSiteSim.h"
#include "Organism.h"

/* Default constructor, just sets pointers.
 */	
SingleSiteSim::SingleSiteSim():
	Simulation()
{
}

/* Constructor from a single original organism (produce homogenous
 * population of clones).
 */
SingleSiteSim::SingleSiteSim(int generations, shared_ptr<Organism> origSample, int sampleSize):
	Simulation(generations, origSample, sampleSize)
{
	/* Any need for anything very specific here?
	 */
}

/* Constructor from a pre-produced set of Organisms.
 */
SingleSiteSim::SingleSiteSim(int generations, vector< shared_ptr<Organism> > inputPopulation):
	Simulation(generations, inputPopulation)
{
	/* Do we need anything specific?
	 */
}

/* Calculate [RELEVANT STATISTICS] for the whole popluation.
 */
void SingleSiteSim::calcPopulationStats(){
	// Fitness? Divergence from starting point?
}

/* Accessor method for population DNA sequences.
 */
vector<string> SingleSiteSim::getPopDna() const{
	vector<string> seqs;
	vector<string> tempSeqs;
	unsigned int i;

	/* Since each organism may have multiple sequences,
	 * this requires a little bit of trickery (vector
	 * concatenation).
	 */
	for(i = 0; i < currentPopulation.size(); i++){
		tempSeqs = currentPopulation[i]->getDnaSequences();
		seqs.insert(seqs.end(), tempSeqs.begin(), tempSeqs.end());
	}

	return seqs;
}

/* Accessor method for population protein sequences.
 */
vector<string> SingleSiteSim::getPopProtein() const{
	vector<string> seqs;
	vector<string> tempSeqs;
	unsigned int i;

	/* Since each organism may have multiple sequences,
	 * this requires a little bit of trickery (vector
	 * concatenation).
	 */
	for(i = 0; i < currentPopulation.size(); i++){
		tempSeqs = currentPopulation[i]->getProteinSequences();
		seqs.insert(seqs.end(), tempSeqs.begin(), tempSeqs.end());
	}

	return seqs;
}
		
/* Accessor method for phenotypes in the population.
 * Assumes that any Organism only has one phenotype.
 */
vector<string> SingleSiteSim::getPopPhenotypes() const{
	
	vector<string> types;
	unsigned int i;

	/* Assume that each organism has only one phenotype.
	 */
	for(i = 0; i < currentPopulation.size(); i++){
		types.push_back(currentPopulation[i]->getPhenotype());
	}

	return types;
}
