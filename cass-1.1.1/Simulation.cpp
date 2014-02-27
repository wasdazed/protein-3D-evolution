#include "Simulation.h"
#include "Organism.h"
#include "randomgen/randomc.h"
#include <typeinfo>
#include <stdexcept>
#include <iostream>

/* Shared static PRNG.
 */
shared_ptr<CRandomMersenne> Simulation::pMersenneGen;

/* Default constructor, doesn't do much.
 */
Simulation::Simulation(){
}

/* Constructor from a single original organism (produce homogenous
 * population of clones).
 */
Simulation::Simulation(int generations, shared_ptr<Organism> origSample, int sampleSize){

	unsigned int i;

	/* Make sure the class has its PRNGs handy.
	 */
	if(Simulation::pMersenneGen.get() == 0){
		throw runtime_error("You must supply a random number generator via Simulation::setPRNG()!");
	}


	/* Initialize the simple parameters
	 */
	numGenerations = generations;
	pStartingSample = origSample;
	numSamples = sampleSize;
	currentStep = 0;

	/* Set up a population of identical Organisms.
	 */
	for(i = 0; i < numSamples; i++){
		try{
			currentPopulation.push_back(origSample->Clone(false)); 	
		}
		catch(bad_alloc&){
			cerr << "Couldn't allocate enough memory for that many Organisms!" << endl;
			throw;
		}
	}
}

/* Constructor from a pre-produced set of Organisms.
 */
Simulation::Simulation(int generations, vector< shared_ptr<Organism> > inputPopulation){
	
	/* Set some atomic data members.
	 */
	numGenerations = generations;
	pStartingSample = shared_ptr<Organism>(); // No one starting sample here
	currentStep = 0;

	/* Just set the current population to whatever the input
	 * population is.
	 */
	currentPopulation = inputPopulation;

	/* Number of samples determined post-hoc from number of
	 * restored organisms.
	 */
	numSamples = currentPopulation.size();
}


/* Accessor method for number of generations
 */
int Simulation::getNumGenerations() const{
	return numGenerations;
}

/* Accessor method for number of samples
 */
int Simulation::getNumSamples() const{
	return numSamples;
}

/* Accessor method for how many steps the simulation
 * has gone.
 */
int Simulation::getNumSteps() const{
	return currentStep;
}

/* Performs one generation/time step of the simulation
 */
void Simulation::doTimeStep(){

	/* Calculate some statistics for this generation.
	 */
	calcPopulationStats();

	/* Select the members of the next generation.
	 */
	selectNewPopulation();

	/* DEBUG
	 *
	cout << "One timestep finished." << endl << "--------" << endl;
	*/

	currentStep++;
}

/* Saves some properties of the population of Organisms to file or
 * other stream.
 *
 * Returns true on success, false if the file stream has some problem
 * (i.e. a "badbit" is set).
 */
bool Simulation::savePopulation(ostream &outFile){

	unsigned int i;

	/* Check if the stream is okay.
	 */
	if(outFile.bad() == true){
		return false;
	}

	/* Write a representation of each Organism to file.
	 */
	for(i = 0; i < currentPopulation.size(); i++){
		outFile << currentStep << "\t" << currentPopulation[i]->toString() << endl;
	}

	/* Everything succeeded.
	 */
	return true;
}

/* Returns a copy of the current set of Organisms.
 */
vector< shared_ptr<Organism> > Simulation::getPop() const{
	return currentPopulation;
}

/* Calculates the fitness of all organisms in the population.
 */
void Simulation::calcPopulationFitness(){
	unsigned int i;

	/* Clear the old fitness values.
	 */
	populationFitness.clear();

	/* For all Organisms, get the fitness.
	 */
	for(i = 0; i < numSamples; i++){ 
		populationFitness.push_back(currentPopulation[i]->getFitness()); 
	}
}

/* Selects a new generation of hosts and parasites from the current
 * generation, based on a modified Wright-Fisher model of population 
 * genetics.
 *
 * See [some paper] for a description of the method.
 */
void Simulation::selectNewPopulation(){

	unsigned int randomOrganismNum;
	float randomFloat;

	/* DEBUG
	 *
	cout << "Selecting the next generation..." << endl;
	*/

	/* Select new generation of organisms as for Wright-Fisher model, i.e. random 
	 * drawings with replacement, depending on fitness (if a random number is 
	 * below fitness, it lives). The selected organisms may mutate during the
	 * transfer.
	 */
	vector<float> zeroes(currentPopulation.size(), 0.0);
	calcPopulationFitness();
	if(populationFitness != zeroes){

		while(newPopulation.size() < numSamples){

			/* DEBUG
			 *
			 cout << "Next generation now has " << newPopulation.size() << " individuals." << endl;
			 */

			/* Draw a random organism, with replacement, and draw a 
			 * random number on (0,1) (need to exclude 0, which the generator
			 * does not).
			 */
			randomOrganismNum = Simulation::pMersenneGen->IRandom(0,numSamples-1);

			/* DEBUG
			 *
			 cout << "Testing organism #" << randomOrganismNum << endl;
			 cout << "Fitness: " << currentPopulation[randomOrganismNum]->getFitness() << endl;
			 cout << "Phenotype: " << currentPopulation[randomOrganismNum]->getPhenotype() << endl;
			 */

			randomFloat = Simulation::pMersenneGen->Random();
			while(randomFloat == 0){
				randomFloat = Simulation::pMersenneGen->Random();
			}

			/* Keep the host if its fitness is above the random number.
			 * Allow possibility of mutation as it is copied into the
			 * new population.
			 */
			if(randomFloat < currentPopulation[randomOrganismNum]->getFitness()){

				/* DEBUG
				 *
				 cout << randomFloat << " < " << currentPopulation[randomOrganismNum]->getFitness() << ": selected." << endl;
				 */

				try{
					/* DEBUG
					 *
					cout << "Attempting to copy this organism with mutations..." << endl;
					*/

					newPopulation.push_back( currentPopulation[randomOrganismNum]->Clone(true) );
				}
				catch(bad_alloc&){
					cerr << "Couldn't copy that many Organisms to the next generation!" << endl;
					throw;
				}
			}
		}
	}
	/* Die in an informative fashion.
	 */ 	
	else{
		throw runtime_error("FATAL ERROR: All Organisms dead, simulation cannot continue!");
	}

	/* Replace the old population -- delete it
	 * and save the new pointers in its stead.
	 */
	currentPopulation = newPopulation;
	newPopulation.clear();
}

void Simulation::setPRNG(CRandomMersenne &rangeGenerator){
	Simulation::pMersenneGen = shared_ptr<CRandomMersenne>(new CRandomMersenne(rangeGenerator));
}
