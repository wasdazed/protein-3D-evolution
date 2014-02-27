/*
 * Program to test the sampling algorithm utilized by CASS for the
 * population genetics portion of the simulations.
 */

/* Includes.
 */
#include "randomgen/randomc.h"
#include <stdexcept>
#include <iostream>
#include <vector>
#include <memory>

using namespace std;

/* Mock Organism class for testing purposes.
 *
 * Only supports initalization, fitness/allele retrieval and
 * a Clone() method that produces an exact copy of the object.
 */
class Organism{
	public:
		Organism(float init_fitness, char allelic_state){
			this->allele = allelic_state;
			this->fitness = init_fitness;
		}
		
		float getFitness(){
			return this->fitness;
		}

		char getAllele(){
			return this->allele;
		}

		shared_ptr<Organism> Clone(bool mutate){
			return shared_ptr<Organism>(new Organism(this->fitness, this->allele));
		}
		
	protected:
		float fitness;
		char allele;
};

/* Mock Simulation class for testing purposes.
 *
 * Only supports statically setting and accessing a pointer to a MersenneRandom object.
 */
class Simulation{
	public:
		static shared_ptr<CRandomMersenne> pMersenneGen;
		Simulation(CRandomMersenne &randGen){
			Simulation::pMersenneGen = shared_ptr<CRandomMersenne>(new CRandomMersenne(randGen));
		}
};

/* Global variables.
 *
 * Stand-ins for the state variables provided by the Simulation class in the full
 * code.
 */
vector< shared_ptr<Organism> > currentPopulation;
vector< shared_ptr<Organism> > newPopulation;
vector<float> populationFitness;
unsigned int numSamples = 1000; // Population size
shared_ptr<CRandomMersenne> Simulation::pMersenneGen;

/* Mock method to enable running selectNewPopulation(). Does nothing.
 */
void calcPopulationFitness(){
}

/* START TESTED CODE SNIPPET.
 *
 * This method was copied verbatim from Simulation::selectNewPopulation() on
 * 2013-12-02 for testing purposes.
 */
void selectNewPopulation(){
	unsigned int randomOrganismNum;
	float randomFloat;

	vector<float> zeroes(currentPopulation.size(), 0.0);
	calcPopulationFitness();
	if(populationFitness != zeroes){
		while(newPopulation.size() < numSamples){
			randomOrganismNum = Simulation::pMersenneGen->IRandom(0,numSamples-1);
			randomFloat = Simulation::pMersenneGen->Random();
			while(randomFloat == 0){
				randomFloat = Simulation::pMersenneGen->Random();
			}

			if(randomFloat < currentPopulation[randomOrganismNum]->getFitness()){
				try{
					newPopulation.push_back( currentPopulation[randomOrganismNum]->Clone(true) );
				}
				catch(bad_alloc&){
					cerr << "Couldn't copy that many Organisms to the next generation!" << endl;
					throw;
				}
			}
		}
	}
	else{
		throw runtime_error("FATAL ERROR: All Organisms dead, simulation cannot continue!");
	}

	currentPopulation = newPopulation;
	newPopulation.clear();
}
/* END TESTED CODE SNIPPET.
 */

/* Entry point.
 */
int main(int argc, char* argv[]){
	/* Start the random number generator.
	 */
	unsigned int seed = 123;
	CRandomMersenne randGen(seed);
	Simulation fakeSim(randGen);

	/* Set the number of attempts per test type.
	 */
	unsigned int num_tests = 1000;

	/* Test the sampling method for no selection case.
	*/
	cout << "Running test #1: No selection." << endl;
	float mean_num_A = 0.0;
	float baseline_fitness = 0.5; // As per all the Organism classes
	float s = 0.00;
	float testing_fitness = (1.0+s)*0.5;
	unsigned int pop_size = numSamples;
	for(unsigned int n = 0; n < num_tests; n++){
		currentPopulation.clear();
		for(unsigned int i = 0; i < pop_size/2; i++){
			currentPopulation.push_back(shared_ptr<Organism>(new Organism(baseline_fitness, 'a')));
		}
		for(unsigned int i = 0; i < pop_size/2; i++){
			currentPopulation.push_back(shared_ptr<Organism>(new Organism(testing_fitness, 'A')));
		}

		selectNewPopulation();

		unsigned int num_A = 0;
		for(unsigned int i = 0; i < pop_size; i++){
			if(currentPopulation[i]->getAllele() == 'A'){
				num_A++;
			}
		}
		mean_num_A += num_A;
	}
	mean_num_A /= num_tests;
	cout << "\tPopulation size: " << pop_size << endl;
	cout << "\tNumber of test runs: " << num_tests << endl;
	cout << "\tAllele 'a' fitness: " << baseline_fitness << endl;
	cout << "\tAllele 'A' fitness: " << testing_fitness << endl;
	cout << "\tExpected frequency of 'A' : " << (baseline_fitness + ((testing_fitness-baseline_fitness)/2))*pop_size << endl;
	cout << "\tMean observed frequency of 'A': " << mean_num_A << endl;

	/* Test the sampling method for negative selection case (5% drop).
	 */
	cout << "Running test #2: Negative selection (5%)." << endl;
	mean_num_A = 0.0;
	s = -0.05;
	testing_fitness = (1.0+s)*baseline_fitness;
	for(unsigned int n = 0; n < num_tests; n++){
		currentPopulation.clear();
		for(unsigned int i = 0; i < pop_size/2; i++){
			currentPopulation.push_back(shared_ptr<Organism>(new Organism(baseline_fitness, 'a')));
		}
		for(unsigned int i = 0; i < pop_size/2; i++){
			currentPopulation.push_back(shared_ptr<Organism>(new Organism(testing_fitness, 'A')));
		}

		selectNewPopulation();

		unsigned int num_A = 0;
		for(unsigned int i = 0; i < pop_size; i++){
			if(currentPopulation[i]->getAllele() == 'A'){
				num_A++;
			}
		}
		mean_num_A += num_A;
	}
	mean_num_A /= num_tests;
	cout << "\tPopulation size: " << pop_size << endl;
	cout << "\tNumber of test runs: " << num_tests << endl;
	cout << "\tAllele 'a' fitness: " << baseline_fitness << endl;
	cout << "\tAllele 'A' fitness: " << testing_fitness << endl;
	cout << "\tExpected frequency of 'A': " << (baseline_fitness + ((testing_fitness-baseline_fitness)/2))*pop_size << endl;
	cout << "\tMean observed frequency of 'A': " << mean_num_A << endl;

	/* Test the sampling method for positive selection case (5% increase).
	 */	
	cout << "Running test #3: Positive selection (5%)." << endl;
	mean_num_A = 0.0;
	s = 0.05;
	testing_fitness = (1.0+s)*baseline_fitness;
	for(unsigned int n = 0; n < num_tests; n++){
		currentPopulation.clear();
		for(unsigned int i = 0; i < pop_size/2; i++){
			currentPopulation.push_back(shared_ptr<Organism>(new Organism(baseline_fitness, 'a')));
		}
		for(unsigned int i = 0; i < pop_size/2; i++){
			currentPopulation.push_back(shared_ptr<Organism>(new Organism(testing_fitness, 'A')));
		}

		selectNewPopulation();

		unsigned int num_A = 0;
		for(unsigned int i = 0; i < pop_size; i++){
			if(currentPopulation[i]->getAllele() == 'A'){
				num_A++;
			}
		}
		mean_num_A += num_A;
	}
	mean_num_A /= num_tests;
	cout << "\tPopulation size: " << pop_size << endl;
	cout << "\tNumber of test runs: " << num_tests << endl;
	cout << "\tAllele 'a' fitness: " << baseline_fitness << endl;
	cout << "\tAllele 'A' fitness: " << testing_fitness << endl;
	cout << "\tExpected frequency of 'A': " << (baseline_fitness + ((testing_fitness-baseline_fitness)/2))*pop_size << endl;
	cout << "\tMean observed frequency of 'A': " << mean_num_A << endl;

	/* Finish.
	 */
	return 0;
}
