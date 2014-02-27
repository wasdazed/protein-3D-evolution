/* Abstract class that defines some common properties of all simluations.
 */

/* Ensure idempotency (avoid multiple includes).
 */
#ifndef SIMULATION_H
#define SIMULATION_H

/* Forward declarations.
 */
class Organism;
class CRandomMersenne;

/* Includes.
 */
#include <vector>
#include <fstream>
#include <tr1/memory>

using namespace std;

class Simulation {
	
	public:
		Simulation();
		Simulation(int generations, shared_ptr<Organism> origSample, int sampleSize); 
		Simulation(int generations, vector< shared_ptr<Organism> > inputPopulation);

		int getNumGenerations() const;
		int getNumSamples() const;
		int getNumSteps() const;
		void doTimeStep();

		/* For dumping simulation data to file or elsewhere.
		 */
		bool savePopulation(ostream &outFile);
		vector< shared_ptr<Organism> > getPop() const;

		// EXPERIMENTING W/ STATIC PRNGS...
		static void setPRNG(CRandomMersenne &rangeGenerator);

	protected:
		unsigned int numGenerations;
		unsigned int numSamples;
		unsigned int currentStep;

		// EXPERIMENTING W/ STATIC PRNGS...
		static shared_ptr<CRandomMersenne> pMersenneGen;

		/* Pointers to organisms being simulated.
		 */
		shared_ptr<Organism> pStartingSample;
		vector< shared_ptr<Organism> > currentPopulation;
		vector< shared_ptr<Organism> > newPopulation;

		vector<float> populationFitness; // ONLY USEFUL TO DO FITNESS MEASUREMENTS...
		
		void calcPopulationFitness();
		void selectNewPopulation();
		virtual void calcPopulationStats() =0;
};

/* Idempotency end.
 */
#endif
