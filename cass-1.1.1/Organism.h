/* An Organism class to hold a bunch of methods and properties that
 * the all simulated organisms have in common.
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef ORGANISM_H
#define ORGANISM_H

/* Forward declarations.
 */
class Structure;
class EnergyModel;
class CRandomMersenne;
class StochasticLib1;

/* Includes
 */
#include <sys/time.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <tr1/memory>

using namespace std;

/* Definitions
 */
typedef unordered_map<string, float> string_float_map;
typedef unordered_map<string, shared_ptr<Structure> > string_structure_map;
typedef string_float_map::value_type stringAndFloat;
typedef string_structure_map::value_type stringAndStructure;

/* Interface.
 */
class Organism{
	public:

		/* Constructors and destructors.
		 */
		Organism();
		Organism(float mutationRate, vector<string> codingSequences, vector< shared_ptr<Structure> > structures, shared_ptr<EnergyModel> engModel, float foldingThreshold);
		virtual ~Organism();

		/* Copy constructors.
		 */
		Organism(const Organism& org);
		Organism(const Organism& org, bool evolve);
		virtual shared_ptr<Organism> Clone(bool evolve);

		/* Acessor methods.
		 *
		 * FUNCTION DEFINITION ORDER DOESN'T MAP 1-TO-1 WITH THE ACTUAL IMPLEMENTATION...
		 */
		float getMutationRate() const;
		float getFitness();
		int getNumberOfMutations() const;
		string getPhenotype();

		vector<string> getDnaSequences() const;
		vector<string> getProteinSequences() const;
		vector< shared_ptr<Structure> > getStructures() const;

		/* Properties will vary among subclasses, so string-ifying
		 * might be best done there (depending on what you're looking for).
		 */
		virtual string toString() const;

		// EXPERIMENTING W/ STATIC PRNGS...
		static void setPRNGs(CRandomMersenne &rangeGenerator, StochasticLib1 &stochasticGenerator);

	protected:

		/* Some data members.
		 */
		float mutRate;
		int cumNumMutations;
		bool changeHappened;

		float fitness;
		string phenotype;

		float foldThreshold;

		int id;
		int parentId;
		static unsigned long int totalCount;
		
		vector<string> dnaSeqs;
		vector< shared_ptr<Structure> > structs; 

		// EXPERIMENTING W/ STATIC PRNGS...
		static shared_ptr<CRandomMersenne> pMersenneGen;
		static shared_ptr<StochasticLib1> pStochastGen;

		/* Pointer to whatever energy model you're using for
		 * the calculations.
		 */
		shared_ptr<EnergyModel> pEngModel;

		/* Fitness and energy functions will depend on what phenotypic
		 * model is used, so they're virtual.
		 */
		virtual void calcFitnessAndPhenotype();
		
		/* Functions to calculate folding energy for all protein structures 
		 * in the Organism.
		 */
		virtual vector<float> calcFoldingEnergies();

		/* Function to calculate binding energy of two Structures.
		 *
		 * TRANSFER RESPONSIBILITY TO EnergyModel? OR JUST DO
		 * THAT DOWNSTREAM? YOU MAY NEED SOME WRAPPING FUNCTIONALITY,
		 * FOR INSTANCE.
		 */
		virtual float calcBindingEnergy(shared_ptr<Structure> bindingStructure, shared_ptr<Structure> proteinStructure); 

		/* Keeps track of what sequences in what structures
		 * correspond to what energies (saves computational
		 * time in the long run).
		 *
		 * OF SOMEWHAT DUBIOUS USE...
		 */
		static vector<string_float_map> knownFoldEnergies;
		static vector<string_structure_map> knownStructures;
};

/* Idempotency end.
 */
#endif
