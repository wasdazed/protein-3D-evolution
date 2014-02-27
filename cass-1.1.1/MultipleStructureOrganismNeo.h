/* Class to represent an organism with one or more proteins, one or more
 * decoy conformations, one or more original ligands, one or more
 * decoy ligands, and one or more novel ligands.
 *
 * Fitnesses are calculated as:
 *
 * - Organism which loses folding/binding such that it cannot bind
 *   the original ligand(s) has a fitness of 0.
 *
 * - Organism which binds the native ligand(s), but nothing else,
 *   has a fitness of 0.5.
 *
 * - Organism that binds any decoy ligand(s) receives a X% (10%) 
 *   multiplicative cumulative penalty to fitness, e.g. fitness of 
 *   base_fitness*(0.9X)^N.
 *
 * - Organism that binds any novel ligand(s) receives a Y% (5%)
 *   multiplicative cumulative increase in fitness, e.g. fitness of
 *   adjusted_fitness*(1.0Y)^N.
 */

/* Idempotency.
 */
#ifndef MULTIPLESTRUCTUREORGANISMNEO_H
#define MULTIPLESTRUCTUREORGANISMNEO_H

/* Forward declarations.
 */

/* Includes.
 */
#include "MultipleStructureOrganism.h"

class MultipleStructureOrganismNeo : public MultipleStructureOrganism {
	public:

		MultipleStructureOrganismNeo();
		MultipleStructureOrganismNeo(float mutationRate, vector<string> codingSequences, vector< shared_ptr<Structure> > nativeStructures, vector< shared_ptr<Structure> > nativeLigands, vector< shared_ptr<Structure> > decoyLigands, vector< shared_ptr<Structure> > novelLigands, shared_ptr<EnergyModel> engModel, float foldEnergyThreshold, float bindEnergyThreshold, shared_ptr<Structure> geometrySample);
		virtual ~MultipleStructureOrganismNeo();

		MultipleStructureOrganismNeo(const MultipleStructureOrganismNeo &org, bool evolve);

		virtual shared_ptr<Organism> Clone(bool evolve);

	protected:

		/* Ligands.
		 */
		vector< shared_ptr<Structure> > neoLigs;

		/* Folding, binding and fitness calculations.
		 */
		virtual vector<float> calcBindingEnergies();
		virtual void calcFitnessAndPhenotype();
};

/* Idempotency end.
 */
#endif
