/* Class to represent an organism with one or more proteins, one or more
 * decoy conformations, one or more original ligands, and one or more
 * decoy ligands.
 *
 * Fitnesses are calculated as:
 *
 * - Organism which loses folding/binding such that it cannot bind
 *   the original ligand(s) has a fitness of 0.
 *
 * - Organism which binds the native ligand(s), but nothing else,
 *   has a fitness of 0.5.
 *
 * - Organism that binds any decoy ligand(s) receives a X% (HOW MUCH?) 
 *   multiplicative cumulative penalty to fitness, e.g. fitness of 
 *   0.5*(0.9X)^N.
 */

/* Idempotency.
 */
#ifndef MULTIPLESTRUCTUREORGANISM_H
#define MULTIPLESTRUCTUREORGANISM_H

/* Forward declarations.
 */

/* Includes.
 */
#include "Organism.h"

class MultipleStructureOrganism : public Organism {

	public:

		MultipleStructureOrganism();
		MultipleStructureOrganism(float mutationRate, vector<string> codingSequences, vector< shared_ptr<Structure> > nativeStructures, vector< shared_ptr<Structure> > nativeLigands, vector< shared_ptr<Structure> > decoyLigands, shared_ptr<EnergyModel> engModel, float foldEnergyThreshold, float bindEnergyThreshold, shared_ptr<Structure> geometrySample);
		virtual ~MultipleStructureOrganism();

		MultipleStructureOrganism(const MultipleStructureOrganism &org, bool evolve);

		virtual shared_ptr<Organism> Clone(bool evolve);

	protected:

		/* The binding threshold.
		 */
		float bindThreshold;

		/* Ligands.
		 */
		vector< shared_ptr<Structure> > nativeLigs;
		vector< shared_ptr<Structure> > decoyLigs;

		/* Folding, binding and fitness calculations.
		 */
		virtual vector<float> calcFoldingEnergies();
		virtual vector<float> calcBindingEnergies();
		virtual void calcFitnessAndPhenotype();

		/* Storage for a structure to sample geometrical
		 * measurements from (for folding energy
		 * calculations.
		 */
		shared_ptr<Structure> invariant;
};

/* Idempotency end.
 */
#endif
