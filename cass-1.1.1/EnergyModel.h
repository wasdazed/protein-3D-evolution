/* More or less abstract class to define a common interface for
 * all types of energy (or scoring) models.
 *
 * Doesn't do any calculations on its own, strictly for deriving
 * from.
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef ENERGYMODEL_H
#define ENERGYMODEL_H

/* Foward declarations.
 */
class Structure;

/* Includes.
 */
#include <sys/time.h> // Used for timing calculations in subclasses
#include <tr1/memory>

/* Definitions.
 */
#define SEQ_SEPARATION 4 // Minimum number of residues between pairwise contacts
#define NUM_RANDOM_STRUCTURES 100 // Number of samples from distributions in Random Energy Model.


using namespace std;

class EnergyModel{
	public:

		/* Constructors and destructors.
		 */
		EnergyModel();
		virtual ~EnergyModel();

		/* Calculating folding and interaction scores/energies.
		 */
		virtual float calcFoldScore(shared_ptr<Structure> protein) =0;
		virtual float calcInteractionScore(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo) =0;

		// EXPERIMENTAL
		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein) =0;
		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation) =0;

	protected:

		/* Data members.
		 *
		 * DO THEY ALL HAVE ANY IN COMMON?
		 */
};

/* Idempotency end.
 */
#endif
