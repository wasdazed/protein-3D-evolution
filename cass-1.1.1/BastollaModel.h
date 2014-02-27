/* Energy model due to Bastolla et al, 2001, Proteins 44:2. Knowledge-
 * based residue-residue contact model, with a quasi-chemical 
 * approximation parameterized by optimizing contact potentials over
 * a large number of protein structures.
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef BASTOLLAMODEL_H
#define BASTOLLAMODEL_H

/* Foward declarations.
 */
class Structure;

/* Includes.
 */
#include <string>
#include "EnergyModel.h"

class BastollaModel : public EnergyModel{
	public:

		/* Constructors and destructors.
		 */
		BastollaModel();
		virtual ~BastollaModel();

		/* Calculating folding and interaction scores/energies.
		 */
		virtual float calcFoldScore(shared_ptr<Structure> protein);
		virtual float calcInteractionScore(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo);

		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein);
		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation);

	protected:

		/* Data members, representing the interaction energy
		 * matrix of Bastolla.
		 */
		static float bastollaMatrix[20][20];
		static string bastollaKey;
};

/* End idempotency.
 */
#endif
