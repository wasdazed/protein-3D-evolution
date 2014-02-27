/* Energy model expanding on that due to Bastolla et al, 2001, Proteins 44:2. 
 * Knowledge- based residue-residue contact model, with a quasi-chemical 
 * approximation parameterized by optimizing contact potentials over
 * a large number of protein structures. 
 *
 * ADD TWO-BEAD SUPPORT, SOLVATION.
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef BASTOLLAAUGMENTEDMODEL_H
#define BASTOLLAAUGMENTEDMODEL_H

/* Foward declarations.
 */

/* Includes.
 */
#include "BastollaModel.h"
#include <vector>

class BastollaAugmentedModel : public BastollaModel{
	public:

		/* Constructors and destructors.
		 */
		BastollaAugmentedModel();
		// INCLUDE CONSTRUCTOR WITH WEIGHTS
		virtual ~BastollaAugmentedModel();

		/* Calculating folding and interaction scores/energies.
		 */
		virtual float calcFoldScore(shared_ptr<Structure> protein);
		virtual float calcInteractionScore(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo);

		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein);
		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation);

		// EXPERIMENTAL, NOT VERY PHYSICAL
		vector<float> residueFoldContributions(shared_ptr<Structure> protein);
		vector<float> residueBindContributions(shared_ptr<Structure> protein, shared_ptr<Structure> interactor);

		virtual float singleContact(char residueOne, char residueTwo);

	protected:

		/* WEIGHTING SCHEME GOES HERE.
		 */

};

/* End idempotency.
 */
#endif
