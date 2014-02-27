/* Semi-empirical force field model due to Rastogi et al, 2006, 124:2. Extends
 * the model of Mukherjee and Bagchi, 2003, J Chem Phys 118:10 to consider
 * beta sheet formation.
 *
 * See MukherjeeModel for further details.
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef RASTOGIMODEL_H
#define RASTOGIMODEL_H

/* Forward declarations.
 */
//class Structure;

/* Includes.
 */
#include "MukherjeeModel.h"

class RastogiModel : public MukherjeeModel{
	public:

		/* Constructors and destructors.
		 */
		RastogiModel();
		virtual ~RastogiModel();

		/* Calculating folding and interaction scores/energies.
		 */
		virtual float calcFoldScore(shared_ptr<Structure> protein);
		virtual float calcInteractionScore(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo);

		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein);
		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation);

	protected:

		/* Helper functions.
		 */
		float calcBetaPotential(shared_ptr<Structure> protStruct);
		vector<float> calcBetaPotential_pos(shared_ptr<Structure> protStruct);
		vector<float> calcBetaPotential_pos2(shared_ptr<Structure> protStruct);
};

/* Idempotency end.
 */
#endif
