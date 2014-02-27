/* Model of potential of mean force of protein folding due to
 * Mukherjee and Bagchi, 2003, J Chem Phys 118:10. Based on
 * a coarse-grained description of protein structure, where each
 * residue is represented by a C-alpha and C-beta bead. Originally
 * used as a force field in molecular dynamics studies of protein folding.
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef MUKHERJEEMODEL_H
#define MUKHERJEEMODEL_H

/* Forward declarations.
 */

/* Includes.
 */
#include "EnergyModel.h"
#include <string>
#include <vector>

class MukherjeeModel : public EnergyModel{
	public:

		/* Constructors and destructors.
		 */
		MukherjeeModel();
		virtual ~MukherjeeModel();

		/* Calculating folding and interaction scores/energies.
		 */
		virtual float calcFoldScore(shared_ptr<Structure> protein);
		virtual float calcInteractionScore(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo);

		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein);
		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation);

	protected:

		/* Helper functions.
		 */
		float calcBondPotential(shared_ptr<Structure> protStruct);
		float calcBendPotential(shared_ptr<Structure> protStruct);
		float calcTorsionPotential(shared_ptr<Structure> protStruct);
		float calcNonbondPotential(shared_ptr<Structure> protStruct);
		vector<float> calcNonbondPotential_pos(shared_ptr<Structure> protStruct);
		float calcHelixPotential(shared_ptr<Structure> protStruct);
		vector<float> calcHelixPotential_pos(shared_ptr<Structure> protStruct);
		vector<float> calcHelixPotential_pos2(shared_ptr<Structure> protStruct);

		/* Data members.
		 */
		static float forceFieldMatrix[20][8];
		static string forceFieldKey;
};

/* Idempotency end.
 */
#endif
