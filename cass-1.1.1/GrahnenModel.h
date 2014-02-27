/* Semi-empirical force field model due to Grahnen et al, 20XX, XXXXX.
 * Extends and modifies the model of Rastogi et al, 2006, 124:2 in order to
 * deal with a static backbone, charge-charge interactions, solvation
 * effects and disulphide bridges.
 *
 * ALSO CAPABLE OF CALCULATING GAPS ACCORDING TO RANDOM ENERGY MODEL
 * OF GOLDSTEIN, NOW.
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef GRAHNENMODEL_H
#define GRAHNENMODEL_H

/* Forward declarations.
 */
class CRandomMersenne;

/* Includes.
 */
#include "RastogiModel.h"
#include <vector>

/* Definitions.
 */

class GrahnenModel : public RastogiModel{
	public:

		/* Constructors and destructors.
		 */
		GrahnenModel();
		GrahnenModel(float bendWeight, float LJWeight, float helixWeight, float sheetWeight, float chargeWeight, float solvWeight, float cysWeight, float LJWeightBind, float chargeWeightBind, float solvWeightBind);
		virtual ~GrahnenModel();

		/* Calculating folding and interaction scores/energies
		 * with the default parameters.
		 */
		virtual float calcFoldScore(shared_ptr<Structure> protein);
		virtual float calcInteractionScore(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo);

		/* Calculating individual terms in the energy/scoring functions.
		 */
		virtual vector<float> getFoldScoreTerms(shared_ptr<Structure> protein);
		virtual vector<vector<float>> getFoldScoreTerms_pos(shared_ptr<Structure> protein);
		virtual vector<float> getBindScoreTerms(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo);

		/* Calculating the energy/scoring functions with arbitrary
		 * parameters.
		 */
		virtual float getFoldScoreByParams(shared_ptr<Structure> protein, float bendWeight, float LJWeight, float helixWeight, float sheetWeight, float chargeWeight, float solvWeight, float cysWeight);
		virtual float getBindScoreByParams(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo, float LJWeight, float chargeWeight, float solvWeight);

		/* Calculating the fold score relative to the unfolded
		 * state.
		 */
		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein);
		virtual float calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation);
		virtual vector<float> sampleFoldGapTerms(string sequence, vector< pair<float,float> > bondAngles, vector<float> cBetaCAlphaDists, vector<float> cBetaDists, vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, vector<string> betaFourMers, vector<float> backboneTorsions, vector<float> exposures, CRandomMersenne& randGen);
		virtual vector<vector<float>> sampleFoldGapTerms_pos(string sequence, vector< pair<float,float> > bondAngles, vector<float> cBetaCAlphaDists, vector<float> cBetaDists, vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, vector<string> betaFourMers, vector<float> backboneTorsions, vector<float> exposures, CRandomMersenne& randGen);

		/* Helper functions to calculate special versions of
		 * energy function under his model.
		 */
		static float scwrlRepulsionEnergy(shared_ptr<Structure> protein, vector<bool> residueCanMove);
		static float LJRepulsionEnergy(shared_ptr<Structure> protein, vector<bool> residueCanMove);
		static double getDefaultBondLength(string resType);

	protected:

		/* Helper functions.
		 */
		float calcSidechainBendPotential(shared_ptr<Structure> protStruct);
		vector<float> calcSidechainBendPotential_pos(shared_ptr<Structure> protStruct);
		float calcChargePotential(shared_ptr<Structure> protStruct);
		vector<float> calcChargePotential_pos(shared_ptr<Structure> protStruct);
		float calcSolvationPotential(shared_ptr<Structure> protStruct, bool isAComplex);
		vector<float> calcSolvationPotential_pos(shared_ptr<Structure> protStruct, bool isAComplex);
		float calcCysteinePotential(shared_ptr<Structure> protStruct);
		vector<float> calcCysteinePotential_pos(shared_ptr<Structure> protStruct);
		float calcSolvationChange(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo);
		float dielectricConstant(float exposureOne, float exposureTwo);

		float sampleBend(string sequence, vector< pair<float,float> > bondAngles, CRandomMersenne& randGen);
		vector<float> sampleBend_pos(string sequence, vector< pair<float,float> > bondAngles, CRandomMersenne& randGen);
		float sampleLJ(string sequence, vector<float> cBetaCAlphaDists, vector<float> cBetaDists, CRandomMersenne& randGen);
		vector<float> sampleLJ_pos(string sequence, vector<float> cBetaCAlphaDists, vector<float> cBetaDists, CRandomMersenne& randGen);
		float sampleHelix(vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, CRandomMersenne& randGen);
		vector<float> sampleHelix_pos(vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, CRandomMersenne& randGen);
		vector<float> sampleHelix_pos2(vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, CRandomMersenne& randGen);
		float sampleBeta(vector<string> betaFourMers, vector<float> backboneTorsions, CRandomMersenne& randGen);
		vector<float> sampleBeta_pos(vector<string> betaFourMers, vector<float> backboneTorsions, CRandomMersenne& randGen);
		vector<float> sampleBeta_pos2(vector<string> betaFourMers, vector<float> backboneTorsions, CRandomMersenne& randGen);
		float sampleCharge(string sequence, vector<float> cBetaDists, vector<float> exposures, CRandomMersenne& randGen);
		vector<float> sampleCharge_pos(string sequence, vector<float> cBetaDists, vector<float> exposures, CRandomMersenne& randGen);
		float sampleSolv(string sequence, vector<float> exposures, CRandomMersenne& randGen);
		vector<float> sampleSolv_pos(string sequence, vector<float> exposures, CRandomMersenne& randGen);
		float sampleSSBond(string sequence, vector<float> cBetaDists, CRandomMersenne& randGen);
		vector<float> sampleSSBond_pos(string sequence, vector<float> cBetaDists, CRandomMersenne& randGen);

		vector<float> sampleHelixPosition(int protSize);

		/* Data members, mainly parameters.
		 */
		float w_bend;
		float w_LJ; 
		float w_helix;
		float w_beta;
		float w_charge;
	       	float w_solv; 
		float w_cys;

		float w_bind_LJ;
		float w_bind_charge;
		float w_bind_solv;
};

/* End idempotency.
 */
#endif
