#include "RastogiModel.h"
#include "Structure.h"
#include "TwoBeadStructure.h"
#include <typeinfo>
#include <stdexcept>

/* Default constructor, doesn't do much.
 */
RastogiModel::RastogiModel():
	MukherjeeModel()
{
}

/* Destructor, required due to gcc being
 * finicky.
 */
RastogiModel::~RastogiModel(){
}

/* Calculating folding score/energy. Simply adds one
 * more term to Mukherjee's model.
 *
 * REMOVE TIMING WHEN DONE WITH EXPERIMENTS.
 *
 * REQUIRE BEAD MODEL STRUCTURE HERE?
 */
float RastogiModel::calcFoldScore(shared_ptr<Structure> protein){
	
	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only TwoBeadStructure is valid for RastogiModel.");
	}

	/* Calculate Mukherjee's terms.
	 */
	float eng = MukherjeeModel::calcFoldScore(protein);
	
	/* Start the clock.
	 *
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);
	*/

	/* Calculate the beta sheet term.
	 */
	float VBeta = calcBetaPotential(protein);

	/* Stop the clock and print.
	 *
	gettimeofday(&tim, NULL);
	double timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
	cout << (timeStop - timeStart) << endl;
	*/

	return eng + VBeta;
}

// DUMMY IMPLEMENTATION FOR NOW...
float RastogiModel::calcFoldEnergyGap(shared_ptr<Structure> protein){
	return calcFoldScore(protein);
}

// DUMMY IMPLEMENTATION FOR NOW...
float RastogiModel::calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation){
	return calcFoldScore(protein);
}

/* Helper function to calculate beta sheet potential in this
 * model.
 *
 * STILL NOT SURE ON EXACT FORM OF THIS POTENTIAL...
 */
float RastogiModel::calcBetaPotential(shared_ptr<Structure> protStruct){

	/* Beta sheet potential is calculated as
	 *
	 * V_beta = sum_2toN-2[ K_i_1-4 * (C_s * (phi_i,i+1 - phi_b))^2 ]
	 * 
	 * where
	 * 
	 * N = # residues
	 * K_i_1-4 = 1/4*(K_i-1, K_i, K_i+1, K_i+2)
	 * K_i = some linearly scaled version of Kim and Berg's beta propensity for i
	 * C_s = scaling constant = 0.01 (scales torsion angle diffs to same magnitude as distance diffs in helical function)
	 * phi_i,i+1 = virtual torsion bond angle between Ca_i and Ca_i+1
	 * phi_b = equilibrium torsional angle for beta sheets = 210 degrees
	 */
	int i, resNum1, resNum2, resNum3, resNum4;
	string resType1, resType2, resType3, resType4;
	float K_1_4;
	float potential = 0.0;
	float phi_b = -150; // Equilibrium Ca-Ca torsion angle (Ca_i-1->Ca_i->Ca_i+1->Ca+i+2 torsion) for beta sheets. Equivalent to 210 positive degrees.
	float C_s = 0.01; // Scaling constant to make similar to calcHelixPotential()

	/* Go through residues 2 to N-2.
	 */
	for(i = 1; i < protStruct->getNumResidues()-3; i++){

		/* Get the residue types.
		 */
		resType1 = protStruct->getSingleResidue(i-1).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType4 = protStruct->getSingleResidue(i+2).getOneLetterType();

		/* And their numbers in the lookup table.
		 */
		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);
		resNum4 = forceFieldKey.find(resType4);

		/* Calculate the force constant.
		 */
		K_1_4 = 0.25*( forceFieldMatrix[resNum1][5] + forceFieldMatrix[resNum2][5] + forceFieldMatrix[resNum3][5] + forceFieldMatrix[resNum4][5]);

		/* Get the Ca coordinates.
		 */
		coord cAlpha1 = protStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");
		coord cAlpha2 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		coord cAlpha3 = protStruct->getSingleResidue(i+1).getAtomCoordsByType("CA");
		coord cAlpha4 = protStruct->getSingleResidue(i+2).getAtomCoordsByType("CA");

		/* Add energy to the potential.
		 */
		potential += K_1_4*pow( C_s*(dihedralAngle(cAlpha1, cAlpha2, cAlpha3, cAlpha4) - phi_b), 2);

		/* DEBUG
		 *
		cout << "Virtual phi angle: " << dihedralAngle(cAlpha1, cAlpha2, cAlpha3, cAlpha4) << endl;
		*/
	}

	/* Return the total potential.
	 */
	return potential;
}

vector<float> RastogiModel::calcBetaPotential_pos(shared_ptr<Structure> protStruct){

	/* Beta sheet potential is calculated as
	 *
	 * V_beta = sum_2toN-2[ K_i_1-4 * (C_s * (phi_i,i+1 - phi_b))^2 ]
	 *
	 * where
	 *
	 * N = # residues
	 * K_i_1-4 = 1/4*(K_i-1, K_i, K_i+1, K_i+2)
	 * K_i = some linearly scaled version of Kim and Berg's beta propensity for i
	 * C_s = scaling constant = 0.01 (scales torsion angle diffs to same magnitude as distance diffs in helical function)
	 * phi_i,i+1 = virtual torsion bond angle between Ca_i and Ca_i+1
	 * phi_b = equilibrium torsional angle for beta sheets = 210 degrees
	 */

	int i, resNum1, resNum2, resNum3, resNum4;
	string resType1, resType2, resType3, resType4;
	float K_1_4;
	float potential = 0.0;
	float phi_b = -150; // Equilibrium Ca-Ca torsion angle (Ca_i-1->Ca_i->Ca_i+1->Ca+i+2 torsion) for beta sheets. Equivalent to 210 positive degrees.
	float C_s = 0.01; // Scaling constant to make similar to calcHelixPotential()

	int seqLength = protStruct->getNumResidues();
	vector<float>res(seqLength,0.0);

	/* Go through residues 2 to N-2.
	 */
	for(i = 1; i < seqLength-3; i++){

		/* Get the residue types.
		 */
		resType1 = protStruct->getSingleResidue(i-1).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType4 = protStruct->getSingleResidue(i+2).getOneLetterType();

		/* And their numbers in the lookup table.
		 */
		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);
		resNum4 = forceFieldKey.find(resType4);

		/* Calculate the force constant.
		 */
		K_1_4 = 0.25*( forceFieldMatrix[resNum1][5] + forceFieldMatrix[resNum2][5] + forceFieldMatrix[resNum3][5] + forceFieldMatrix[resNum4][5]);

		/* Get the Ca coordinates.
		 */
		coord cAlpha1 = protStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");
		coord cAlpha2 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		coord cAlpha3 = protStruct->getSingleResidue(i+1).getAtomCoordsByType("CA");
		coord cAlpha4 = protStruct->getSingleResidue(i+2).getAtomCoordsByType("CA");

		/* Add energy to the potential.
		 */
		res[i] += K_1_4*pow( C_s*(dihedralAngle(cAlpha1, cAlpha2, cAlpha3, cAlpha4) - phi_b), 2);

		/* DEBUG
		 *
		cout << "Virtual phi angle: " << dihedralAngle(cAlpha1, cAlpha2, cAlpha3, cAlpha4) << endl;
		*/
	}

	/* Return the total potential.
	 */
	return res;
}

vector<float> RastogiModel::calcBetaPotential_pos2(shared_ptr<Structure> protStruct){

	/* Beta sheet potential is calculated as
	 *
	 * V_beta = sum_2toN-2[ K_i_1-4 * (C_s * (phi_i,i+1 - phi_b))^2 ]
	 *
	 * where
	 *
	 * N = # residues
	 * K_i_1-4 = 1/4*(K_i-1, K_i, K_i+1, K_i+2)
	 * K_i = some linearly scaled version of Kim and Berg's beta propensity for i
	 * C_s = scaling constant = 0.01 (scales torsion angle diffs to same magnitude as distance diffs in helical function)
	 * phi_i,i+1 = virtual torsion bond angle between Ca_i and Ca_i+1
	 * phi_b = equilibrium torsional angle for beta sheets = 210 degrees
	 */

	int i, resNum1, resNum2, resNum3, resNum4;
	string resType1, resType2, resType3, resType4;
	float K_1_4;
	float potential = 0.0;
	float phi_b = -150; // Equilibrium Ca-Ca torsion angle (Ca_i-1->Ca_i->Ca_i+1->Ca+i+2 torsion) for beta sheets. Equivalent to 210 positive degrees.
	float C_s = 0.01; // Scaling constant to make similar to calcHelixPotential()

	int seqLength = protStruct->getNumResidues();
	vector<float>res(seqLength,0.0);

	/* Go through residues 2 to N-2.
	 */
	for(i = 1; i < seqLength-3; i++){

		/* Get the residue types.
		 */
		resType1 = protStruct->getSingleResidue(i-1).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType4 = protStruct->getSingleResidue(i+2).getOneLetterType();

		/* And their numbers in the lookup table.
		 */
		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);
		resNum4 = forceFieldKey.find(resType4);

		/* Calculate the force constant.
		 */
//		K_1_4 = 0.25*( forceFieldMatrix[resNum1][5] + forceFieldMatrix[resNum2][5] + forceFieldMatrix[resNum3][5] + forceFieldMatrix[resNum4][5]);

		/* Get the Ca coordinates.
		 */
		coord cAlpha1 = protStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");
		coord cAlpha2 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		coord cAlpha3 = protStruct->getSingleResidue(i+1).getAtomCoordsByType("CA");
		coord cAlpha4 = protStruct->getSingleResidue(i+2).getAtomCoordsByType("CA");

		/* Add energy to the potential.
		 */
		float k = 0.25*pow( C_s*(dihedralAngle(cAlpha1, cAlpha2, cAlpha3, cAlpha4) - phi_b), 2);

		res[i-1]+=forceFieldMatrix[resNum1][5]*k;
		res[i]+=forceFieldMatrix[resNum2][5]*k;
		res[i+1]+=forceFieldMatrix[resNum3][5]*k;
		res[i+2]+=forceFieldMatrix[resNum4][5]*k;

		/* DEBUG
		 *
		cout << "Virtual phi angle: " << dihedralAngle(cAlpha1, cAlpha2, cAlpha3, cAlpha4) << endl;
		*/
	}

	/* Return the total potential.
	 */
	return res;
}


/* Calculating interaction score/energy.
 *
 * Dummy fuction, always returns "1e10". Override in subclasses.
 *
 * SHOULD THERE EVEN BE AN IMPLEMENTATION? PROBABLY, OTHERWISE
 * THIS BECOMES AN ABSTRACT CLASS...
 */
float RastogiModel::calcInteractionScore(shared_ptr<Structure> interactorOne, shared_ptr<Structure> interactorTwo){
	return 1000000000.0;
}
