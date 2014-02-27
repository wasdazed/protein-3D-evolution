#include "MukherjeeModel.h"
#include "Structure.h"
#include "TwoBeadStructure.h"
#include <typeinfo>
#include <stdexcept>

/* Static member: Table for looking up residue-specific
 * properties.
 *
 * Defines Cb bead equilibrium bond lengths, Cb bead radii, 
 * equilibrium Cb bond angles (degrees), epsilon_ii values, 
 * modified helix propensities, modified beta propensities,
 * sidechains charges and a "polarity index" (inverse of
 * epsilon_ii values, see GrahnenModel::calcSolvationPotential()).
 *
 * Only columns 1-5 are used by this model, the
 * remainder refer to properties used by subclasses.
 */
float MukherjeeModel::forceFieldMatrix[20][8] = 
{
// bond, radius, angle, epsilon, helix, beta, charge, polarity
{0.77, 0.77, 121.9,  7.76,  0.00, 10.75,  0, 3.44}, // A
{1.49, 1.29, 121.7, 10.64,   6.1, 16.28,  0, 0.56}, // V
{2.08, 1.54, 118.1, 10.16,  2.1, 14.74,  0, 1.04}, // L
{1.83, 1.56, 118.9, 11.00,  4.1, 17.20,  0, 0.2}, // I
{1.38, 1.22, 113.7,  8.60,   6.8, 14.44,  0, 2.6}, // C
{2.34, 1.80, 113.1,  7.88,  2.4, 14.13,  0, 3.32}, // M
{1.42, 1.25,  81.9,  7.52, 31.6,  7.06,  0, 3.68}, // P
{2.97, 1.90, 118.2,  8.96,   5.4, 16.89,  0, 2.24}, // F
{3.36, 2.13, 110.0,  4.04,   5.3, 15.36,  0, 7.16}, // Y
{3.58, 2.21, 118.4,  4.52,   4.9, 14.74,  0, 6.68}, // W
{1.99, 1.43, 121.2,  1.40,   6.9, 12.59, -1, 9.8}, // D
{1.98, 1.45, 117.9,  1.40,   6.5, 11.67,  0, 9.8}, // N
{2.58, 1.75, 118.0,  1.40,  3.9, 12.29,  0, 9.8}, // Q
{2.76, 1.78, 118.2,  1.76,   5.6, 14.13,  1, 9.44}, // H
{2.63, 1.77, 118.2,  1.40,  4.0, 12.59, -1, 9.8}, // E
{1.28, 1.08, 117.9,  4.64,   5.0, 11.98,  0, 6.56}, // S
{1.43, 1.24, 117.1,  4.76,   6.6, 14.74,  0, 6.44}, // T
{3.72, 2.38, 121.4,  0.20,  2.1, 13.51,  1, 11}, // R
{2.94, 2.08, 122.0,  0.92,  2.6, 12.59,  1, 10.28}, // K
{ 0.0,  0.0,   0.0,  5.12,    10.0,   0.0,  0, 6.08} // G
};               
string MukherjeeModel::forceFieldKey = "AVLICMPFYWDNQHESTRKG";

/* Default constructor, doesn't do much.
 */
MukherjeeModel::MukherjeeModel():
	EnergyModel()
{
}

/* Destructor, required due to gcc being
 * finicky.
 */
MukherjeeModel::~MukherjeeModel(){
}

/* Calculating folding score/energy.
 *
 * REMOVE TIMER WHEN EXPERIMENTS DONE.
 *
 * REQUIRE BEAD MODEL STRUCTURE AS INPUT INSTEAD?
 */
float MukherjeeModel::calcFoldScore(shared_ptr<Structure> protein){
	
	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only TwoBeadStructure is valid for MukherjeeModel.");
	}

	/* Start the clock.
	 */
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);

	/* Calculate the bonding potential.
	 */
	float VBond = calcBondPotential(protein);

	/* Calculate the bending potential.
	 */
	float VBend = calcBendPotential(protein);

	/* Calculate the torsional potential.
	 */
	float VTors = calcTorsionPotential(protein);

	/* Calculate the non-bonding (Lennard-Jones) potential.
	 */
	float VNonbond = calcNonbondPotential(protein);

	/* Finally, calculate the helix potential.
	 */
	float VHelix = calcHelixPotential(protein);

	/* Stop the clock and print.
	 */
	gettimeofday(&tim, NULL);
	double timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
	cout << (timeStop - timeStart) << endl;

	/* DEBUG
	 *
	cout << "VBond\tVBend\tVTors\tVNonbond\tVHelix" << endl;
	cout << VBond << "\t" << VBend << "\t" << VTors << "\t" << VNonbond << "\t" << VHelix << endl;
	*/

	/* Add up to the total potential and return.
	 */
	return VBond + VBend + VTors + VNonbond + VHelix;
}

// DUMMY IMPLEMENTATION FOR NOW...
float MukherjeeModel::calcFoldEnergyGap(shared_ptr<Structure> protein){
	return calcFoldScore(protein);
}

// DUMMY IMPLEMENTATION FOR NOW...
float MukherjeeModel::calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation){
	return calcFoldScore(protein);
}

/* Helper function to calculate the bonding potential under
 * Mukherjee's force field model.  
 * 
 * SOMEWHAT BOGUS: DUE TO NON-MODIFICATION
 * OF STRUCTURE BACKBONE DURING THREADING, CA-CA
 * POTENTIALS NEVER CHANGE. AND WITH SARA AS A REPLACEMENT
 * METHOD, CA-CB POTENTIALS DON'T CHANGE EITHER.
 *
 * SO TURNING IT OFF.
 */
float MukherjeeModel::calcBondPotential(shared_ptr<Structure> protStruct){

	/* Vbond = 0.5*K_r*sum_i=2toN[(r_i,i-1 - r_0)^2] + 0.5*K_rs*sum_i=1toN[(r_is,i - r_s0(i))^2]
	 * 
	 * where 
	 * 
	 * N = total number of residues/Ca's
	 * r_i,i-1 = dist(r_i,r_i-1) in R^3 space
	 * r_si,i = dist(r_si,r_i) in R^3 space
	 * r_i = position of Ca bead i
	 * r_si = position of side-chain bead i
	 * r_0 = equilibrium bond length bewteen Ca's = 3.81 A
	 * r_s0(i) = equilibrium bond length between Ca i and side-chain bead i, depends on size of side-chain (but how? See Leavitt tables 1 and 2?)
	 * K_r = force constant for Ca bonds = 43.0 kJ/mol/A^2
	 * K_rs = force constant for Ca-to-side-chain bonds = 8.6 kJ/mol/A^2
	 */

//	int i, resNum;
//
//	float cAlphaPot = 0;
//	float cBetaPot = 0;
//	
//	float r_0 = 3.81; 
//	float K_r = 43.0; 
//	float K_rs = 8.6; 
//
//	string resType;
//
//	int seqLength = protStruct->getNumResidues();
//
//	
//	/* Sum up the Ca-Ca bond length differences. Since indexing
//	 * is 0-based, 2->N becomes 1->N-1.
//	 */
//	for(i = 1; i < seqLength; i++){
//		coord cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
//		coord cAlpha2 = protStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");
//
//		/* THIS TERM WILL BE CONSTANT THROUGHOUT THE GENERATIONS,
//		 * AS THE SECONDARY STRUCTURE DOESN'T CHANGE...
//		 */
//		cAlphaPot += pow( ( NDDist(cAlpha1, cAlpha2) - r_0 ) , 2);
//	}
//
//	/* Sum up the Ca-Cb bond length differences. Indexing as
//	 * above. Don't do anything with glycines (have no Cb bead).
//	 */
//	for(i = 0; i < seqLength; i++){
//		resType = protStruct->getSingleResidue(i).getOneLetterType();
//
//		coord cAlpha = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
//
//		if(resType.compare("G") != 0){
//			coord cBeta = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");
//			resNum = forceFieldKey.find(resType);
//
//			// THIS TERM IS CONSTANT FOR IDENTICAL ROTAMERS IN
//			// THE SCWRL LIBRARY.
//			cBetaPot += pow( ( NDDist(cAlpha, cBeta) - forceFieldMatrix[resNum][0] ), 2);
//		}
//	}
//
//	/* DEBUG
//	 *
//	cout << "CAlphaBond = " << cAlphaPot << ", CBetaBond = " << cBetaPot << endl;
//	*/
//
//
//	/* Convert to potentials and return the total
//	 * potential.
//	 */
//	return 0.5*(K_r*cAlphaPot + K_rs*cBetaPot);
	
	return 0.0;
}	

/* Function to calculate the bending potential term in
 * Mukherjee's force field model.  
 */
float MukherjeeModel::calcBendPotential(shared_ptr<Structure> protStruct){
	/* V_theta, the bending potential:
	*
	* The sum of 3 bending potential terms around a central Ca bead, involving the 
	* 2 neighboring Ca's and a side-chain bead.
	* 
	* V_theta = 0.5*K_theta*sum_i=2toN-1[(theta_i-1,i,i+1 - theta_0)^2]
	* 	  + 0.5*K_theta*sum_i=2toN[(thetaS_i-1,i,i - thetaS_0(i))^2]
	* 	  + 0.5*K_theta*sum_i=1toN-1[(thetaS_i,i,i+1 - thetaS_0(i))^2]
	* where
	* 
	* N = all the Ca beads
	* theta_i-1,i,i+1 = angle between Ca's i-1, i, i+1
	* thetaS_i-1,i,i = angle between Ca i-1, Ca i and side-chain bead i
	* thetaS_i,i,i_1 = angle between side chain bead i, Ca i and Ca i+1
	* theta_0 = equilibrium bond angle between Ca's = 96 degrees
	* thetaS_0(i) = equilibirum bond angle between Ca i and side chain bead i, see 
	*               Table 1
	* K_theta = force constant for harmonic bending potential = 10.0 kJ/mol/rad^2 (Note: radians, not degrees!)
	*/
	float K_theta = 10.0;
	float theta_0 = 96;  // In degrees.
	float backboneTerm = 0.0;
	float sidechainTerm1 = 0.0;
	float sidechainTerm2 = 0.0;

	int seqLength = protStruct->getNumResidues();

	int i, resNum;
	string resType;
	coord cAlpha1, cAlpha2, cAlpha3, cBeta;

	/* Term 1: bending of backbone around Ca_i
	 */
	for(i = 1; i < seqLength-1; i++){
		cAlpha1 = protStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");
		cAlpha2 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		cAlpha3 = protStruct->getSingleResidue(i+1).getAtomCoordsByType("CA");

		/* DEBUG
		 *
		cout << "i: " << i << ", Ca_i: " << cAlpha2[0] << "," << cAlpha2[1] << "," << cAlpha2[2] << endl;
		*/

		backboneTerm += pow( (NDAngle(cAlpha1, cAlpha2, cAlpha3) - theta_0) * (PI/180.0), 2); // Conversion to radians important
	}

	/* Term 2: sidechain bending w r t to Ca_i-1
	 */
	for(i = 1; i < seqLength; i++){
		
		resType = protStruct->getSingleResidue(i).getOneLetterType();

		if(resType.compare("G") != 0){
			cAlpha1 = protStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");
			cAlpha2 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
			cBeta = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");
			resNum = forceFieldKey.find(resType);

			/* DEBUG
			*
			cout << "Residue type: " << resType << ", Residue #: " << resNum << ", Cb coords: " << cBeta[0] << "," << cBeta[1] << "," << cBeta[2] << endl;
			*/

			sidechainTerm1 += pow( (NDAngle(cAlpha1, cAlpha2, cBeta) - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
	}

	/* Term 3: sidechain bending w r t to Ca_i+1
	 */
	for(i = 0; i < seqLength-1; i++){

		resType = protStruct->getSingleResidue(i).getOneLetterType();

		if(resType.compare("G") != 0){
			cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
			cAlpha2 = protStruct->getSingleResidue(i+1).getAtomCoordsByType("CA");
			cBeta = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");
			resNum = forceFieldKey.find(resType);

			sidechainTerm2 += pow( (NDAngle(cBeta, cAlpha1, cAlpha2) - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
	}

	/* DEBUG
	 *
	cout << "Backbone bend: " << backboneTerm << ", sidechain bend 1: " << sidechainTerm1 << ", sidechain bend 2: " << sidechainTerm2 << endl;
	*/

	/* Sum up terms, convert to a potential and return.
	 */
	return 0.5*K_theta*(backboneTerm + sidechainTerm1 + sidechainTerm2);
}

/* Function to calculate the torsional potential under Mukherjee's
 * force field model. 
 *
 * ON THE BOGUS SIDE: OUR CURRENT METHOD DOESN'T CHANGE THE BACKBONE
 * OF THE STRUCTURE.
 */
float MukherjeeModel::calcTorsionPotential(shared_ptr<Structure> protStruct){
	/* V_t, the torsional potential:
	 * 
	 * The sum of the torsional angle in the pseudo-peptide bond 
	 * (not there for terminal bonds), i.e. phi.
	 *
	 * V_t = epsilon_t*sum_angles[0.5*(1+cos(3*angle))]
	 * 
	 * where
	 * 
	 * epsilon_t = 1 kJ/mol
	 * angle = pseudo-dihedral angle for the Ca_i-1-to-Ca_i bond, i.e.
	 *         the dihedral angle over Ca_i-2, Ca_i-1, Ca_i, Ca_i+1
	 */
	float epsilon_t = 1.0;
	float torsPot = 0.0;
	int i, seqLength;
	float dihedral;
	coord cAlpha1, cAlpha2, cAlpha3, cAlpha4;

	/* Calculate the dihedral angle for bonds where it exists
	 * (i.e. from 3 to N-1) and add up the potentials.
	 */
	seqLength = protStruct->getNumResidues();
	for(i = 2; i < seqLength-1; i++){
		cAlpha1 = protStruct->getSingleResidue(i-2).getAtomCoordsByType("CA");
		cAlpha2 = protStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");
		cAlpha3 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		cAlpha4 = protStruct->getSingleResidue(i+1).getAtomCoordsByType("CA");

		dihedral = dihedralAngle(cAlpha1, cAlpha2, cAlpha3, cAlpha4);
		torsPot += 1 + cos(3*dihedral*(PI/180)); // Conversion to radians important

		/* DEBUG
		 *
		cout << "Dihedral of residues " << i << "-" << i + 1 << " = " << dihedral << endl;
		*/
	}

	/* Scale the potential and return.
	 */
	return 0.5*epsilon_t*torsPot;
}

/* Function to calculate the non-bonding (van der Waals) potential
 * under Mukherjee's force field model. See mukherjeeFoldEnergy()
 * for more details.
 *
 * TURNING ON PROCESSING OF GLYCINES, LET'S SEE HOW THIS GOES...
 */
float MukherjeeModel::calcNonbondPotential(shared_ptr<Structure> protStruct){
	/* V_lj, non-bonding potential:
	 * 
	 * The Lennard-Jones potential between two atoms.
	 * 
	 * V_lj = 4*sum_i,j[epsilon_ij( (sigma_ij/r_ij)^12 - (sigma_ij/r_ij)^6 )]
	 *
	 * where
	 *
	 * i,j = every possible bead combination, but no self-interactions or
	 *       interactions with bonded partners.
	 *
	 * epsilon_ij = sqrt(epsilon_ii*epsilon_jj), epsilon_ii in Table 1 or 0.05 kJ/mol for Ca beads
	 * 
	 * sigma_ij = collision diameter for i,j, i.e. the twice the r_gyr (see Table 1 or Ca = 1.8 A)
	 *
	 * r_ij = dist(i,j)
	 */
	float epsilon_ca = 0.05;
	float epsilon_ij_ca = epsilon_ca; // I.e. sqrt(epsilon_ca*epsilon_ca)
	float r_gyr_ca = 1.8;
	float sigma_ij_ca = 2.0*r_gyr_ca;
	float ljPot = 0.0;
	int seqLength = protStruct->getNumResidues();

	int i, j, resNum1, resNum2;
	float epsilon_ij, sigma_ij;
	string resType2;
	coord cAlpha1, cAlpha2, cBeta1, cBeta2;

	/* Calculate the non-bonded Ca-Ca potentials
	 *
	 * STATIC FOR STATIC BACKBONE, IGNORE FOR NOW.
	*/
//	for(i = 0; i < seqLength; i++){
//
//		cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
//
//		for(j = i + SEQ_SEPARATION; j < seqLength; j++){
//
//			/* DEBUG
//			 *
//			 cout << "Ca " << i+1 << "-Ca " << j+1 << endl;
//			 */
//
//			/* Interaction with other Ca bead.
//			*/
//			cAlpha2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CA");
//			ljPot += epsilon_ij_ca*( pow( sigma_ij_ca/NDDist(cAlpha1, cAlpha2) ,12) - pow( sigma_ij_ca/NDDist(cAlpha1, cAlpha2) ,6) );
//		}
//	}

	/* DEBUG
	 *
	 cout << "Ca-Ca LJ: " << ljPot;
	 */

	/* Calculate the non-bonded Ca-Cb potentials.
	*/
	for(i = 0; i < seqLength; i++){

		cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		resNum1 = forceFieldKey.find(protStruct->getSingleResidue(i).getOneLetterType());

		for(j = i + SEQ_SEPARATION; j < seqLength; j++){

			/* Interaction with the other Cb bead, unless it's a
			 * glycine (avoids double-counting in the next step). 
			 */
			resType2 = protStruct->getSingleResidue(j).getOneLetterType();

			if(resType2.compare("G") == 0){
				/* DEBUG
				 *
				cout << "Cb is Gly, skipping..." << endl;
				*/

				continue;
			}
			else{
				cBeta2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CB");
			}

			/* DEBUG
			 *
			 cout << "Ca " << i+1 << "-Cb " << j+1 << endl;
			 */

			resNum2 = forceFieldKey.find(resType2);

			/* Compute epsilon_ij and sigma_ij
			*/
			epsilon_ij = sqrt(epsilon_ca * forceFieldMatrix[resNum2][3]);
			sigma_ij = r_gyr_ca + forceFieldMatrix[resNum2][1];

			/* Add to potential.
			*/
			ljPot += epsilon_ij*( pow( sigma_ij/NDDist(cAlpha1, cBeta2) ,12) - pow( sigma_ij/NDDist(cAlpha1, cBeta2) ,6) );
		}
	}

	/* DEBUG
	 *
	 cout << ", +Ca-Cb: " << ljPot;
	 */

	/* Calculate the non-bonded Cb-Cb potentials.
	*/
	for(i = 0; i < seqLength; i++){

		/* Use Ca bead for glycines (no Cb bead
		 * exists).
		 */
		if(protStruct->getSingleResidue(i).getOneLetterType().compare("G") == 0){

			/* DEBUG
			 *
			cout << "#1 is Glycine bead, special treatment." << endl;
			*/

			cBeta1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		}
		else{
			cBeta1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");
		}

		resNum1 = forceFieldKey.find(protStruct->getSingleResidue(i).getOneLetterType());

		for(j = i + SEQ_SEPARATION; j < seqLength; j++){

			/* Interaction with the other Cb bead, unless it's a
			 * glycine (use Ca bead as above)
			 */
			resType2 = protStruct->getSingleResidue(j).getOneLetterType();
			if(resType2.compare("G") == 0){

				/* DEBUG
				 *
				cout << "#2 is Glycine bead, special treatment." << endl;
				*/

				cBeta2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CA");
			}
			else{
				cBeta2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CB");
			}

			/* DEBUG
			 *
			 cout << "Cb " << i+1 << "-Cb " << j+1 << endl;
			 */

			resNum2 = forceFieldKey.find(resType2);

			/* Compute sigma_ij and epsilon_ij
			*/
			sigma_ij = forceFieldMatrix[resNum1][1] + forceFieldMatrix[resNum2][1];
			epsilon_ij = sqrt(forceFieldMatrix[resNum1][3] * forceFieldMatrix[resNum2][3]);

			/* Add to potential.
			*/
			ljPot += epsilon_ij*( pow( sigma_ij/NDDist(cBeta1, cBeta2) ,12) - pow( sigma_ij/NDDist(cBeta1, cBeta2) ,6) );
		}
	}

	/* DEBUG
	 *
	 cout << ", + Cb-Cb: " << ljPot << endl;
	 */

	/* Scale the potential and return.
	*/
	return 4*ljPot;
}

/* Function to calculate the helical potential under
 * Mukherjee's force field model. See mukherjeeFoldEnergy()
 * for more details.
 */
float MukherjeeModel::calcHelixPotential(shared_ptr<Structure> protStruct){

	/* V_helix = sum_i=3toN-3[ 0.5*K_1-3_i*(r_i,i+2-r_h)^2 + 0.5*K_1-4_i*(r_i,i+3-r_h)^2 ]
	 * 
	 * where
	 * 
	 * r_h = equilibrium distance in helix = 5.5 A
	 * K_1-3_i = (1/3)*(K_i + K_i+1 + K_i+2)
	 * K_1-4_i = (1/4)*(K_i + K_i+1 + K_i+2 + K_i+3)
	 * K_i = scaled helix propensity in kJ/mol/A^2 for residue i
	 *
	 * K's here are actually big Kappas.
	 */
	float r_h = 5.5;
	float helixPot = 0.0;
	float K_1_3, K_1_4;
	int i;

	string resType1, resType2, resType3, resType4;
	int resNum1, resNum2, resNum3, resNum4;

	int seqLength = protStruct->getNumResidues();

	/* Sum from 3 to N-3.
	 *
	 * REALLY START THIS AT POSITION 3? WHY?
	 * NOT IMMEDIATELY APPARENT IN PAPER. 	
	 */
	for(i = 2; i < seqLength-3; i++){
		/* Calculate the first term.
		 */

		resType1 = protStruct->getSingleResidue(i).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+2).getOneLetterType();

		/* DEBUG
		 *
		cout << "Residue type: " << resType1 << endl;
		*/

		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);

		K_1_3 = (1.0/3.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4]);

		coord cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		coord cAlpha2 = protStruct->getSingleResidue(i+2).getAtomCoordsByType("CA");


//Added to remove term values that are not under helices and distance
///		if ((i >= 12 && i < 22) || (i >= 79 && i < 86))
///		{
			helixPot += K_1_3;
///		}


//orig		helixPot += K_1_3*pow( NDDist(cAlpha1, cAlpha2) - r_h ,2);

		/* Calculate the second term.
		 */
		resType1 = protStruct->getSingleResidue(i).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+2).getOneLetterType();
		resType4 = protStruct->getSingleResidue(i+3).getOneLetterType();

		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);
		resNum4 = forceFieldKey.find(resType4);

		K_1_4 = (1.0/4.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4] + forceFieldMatrix[resNum4][4]);

		cAlpha2 = protStruct->getSingleResidue(i+3).getAtomCoordsByType("CA");

//Added to remove term values that are not under helices and distance
///		if ((i >= 12 && i < 22) || (i >= 79 && i < 86))
///		{
			helixPot += K_1_4;
//	cout << "MUKHERJEE Checkpoint" << endl;
///		}



		
//orig		helixPot += K_1_4*pow( NDDist(cAlpha1, cAlpha2) - r_h ,2);
	}

	/* Return the final potential.
	 */
	return 0.5*helixPot;
}


vector<float> MukherjeeModel::calcNonbondPotential_pos(shared_ptr<Structure> protStruct){
	/* V_lj, non-bonding potential:
	 *
	 * The Lennard-Jones potential between two atoms.
	 *
	 * V_lj = 4*sum_i,j[epsilon_ij( (sigma_ij/r_ij)^12 - (sigma_ij/r_ij)^6 )]
	 *
	 * where
	 *
	 * i,j = every possible bead combination, but no self-interactions or
	 *       interactions with bonded partners.
	 *
	 * epsilon_ij = sqrt(epsilon_ii*epsilon_jj), epsilon_ii in Table 1 or 0.05 kJ/mol for Ca beads
	 *
	 * sigma_ij = collision diameter for i,j, i.e. the twice the r_gyr (see Table 1 or Ca = 1.8 A)
	 *
	 * r_ij = dist(i,j)
	 */
	float epsilon_ca = 0.05;
	float epsilon_ij_ca = epsilon_ca; // I.e. sqrt(epsilon_ca*epsilon_ca)
	float r_gyr_ca = 1.8;
	float sigma_ij_ca = 2.0*r_gyr_ca;
	float ljPot = 0.0;
	int seqLength = protStruct->getNumResidues();

	vector<float>res(seqLength,0.0);

	int i, j, resNum1, resNum2;
	float epsilon_ij, sigma_ij;
	string resType2;
	coord cAlpha1, cAlpha2, cBeta1, cBeta2;

	/* Calculate the non-bonded Ca-Ca potentials
	 *
	 * STATIC FOR STATIC BACKBONE, IGNORE FOR NOW.
	*/
//	for(i = 0; i < seqLength; i++){
//
//		cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
//
//		for(j = i + SEQ_SEPARATION; j < seqLength; j++){
//
//			/* DEBUG
//			 *
//			 cout << "Ca " << i+1 << "-Ca " << j+1 << endl;
//			 */
//
//			/* Interaction with other Ca bead.
//			*/
//			cAlpha2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CA");
//			ljPot += epsilon_ij_ca*( pow( sigma_ij_ca/NDDist(cAlpha1, cAlpha2) ,12) - pow( sigma_ij_ca/NDDist(cAlpha1, cAlpha2) ,6) );
//		}
//	}

	/* DEBUG
	 *
	 cout << "Ca-Ca LJ: " << ljPot;
	 */

	/* Calculate the non-bonded Ca-Cb potentials.
	*/
	for(i = 0; i < seqLength; i++){

		cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		resNum1 = forceFieldKey.find(protStruct->getSingleResidue(i).getOneLetterType());

		for(j = 0; j < seqLength; j++){

			if ((i>=j&&(i-j)<SEQ_SEPARATION)||(i<j&&(j-i)<SEQ_SEPARATION)) {
							continue;
						}

			/* Interaction with the other Cb bead, unless it's a
			 * glycine (avoids double-counting in the next step).
			 */
			resType2 = protStruct->getSingleResidue(j).getOneLetterType();

			if(resType2.compare("G") == 0){
				/* DEBUG
				 *
				cout << "Cb is Gly, skipping..." << endl;
				*/

				continue;
			}
			else{
				cBeta2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CB");
			}

			/* DEBUG
			 *
			 cout << "Ca " << i+1 << "-Cb " << j+1 << endl;
			 */

			resNum2 = forceFieldKey.find(resType2);

			/* Compute epsilon_ij and sigma_ij
			*/
			epsilon_ij = sqrt(epsilon_ca * forceFieldMatrix[resNum2][3]);
			sigma_ij = r_gyr_ca + forceFieldMatrix[resNum2][1];

			/* Add to potential.
			*/
			res[i] += 4*epsilon_ij*( pow( sigma_ij/NDDist(cAlpha1, cBeta2) ,12) - pow( sigma_ij/NDDist(cAlpha1, cBeta2) ,6) );
		}
	}

	/* DEBUG
	 *
	 cout << ", +Ca-Cb: " << ljPot;
	 */

	/* Calculate the non-bonded Cb-Cb potentials.
	*/
	for(i = 0; i < seqLength; i++){

		/* Use Ca bead for glycines (no Cb bead
		 * exists).
		 */
		if(protStruct->getSingleResidue(i).getOneLetterType().compare("G") == 0){

			/* DEBUG
			 *
			cout << "#1 is Glycine bead, special treatment." << endl;
			*/

			cBeta1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		}
		else{
			cBeta1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");
		}

		resNum1 = forceFieldKey.find(protStruct->getSingleResidue(i).getOneLetterType());

		for(j = i + SEQ_SEPARATION; j < seqLength; j++){

			/* Interaction with the other Cb bead, unless it's a
			 * glycine (use Ca bead as above)
			 */
			resType2 = protStruct->getSingleResidue(j).getOneLetterType();
			if(resType2.compare("G") == 0){

				/* DEBUG
				 *
				cout << "#2 is Glycine bead, special treatment." << endl;
				*/

				cBeta2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CA");
			}
			else{
				cBeta2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CB");
			}

			/* DEBUG
			 *
			 cout << "Cb " << i+1 << "-Cb " << j+1 << endl;
			 */

			resNum2 = forceFieldKey.find(resType2);

			/* Compute sigma_ij and epsilon_ij
			*/
			sigma_ij = forceFieldMatrix[resNum1][1] + forceFieldMatrix[resNum2][1];
			epsilon_ij = sqrt(forceFieldMatrix[resNum1][3] * forceFieldMatrix[resNum2][3]);

			/* Add to potential.
			*/
			res[i] += 2*epsilon_ij*( pow( sigma_ij/NDDist(cBeta1, cBeta2) ,12) - pow( sigma_ij/NDDist(cBeta1, cBeta2) ,6) );
			res[j] += 2*epsilon_ij*( pow( sigma_ij/NDDist(cBeta1, cBeta2) ,12) - pow( sigma_ij/NDDist(cBeta1, cBeta2) ,6) );
		}
	}

	/* DEBUG
	 *
	 cout << ", + Cb-Cb: " << ljPot << endl;
	 */

	/* Scale the potential and return.
	*/
	return res;
}

vector<float> MukherjeeModel::calcHelixPotential_pos(shared_ptr<Structure> protStruct){

	/* V_helix = sum_i=3toN-3[ 0.5*K_1-3_i*(r_i,i+2-r_h)^2 + 0.5*K_1-4_i*(r_i,i+3-r_h)^2 ]
	 *
	 * where
	 *
	 * r_h = equilibrium distance in helix = 5.5 A
	 * K_1-3_i = (1/3)*(K_i + K_i+1 + K_i+2)
	 * K_1-4_i = (1/4)*(K_i + K_i+1 + K_i+2 + K_i+3)
	 * K_i = scaled helix propensity in kJ/mol/A^2 for residue i
	 *
	 * K's here are actually big Kappas.
	 */
	float r_h = 5.5;
	float helixPot = 0.0;
	float K_1_3, K_1_4;
	int i;

	string resType1, resType2, resType3, resType4;
	int resNum1, resNum2, resNum3, resNum4;

	int seqLength = protStruct->getNumResidues();
	vector<float>res(seqLength,0.0);

	/* Sum from 3 to N-3.
	 *
	 * REALLY START THIS AT POSITION 3? WHY?
	 * NOT IMMEDIATELY APPARENT IN PAPER.
	 */
	for(i = 2; i < seqLength-3; i++){
		/* Calculate the first term.
		 */
		resType1 = protStruct->getSingleResidue(i).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+2).getOneLetterType();

		/* DEBUG
		 *
		cout << "Residue type: " << resType1 << endl;
		*/

		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);

		K_1_3 = (1.0/3.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4]);

		coord cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		coord cAlpha2 = protStruct->getSingleResidue(i+2).getAtomCoordsByType("CA");

//added for the new scale without distance
		res[i] += 0.5*K_1_3;

//orig		res[i] += 0.5*K_1_3*pow( NDDist(cAlpha1, cAlpha2) - r_h ,2);

		/* Calculate the second term.
		 */
		resType1 = protStruct->getSingleResidue(i).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+2).getOneLetterType();
		resType4 = protStruct->getSingleResidue(i+3).getOneLetterType();

		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);
		resNum4 = forceFieldKey.find(resType4);

		K_1_4 = (1.0/4.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4] + forceFieldMatrix[resNum4][4]);

		cAlpha2 = protStruct->getSingleResidue(i+3).getAtomCoordsByType("CA");
//added for the new scale without distance
		res[i] += 0.5*K_1_4;

//orig		res[i] += 0.5*K_1_4*pow( NDDist(cAlpha1, cAlpha2) - r_h ,2);
	}

	/* Return the final potential.
	 */
	return res;
}

vector<float> MukherjeeModel::calcHelixPotential_pos2(shared_ptr<Structure> protStruct){

	/* V_helix = sum_i=3toN-3[ 0.5*K_1-3_i*(r_i,i+2-r_h)^2 + 0.5*K_1-4_i*(r_i,i+3-r_h)^2 ]
	 *
	 * where
	 *
	 * r_h = equilibrium distance in helix = 5.5 A
	 * K_1-3_i = (1/3)*(K_i + K_i+1 + K_i+2)
	 * K_1-4_i = (1/4)*(K_i + K_i+1 + K_i+2 + K_i+3)
	 * K_i = scaled helix propensity in kJ/mol/A^2 for residue i
	 *
	 * K's here are actually big Kappas.
	 */
	float r_h = 5.5;
	float helixPot = 0.0;
	float K_1_3, K_1_4;
	int i;

	string resType1, resType2, resType3, resType4;
	int resNum1, resNum2, resNum3, resNum4;

	int seqLength = protStruct->getNumResidues();
	vector<float>res(seqLength,0.0);

	/* Sum from 3 to N-3.
	 *
	 * REALLY START THIS AT POSITION 3? WHY?
	 * NOT IMMEDIATELY APPARENT IN PAPER.
	 */
	for(i = 2; i < seqLength-3; i++){
		/* Calculate the first term.
		 */
		resType1 = protStruct->getSingleResidue(i).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+2).getOneLetterType();

		/* DEBUG
		 *
		cout << "Residue type: " << resType1 << endl;
		*/

		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);

//		K_1_3 = (1.0/3.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4]);

		coord cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
		coord cAlpha2 = protStruct->getSingleResidue(i+2).getAtomCoordsByType("CA");

		float k = 0.5*(1.0/3.0)*pow( NDDist(cAlpha1, cAlpha2) - r_h ,2);

		res[i]+=forceFieldMatrix[resNum1][4]*k;
		res[i+1]+=forceFieldMatrix[resNum2][4]*k;
		res[i+2]+=forceFieldMatrix[resNum3][4]*k;

		/* Calculate the second term.
		 */
		resType1 = protStruct->getSingleResidue(i).getOneLetterType();
		resType2 = protStruct->getSingleResidue(i+1).getOneLetterType();
		resType3 = protStruct->getSingleResidue(i+2).getOneLetterType();
		resType4 = protStruct->getSingleResidue(i+3).getOneLetterType();

		resNum1 = forceFieldKey.find(resType1);
		resNum2 = forceFieldKey.find(resType2);
		resNum3 = forceFieldKey.find(resType3);
		resNum4 = forceFieldKey.find(resType4);

//		K_1_4 = (1.0/4.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4] + forceFieldMatrix[resNum4][4]);

		cAlpha2 = protStruct->getSingleResidue(i+3).getAtomCoordsByType("CA");

		k = 0.5*(1.0/4.0)*pow( NDDist(cAlpha1, cAlpha2) - r_h ,2);

		res[i]+=forceFieldMatrix[resNum1][4]*k;
		res[i+1]+=forceFieldMatrix[resNum2][4]*k;
		res[i+2]+=forceFieldMatrix[resNum3][4]*k;
		res[i+3]+=forceFieldMatrix[resNum4][4]*k;
	}

	/* Return the final potential.
	 */
	return res;
}




/* Calculating interaction score/energy.
 *
 * Dummy fuction, always returns "1e10". Override in subclasses.
 */
float MukherjeeModel::calcInteractionScore(shared_ptr<Structure> interactorOne, shared_ptr<Structure> interactorTwo){
	return 10000000000.0;
}
