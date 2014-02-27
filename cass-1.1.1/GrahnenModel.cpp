#include "GrahnenModel.h"
#include "Structure.h"
#include "TwoBeadStructure.h"
#include <typeinfo>
#include <stdexcept>
#include <math.h>
#include <limits>
#include <functional>

/* Helper functor object to provide random numbers for
 * contact map shuffling. See docs for random_shuffle()
 * for more details.
 * 
 * Adapted from 
 * http://stackoverflow.com/questions/147391/using-boostrandom-as-the-rng-for-stdrandom-shuffle
 *
 * ABSTRACT THIS OUT AND INCLUDE IT?
 */
struct myRng : std::unary_function<unsigned, unsigned> {
	CRandomMersenne& _state;
	unsigned operator()(unsigned i) {
		return _state.IRandom(0,i-1);
	}
	myRng(CRandomMersenne& state) : _state(state) {}
};

/* Default constructor, initializes all weights
 * to 1.0.
 */
GrahnenModel::GrahnenModel():
	RastogiModel()
{
	w_bend = 1.0;
	w_LJ = 1.0;
	w_helix = 1.0;
	w_beta = 1.0;
	w_charge = 1.0;
	w_solv = 1.0;
	w_cys = 1.0;

	w_bind_LJ = 1.0;
	w_bind_charge = 1.0;
	w_bind_solv = 1.0;
}

/* Constructor to specify default weightings of individual terms in the energy
 * function.
 */
GrahnenModel::GrahnenModel(float bendWeight, float LJWeight, float helixWeight, float sheetWeight, float chargeWeight, float solvWeight, float cysWeight, float LJWeightBind, float chargeWeightBind, float solvWeightBind):
	RastogiModel()
{
	w_bend = bendWeight;
	w_LJ = LJWeight;
	w_helix = helixWeight;
	w_beta = sheetWeight;
	w_charge = chargeWeight;
	w_solv = solvWeight;
	w_cys = cysWeight;
	
	w_bind_LJ = LJWeightBind;
	w_bind_charge = chargeWeightBind;
	w_bind_solv = solvWeightBind;
}

/* Destructor, required due to gcc being
 * finicky.
 */
GrahnenModel::~GrahnenModel(){
}

/* Calculate folding score/energy with the current parameter weights.
 */
float GrahnenModel::calcFoldScore(shared_ptr<Structure> protein){

	/* DEBUG
	 *
	cout << "Checking for correct type..." << endl;
	*/

	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only TwoBeadStructure is valid for GrahnenModel.");
	}

	/* DEBUG
	 *
	cout << "All's well with the type, proceeding to calculate folding score." << endl;
	*/

	return getFoldScoreByParams(protein, w_bend, w_LJ, w_helix, w_beta, w_charge, w_solv, w_cys);
}

/* Calculate interaction score/energy with the current parameter weights.
 */
float GrahnenModel::calcInteractionScore(shared_ptr<Structure> interactorOne, shared_ptr<Structure> interactorTwo){

	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*interactorOne) || typeid(test) != typeid(*interactorTwo)){
		throw runtime_error("Only TwoBeadStructure is valid for GrahnenModel.");
	}

	return getBindScoreByParams(interactorOne, interactorTwo, w_bind_LJ, w_bind_charge, w_bind_solv);
}

/* Calculate individual terms in the folding energy/scoring function.
 *
 * NOT DOING A TYPE CHECK HERE IS A LITTLE DANGEROUS...
 */
vector<float> GrahnenModel::getFoldScoreTerms(shared_ptr<Structure> protein){
	vector<float> values;
	int i;

	/* Re-cast the pointer.
	 */
	shared_ptr<TwoBeadStructure> prot = dynamic_pointer_cast<TwoBeadStructure>(protein);

	/* Check that exposure values exist, otherwise throw an exception.
	 */
	int protSize = prot->getNumResidues();
	vector<float> noCalcs(protSize, 0.0);
	vector<float> exposures;
	for(i = 0; i < protSize; i++){
		exposures.push_back(prot->getSingleResidue(i).getApoExposure());
	}
	if(noCalcs == exposures){
		throw runtime_error("No apo-exposure values for residues! Pre-calculate with TwoBeadStructure::calcExposure().");
	}
	
	/* Start the clock.
	 *
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);
	double timeStop;
	bool clockRunning = true;
	*/

	/* Calculate the bending potential.
	 */
	values.push_back(calcSidechainBendPotential(prot));
	
	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_LJ == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(bend)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the non-bonding (Lennard-Jones) potential.
	 */
	values.push_back(calcNonbondPotential(prot));
	
	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_helix == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(LJ)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the helix potential.
	 */
	values.push_back(calcHelixPotential(prot));
	
	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_beta == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(helix)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the beta-sheet potential.
	 */
	values.push_back(calcBetaPotential(prot));
	
	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_charge == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(beta)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the electrostatic potential 	 
	 */
	values.push_back(calcChargePotential(prot));
	
	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_solv == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(charge)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate solvation potential.
	 */
	values.push_back(calcSolvationPotential(prot, false));
	
	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_cys == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(solv)" << endl;
		clockRunning = false;
	}
	*/

	// AND THE CYSTEINE ONE AS THAT IS DEVELOPED.
	values.push_back(calcCysteinePotential(prot));
	
	/* Alway stop the clock here.
	 *
	if(clockRunning){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(cys)" << endl;
		clockRunning = false;
	}
	*/

	return values;
}

/* Calculate individual terms in the interaction energy/scoring function.
 *
 * NEW AND IMPROVED WITH 6 A CAP.
 *
 * NOT DOING A TYPE CHECK BEFORE CASTING THE POINTERS IS A TOUCH
 * DANGEROUS...
 */
vector<float> GrahnenModel::getBindScoreTerms(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo){

	/* V_lj, non-bonding potential:
	 * 
	 * The Lennard-Jones potential between two atoms.
	 * 
	 * V_lj = 4*sum_i,j[epsilon_ij( (sigma_ij/r_ij)^12 - (sigma_ij/r_ij)^6 )]
	 *
	 * where
	 *
	 * i,j = every possible Ca-Cb and Cb-Cb combination between proteins.
	 *
	 * epsilon_ij = sqrt(epsilon_ii*epsilon_jj), epsilon_ii in Table 1 or 0.05 kJ/mol for Ca beads
	 * 
	 * sigma_ij = collision diameter for i,j, i.e. the twice the r_gyr (see Table 1 or Ca = 1.8 A)
	 *
	 * r_ij = dist(i,j) (WITH A CAP OF 6 A)
	 */
	float epsilon_ca = 0.05;
	float r_gyr_ca = 1.8;
	float r_max = 6.0;
	float bindEng = 0.0;
	
	/* Sum of force due to pairwise charge interactions is calculated as
	 *
	 * V_charge = sum_i,j[k_e * (q_i * q_j / r_ij)]
	 *
	 * where 
	 *
	 * i,j = every unique pair of charged residues
	 * k_e = 1 (Coloumbs constant, provisionally set to 1 instead of 1/(4*pi*epsilon_0)
	 * q_i/q_j = charge of sidechain i/j
	 * r_ij = distance between centers of Cb beads for residues i and j
	 */
	float k_e = 1.0; 
	float chargeScaleFactor = 10000.0; // EXPERIMENTAL

	float bsaScaleFactor = 10.0; // EXPERIMENTAL

	int seqLength = proteinTwo->getNumResidues();
	int protSizeOne = proteinOne->getNumResidues();
	int protSizeTwo = proteinTwo->getNumResidues();

	int i, j, resNum1, resNum2;
	float epsilon_ij, sigma_ij;
	string resType1, resType2;
	coord cAlpha1, cAlpha2, cBeta1, cBeta2;

	float ljTerm = 0.0;
	float chargeTerm = 0.0;
	float buried = 0.0;
	float distance;
	vector<float> values;

	/* Re-cast the pointers.
	 */
	shared_ptr<TwoBeadStructure> protOne = dynamic_pointer_cast<TwoBeadStructure>(proteinOne);
	shared_ptr<TwoBeadStructure> protTwo = dynamic_pointer_cast<TwoBeadStructure>(proteinTwo);

	/* Check that exposure values exist, otherwise throw an exception.
	 */
	vector<float> noCalcsOne(protSizeOne, 0.0);
	vector<float> noCalcsTwo(protSizeTwo, 0.0);
	vector<float> apoExposures;
	vector<float> holoExposuresOne;
	vector<float> holoExposuresTwo;
	shared_ptr<TwoBeadStructure> newProtOne;
	shared_ptr<TwoBeadStructure> newProtTwo;
	vector< shared_ptr<TwoBeadStructure> > wholeComplex;

	/* Check and possibly calculate exposure of the parts of the
	 * complex.
	 */
	for(i = 0; i < protSizeOne; i++){
		apoExposures.push_back(protOne->getSingleResidue(i).getApoExposure());
		holoExposuresOne.push_back(protOne->getSingleResidue(i).getHoloExposure());
	}
	if(noCalcsOne == apoExposures){
		/* DEBUG
		 *
		cout << "Protein #1 doesn't have apo-exposure values, calculating..." << endl;
		*/
		throw runtime_error("No apo-exposure values for protein 1! Pre-calculate with TwoBeadStructure::calcExposure().");
	}
	apoExposures.clear();
	for(i = 0; i < protSizeTwo; i++){
		apoExposures.push_back(protTwo->getSingleResidue(i).getApoExposure());
		holoExposuresTwo.push_back(protTwo->getSingleResidue(i).getHoloExposure());
	}
	if(noCalcsTwo == apoExposures){
		/* DEBUG
		 *
		cout << "Protein #2 doesn't have apo-exposure values, calculating..." << endl;
		*/
		throw runtime_error("No apo-exposure values for protein 2! Pre-calculate with TwoBeadStructure::calcExposure().");
	}

	/* Check for and possibly calculate exposure in the assembled complex.
	 */
	if(noCalcsOne == holoExposuresOne || noCalcsTwo == holoExposuresTwo){

		/* DEBUG
		 *
		cout << "One or both of the proteins lacks holo-exposure values. ";
		if(noCalcsOne == holoExposuresOne){
			cout << "#1 is missing them. ";
		}
		if(noCalcsTwo == holoExposuresTwo){
			cout << "#2 is missing them. ";
		}
		*/
		throw runtime_error("No holo-exposure values for one or both of the proteins! Pre-calculate with TwoBeadStructure::calcComplexExposure().");
	}
	
	/* Start the clock.
	 *
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);
	double timeStop;
	bool clockRunning = true;
	*/

	/* Calculate the non-bonded Ca-Cb potentials.
	*/
	for(i = 0; i < seqLength; i++){

		cAlpha1 = protTwo->getSingleResidue(i).getAtomCoordsByType("CA");
		resNum1 = forceFieldKey.find(protTwo->getSingleResidue(i).getOneLetterType());

		for(j = 0; j < protOne->getNumResidues(); j++){

			/* Interaction with the other Cb bead, unless it's a
			 * glycine (avoids double-counting in the next step).
			 */
			resType2 = protOne->getSingleResidue(j).getOneLetterType();

			if(resType2.compare("G") == 0){
				/* DEBUG
				 *
				cout << "#2 is a glycine, skip it." << endl;
				*/

				continue;
			}
			else{
				cBeta2 = protOne->getSingleResidue(j).getAtomCoordsByType("CB");
			}

			/* DEBUG
			 *
			cout << "Ca " << i+1 << "-Cb " << j+1 << endl;
			*/

			resNum2 = forceFieldKey.find(resType2);

			/* Compute epsilon_ij and sigma_ij and the distance
			*/
			epsilon_ij = sqrt(epsilon_ca * forceFieldMatrix[resNum2][3]);
			sigma_ij = r_gyr_ca + forceFieldMatrix[resNum2][1];
			distance = NDDist(cAlpha1, cBeta2);

			/* Add to potential if they close enough.
			*/
			if(distance < r_max){
				bindEng += 4*epsilon_ij*( pow( sigma_ij/distance ,12) - pow( sigma_ij/distance ,6) );

				ljTerm += 4*epsilon_ij*( pow( sigma_ij/distance ,12) - pow( sigma_ij/distance ,6) );
			}
		}
	}

	/* DEBUG
	 *
	cout << ", +Ca-Cb: " << bindEng;
	*/

	/* Calculate the non-bonded Cb-Cb potentials and any charged interactions.
	 */
	for(i = 0; i < seqLength; i++){

		resType1 = protTwo->getSingleResidue(i).getOneLetterType();

		/* Use Ca beads for glycines (no Cb bead
		 * exists).
		 */
		if(resType1.compare("G") == 0){
			/* DEBUG
			 *
			cout << "#1 is Glycine, special case." << endl;
			*/

			cBeta1 = protTwo->getSingleResidue(i).getAtomCoordsByType("CA");
		}
		else{
			cBeta1 = protTwo->getSingleResidue(i).getAtomCoordsByType("CB");
		}

		resNum1 = forceFieldKey.find(resType1);

		for(j = 0; j < protOne->getNumResidues(); j++){

			/* Interaction with the other Cb bead, unless it's a
			 * glycine (use Ca bead instead)
			 */
			resType2 = protOne->getSingleResidue(j).getOneLetterType();

			if(resType2.compare("G") == 0){
				/* DEBUG
				 *
				cout << "#2 is glycine, special case..." << endl;
				*/

				cBeta2 = protOne->getSingleResidue(j).getAtomCoordsByType("CA");
			}
			else{
				cBeta2 = protOne->getSingleResidue(j).getAtomCoordsByType("CB");
			}

			/* DEBUG
			 *
			cout << "Cb " << i+1 << "-Cb " << j+1 << endl;
			*/

			resNum2 = forceFieldKey.find(resType2);

			/* Compute sigma_ij and epsilon_ij and distance
			*/
			sigma_ij = forceFieldMatrix[resNum1][1] + forceFieldMatrix[resNum2][1];
			epsilon_ij = sqrt(forceFieldMatrix[resNum1][3] * forceFieldMatrix[resNum2][3]);
			distance = NDDist(cBeta1, cBeta2);

			/* Add to potential if they're close enough.
			*/
			if(distance < r_max){
				bindEng += 4*epsilon_ij*( pow( sigma_ij/distance ,12) - pow( sigma_ij/distance ,6) );

				ljTerm += 4*epsilon_ij*( pow( sigma_ij/distance ,12) - pow( sigma_ij/distance ,6) );

				/* Are both residues charged (and this term in use)? If so, calculate that force.
				*/
				if(w_bind_charge > 0.0 && forceFieldMatrix[resNum1][6] != 0 && forceFieldMatrix[resNum2][6] != 0){

					k_e = chargeScaleFactor*(1/dielectricConstant(protOne->getSingleResidue(j).getHoloExposure(), protTwo->getSingleResidue(i).getHoloExposure()));

					bindEng += k_e * ( forceFieldMatrix[resNum1][6] * forceFieldMatrix[resNum2][6] ) / distance;

					chargeTerm += k_e * ( forceFieldMatrix[resNum1][6] * forceFieldMatrix[resNum2][6] ) / distance;
					
					/* DEBUG 
					 *
					cout << (i+1) << resType1 << "-" << (j+1) << resType2 << " (" << protOne->getSingleResidue(j).getHoloExposure() << "," << protTwo->getSingleResidue(i).getHoloExposure() << ") -> k_e = " << k_e << endl;
					*/
				}
			}
		}
	}
	
	/* Stop the clock and print, if necessary.
	 *
	if(clockRunning && (w_bind_LJ > 0.0 && w_bind_charge == 0.0 && w_bind_solv == 0.0) || (w_bind_LJ > 0.0 && w_bind_charge > 0.0 && w_bind_solv == 0.0)){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(LJ or charge)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the change in solvation potential on complex
	 * formation.
	 */
	buried = bsaScaleFactor*calcSolvationChange(protOne, protTwo);
	
	/* Always stop the clock here.
	 *
	if(clockRunning){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(dSolv)" << endl;
	}
	*/

	/* DEBUG
	 *
	cout << "-1*(BSA): " << buried << endl;
	*/
	
	values.push_back(ljTerm);
	values.push_back(chargeTerm);
	values.push_back(buried);

	return values;
}

/* Calculate the energy/scoring functions with arbitrary
 * parameters.
 */
float GrahnenModel::getFoldScoreByParams(shared_ptr<Structure> protein, float bendWeight, float LJWeight, float helixWeight, float sheetWeight, float chargeWeight, float solvWeight, float cysWeight){

	/* Calculate the raw values of each potential term.
	 */
	vector<float> terms = getFoldScoreTerms(protein);

	/* .. and scale them.
	 */
	float VBend = bendWeight*terms[0];

	float VNonbond = LJWeight*terms[1];

	// 0.1 ISN'T A BAD DEFAULT
	float VHelix = helixWeight*terms[2];

	// 0.01 IS AN OKAY DEFAULT
	float VBeta = sheetWeight*terms[3];

	// 1000 SEEMS LIKE AN OKAY DEFAULT
	float VCharge = chargeWeight*terms[4];
	
	float VSolv = solvWeight*terms[5];

	// DO CYSTEINE CALCULATIONS HERE AS NECESSARY
	float VCys = cysWeight*terms[6];

	/* DEBUG
	 *
	cerr << "VBend\tVNonbond\tVHelix\tVBeta\tVCharge\tVSolv\tVCys" << endl;
	cerr << VBend << "\t" << VNonbond << "\t" << VHelix << "\t" << VBeta << "\t" << VCharge << "\t" << VSolv << "\t" << VCys << endl;
	*/

	/* Add up to the total potential and return.
	 */
	return VBend + VNonbond + VHelix + VBeta + VCharge + VSolv + VCys;
}

/* Calculate the interaction score with arbitrary parameters.
 */
float GrahnenModel::getBindScoreByParams(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo, float LJWeight, float chargeWeight, float solvWeight){
	

	/* Calculate the raw terms of the energy function.
	 */
	vector<float> terms = getBindScoreTerms(proteinOne, proteinTwo);

	/* ... and scale them.
	 */
	float VLJ = LJWeight*terms[0];

	float VCharge = chargeWeight*terms[1];

	float VSolv = solvWeight*terms[2];

	/* DEBUG
	 *
	cout << "VLJ\tVCharge\tVSolv" << endl;
	cout << VLJ << "\t" << VCharge << "\t" << VSolv << endl;
	*/

	return VLJ + VCharge + VSolv;
}

/* Calculating the fold score relative to the unfolded
 * state.
 *
 * Specifically, calculate the quantity
 *
 * dG = G_NS + (var(G) - 2*mean(G))/2
 *
 * where G_NS is the folding score in the supplied
 * sequence+conformation, and G is the distribution of
 * scores in the misfolded compact non-native conformations
 * (as appximated by the Random Energy Model).
 *
 * See Goldstein (2011, Proteins 79) for the formula
 * itself, and Bryngelson and Wolynes (1987, PNAS 84)
 * as well as Goldstein et al (1992, PNAS 89) for the
 * approximation of the non-native conformations.
 *
 * STILL HIGHLY EXPERIMENTAL...
 *
 * PROBABLY WON'T END UP USING GAUSSIAN APPROX.
 */
float GrahnenModel::calcFoldEnergyGap(shared_ptr<Structure> protein){
	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only TwoBeadStructure is valid for GrahnenModel.");
	}
	shared_ptr<TwoBeadStructure> prot = dynamic_pointer_cast<TwoBeadStructure>(protein);

	/* Useful variables.
	 */
	int i, j, resNum, resNum1, resNum2, resNum3, resNum4;
	string resType, resType1, resType2, resType3, resType4;
	float bbDist, cBDistance, cAlphaDist, q_1, q_2;
	float torsion, exposure, randBend, randLj, randHelix, randBeta;
	float randCharge, randSolv, randSS, natScore, gapScore;

	/* Some pertinent information.
	 */
	int seqLength = protein->getNumResidues();
	string sequence = protein->getSequence();

	/* Distributions of native structure geometry.
	 */
	vector< pair<float, float> > bondAngles = prot->bondAngles();
	vector<float> resDists = prot->cbCbDistances();
	vector<float> resBBDists = prot->caCbDistances();
	vector< pair<float, float> > backboneDists = prot->backboneDistances();
	vector<float> backboneTorsions = prot->backboneTorsions();
	vector<float> exposures = prot->residueExposures();
	vector<string> threeMers = prot->threeMers();
	vector<string> helixFourMers = prot->helixFourMers();
	vector<string> betaFourMers = prot->betaFourMers();

	vector<float> randEngDist;
	vector<float> quarts;
	vector<float> sampEngTerms;
	
	/* Bond angle term.
	 */	
	float K_theta = 10.0;
	float bendPot = 0.0;
	for(i = 0; i < bondAngles.size(); i++){
		resType = sequence.substr(i,1);
		resNum = forceFieldKey.find(resType);

		/* No defined angles for glycines.
		 */
		if(resType.compare("G") == 0){
			continue;
		}

		/* Most residues have two defined angles.
		 */
		if(i != 0 && i != sequence.length()-1){
			bendPot += pow((bondAngles[i].first - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
			bendPot += pow((bondAngles[i].second - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
		/* N-terminal residues lack the first angle.
		 */
		else if(i == 0){
			bendPot += pow((bondAngles[i].second - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
		/* C-terminal residues lack the second angle.
		 */
		else if(i == sequence.length()-1){
			bendPot += pow((bondAngles[i].first - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
	}
	/* And scale it.
	 */
	bendPot *= 0.5*K_theta*w_bend;

	/* DEBUG
	 *
	cerr << "Bending potential: " << bendPot << endl;
	*/

	/* Pair-wise terms: Lennard-Jones, and charge and disulfide bond prediction.
	 */
	float epsilon_ij, sigma_ij, k_e;
	float epsilon_ca = 0.05;
	float r_gyr_ca = 1.8;
	double chargeScaleFactor = 1000.0; // Scales charge term to approx same magnitude as LJ
	float maxSSDistance = 4.5; // Threshold for disulfide bond distance
	float ljPot = 0.0;
	float chargePot = 0.0;
	float ssbondPot = 0.0;
	unsigned int counter = 0;
	/* Non-bonded Ca-Cb potentials.
	 */
	for(i = 0; i < seqLength; i++){
		resNum1 = forceFieldKey.find(sequence[i]);

		/* Only count Ca-Cb distances that are separated by at least
		 * a certain number of amino acids.
		 * 
		 * 2013-05-22: Fixed a bug in not counting all the distances
		 * (Nadia Bykova).
		 */
		for(j = 0; j < sequence.length(); j++){

			if ((i >= j && (i - j) < SEQ_SEPARATION) || (i < j && (j - i) < SEQ_SEPARATION)){
				continue;
			}
			/* End of edit
			 */

			resType2 = sequence[j];
			
			/* Don't treat glycines (avoids double-counting in the next step).
			 */
			if(resType2.compare("G") == 0){
				counter++;
				continue;
			}

			resNum2 = forceFieldKey.find(resType2);

			/* Compute epsilon_ij and sigma_ij and distance.
			 */
			epsilon_ij = sqrt(epsilon_ca * forceFieldMatrix[resNum2][3]);
			sigma_ij = r_gyr_ca + forceFieldMatrix[resNum2][1];
			bbDist = resBBDists[counter];

			/* Add to potential.
			 */
			ljPot += epsilon_ij*( pow( sigma_ij/bbDist ,12) - pow( sigma_ij/bbDist ,6) );

			counter++;
		}
	}
	/* Non-bonded Cb-Cb potentials (all three).
	 */
	counter = 0;
	for(i = 0; i < seqLength; i++){
		resType1 = sequence.substr(i,1);
		resNum1 = forceFieldKey.find(resType1);

		for(j = i + SEQ_SEPARATION; j < seqLength; j++){
			resType2 = sequence[j];
			resNum2 = forceFieldKey.find(resType2);

			/* Compute sigma_ij, epsilon_ij, charges and the distance.
			*/
			sigma_ij = forceFieldMatrix[resNum1][1] + forceFieldMatrix[resNum2][1];
			epsilon_ij = sqrt(forceFieldMatrix[resNum1][3] * forceFieldMatrix[resNum2][3]);
			q_1 = forceFieldMatrix[resNum1][6];
			q_2 = forceFieldMatrix[resNum2][6];
			cBDistance = resDists[counter];

			ljPot += epsilon_ij*( pow( sigma_ij/cBDistance ,12) - pow( sigma_ij/cBDistance ,6) );

			/* Are both residues cysteines? If so, might have a disulfide bond.
			 */
			if(resType1.compare("C") == 0 && resType2.compare("C") == 0){
				if(cBDistance < maxSSDistance){
					ssbondPot -= 1.0;
				}
			}
			/* If both are charged, go head and calculate that interaction too.
			 */
			else if(q_1 != 0 && q_2 != 0){
				k_e = chargeScaleFactor*(1/dielectricConstant(exposures[i], exposures[j]));
				chargePot += k_e*((q_1*q_2)/cBDistance);
			}

			counter++;
		}
	}
	/* Scale the terms.
	 */
	ljPot *= 4*w_LJ;
	chargePot *= w_charge;
	ssbondPot *= w_cys;

	/* DEBUG
	 *
	cerr << "LJ potential: " << ljPot << endl;
	*/

	/* Helix term.
	 */
	float K_1_3, K_1_4;
	float r_h = 5.5;
	float helixPot = 0.0;
	for(i = 0; i < threeMers.size(); i++){
		/* Calculate the first half.
		 */
		resNum1 = forceFieldKey.find(threeMers[i][0]);
		resNum2 = forceFieldKey.find(threeMers[i][1]);
		resNum3 = forceFieldKey.find(threeMers[i][2]);

		K_1_3 = (1.0/3.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4]);
		cAlphaDist = backboneDists[i].first;


//Added to remove term values that are not under helices and distance
///		if ((i >= 12 && i < 22) || (i >= 79 && i < 86))
///		{
			helixPot += K_1_3;
///	cout << "calcFoldEnergyGap::GrahnenModel checkpoint" << endl;
///		}


//orig		helixPot += K_1_3*pow(cAlphaDist - r_h ,2);

		/* And the second half.
		 */
		resNum1 = forceFieldKey.find(helixFourMers[i][0]);
		resNum2 = forceFieldKey.find(helixFourMers[i][1]);
		resNum3 = forceFieldKey.find(helixFourMers[i][2]);
		resNum4 = forceFieldKey.find(helixFourMers[i][3]);

		K_1_4 = (1.0/4.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4] + forceFieldMatrix[resNum4][4]);
		cAlphaDist = backboneDists[i].second;

//Added to remove term values that are not under helices and distance
///		if ((i >= 12 && i < 22) || (i >= 79 && i < 86))
///		{
			helixPot += K_1_4;
///		}


//orig		helixPot += K_1_4*pow( cAlphaDist - r_h ,2);
	}
	helixPot *= 0.5*w_helix;

	/* DEBUG
	 *
	cerr << "Helix potential: " << helixPot << endl;
	*/

	/* Beta term.
	 */
	float betaPot = 0.0;
	float phi_b = -150; // Equilibrium Ca-Ca torsion angle (Ca_i-1->Ca_i->Ca_i+1->Ca+i+2 torsion) for beta sheets. Equivalent to 210 positive degrees.
	float C_s = 0.01; // Scaling constant to make similar to calcHelixPotential()
	for(i = 0; i < betaFourMers.size(); i++){
		resNum1 = forceFieldKey.find(betaFourMers[i][0]);
		resNum2 = forceFieldKey.find(betaFourMers[i][1]);
		resNum3 = forceFieldKey.find(betaFourMers[i][2]);
		resNum4 = forceFieldKey.find(betaFourMers[i][3]);
		
		K_1_4 = 0.25*( forceFieldMatrix[resNum1][5] + forceFieldMatrix[resNum2][5] + forceFieldMatrix[resNum3][5] + forceFieldMatrix[resNum4][5]);
		torsion = backboneTorsions[i];
		
		betaPot += K_1_4*pow( C_s*(torsion - phi_b), 2);
	}
	betaPot *= w_beta;

	/* DEBUG
	 *
	cerr << "Beta potential: " << betaPot << endl;
	cerr << "Charge potential: " << chargePot << endl;
	*/

	/* Solvation term.
	 */
	float solvPot = 0.0;
	for(i = 0; i < exposures.size(); i++){
		exposure = exposures[i];
		resNum = forceFieldKey.find(sequence[i]);

		/* Weight by hydrophobicity and polarity.
		*/
		solvPot += forceFieldMatrix[resNum][3]*exposure + forceFieldMatrix[resNum][7]*(1-exposure);
	}
	solvPot *= w_solv;

	/* DEBUG
	 *
	cerr << "Native:";
	cerr << "\t" << bendPot;
	cerr << "\t" << ljPot;
	cerr << "\t" << helixPot;
	cerr << "\t" << betaPot;
	cerr << "\t" << chargePot;
	cerr << "\t" << solvPot;
	cerr << "\t" << ssbondPot;
	cerr << endl;
	cerr << "\tBend\tLJ\tHelix\tBeta\tCharge\tSolv\tSSBond" << endl;
	*/

	/* Sample a large number of randomly perturbed structural
	 * contexts with the same sequence. And make sure this
	 * sample is unique and constant for each unique sequence.
	*/
	std::hash<string> h;
	int seed = h(sequence);
	CRandomMersenne newRNG(seed);
	for(i = 0; i < NUM_RANDOM_STRUCTURES; i++){
		sampEngTerms = sampleFoldGapTerms(sequence, bondAngles, resBBDists, resDists, threeMers, helixFourMers, backboneDists, betaFourMers, backboneTorsions, exposures, newRNG);

		randBend = w_bend * sampEngTerms[0];

		//cerr << "\t" << randBend;

		randLj = w_LJ * sampEngTerms[1];
		
		//cerr << "\t" << randLj;

		randHelix = w_helix * sampEngTerms[2];
	
		//cerr << "\t" << randHelix;

		randBeta = w_beta * sampEngTerms[3];

		//cerr << "\t" << randBeta;

		randCharge = w_charge * sampEngTerms[4];

		//cerr << "\t" << randCharge;

		randSolv = w_solv * sampEngTerms[5];

		//cerr << "\t" << randSolv;

		randSS = w_cys * sampEngTerms[6];

		//cerr << "\t" << randSS;

		//cerr << endl;

		randEngDist.push_back(randBend + randLj + randHelix + randBeta + randCharge + randSolv + randSS);
	}

	/* Calculate the gap between the distribution and
	 * the native value, as per Goldstein.
	 *
	 * OR SOME SORT OF OTHER Z-SCORE? HOW ABOUT
	 * THE HALF-IQR SCORE?
	 */
	natScore = bendPot + ljPot + helixPot + betaPot + chargePot + solvPot + ssbondPot;
	quarts = quartiles(randEngDist);
	gapScore = natScore - quarts[2];

	/* DEBUG
	 *
	cerr << "dG = " << natScore << " - " << quarts[2] << " = " << gapScore << endl;
	*/

	return gapScore;
}

// SEE THE OTHER METHOD SIGNATURE. THIS JUST ALLOWS FOR A
// CONSTANT NON-NATIVE SOURCE OF SAMPLED GEOMETRY.
float GrahnenModel::calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation){
	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein) || typeid(test) != typeid(*compactConformation)){
		throw runtime_error("Only TwoBeadStructure is valid for GrahnenModel.");
	}
	shared_ptr<TwoBeadStructure> prot = dynamic_pointer_cast<TwoBeadStructure>(protein);
	shared_ptr<TwoBeadStructure> compact = dynamic_pointer_cast<TwoBeadStructure>(compactConformation);

	/* Useful variables.
	 */
	int i, j, resNum, resNum1, resNum2, resNum3, resNum4;
	string resType, resType1, resType2, resType3, resType4;
	float bbDist, cBDistance, cAlphaDist, q_1, q_2;
	float torsion, exposure, randBend, randLj, randHelix, randBeta;
	float randCharge, randSolv, randSS, natScore, gapScore;

	/* Some pertinent information.
	 */
	int seqLength = prot->getNumResidues();
	string sequence = prot->getSequence();

	/* Distributions of native structure geometry.
	 */
	vector< pair<float, float> > bondAngles = prot->bondAngles();
	vector<float> resDists = prot->cbCbDistances();
	vector<float> resBBDists = prot->caCbDistances();
	vector< pair<float, float> > backboneDists = prot->backboneDistances();
	vector<float> backboneTorsions = prot->backboneTorsions();
	vector<float> exposures = prot->residueExposures();
	vector<string> threeMers = prot->threeMers();
	vector<string> helixFourMers = prot->helixFourMers();
	vector<string> betaFourMers = prot->betaFourMers();

	vector<float> randEngDist;
	vector<float> quarts;
	vector<float> sampEngTerms;
	
	/* Bond angle term.
	 */	
	float K_theta = 10.0;
	float bendPot = 0.0;
	for(i = 0; i < bondAngles.size(); i++){
		resType = sequence.substr(i,1);
		resNum = forceFieldKey.find(resType);

		/* No defined angles for glycines.
		 */
		if(resType.compare("G") == 0){
			/* DEBUG
			 *
			cerr << "Residue " << (i+1) << " is a glycine, has no bond angle." << endl;
			*/

			continue;
		}

		/* DEBUG
		 *
		cerr << (i+1) << resType << ": " << bondAngles[i].first << "," << bondAngles[i].second << endl;
		*/

		/* Most residues have two defined angles.
		 */
		if(i != 0 && i != sequence.length()-1){
			bendPot += pow((bondAngles[i].first - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
			bendPot += pow((bondAngles[i].second - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
		/* N-terminal residues lack the first angle.
		 */
		else if(i == 0){
			bendPot += pow((bondAngles[i].second - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
		/* C-terminal residues lack the second angle.
		 */
		else if(i == sequence.length()-1){
			bendPot += pow((bondAngles[i].first - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
	}
	/* And scale it.
	 */
	bendPot *= 0.5*K_theta*w_bend;

	/* DEBUG
	 *
	cerr << "Bending potential: " << bendPot << endl;
	*/

	/* Pair-wise terms: Lennard-Jones, and charge and disulfide bond prediction.
	 */
	float epsilon_ij, sigma_ij, k_e;
	float epsilon_ca = 0.05;
	float r_gyr_ca = 1.8;
	double chargeScaleFactor = 1000.0; // Scales charge term to approx same magnitude as LJ
	float maxSSDistance = 4.5; // Threshold for disulfide bond distance
	float ljPot = 0.0;
	float chargePot = 0.0;
	float ssbondPot = 0.0;
	unsigned int counter = 0;
	/* Non-bonded Ca-Cb potentials.
	 */
	for(i = 0; i < seqLength; i++){
		resNum1 = forceFieldKey.find(sequence[i]);
		
		/* Use only Ca-Cb distances separated by a certain number of
		 * residues in the sequence.
		 *
		 * 2013-05-22: Fixed bug that missed counting some distances
		 * (Nadia Bykova).
		 */
		for(j = 0; j < sequence.length(); j++){
			if ((i >= j && (i - j) < SEQ_SEPARATION) || (i < j && (j - i) < SEQ_SEPARATION)){
				continue;
			}
			/* End of edit
			 */
		
			resType2 = sequence[j];

			/* Don't treat glycines (avoids double-counting in the next step).
			 */
			if(resType2.compare("G") == 0){
				counter++;
				continue;
			}

			resNum2 = forceFieldKey.find(resType2);

			/* Compute epsilon_ij and sigma_ij and distance.
			*/
			epsilon_ij = sqrt(epsilon_ca * forceFieldMatrix[resNum2][3]);
			sigma_ij = r_gyr_ca + forceFieldMatrix[resNum2][1];
			bbDist = resBBDists[counter];

			/* Add to potential.
			*/
			ljPot += epsilon_ij*( pow( sigma_ij/bbDist ,12) - pow( sigma_ij/bbDist ,6) );

			counter++;
		}
	}
	/* Non-bonded Cb-Cb potentials (all three).
	 */
	counter = 0;
	for(i = 0; i < seqLength; i++){
		resType1 = sequence.substr(i,1);
		resNum1 = forceFieldKey.find(resType1);

		for(j = i + SEQ_SEPARATION; j < seqLength; j++){
			resType2 = sequence[j];
			resNum2 = forceFieldKey.find(resType2);

			/* Compute sigma_ij, epsilon_ij, charges and the distance.
			 */
			sigma_ij = forceFieldMatrix[resNum1][1] + forceFieldMatrix[resNum2][1];
			epsilon_ij = sqrt(forceFieldMatrix[resNum1][3] * forceFieldMatrix[resNum2][3]);
			q_1 = forceFieldMatrix[resNum1][6];
			q_2 = forceFieldMatrix[resNum2][6];
			cBDistance = resDists[counter];

			ljPot += epsilon_ij*( pow( sigma_ij/cBDistance ,12) - pow( sigma_ij/cBDistance ,6) );

			/* Are both residues cysteines? If so, might have a disulfide bond.
			 */
			if(resType1.compare("C") == 0 && resType2.compare("C") == 0){
				if(cBDistance < maxSSDistance){
					ssbondPot -= 1.0;
				}
			}
			/* If both are charged, go head and calculate that interaction too.
			 */
			else if(q_1 != 0 && q_2 != 0){
				k_e = chargeScaleFactor*(1/dielectricConstant(exposures[i], exposures[j]));
				chargePot += k_e*((q_1*q_2)/cBDistance);
			}

			counter++;
		}
	}
	/* Scale the terms.
	 */
	ljPot *= 4*w_LJ;
	chargePot *= w_charge;
	ssbondPot *= w_cys;

	/* DEBUG
	 *
	cerr << "LJ potential: " << ljPot << endl;
	 */

	/* Helix term.
	 */
	float K_1_3, K_1_4;
	float r_h = 5.5;
	float helixPot = 0.0;
	for(i = 0; i < threeMers.size(); i++){
		/* Calculate the first half.
		 */
		resNum1 = forceFieldKey.find(threeMers[i][0]);
		resNum2 = forceFieldKey.find(threeMers[i][1]);
		resNum3 = forceFieldKey.find(threeMers[i][2]);

		K_1_3 = (1.0/3.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4]);
		cAlphaDist = backboneDists[i].first;

//Added to remove term values that are not under helices and distance
///		if ((i >= 12 && i < 22) || (i >= 79 && i < 86))
///		{
			helixPot += K_1_3;
//	cout << "GrahnenModel::calcFoldEnergyGap checkpoint" << endl;

///		}


//orig		helixPot += K_1_3*pow(cAlphaDist - r_h ,2);

		/* And the second half.
		 */
		resNum1 = forceFieldKey.find(helixFourMers[i][0]);
		resNum2 = forceFieldKey.find(helixFourMers[i][1]);
		resNum3 = forceFieldKey.find(helixFourMers[i][2]);
		resNum4 = forceFieldKey.find(helixFourMers[i][3]);

		K_1_4 = (1.0/4.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4] + forceFieldMatrix[resNum4][4]);
		cAlphaDist = backboneDists[i].second;

//Added to remove term values that are not under helices and distance
///		if ((i >= 12 && i < 22) || (i >= 79 && i < 86))
///		{
			helixPot += K_1_4;
///		}


//orig		helixPot += K_1_4*pow( cAlphaDist - r_h ,2);
	}
	helixPot *= 0.5*w_helix;

	/* DEBUG
	 *
	cerr << "Helix potential: " << helixPot << endl;
	 */

	/* Beta term.
	 */
	float betaPot = 0.0;
	float phi_b = -150; // Equilibrium Ca-Ca torsion angle (Ca_i-1->Ca_i->Ca_i+1->Ca+i+2 torsion) for beta sheets. Equivalent to 210 positive degrees.
	float C_s = 0.01; // Scaling constant to make similar to calcHelixPotential()
	for(i = 0; i < betaFourMers.size(); i++){
		resNum1 = forceFieldKey.find(betaFourMers[i][0]);
		resNum2 = forceFieldKey.find(betaFourMers[i][1]);
		resNum3 = forceFieldKey.find(betaFourMers[i][2]);
		resNum4 = forceFieldKey.find(betaFourMers[i][3]);
		
		K_1_4 = 0.25*( forceFieldMatrix[resNum1][5] + forceFieldMatrix[resNum2][5] + forceFieldMatrix[resNum3][5] + forceFieldMatrix[resNum4][5]);
		torsion = backboneTorsions[i];
		
		betaPot += K_1_4*pow( C_s*(torsion - phi_b), 2);
	}
	betaPot *= w_beta;

	/* DEBUG
	 *
	cerr << "Beta potential: " << betaPot << endl;
	cerr << "Charge potential: " << chargePot << endl;
	 */

	/* Solvation term.
	 */
	float solvPot = 0.0;
	for(i = 0; i < exposures.size(); i++){
		exposure = exposures[i];
		resNum = forceFieldKey.find(sequence[i]);

		/* Weight by hydrophobicity and polarity.
		*/
		solvPot += forceFieldMatrix[resNum][3]*exposure + forceFieldMatrix[resNum][7]*(1-exposure);
	}
	solvPot *= w_solv;

	/* DEBUG
	 *
	cerr << "Native:";
	cerr << "\t" << bendPot;
	cerr << "\t" << ljPot;
	cerr << "\t" << helixPot;
	cerr << "\t" << betaPot;
	cerr << "\t" << chargePot;
	cerr << "\t" << solvPot;
	cerr << "\t" << ssbondPot;
	cerr << endl;
	cerr << "\tBend\tLJ\tHelix\tBeta\tCharge\tSolv\tSSBond" << endl;
	 */

	/* Sample a large number of randomly perturbed structural
	 * contexts, conformation distinct from the native but 
	 * with the same sequence. And make sure this  sample is 
	 * unique and constant for each unique sequence.
	 */
	std::hash<string> h;
	int seed = h(sequence);
	CRandomMersenne newRNG(seed);

	bondAngles = compact->bondAngles();
	resDists = compact->cbCbDistances();
	resBBDists = compact->caCbDistances();
	backboneDists = compact->backboneDistances();
	backboneTorsions = compact->backboneTorsions();
	exposures = compact->residueExposures();

	for(i = 0; i < NUM_RANDOM_STRUCTURES; i++){
		sampEngTerms = sampleFoldGapTerms(sequence, bondAngles, resBBDists, resDists, threeMers, helixFourMers, backboneDists, betaFourMers, backboneTorsions, exposures, newRNG);

		randBend = w_bend * sampEngTerms[0];

		//cerr << "\t" << randBend;

		randLj = w_LJ * sampEngTerms[1];
		
		//cerr << "\t" << randLj;

		randHelix = w_helix * sampEngTerms[2];
	
		//cerr << "\t" << randHelix;

		randBeta = w_beta * sampEngTerms[3];

		//cerr << "\t" << randBeta;

		randCharge = w_charge * sampEngTerms[4];

		//cerr << "\t" << randCharge;

		randSolv = w_solv * sampEngTerms[5];

		//cerr << "\t" << randSolv;

		randSS = w_cys * sampEngTerms[6];

		//cerr << "\t" << randSS;

		//cerr << endl;

		randEngDist.push_back(randBend + randLj + randHelix + randBeta + randCharge + randSolv + randSS);
	}

	/* Calculate the gap between the distribution and
	 * the native value, as per Goldstein.
	 *
	 * OR SOME SORT OF OTHER Z-SCORE? HOW ABOUT
	 * THE HALF-IQR SCORE?
	 */
	natScore = bendPot + ljPot + helixPot + betaPot + chargePot + solvPot + ssbondPot;
	quarts = quartiles(randEngDist);
	gapScore = natScore - quarts[2];

	/* DEBUG
	 *
	cerr << "dG = " << natScore << " - " << quarts[2] << " = " << gapScore << endl;
	 */

	return gapScore;
}

// EXTRACTION OF INDIVIDUAL SAMPLED TERMS UNDER THE RANDOM
// ENERGY MODEL.
vector<float> GrahnenModel::sampleFoldGapTerms(string sequence, vector< pair<float,float> > bondAngles, vector<float> cBetaCAlphaDists, vector<float> cBetaDists, vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, vector<string> betaFourMers, vector<float> backboneTorsions, vector<float> exposures, CRandomMersenne& randGen){
	vector<float> terms;

	terms.push_back(sampleBend(sequence, bondAngles, randGen));
	terms.push_back(sampleLJ(sequence, cBetaCAlphaDists, cBetaDists, randGen));
	terms.push_back(sampleHelix(threeMers, helixFourMers, backboneDists, randGen));
	terms.push_back(sampleBeta(betaFourMers, backboneTorsions, randGen));
	terms.push_back(sampleCharge(sequence, cBetaDists, exposures, randGen));
	terms.push_back(sampleSolv(sequence, exposures, randGen));
	terms.push_back(sampleSSBond(sequence, cBetaDists, randGen));

	return terms;
}

/* Function to calculate the bending potential term in
 * Mukherjee's force field model, but with only the
 * bending of the sidechains included.
 */
float GrahnenModel::calcSidechainBendPotential(shared_ptr<Structure> protStruct){

	float K_theta = 10.0;
	float sidechainTerm1 = 0.0;
	float sidechainTerm2 = 0.0;

	int seqLength = protStruct->getNumResidues();

	int i, resNum;
	string resType;
	coord cAlpha1, cAlpha2, cAlpha3, cBeta;

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
	return 0.5*K_theta*(sidechainTerm1 + sidechainTerm2);
}

/* Function to calculate the sum of pairwise potentials
 * due to charges in this model. Uses the scalar form of 
 * Coloumb's Law with k_e = 1 and distances in Angstroms
 * (NOT ANYMORE, IT DOESN'T)
 *
 * Originally written by Jason Lai (jklai@asu.edu),
 * modified by Johan Grahnen.
 */
float GrahnenModel::calcChargePotential(shared_ptr<Structure> protStruct)
{
	/* Sum of potential due to pairwise charge interactions is calculated as
	 *
	 * V_charge = sum_i,j[k_e * (q_i * q_j / r_ij)]
	 *
	 * where 
	 *
	 * i,j = every unique pair of charged residues
	 * k_e = 1 (Coloumbs constant, provisionally set to 1 instead of 1/(4*pi*epsilon_0)
	 * q_i/q_j = charge of sidechain i/j
	 * r_ij = distance between centers of Cb beads for residues i and j
	 */
	double q_1,q_2,r_ij;
	string resType1,resType2;
	int resNum1,resNum2,i,j;

	double coulombicPot = 0;
	int seqLength = protStruct->getNumResidues();

	/* Adjustable parameter, but scaling is
	 * best done outside this function.
	 *
	 * Standard value in vacuum is 1/4*pi*epsilon_0 =~ 9*10^9.
	 *
	 * GETS MODIFIED DUE TO EXPOSURES.
	 */
	double k_e = 1;
	double chargeScaleFactor = 1000.0;

	/* Loop through all residues.
         */
	for(i = 0; i < seqLength; i++)
	{
		/* Convert residue type to single letter.
		 */
		resType1 = protStruct->getSingleResidue(i).getOneLetterType();

		/* Search matrix for electrostatic charge of residue i.
		 */
		resNum1 = forceFieldKey.find(resType1);
		q_1 = forceFieldMatrix[resNum1][6];

		/* Only calculate force if there is a charge on this
		 * residue.
		 */
		if(q_1 != 0){

			/* Grab CB bead coordinates of residue i.
			 */
			coord cResidue1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");

			/* Compare with all other residues after this in sequence (ensures
			 * all pairs are calculated only once), that are sufficiently
			 * distant in primary sequence.
			 */  
			for(j = i + SEQ_SEPARATION; j < seqLength; j++){

				/* Convert residue type to single letter
				 */
				resType2 = protStruct->getSingleResidue(j).getOneLetterType();

				/* Search matrix for electrostatic charge of residue j
				 */
				resNum2 = forceFieldKey.find(resType2);
				q_2 = forceFieldMatrix[resNum2][6];

				/* Only calculate potential if there is a charge on this
				 * residue.
				 */
				if(q_2 != 0){

					/* Grab coordinates of residue j
					 */
					coord cResidue2=protStruct->getSingleResidue(j).getAtomCoordsByType("CB");

					/* Calculate distance between residue i and j.
					 *
					 * Assume that we're dealing with point charges
					 * located at the center of the Cb bead.
					 */
					r_ij = NDDist(cResidue1,cResidue2);

					k_e = chargeScaleFactor*(1/dielectricConstant(protStruct->getSingleResidue(j).getApoExposure(), protStruct->getSingleResidue(i).getApoExposure()));

					/* Calculate the potential.
					 */
					coulombicPot += k_e*((q_1*q_2)/r_ij);

					/* DEBUG 
					 *
					cout << (i+1) << resType1 << "-" << (j+1) << resType2 << " (" << protStruct->getSingleResidue(j).getApoExposure() << "," << protStruct->getSingleResidue(i).getApoExposure() << ") -> k_e = " << k_e << ", V_C = " << coulombicPot << endl;
					 */
				}
			}
		}
	}

	/* Return the summed potential.
	 */
	return coulombicPot; 
}

/* Function to calculate a solvation potential (score) 
 * by approximating the exposed hydrophobic and buried polar
 * surface areas (by residue). Surface areas are approximated
 * by the NeighborVector algorithm described by Durham et al,
 * 2009, J Mol Model 15. We then weight the individual SASA's by 
 * hydrophobicity (same scale as the non-bonded potential in 
 * Model8) and polarity (based on hydrophobicity) such that
 * the potential becomes
 *
 * V_solv = sum_1_to_N[ h_i*SASA(i) + p_i*(1-SASA(i)) ]
 *
 * where
 * 
 * h_i = epsilon_ii from calcNonBondPotential() (proportional to the K-T hydrophobicity index)
 * p_i = (h_max + h_min) - h_i ("polarity index", the reverse of h_i)
 * SASA(i) = fractional exposure of residue i, calculated by Durham's 
 *           NeighborVector method
 *
 * NOTE: the use of weights assumes that the input
 * structure has been reduced to the two-bead representation
 * used by MukherjeeModel, RastogiModel and GrahnenModel.
 * Results may not be sensible for e.g. all-atom representations.
 */
float GrahnenModel::calcSolvationPotential(shared_ptr<Structure> protStruct, bool isAComplex){

	/* DEBUG
	 *
	cout << "Calculating a single solvation value..." << endl;
	 */

	/* Variables
	 */
	int i, resNum;
	float exposure;
	string resType1;
	Residue thisResidue;

	int protSize = protStruct->getNumResidues();
	float VSolv = 0.0;

	/* Re-cast the pointer.
	 *
	 * THIS IS LIKELY A LITTLE REDUNDANT BUT SHOULDN'T BE DANGEROUS...
	 */
	shared_ptr<TwoBeadStructure> protein = dynamic_pointer_cast<TwoBeadStructure>(protStruct);

	/* Go over all beads, weight the exposure by hydrophobicity and polarity.
	 */
	for(i = 0; i < protSize; i++){

		thisResidue = protein->getSingleResidue(i);
		resType1 = thisResidue.getOneLetterType();

		/* Search matrix for residue i.
		 */
		resNum = forceFieldKey.find(resType1);
		exposure = isAComplex ? thisResidue.getHoloExposure() : thisResidue.getApoExposure();

		/* Weight by hydrophobicity and polarity.
		 */
		VSolv += forceFieldMatrix[resNum][3]*exposure + forceFieldMatrix[resNum][7]*(1-exposure);

		/* DEBUG
		 *
		cout << "Contribution to VSolv: " << forceFieldMatrix[resNum][3]*exposure + forceFieldMatrix[resNum][7]*(1-exposure) << endl << endl;
		 */
	}

	return VSolv;
}

// COURTESY OF PRIYANKA NANDAKUMAR (pnandaku@andrew.cmu.edu), DOCUMENTATION
// FORTHCOMING.
//
// CURRENTLY JUST RETURNS THE NEGATIVE COUNT OF PREDICTED BONDS.
//
// DROPPED TORSION ANGLE THINGY FOR NOW.
float GrahnenModel::calcCysteinePotential(shared_ptr<Structure> protStruct){

	/* NEVERMIND THE TORSIONS...
	 *
	float angle1 = 44.0, angle2 = 46.0, angle3 = 62.0, angle4 = 64.0; // Thresholds for torsion angles
	 */
	float maxDistance = 4.5; // Threshold for bead distance
	int seqLength = protStruct->getNumResidues();

	float VCys = 0;
	int i, j;
	//float cbDistance, chiTB;
	float cbDistance;

	/* Search for first cysteine in sequence.
	 */
	for(i = 0; i < seqLength; i++){

		/* If a cysteine is found, get coordinates of alpha and beta beads.
		 */
		if(protStruct->getSingleResidue(i).getOneLetterType().compare("C") == 0){
			//coord cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
			coord cBeta1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");

			/* Search for a second cysteine in sequence after first cysteine. 
			 * This ensures all cysteine pairs are accounted for once.
			 */
			for(j = i + SEQ_SEPARATION; j < protStruct->getNumResidues(); j++){
				if(protStruct->getSingleResidue(j).getOneLetterType().compare("C")==0){

					/* If second cysteine is found, get coordinates of its alpha and beta beads.
					 */
					//coord cAlpha2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CA");
					coord cBeta2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CB");

					/* Distance between beta beads of first and second cysteines.
					 */
					cbDistance = NDDist(cBeta1, cBeta2);

					/* Torsion (dihedral) angle over the possible SS-bond.
					 */
					//chiTB = dihedralAngle(cAlpha1, cBeta1, cBeta2, cAlpha2);

					/* If distance is less than 4.5 A and angle is outside both ranges, 
					 * consider it a disulphide bond.
					 */
					//if(cbDistance < maxDistance && ((chiTB < angle1 || chiTB > angle2) && (chiTB < angle3 || chiTB > angle4))){
					if(cbDistance < maxDistance){
						VCys -= 1.0;

						/* DEBUG
						 *
						cout << "SS bond: " << (i+1) << "-" << (j+1) << endl;
						 */
					}
				}
			}
		}
	}

	return VCys;
}

/* Calculate the change in solvation potential (analogous to
 * the amount of buried surface area) on complex formation
 * in a protein-protein interaction, e.g.
 *
 * deltaSolv = V_solv(complex) - (V_solv(proteinOne) + V_solv(proteinTwo))
 *
 * where smaller/negative values are better.
 *
 * See calcSolvationPotential() for more details.
 */
float GrahnenModel::calcSolvationChange(shared_ptr<Structure> proteinOne, shared_ptr<Structure> proteinTwo){
	
	/* DEBUG
	 *
	cout << "Calculating solvation change..." << endl;
	 */

	/* Define some variables.
	 */
	float proteinOnePotential, proteinTwoPotential, complexPotential;
	vector<Residue> combinedResidues;
	vector<Residue> bindResidues;
	shared_ptr<TwoBeadStructure> newProtOne;
	shared_ptr<TwoBeadStructure> newProtTwo;
	vector< shared_ptr<TwoBeadStructure> > wholeComplex;

	/* Re-cast the pointers.
	 *
	 * PROBABLY REDUNDANT BUT SHOULDN'T BE DANGEROUS...
	 */
	shared_ptr<TwoBeadStructure> protOne = dynamic_pointer_cast<TwoBeadStructure>(proteinOne);
	shared_ptr<TwoBeadStructure> protTwo = dynamic_pointer_cast<TwoBeadStructure>(proteinTwo);

	/* Construct the complex from the sum of the parts.
	 */
	combinedResidues = protOne->getAllResidues();
	bindResidues = protTwo->getAllResidues();
	combinedResidues.insert(combinedResidues.end(), bindResidues.begin(), bindResidues.end());
	shared_ptr<TwoBeadStructure> combinedStruct(new TwoBeadStructure(combinedResidues));

	/* Calculate solvation potential of the complex.
	 *
	 * IF APO- AND HOLO-EXPOSURES ARE ALREADY CALCULATED, WHAT DOES THIS EVEN DO?
	 * DOES IT USE THE HOLO-EXPOSURES? IF NOT, WHY NOT? SHOULD WE HAVE A
	 * "IS IN COMPLEX" SWITCH HERE TO DIFFERENTIATE?
	 */
	complexPotential = calcSolvationPotential(combinedStruct, true);

	/* Calculate potential of the separate parts.
	 */
	proteinOnePotential = calcSolvationPotential(protOne, false);
	proteinTwoPotential = calcSolvationPotential(protTwo, false);

	/* DEBUG
	 *
	cout << "Protein surface area: " << proteinArea << endl;
	cout << "Ligand surface area: " << ligandArea << endl;
	cout << "Complex surface area: " << complexArea << endl;
	 */

	/* And finally calculate the difference: smaller is better.
	 */
	return complexPotential - (proteinOnePotential + proteinTwoPotential);
}

/* Function to calculate an exposure-dependent dielectric constant
 * for two interacting residues. Considers two classes of
 * interaction:
 *
 * 1) One or both exposed: infinite dielectric (charged interactions
 *    completely screened)
 * 2) Both buried: low dielectric (3)
 *
 * "Exposed" is defined as a fractional exposure of more 0.25.
 * See Finkelstein, "Protein physics", ch 6, (YEAR) for 
 * more depth on this phenomenon, and (REF TO ME) for a rationale of 
 * the categorization.
 *
 * NOTE NOTE NOTE: You'll get some pretty funny return
 * values unless your compiler supports IEEE 754 floating 
 * point numbers...
 */
float GrahnenModel::dielectricConstant(float exposureOne, float exposureTwo){
	float lowEpsilon = 3;
	float infEpsilon = INFINITY; // Hopefully your compiler supports this
	float threshold = 0.25;
	float dielectric;

	/* Decide on the exposure situation.
	 */
	if(exposureOne > threshold || exposureTwo > threshold){
		dielectric = infEpsilon;
	}
	else{
		dielectric = lowEpsilon;
	}

	return dielectric;
}

float GrahnenModel::sampleBend(string sequence, vector< pair<float,float> > bondAngles, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(sequence.length() != bondAngles.size()){
		string message = "ERROR: number of residues not equivalent to number of supplied bond angle pairs!";
		throw runtime_error(message);
	}

	unsigned int i, sampNum;
	int resNum;
	string resType;
	float K_theta = 10.0;
	float sample = 0.0;
	
	myRng random(randGen);
	random_shuffle(bondAngles.begin(), bondAngles.end(), random);

	for(i = 0; i < sequence.length(); i++){
		resType = sequence.substr(i,1);
		resNum = forceFieldKey.find(resType);

		/* No defined angles for glycines.
		 */
		if(resType.compare("G") == 0){
			/* DEBUG
			 *
			cerr << "Residue " << (i+1) << " is a glycine, has no bond angle." << endl;
			 */

			continue;
		}

		/* Make sure we're sampling an angle that exists: glycines
		 * in the structure produce NaN angles.
		 */
		sampNum = i;
		while(isnan(bondAngles[sampNum].first) || isnan(bondAngles[sampNum].second)){
			/* DEBUG
			 *
			cerr << "Angle #" << sampNum << " doesn't exist, picking another one..." << endl;
			 */

			//sampNum = rand() % bondAngles.size();
			sampNum++;
			sampNum = sampNum % bondAngles.size();
		}

		/* Most residues have two defined angles.
		 */
		if(i != 0 && i != sequence.length()-1){
			sample += pow((bondAngles[sampNum].first - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
			sample += pow((bondAngles[sampNum].second - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
		/* N-terminal residues lack the first angle.
		 */
		else if(i == 0){
			sample += pow((bondAngles[sampNum].second - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
		/* C-terminal residues lack the second angle.
		 */
		else if(i == sequence.length()-1){
			sample += pow((bondAngles[sampNum].first - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
	}
	sample *= 0.5*K_theta;

	return sample;
}

float GrahnenModel::sampleLJ(string sequence, vector<float> cBetaCAlphaDists, vector<float> cBetaDists, CRandomMersenne& randGen){
	/* Validate input.
	 *
	 * In an NxN matrix of pairwise distances, excluding those separated by
	 * less than a distance s in primary sequence, we should find
	 *
	 * (N-s)*((N-s) - 1) / 2 + (N-s)
	 *
	 * such distances measured in total.
	 */
	unsigned int expectedCount = (sequence.length()-SEQ_SEPARATION)*( (sequence.length()-SEQ_SEPARATION) - 1 ) / 2 + (sequence.length() - SEQ_SEPARATION);
	unsigned int i, j;
	int sampNum, resNum1, resNum2;
	float epsilon_ij, sigma_ij;

	unsigned int counter = 0;
	float epsilon_ca = 0.05; // C-alpha bead interaction potential
	float r_gyr_ca = 1.8; // C-alpha bead radius
	float sample = 0.0;
	
	myRng random(randGen);
	random_shuffle(cBetaCAlphaDists.begin(), cBetaCAlphaDists.end(), random);
	random_shuffle(cBetaDists.begin(), cBetaDists.end(), random);

	/* Ca-Cb bead interactions.
	 */
	for(i = 0; i < sequence.length(); i++){
		/* Count all Ca-Cb interactions that are separated by at least
		 * a certain number of residues.
		 *
		 * 2013-05-21: Fixed a bug where all residues were not being
		 * examined (Nadia Bykova).
		 */
                for(j = 0; j < sequence.length(); j++){

                        if ((i >= j && (i - j) < SEQ_SEPARATION) || (i < j && (j - i) < SEQ_SEPARATION)){
                                continue;
                        }
			/* End of edit
			 */
 
			/* Don't double-count Gly interactions.
			 */
			if(sequence[j] == 'G'){
				counter++;
				continue;
			}
			else{
				resNum2 = forceFieldKey.find(sequence[j]);
				epsilon_ij = sqrt(epsilon_ca * forceFieldMatrix[resNum2][3]);
				sigma_ij = r_gyr_ca + forceFieldMatrix[resNum2][1];

				/* Gly residues in original structure have distances
				 * of "NaN", but we need a real one here.
				 */
				sampNum = counter;
				while(isnan(cBetaCAlphaDists[sampNum])){
					/* DEBUG
					 *
					cerr << "Distance #" << sampNum << " doesn't exist, picking another one..." << endl;
					 */

					//sampNum = rand() % cBetaCAlphaDists.size();
					sampNum++;
					sampNum = sampNum % cBetaCAlphaDists.size();
				}

				sample += epsilon_ij*( pow( sigma_ij/cBetaCAlphaDists[sampNum] ,12) - pow( sigma_ij/cBetaCAlphaDists[sampNum] ,6) );

				counter++;
			}
		}
	}

	/* DEBUG
	 *
	cerr << "Due to Ca-Cb interactions: " << sample << endl;
	 */

	/* Cb-Cb bead interactions.
	 */
	counter = 0;
	for(i = 0; i < sequence.length(); i++){
		resNum1 = forceFieldKey.find(sequence[i]);

		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){
			resNum2 = forceFieldKey.find(sequence[j]);

			sigma_ij = forceFieldMatrix[resNum1][1] + forceFieldMatrix[resNum2][1];
			epsilon_ij = sqrt(forceFieldMatrix[resNum1][3] * forceFieldMatrix[resNum2][3]);

			sample += epsilon_ij*( pow( sigma_ij/cBetaDists[counter] ,12) - pow( sigma_ij/cBetaDists[counter] ,6) );

			counter++;
		}
	}

	/* DEBUG
	 *
	cerr << "Due to all interactions: " << sample << endl;
	 */

	sample *= 4;

	return sample;
}

float GrahnenModel::sampleHelix(vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(threeMers.size() != helixFourMers.size() || threeMers.size() != backboneDists.size()){
		string message = "ERROR: number of C-alpha distances does not match number of k-mers!";
		throw runtime_error(message);
	}

	unsigned int i;
	int resNum1, resNum2, resNum3, resNum4;
	float K_1_3, K_1_4;
	float r_h = 5.5;
	float sample = 0.0;
	
	myRng random(randGen);
	random_shuffle(backboneDists.begin(), backboneDists.end(), random);

///Added to reflect the distribution of the helices from the PDB
///	srand(time(NULL));
	vector<float> positions = sampleHelixPosition(threeMers.size());



	for(i = 0; i < threeMers.size(); i++){
		resNum1 = forceFieldKey.find(threeMers[i][0]);
		resNum2 = forceFieldKey.find(threeMers[i][1]);
		resNum3 = forceFieldKey.find(threeMers[i][2]);
		
		K_1_3 = (1.0/3.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4]);

///Added to reflect the distribution of the helices from the PDB
			sample += K_1_3*positions[i];

//orig		sample += K_1_3*pow(backboneDists[i].first - r_h ,2);

		resNum1 = forceFieldKey.find(helixFourMers[i][0]);
		resNum2 = forceFieldKey.find(helixFourMers[i][1]);
		resNum3 = forceFieldKey.find(helixFourMers[i][2]);
		resNum4 = forceFieldKey.find(helixFourMers[i][3]);
			
		K_1_4 = (1.0/4.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4] + forceFieldMatrix[resNum4][4]);

///Added to reflect the distribution of the helices from the PDB
		sample += K_1_4*positions[i];

//orig		sample += K_1_4*pow(backboneDists[i].second - r_h ,2);
	}
	sample *= 0.5;

	return sample;
}

///Added to get random positions of helices
///Feb 11th, 2014 by Dohyup Kim
///Return float vectors for the calculation
vector<float> GrahnenModel::sampleHelixPosition(int protSize){
	/*Distribution of helix sizes from PDB
	 */
//Distribution of helix Sizes of 2424 protein structures from the PDB
	int helixSizes[110] = 
	{0, 1, 6, 401, 227, 631, 534, 332, 337, 381, 406, 415, 393, 437, 336, 289, 252, 223, 151, 147, 110,
	 73, 73, 56, 29, 20, 23, 26, 17, 21, 20, 15, 14, 14, 7, 3, 5, 7, 7, 2, 1,
	 1, 1, 0, 0, 2, 2, 3, 1, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
	 1, 0, 0, 0, 0, 0, 0, 0, 1};

//Distribution of helix numbers of 2424 proteins structures from the PDB
	int helixNumbers[11] = 
	{549, 285, 397, 354, 322, 250, 144, 74, 35, 13, 1};

//Seed the random number generator
///	srand(time(NULL));

//Get the random number between 1 to 2424 to determine the number of helices in a protein
	int randNum = rand() % 2424 + 1;
	int helixNum = 0;

//The size of the native protein.
	int proteinSize = protSize;
	
//For loop to select the number of helices in a protein.
//The distribution of helix numbers are added to position the random number with in the range.
	int count = 0;
	for(int i = 0; i < 11; i++){
		count += helixNumbers[i];
		if (randNum <= count)
		{
			helixNum = i;
			cout << i << endl;
			break;
		}
	}
	cout << randNum << endl;

//vectors to prepare the output.
	vector<int> sizes;
	vector<float> positions;
	float positionArray[proteinSize];	
//Initialize the helical positions.
//0 indicates no helix, 1 indicates helix, and numbers more than 1 indicates overlapping helices, which will be recalculated.
	for(int i = 0; i < proteinSize; i++)
	{
		positionArray[i] = 0;
	}
	for(int i = 0; i < proteinSize; i++)
	{
		positions.push_back(0);
	}	

	while(1)
	{

//for loop to get the sizes of helices according to the number of helices generated above.
//Same weighted sampling is done.
		int sizeSum = 0;
		for(int i = 0; i < helixNum; i++)
		{
			randNum = rand() % 6467;
			count = 0;
			for(int j = 0; j < 110; j++)
			{
				count += helixSizes[j];
				if (randNum <= count)
				{
					cout << j << endl;
					sizes.push_back(j);
					sizeSum += j;
					break;
				}
			}
			cout << randNum << endl;
		}

//We cannot have the sum of helix sizes over the length of the protein.
		if(sizeSum < proteinSize)
		{
			break;
		}
	}

	bool stopCase = true;
	while(stopCase)
	{
		for(int i = 0; i < helixNum; i++)
		{

//Set the position for a helix onto the native protein and mark it with 1.
//If the helices are overlapping, the number will be bigger than 1 so that it will regenerate the position.
			while(1)
			{
				randNum = rand()%proteinSize;
				if(randNum+sizes[i] < proteinSize)
				{
					for(int j = randNum; j < (randNum+sizes[i]); j++)
					{
						positionArray[j] += 1.0;
					}
					break;
				}
			}
		}
		int count = 0;
		for(int k = 0; k < proteinSize; k++)
		{
			count = k;
			if(positionArray[k] > 1)
			{
				for(int l = 0; l < proteinSize; l++)
				{
					positionArray[l] = 0;
				}
				break;
			}
		}

		if(count >= proteinSize-1)
		{
			stopCase = false;
		}
	}

//Dump array into a vector so that the return output can be set up.	
	for(int i = 0; i < proteinSize; i++)
	{
		positions[i] = positionArray[i];

	}
	return positions;
}


float GrahnenModel::sampleBeta(vector<string> betaFourMers, vector<float> backboneTorsions, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(betaFourMers.size() != backboneTorsions.size()){
		string message = "ERROR: number of k-mers do not match number of torsions!";
		throw runtime_error(message);
	}
	unsigned int i;
	int resNum1, resNum2, resNum3, resNum4;
	float K_1_4;
	float phi_b = -150; // Equilibrium Ca-Ca torsion angle (Ca_i-1->Ca_i->Ca_i+1->Ca+i+2 torsion) for beta sheets. Equivalent to 210 positive degrees.
	float C_s = 0.01; // Scaling constant to make similar to calcHelixPotential()
	float sample = 0.0;
	
	myRng random(randGen);
	random_shuffle(backboneTorsions.begin(), backboneTorsions.end(), random);

	for(i = 0; i < betaFourMers.size(); i++){
		resNum1 = forceFieldKey.find(betaFourMers[i][0]);
		resNum2 = forceFieldKey.find(betaFourMers[i][1]);
		resNum3 = forceFieldKey.find(betaFourMers[i][2]);
		resNum4 = forceFieldKey.find(betaFourMers[i][3]);

		K_1_4 = 0.25*( forceFieldMatrix[resNum1][5] + forceFieldMatrix[resNum2][5] + forceFieldMatrix[resNum3][5] + forceFieldMatrix[resNum4][5]);

		sample += K_1_4*pow( C_s*(backboneTorsions[i] - phi_b), 2);
	}

	return sample;
}

float GrahnenModel::sampleCharge(string sequence, vector<float> cBetaDists, vector<float> exposures, CRandomMersenne& randGen){
	/* Validate input.
	 *
	 * In an NxN matrix of pairwise distances, excluding those separated by
	 * less than a distance s in primary sequence, we should find
	 *
	 * (N-s)*((N-s) - 1) / 2 + (N-s)
	 *
	 * such distances measured in total.
	 */
	unsigned int expectedCount = (sequence.length()-SEQ_SEPARATION)*( (sequence.length()-SEQ_SEPARATION) - 1 ) / 2 + (sequence.length() - SEQ_SEPARATION);
	if(expectedCount != cBetaDists.size() || sequence.length() != exposures.size()){
		string message = "ERROR: number of bead distances or exposures do not match sequence length!";
		throw runtime_error(message);
	}

	unsigned int i, j;
	int resNum1, resNum2;
	float q_1, q_2, k_e;
	
	double chargeScaleFactor = 1000.0; // Scales charge term to approx same magnitude as LJ
	float sample = 0.0;
	unsigned int counter = 0;
	
	myRng random(randGen);
	random_shuffle(cBetaDists.begin(), cBetaDists.end(), random);
	random_shuffle(exposures.begin(), exposures.end(), random);
	
	for(i = 0; i < sequence.length(); i++){
		resNum1 = forceFieldKey.find(sequence[i]);
		q_1 = forceFieldMatrix[resNum1][6];
		
		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){
			resNum2 = forceFieldKey.find(sequence[j]);
			q_2 = forceFieldMatrix[resNum2][6];
			
			if(q_1 != 0 && q_2 != 0){
				k_e = chargeScaleFactor*(1/dielectricConstant(exposures[i], exposures[j]));
				sample += k_e*((q_1*q_2)/cBetaDists[counter]);
			}

			counter++;
		}
	}
	
	return sample;
}

float GrahnenModel::sampleSolv(string sequence, vector<float> exposures, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(sequence.length() != exposures.size()){
		string message = "ERROR: number of exposure values does not match sequence length!";
		throw runtime_error(message);
	}

	unsigned int i;
	int resNum;
	float sample = 0.0;
	
	myRng random(randGen);
	random_shuffle(exposures.begin(), exposures.end(), random);

	for(i = 0; i < sequence.length(); i++){
		resNum = forceFieldKey.find(sequence[i]);
		sample += forceFieldMatrix[resNum][3]*exposures[i] + forceFieldMatrix[resNum][7]*(1-exposures[i]);
	}

	return sample;
}

float GrahnenModel::sampleSSBond(string sequence, vector<float> cBetaDists, CRandomMersenne& randGen){
	/* Validate input.
	 *
	 * In an NxN matrix of pairwise distances, excluding those separated by
	 * less than a distance s in primary sequence, we should find
	 *
	 * (N-s)*((N-s) - 1) / 2 + (N-s)
	 *
	 * such distances measured in total.
	 */
	unsigned int expectedCount = (sequence.length()-SEQ_SEPARATION)*( (sequence.length()-SEQ_SEPARATION) - 1 ) / 2 + (sequence.length() - SEQ_SEPARATION);
	if(expectedCount != cBetaDists.size()){
		string message = "ERROR: number of bead distances does not match sequence length!";
		throw runtime_error(message);
	}

	unsigned int i, j;

	float maxSSDistance = 4.5; // Threshold for disulfide bond distance
	float sample = 0.0;
	unsigned int counter = 0;
		
	myRng random(randGen);
	random_shuffle(cBetaDists.begin(), cBetaDists.end(), random);
	
	for(i = 0; i < sequence.length(); i++){
		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){
			if(sequence[i] == 'C' && sequence[j] == 'C' && cBetaDists[counter] < maxSSDistance){
					sample -= 1.0;
			}
			counter++;
		}
	}
	
	return sample;
}


vector<vector<float>> GrahnenModel::getFoldScoreTerms_pos(shared_ptr<Structure> protein){
	vector<vector<float>> values;
	int i;

	/* Re-cast the pointer.
	 */
	shared_ptr<TwoBeadStructure> prot = dynamic_pointer_cast<TwoBeadStructure>(protein);

	/* Check that exposure values exist, otherwise throw an exception.
	 */
	int protSize = prot->getNumResidues();
	vector<float> noCalcs(protSize, 0.0);
	vector<float> exposures;
	for(i = 0; i < protSize; i++){
		exposures.push_back(prot->getSingleResidue(i).getApoExposure());
	}
	if(noCalcs == exposures){
		throw runtime_error("No apo-exposure values for residues! Pre-calculate with TwoBeadStructure::calcExposure().");
	}

	/* Start the clock.
	 *
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);
	double timeStop;
	bool clockRunning = true;
	*/

	/* Calculate the bending potential.
	 */
	values.push_back(calcSidechainBendPotential_pos(prot));

	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_LJ == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(bend)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the non-bonding (Lennard-Jones) potential.
	 */
	values.push_back(calcNonbondPotential_pos(prot));

	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_helix == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(LJ)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the helix potential.
	 */
	values.push_back(calcHelixPotential_pos(prot));
//	values.push_back(calcHelixPotential_pos2(prot));

	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_beta == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(helix)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the beta-sheet potential.
	 */
	values.push_back(calcBetaPotential_pos(prot));
//	values.push_back(calcBetaPotential_pos2(prot));

	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_charge == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(beta)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate the electrostatic potential
	 */
	values.push_back(calcChargePotential_pos(prot));

	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_solv == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(charge)" << endl;
		clockRunning = false;
	}
	*/

	/* Calculate solvation potential.
	 */
	values.push_back(calcSolvationPotential_pos(prot, false));

	/* Stop the clock and print if
	 * necessary.
	 *
	if(clockRunning && w_cys == 0.0){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(solv)" << endl;
		clockRunning = false;
	}
	*/

	// AND THE CYSTEINE ONE AS THAT IS DEVELOPED.
	values.push_back(calcCysteinePotential_pos(prot));

	/* Alway stop the clock here.
	 *
	if(clockRunning){
		gettimeofday(&tim, NULL);
		timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
		cout << (timeStop - timeStart) << "	(cys)" << endl;
		clockRunning = false;
	}
	*/

	return values;
}

vector<vector<float>> GrahnenModel::sampleFoldGapTerms_pos(string sequence, vector< pair<float,float> > bondAngles, vector<float> cBetaCAlphaDists, vector<float> cBetaDists, vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, vector<string> betaFourMers, vector<float> backboneTorsions, vector<float> exposures, CRandomMersenne& randGen){
	vector<vector<float>> terms;

	terms.push_back(sampleBend_pos(sequence, bondAngles, randGen));
	terms.push_back(sampleLJ_pos(sequence, cBetaCAlphaDists, cBetaDists, randGen));
	terms.push_back(sampleHelix_pos(threeMers, helixFourMers, backboneDists, randGen));
//	terms.push_back(sampleHelix_pos2(threeMers, helixFourMers, backboneDists, randGen));
	terms.push_back(sampleBeta_pos(betaFourMers, backboneTorsions, randGen));
//	terms.push_back(sampleBeta_pos2(betaFourMers, backboneTorsions, randGen));
	terms.push_back(sampleCharge_pos(sequence, cBetaDists, exposures, randGen));
	terms.push_back(sampleSolv_pos(sequence, exposures, randGen));
	terms.push_back(sampleSSBond_pos(sequence, cBetaDists, randGen));

	return terms;
}

vector<float> GrahnenModel::calcSidechainBendPotential_pos(shared_ptr<Structure> protStruct){

	float K_theta = 10.0;
	float sidechainTerm1 = 0.0;
	float sidechainTerm2 = 0.0;

	int seqLength = protStruct->getNumResidues();

	vector<float>res(seqLength,0.0);

	int i, resNum;
	string resType;
	coord cAlpha1, cAlpha2, cAlpha3, cBeta;

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

			res[i] += 0.5*K_theta*pow( (NDAngle(cAlpha1, cAlpha2, cBeta) - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
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

			res[i] += 0.5*K_theta*pow( (NDAngle(cBeta, cAlpha1, cAlpha2) - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
	}

	/* DEBUG
	 *
	cout << "Backbone bend: " << backboneTerm << ", sidechain bend 1: " << sidechainTerm1 << ", sidechain bend 2: " << sidechainTerm2 << endl;
	 */

	/* Sum up terms, convert to a potential and return.
	 */
	return res;
}
vector<float> GrahnenModel::calcChargePotential_pos(shared_ptr<Structure> protStruct)
{
	/* Sum of potential due to pairwise charge interactions is calculated as
	 *
	 * V_charge = sum_i,j[k_e * (q_i * q_j / r_ij)]
	 *
	 * where
	 *
	 * i,j = every unique pair of charged residues
	 * k_e = 1 (Coloumbs constant, provisionally set to 1 instead of 1/(4*pi*epsilon_0)
	 * q_i/q_j = charge of sidechain i/j
	 * r_ij = distance between centers of Cb beads for residues i and j
	 */
	double q_1,q_2,r_ij;
	string resType1,resType2;
	int resNum1,resNum2,i,j;

	double coulombicPot = 0;
	int seqLength = protStruct->getNumResidues();

	vector<float> res(seqLength,0.0);
	/* Adjustable parameter, but scaling is
	 * best done outside this function.
	 *
	 * Standard value in vacuum is 1/4*pi*epsilon_0 =~ 9*10^9.
	 *
	 * GETS MODIFIED DUE TO EXPOSURES.
	 */
	double k_e = 1;
	double chargeScaleFactor = 1000.0;

	/* Loop through all residues.
         */
	for(i = 0; i < seqLength; i++)
	{
		/* Convert residue type to single letter.
		 */
		resType1 = protStruct->getSingleResidue(i).getOneLetterType();

		/* Search matrix for electrostatic charge of residue i.
		 */
		resNum1 = forceFieldKey.find(resType1);
		q_1 = forceFieldMatrix[resNum1][6];

		/* Only calculate force if there is a charge on this
		 * residue.
		 */
		if(q_1 != 0){

			/* Grab CB bead coordinates of residue i.
			 */
			coord cResidue1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");

			/* Compare with all other residues after this in sequence (ensures
			 * all pairs are calculated only once), that are sufficiently
			 * distant in primary sequence.
			 */
			for(j = i + SEQ_SEPARATION; j < seqLength; j++){

				/* Convert residue type to single letter
				 */
				resType2 = protStruct->getSingleResidue(j).getOneLetterType();

				/* Search matrix for electrostatic charge of residue j
				 */
				resNum2 = forceFieldKey.find(resType2);
				q_2 = forceFieldMatrix[resNum2][6];

				/* Only calculate potential if there is a charge on this
				 * residue.
				 */
				if(q_2 != 0){

					/* Grab coordinates of residue j
					 */
					coord cResidue2=protStruct->getSingleResidue(j).getAtomCoordsByType("CB");

					/* Calculate distance between residue i and j.
					 *
					 * Assume that we're dealing with point charges
					 * located at the center of the Cb bead.
					 */
					r_ij = NDDist(cResidue1,cResidue2);

					k_e = chargeScaleFactor*(1/dielectricConstant(protStruct->getSingleResidue(j).getApoExposure(), protStruct->getSingleResidue(i).getApoExposure()));

					/* Calculate the potential.
					 */
					res[i] += 0.5*k_e*((q_1*q_2)/r_ij);
					res[j] += 0.5*k_e*((q_1*q_2)/r_ij);

					/* DEBUG
					 *
					cout << (i+1) << resType1 << "-" << (j+1) << resType2 << " (" << protStruct->getSingleResidue(j).getApoExposure() << "," << protStruct->getSingleResidue(i).getApoExposure() << ") -> k_e = " << k_e << ", V_C = " << coulombicPot << endl;
					 */
				}
			}
		}
	}

	/* Return the summed potential.
	 */
	return res;
}
vector<float> GrahnenModel::calcSolvationPotential_pos(shared_ptr<Structure> protStruct, bool isAComplex){

	/* DEBUG
	 *
	cout << "Calculating a single solvation value..." << endl;
	 */

	/* Variables
	 */
	int i, resNum;
	float exposure;
	string resType1;
	Residue thisResidue;

	int protSize = protStruct->getNumResidues();
	vector<float>res(protSize,0.0);

	float VSolv = 0.0;

	/* Re-cast the pointer.
	 *
	 * THIS IS LIKELY A LITTLE REDUNDANT BUT SHOULDN'T BE DANGEROUS...
	 */
	shared_ptr<TwoBeadStructure> protein = dynamic_pointer_cast<TwoBeadStructure>(protStruct);

	/* Go over all beads, weight the exposure by hydrophobicity and polarity.
	 */
	for(i = 0; i < protSize; i++){

		thisResidue = protein->getSingleResidue(i);
		resType1 = thisResidue.getOneLetterType();

		/* Search matrix for residue i.
		 */
		resNum = forceFieldKey.find(resType1);
		exposure = isAComplex ? thisResidue.getHoloExposure() : thisResidue.getApoExposure();

		/* Weight by hydrophobicity and polarity.
		 */
		res[i] += forceFieldMatrix[resNum][3]*exposure + forceFieldMatrix[resNum][7]*(1-exposure);

		/* DEBUG
		 *
		cout << "Contribution to VSolv: " << forceFieldMatrix[resNum][3]*exposure + forceFieldMatrix[resNum][7]*(1-exposure) << endl << endl;
		 */
	}

	return res;
}
vector<float> GrahnenModel::calcCysteinePotential_pos(shared_ptr<Structure> protStruct){

	/* NEVERMIND THE TORSIONS...
	 *
	float angle1 = 44.0, angle2 = 46.0, angle3 = 62.0, angle4 = 64.0; // Thresholds for torsion angles
	 */
	float maxDistance = 4.5; // Threshold for bead distance
	int seqLength = protStruct->getNumResidues();

	vector<float>res(seqLength,0.0);

	float VCys = 0;
	int i, j;
	//float cbDistance, chiTB;
	float cbDistance;

	/* Search for first cysteine in sequence.
	 */
	for(i = 0; i < seqLength; i++){

		/* If a cysteine is found, get coordinates of alpha and beta beads.
		 */
		if(protStruct->getSingleResidue(i).getOneLetterType().compare("C") == 0){
			//coord cAlpha1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CA");
			coord cBeta1 = protStruct->getSingleResidue(i).getAtomCoordsByType("CB");

			/* Search for a second cysteine in sequence after first cysteine.
			 * This ensures all cysteine pairs are accounted for once.
			 */
			for(j = i + SEQ_SEPARATION; j < protStruct->getNumResidues(); j++){
				if(protStruct->getSingleResidue(j).getOneLetterType().compare("C")==0){

					/* If second cysteine is found, get coordinates of its alpha and beta beads.
					 */
					//coord cAlpha2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CA");
					coord cBeta2 = protStruct->getSingleResidue(j).getAtomCoordsByType("CB");

					/* Distance between beta beads of first and second cysteines.
					 */
					cbDistance = NDDist(cBeta1, cBeta2);

					/* Torsion (dihedral) angle over the possible SS-bond.
					 */
					//chiTB = dihedralAngle(cAlpha1, cBeta1, cBeta2, cAlpha2);

					/* If distance is less than 4.5 A and angle is outside both ranges,
					 * consider it a disulphide bond.
					 */
					//if(cbDistance < maxDistance && ((chiTB < angle1 || chiTB > angle2) && (chiTB < angle3 || chiTB > angle4))){
					if(cbDistance < maxDistance){
						res[i] -= 0.5;
						res[j] -= 0.5;

						/* DEBUG
						 *
						cout << "SS bond: " << (i+1) << "-" << (j+1) << endl;
						 */
					}
				}
			}
		}
	}

	return res;
}
vector<float> GrahnenModel::sampleBend_pos(string sequence, vector< pair<float,float> > bondAngles, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(sequence.length() != bondAngles.size()){
		string message = "ERROR: number of residues not equivalent to number of supplied bond angle pairs!";
		throw runtime_error(message);
	}

	unsigned int i, sampNum;
	int resNum;
	string resType;
	float K_theta = 10.0;
	float sample = 0.0;

	vector<float>res(sequence.length(),0.0);

	myRng random(randGen);
	random_shuffle(bondAngles.begin(), bondAngles.end(), random);

	for(i = 0; i < sequence.length(); i++){
		resType = sequence.substr(i,1);
		resNum = forceFieldKey.find(resType);

		/* No defined angles for glycines.
		 */
		if(resType.compare("G") == 0){
			/* DEBUG
			 *
			cerr << "Residue " << (i+1) << " is a glycine, has no bond angle." << endl;
			 */

			continue;
		}

		/* Make sure we're sampling an angle that exists: glycines
		 * in the structure produce NaN angles.
		 */
		sampNum = i;
		while(isnan(bondAngles[sampNum].first) || isnan(bondAngles[sampNum].second)){
			/* DEBUG
			 *
			cerr << "Angle #" << sampNum << " doesn't exist, picking another one..." << endl;
			 */

			//sampNum = rand() % bondAngles.size();
			sampNum++;
			sampNum = sampNum % bondAngles.size();
		}

		/* Most residues have two defined angles.
		 */
		if(i != 0 && i != sequence.length()-1){
			res[i] += 0.5*K_theta*pow((bondAngles[sampNum].first - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
			res[i] += 0.5*K_theta*pow((bondAngles[sampNum].second - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
		/* N-terminal residues lack the first angle.
		 */
		else if(i == 0){
			res[i] += 0.5*K_theta*pow((bondAngles[sampNum].second - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
		/* C-terminal residues lack the second angle.
		 */
		else if(i == sequence.length()-1){
			res[i] += 0.5*K_theta*pow((bondAngles[sampNum].first - forceFieldMatrix[resNum][2]) * (PI/180.0), 2);
		}
	}

	return res;
}

vector<float> GrahnenModel::sampleLJ_pos(string sequence, vector<float> cBetaCAlphaDists, vector<float> cBetaDists, CRandomMersenne& randGen){
	/* Validate input.
	 *
	 * In an NxN matrix of pairwise distances, excluding those separated by
	 * less than a distance s in primary sequence, we should find
	 *
	 * (N-s)*((N-s) - 1) / 2 + (N-s)
	 *
	 * such distances measured in total.
	 */
	unsigned int expectedCount = (sequence.length()-SEQ_SEPARATION)*( (sequence.length()-SEQ_SEPARATION) - 1 ) / 2 + (sequence.length() - SEQ_SEPARATION);
	unsigned int i, j;
	int sampNum, resNum1, resNum2;
	float epsilon_ij, sigma_ij;

	unsigned int counter = 0;
	float epsilon_ca = 0.05; // C-alpha bead interaction potential
	float r_gyr_ca = 1.8; // C-alpha bead radius
	float sample = 0.0;

	myRng random(randGen);
	random_shuffle(cBetaCAlphaDists.begin(), cBetaCAlphaDists.end(), random);
	random_shuffle(cBetaDists.begin(), cBetaDists.end(), random);

	vector<float>res(sequence.length(),0.0);

	/* Ca-Cb bead interactions.
	 */
	for(i = 0; i < sequence.length(); i++){
		/* Count all Ca-Cb interactions that are separated by at least
		 * a certain number of residues.
		 *
		 * 2013-05-21: Fixed a bug where all residues were not being
		 * examined (Nadia Bykova).
		 */
                for(j = 0; j < sequence.length(); j++){

                        if ((i >= j && (i - j) < SEQ_SEPARATION) || (i < j && (j - i) < SEQ_SEPARATION)){
                                continue;
                        }
			/* End of edit
			 */

			/* Don't double-count Gly interactions.
			 */
			if(sequence[j] == 'G'){
				counter++;
				continue;
			}
			else{
				resNum2 = forceFieldKey.find(sequence[j]);
				epsilon_ij = sqrt(epsilon_ca * forceFieldMatrix[resNum2][3]);
				sigma_ij = r_gyr_ca + forceFieldMatrix[resNum2][1];

				/* Gly residues in original structure have distances
				 * of "NaN", but we need a real one here.
				 */
				sampNum = counter;
				while(isnan(cBetaCAlphaDists[sampNum])){
					/* DEBUG
					 *
					cerr << "Distance #" << sampNum << " doesn't exist, picking another one..." << endl;
					 */

					//sampNum = rand() % cBetaCAlphaDists.size();
					sampNum++;
					sampNum = sampNum % cBetaCAlphaDists.size();
				}

				res[i] += 4*epsilon_ij*( pow( sigma_ij/cBetaCAlphaDists[sampNum] ,12) - pow( sigma_ij/cBetaCAlphaDists[sampNum] ,6) );

				counter++;
			}
		}
	}

	/* DEBUG
	 *
	cerr << "Due to Ca-Cb interactions: " << sample << endl;
	 */

	/* Cb-Cb bead interactions.
	 */
	counter = 0;
	for(i = 0; i < sequence.length(); i++){
		resNum1 = forceFieldKey.find(sequence[i]);

		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){
			resNum2 = forceFieldKey.find(sequence[j]);

			sigma_ij = forceFieldMatrix[resNum1][1] + forceFieldMatrix[resNum2][1];
			epsilon_ij = sqrt(forceFieldMatrix[resNum1][3] * forceFieldMatrix[resNum2][3]);

			res[i] += 2*epsilon_ij*( pow( sigma_ij/cBetaDists[counter] ,12) - pow( sigma_ij/cBetaDists[counter] ,6) );
			res[j] += 2*epsilon_ij*( pow( sigma_ij/cBetaDists[counter] ,12) - pow( sigma_ij/cBetaDists[counter] ,6) );

			counter++;
		}
	}

	/* DEBUG
	 *
	cerr << "Due to all interactions: " << sample << endl;
	 */

	return res;
}
vector<float> GrahnenModel::sampleHelix_pos(vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(threeMers.size() != helixFourMers.size() || threeMers.size() != backboneDists.size()){
		string message = "ERROR: number of C-alpha distances does not match number of k-mers!";
		throw runtime_error(message);
	}

	unsigned int i;
	int resNum1, resNum2, resNum3, resNum4;
	float K_1_3, K_1_4;
	float r_h = 5.5;
	float sample = 0.0;

	vector<float>res(threeMers.size()+5,0.0);

	myRng random(randGen);
	random_shuffle(backboneDists.begin(), backboneDists.end(), random);

	for(i = 0; i < threeMers.size(); i++){
		resNum1 = forceFieldKey.find(threeMers[i][0]);
		resNum2 = forceFieldKey.find(threeMers[i][1]);
		resNum3 = forceFieldKey.find(threeMers[i][2]);

		K_1_3 = (1.0/3.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4]);

//added to remove distance
		res[i+2] += 0.5*K_1_3;

//orig		res[i+2] += 0.5*K_1_3*pow(backboneDists[i].first - r_h ,2);

		resNum1 = forceFieldKey.find(helixFourMers[i][0]);
		resNum2 = forceFieldKey.find(helixFourMers[i][1]);
		resNum3 = forceFieldKey.find(helixFourMers[i][2]);
		resNum4 = forceFieldKey.find(helixFourMers[i][3]);

		K_1_4 = (1.0/4.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4] + forceFieldMatrix[resNum4][4]);

//added to remove distance
		res[i+2] += 0.5*K_1_4;

//orig		res[i+2] += 0.5*K_1_4*pow(backboneDists[i].second - r_h ,2);
	}

	return res;
}

vector<float> GrahnenModel::sampleHelix_pos2(vector<string> threeMers, vector<string> helixFourMers, vector< pair<float,float> > backboneDists, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(threeMers.size() != helixFourMers.size() || threeMers.size() != backboneDists.size()){
		string message = "ERROR: number of C-alpha distances does not match number of k-mers!";
		throw runtime_error(message);
	}

	unsigned int i;
	int resNum1, resNum2, resNum3, resNum4;
	float K_1_3, K_1_4;
	float r_h = 5.5;
	float sample = 0.0;

	vector<float>res(threeMers.size()+5,0.0);

	myRng random(randGen);
	random_shuffle(backboneDists.begin(), backboneDists.end(), random);

	for(i = 0; i < threeMers.size(); i++){
		resNum1 = forceFieldKey.find(threeMers[i][0]);
		resNum2 = forceFieldKey.find(threeMers[i][1]);
		resNum3 = forceFieldKey.find(threeMers[i][2]);

//		K_1_3 = (1.0/3.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4]);

		float k = 0.5*(1.0/3.0)*pow(backboneDists[i].first - r_h ,2);

		res[i+2]+=forceFieldMatrix[resNum1][4]*k;
		res[i+3]+=forceFieldMatrix[resNum2][4]*k;
		res[i+4]+=forceFieldMatrix[resNum3][4]*k;


		resNum1 = forceFieldKey.find(helixFourMers[i][0]);
		resNum2 = forceFieldKey.find(helixFourMers[i][1]);
		resNum3 = forceFieldKey.find(helixFourMers[i][2]);
		resNum4 = forceFieldKey.find(helixFourMers[i][3]);

//		K_1_4 = (1.0/4.0)*(forceFieldMatrix[resNum1][4] + forceFieldMatrix[resNum2][4] + forceFieldMatrix[resNum3][4] + forceFieldMatrix[resNum4][4]);

		k = 0.5*(1.0/4.0)*pow(backboneDists[i].second - r_h ,2);

		res[i+2]+=forceFieldMatrix[resNum1][4]*k;
		res[i+3]+=forceFieldMatrix[resNum2][4]*k;
		res[i+4]+=forceFieldMatrix[resNum3][4]*k;
		res[i+5]+=forceFieldMatrix[resNum4][4]*k;

	}

	return res;
}

vector<float> GrahnenModel::sampleBeta_pos(vector<string> betaFourMers, vector<float> backboneTorsions, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(betaFourMers.size() != backboneTorsions.size()){
		string message = "ERROR: number of k-mers do not match number of torsions!";
		throw runtime_error(message);
	}
	unsigned int i;
	int resNum1, resNum2, resNum3, resNum4;
	float K_1_4;
	float phi_b = -150; // Equilibrium Ca-Ca torsion angle (Ca_i-1->Ca_i->Ca_i+1->Ca+i+2 torsion) for beta sheets. Equivalent to 210 positive degrees.
	float C_s = 0.01; // Scaling constant to make similar to calcHelixPotential()
	float sample = 0.0;

	vector<float>res(betaFourMers.size()+4,0.0);

	myRng random(randGen);
	random_shuffle(backboneTorsions.begin(), backboneTorsions.end(), random);

	for(i = 0; i < betaFourMers.size(); i++){
		resNum1 = forceFieldKey.find(betaFourMers[i][0]);
		resNum2 = forceFieldKey.find(betaFourMers[i][1]);
		resNum3 = forceFieldKey.find(betaFourMers[i][2]);
		resNum4 = forceFieldKey.find(betaFourMers[i][3]);

		K_1_4 = 0.25*( forceFieldMatrix[resNum1][5] + forceFieldMatrix[resNum2][5] + forceFieldMatrix[resNum3][5] + forceFieldMatrix[resNum4][5]);

		res[i+1] += K_1_4*pow( C_s*(backboneTorsions[i] - phi_b), 2);
	}

	return res;
}

vector<float> GrahnenModel::sampleBeta_pos2(vector<string> betaFourMers, vector<float> backboneTorsions, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(betaFourMers.size() != backboneTorsions.size()){
		string message = "ERROR: number of k-mers do not match number of torsions!";
		throw runtime_error(message);
	}
	unsigned int i;
	int resNum1, resNum2, resNum3, resNum4;
	float K_1_4;
	float phi_b = -150; // Equilibrium Ca-Ca torsion angle (Ca_i-1->Ca_i->Ca_i+1->Ca+i+2 torsion) for beta sheets. Equivalent to 210 positive degrees.
	float C_s = 0.01; // Scaling constant to make similar to calcHelixPotential()
	float sample = 0.0;

	vector<float>res(betaFourMers.size()+4,0.0);

	myRng random(randGen);
	random_shuffle(backboneTorsions.begin(), backboneTorsions.end(), random);

	for(i = 0; i < betaFourMers.size(); i++){
		resNum1 = forceFieldKey.find(betaFourMers[i][0]);
		resNum2 = forceFieldKey.find(betaFourMers[i][1]);
		resNum3 = forceFieldKey.find(betaFourMers[i][2]);
		resNum4 = forceFieldKey.find(betaFourMers[i][3]);

//		K_1_4 = 0.25*( forceFieldMatrix[resNum1][5] + forceFieldMatrix[resNum2][5] + forceFieldMatrix[resNum3][5] + forceFieldMatrix[resNum4][5]);

		float k = 0.25*pow( C_s*(backboneTorsions[i] - phi_b), 2);

		res[i+1]+=forceFieldMatrix[resNum1][5]*k;
		res[i+2]+=forceFieldMatrix[resNum2][5]*k;
		res[i+3]+=forceFieldMatrix[resNum3][5]*k;
		res[i+4]+=forceFieldMatrix[resNum4][5]*k;
	}

	return res;
}
vector<float> GrahnenModel::sampleCharge_pos(string sequence, vector<float> cBetaDists, vector<float> exposures, CRandomMersenne& randGen){
	/* Validate input.
	 *
	 * In an NxN matrix of pairwise distances, excluding those separated by
	 * less than a distance s in primary sequence, we should find
	 *
	 * (N-s)*((N-s) - 1) / 2 + (N-s)
	 *
	 * such distances measured in total.
	 */
	unsigned int expectedCount = (sequence.length()-SEQ_SEPARATION)*( (sequence.length()-SEQ_SEPARATION) - 1 ) / 2 + (sequence.length() - SEQ_SEPARATION);
	if(expectedCount != cBetaDists.size() || sequence.length() != exposures.size()){
		string message = "ERROR: number of bead distances or exposures do not match sequence length!";
		throw runtime_error(message);
	}

	unsigned int i, j;
	int resNum1, resNum2;
	float q_1, q_2, k_e;

	double chargeScaleFactor = 1000.0; // Scales charge term to approx same magnitude as LJ
	float sample = 0.0;
	unsigned int counter = 0;

	vector<float>res(sequence.length(),0.0);

	myRng random(randGen);
	random_shuffle(cBetaDists.begin(), cBetaDists.end(), random);
	random_shuffle(exposures.begin(), exposures.end(), random);

	for(i = 0; i < sequence.length(); i++){
		resNum1 = forceFieldKey.find(sequence[i]);
		q_1 = forceFieldMatrix[resNum1][6];

		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){
			resNum2 = forceFieldKey.find(sequence[j]);
			q_2 = forceFieldMatrix[resNum2][6];

			if(q_1 != 0 && q_2 != 0){
				k_e = chargeScaleFactor*(1/dielectricConstant(exposures[i], exposures[j]));
				res[i] += 0.5*k_e*((q_1*q_2)/cBetaDists[counter]);
				res[j] += 0.5*k_e*((q_1*q_2)/cBetaDists[counter]);
			}

			counter++;
		}
	}

	return res;
}
vector<float> GrahnenModel::sampleSolv_pos(string sequence, vector<float> exposures, CRandomMersenne& randGen){
	/* Validate input.
	 */
	if(sequence.length() != exposures.size()){
		string message = "ERROR: number of exposure values does not match sequence length!";
		throw runtime_error(message);
	}

	unsigned int i;
	int resNum;
	float sample = 0.0;
	vector<float>res(sequence.length(),0.0);

	myRng random(randGen);
	random_shuffle(exposures.begin(), exposures.end(), random);

	for(i = 0; i < sequence.length(); i++){
		resNum = forceFieldKey.find(sequence[i]);
		res[i] += forceFieldMatrix[resNum][3]*exposures[i] + forceFieldMatrix[resNum][7]*(1-exposures[i]);
	}

	return res;
}
vector<float> GrahnenModel::sampleSSBond_pos(string sequence, vector<float> cBetaDists, CRandomMersenne& randGen){
	/* Validate input.
	 *
	 * In an NxN matrix of pairwise distances, excluding those separated by
	 * less than a distance s in primary sequence, we should find
	 *
	 * (N-s)*((N-s) - 1) / 2 + (N-s)
	 *
	 * such distances measured in total.
	 */
	unsigned int expectedCount = (sequence.length()-SEQ_SEPARATION)*( (sequence.length()-SEQ_SEPARATION) - 1 ) / 2 + (sequence.length() - SEQ_SEPARATION);
	if(expectedCount != cBetaDists.size()){
		string message = "ERROR: number of bead distances does not match sequence length!";
		throw runtime_error(message);
	}

	unsigned int i, j;

	float maxSSDistance = 4.5; // Threshold for disulfide bond distance
	float sample = 0.0;
	unsigned int counter = 0;

	vector<float>res(sequence.length(),0.0);

	myRng random(randGen);
	random_shuffle(cBetaDists.begin(), cBetaDists.end(), random);

	for(i = 0; i < sequence.length(); i++){
		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){
			if(sequence[i] == 'C' && sequence[j] == 'C' && cBetaDists[counter] < maxSSDistance){
					res[i] -= 0.5;
					res[j] -= 0.5;
			}
			counter++;
		}
	}

	return res;
}

/* Helper function to calculate a SCWRL-derived
 * steric repulsion term under this model, to the
 * tune of
 *
 * E(r) = 0			if r > R_ij
 * 	 10			if r < 0.8254*R_ij
 * 	 57.273*(1-r/R_ij) 	if 0.8254*R_ij <= r <= R_ij
 *
 * where r is the Cb_i-Cb_j bead distance and R_ij is the
 * sum of Cb radii, or r is the Cb_i-Ca_j distance and
 * R_ij is the sum of the Cb_i radius and Ca_j radius (always
 * 0.9 A).
 *
 * See Canutescu et al, 2003, Protein Science 1:2001-2014
 * for further details.
 *
 * NOW WITH ONLY ACTIVE-RESIDUE COMPUTATION.
 *
 * RETURNS 1e10 ON ERROR.
 */
float GrahnenModel::scwrlRepulsionEnergy(shared_ptr<Structure> protein, vector<bool> residueCanMove){
	unsigned int i, j, resNumOne, resNumTwo, proteinLength;
	float radiusOne, radiusTwo, distance, collisionDistance;
	string typeOne, typeTwo;
	Residue residueOne, residueTwo;
	coord positionOne, positionTwo;

	float energy = 0.0;
	float caRadius = 0.9;
	proteinLength = protein->getNumResidues();

	/* Check that we have a long enough list of active residues.
	 */
	if(proteinLength != residueCanMove.size()){
		return 1000000000.0;
	}

	/* Do all pairwise comparisons to detect Cb collisions
	 * with Ca and Cb beads, but only for residues that
	 * are marked as being able to move.
	 */
	for(i = 0; i < proteinLength; i++){

		/* Only calculate energies for moving residues.
		 */
		if(residueCanMove[i]){
			residueOne = protein->getSingleResidue(i);
			typeOne = residueOne.getOneLetterType();

			/* Nevermind about glycines, they don't have a CB bead.
			 */
			if(typeOne.compare("G") != 0){
				resNumOne = forceFieldKey.find(typeOne);
				radiusOne = forceFieldMatrix[resNumOne][1];
				positionOne = residueOne.getAtomCoordsByType("CB");

				/* Check for Cb and Ca collision against all other
				 * residues.
				 */
				for(j = 0; j < proteinLength; j++){

					/* But don't compare with yourself, or with residues
					 * you've already compared with.
					 */
					if(i != j && !(residueCanMove[i] && residueCanMove[j] && i > j) ){
						residueTwo = protein->getSingleResidue(j);
						typeTwo = residueTwo.getOneLetterType();

						/* Check for Cb-Cb collision for everything except
						 * glycines (they don't have a Cb bead).
						 */
						if(typeTwo.compare("G") != 0){
							resNumTwo = forceFieldKey.find(typeTwo);
							radiusTwo = forceFieldMatrix[resNumTwo][1];
							positionTwo = residueTwo.getAtomCoordsByType("CB");
							distance = NDDist(positionOne, positionTwo);
							collisionDistance = radiusOne+radiusTwo;

							/* Spheres are overlapping significantly, maximum energy cost.
							 */
							if(distance < 0.8254*collisionDistance){
								energy += 10.0;
							}
							/* Spheres are overlapping a bit, some energy cost.
							 */
							else if(distance >= 0.8254*collisionDistance && distance <= collisionDistance){
								energy += 57.273*(1-(distance/collisionDistance));
							}
						}

						/* Also check for Cb-Ca collision for all residue types.
						 */
						positionTwo = residueTwo.getAtomCoordsByType("CA");
						distance = NDDist(positionOne, positionTwo); 
						collisionDistance = radiusOne + caRadius;

						/* Spheres are overlapping significantly, maximum energy cost.
						 */
						if(distance < 0.8254*collisionDistance){
							energy += 10.0;
						}
						/* Spheres are overlapping a bit, some energy cost.
						 */
						else if(distance >= 0.8254*collisionDistance && distance <= collisionDistance){
							energy += 57.273*(1-(distance/collisionDistance));
						}

					}
				}
			}
		}
	}

	return energy;
}

// UNDER DEVELOPMENT...
//
// ESSENTIALLY calcNonbondPotential BUT ONLY FOR CERTAIN
// RESIDUES, AND WITH CAPPED BOTTOM VALUE AT 0.
float GrahnenModel::LJRepulsionEnergy(shared_ptr<Structure> protein, vector<bool> residueCanMove){
	unsigned int i, j, resNumOne, resNumTwo, proteinLength;
	float radiusOne, radiusTwo, distance, collisionDistance, epsilonOne, epsilonTwo;
	string typeOne, typeTwo;
	Residue residueOne, residueTwo;
	coord positionOne, positionTwo;

	float energy = 0.0;
	float caRadius = 0.9;
	float caEpsilon = 0.05;
	proteinLength = protein->getNumResidues();

	/* Check that we have a long enough list of active residues.
	 */
	if(proteinLength != residueCanMove.size()){
		return 1000000000.0;
	}

	/* Do all pairwise comparisons to calculate Cb interaction
	 * energies with Ca and Cb beads, but only for residues that
	 * are marked as being able to move.
	 */
	for(i = 0; i < proteinLength; i++){

		/* Only calculate energies for moving residues.
		 */
		if(residueCanMove[i]){
			residueOne = protein->getSingleResidue(i);
			typeOne = residueOne.getOneLetterType();

			/* Nevermind about glycines, they don't have a Cb bead.
			 */
			if(typeOne.compare("G") != 0){
				resNumOne = forceFieldKey.find(typeOne);
				radiusOne = forceFieldMatrix[resNumOne][1];
				epsilonOne = forceFieldMatrix[resNumOne][3];
				positionOne = residueOne.getAtomCoordsByType("CB");

				/* Check for Cb and Ca interaction with all other
				 * residues.
				 */
				for(j = 0; j < proteinLength; j++){

					/* But don't compare with yourself, or with residues
					 * you've already compared with.
					 */
					if(i != j && !(residueCanMove[i] && residueCanMove[j] && i > j) ){
						residueTwo = protein->getSingleResidue(j);
						typeTwo = residueTwo.getOneLetterType();

						/* Calculate for Cb-Cb interactions for everything except
						 * glycines (use the Ca bead position there).
						 */
						if(typeTwo.compare("G") == 0){
							positionTwo = residueTwo.getAtomCoordsByType("CA");
						}
						else{
							positionTwo = residueTwo.getAtomCoordsByType("CB");
						}

						resNumTwo = forceFieldKey.find(typeTwo);
						radiusTwo = forceFieldMatrix[resNumTwo][1];
						epsilonTwo = forceFieldMatrix[resNumTwo][3];
						distance = NDDist(positionOne, positionTwo);
						collisionDistance = radiusOne+radiusTwo;

						/* Calculate capped version of LJ potential between
						 * the Cb beads.
						 */
						energy += 4*sqrt(epsilonOne*epsilonTwo)*(pow(collisionDistance/distance,12) - pow(collisionDistance/distance,6) + 0.25);

						/* Also calculate Cb-Ca interactions for all residue types,
						 * except Gly (already used the Ca bead once).
						 */
						if(typeTwo.compare("G") != 0){
							positionTwo = residueTwo.getAtomCoordsByType("CA");
							distance = NDDist(positionOne, positionTwo); 
							collisionDistance = radiusOne + caRadius;

							/* Calculate capped version of LJ potential between
							 * the Cb and Ca beads.
							 */
							energy += 4*sqrt(epsilonOne*caEpsilon)*(pow(collisionDistance/distance,12) - pow(collisionDistance/distance,6) + 0.25);
						}
					}
				}
			}
		}
	}

	return energy;
}

/* Returns the default bond length (Ca-Cb bead distance in
 * Angstroms) under the Mukherjee/Rastogi/Grahnen 2-bead energy 
 * models.
 *
 * Input should be a single character string.
 *
 * Returns -1.0 if input doesn't correspond to an existing
 * one-letter residue abbreviation under this model.
 */
double GrahnenModel::getDefaultBondLength(string resType){

	/* Poorly specified residue type.
	 */
	if(resType.length() != 1){
		return -1.0;
	}

	unsigned int resNum = forceFieldKey.find(resType);

	/* Residue type does not exist.
	 */
	if(resNum == string::npos){
		return -1.0;
	}
	
	return (double) forceFieldMatrix[resNum][0];
}
