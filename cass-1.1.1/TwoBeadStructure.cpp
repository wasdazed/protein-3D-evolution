#include "TwoBeadStructure.h"
#include "AllAtomStructure.h"
#include "randomgen/randomc.h"
#include "randomgen/stocc.h"
#include "MHMCMC.h"
#include "StructureProposal.h"
#include <stdexcept>
#include <limits>
class Proposal;

/* List of bead radii under this model.
 */
float TwoBeadStructure::beadRadii[20] = {0.77, 1.29, 1.54, 1.56, 1.22, 1.8, 1.25, 1.9, 2.13, 2.21, 1.43, 1.45, 1.75, 1.78,1.77, 1.08, 1.24, 2.38, 2.08, 0.0};
string TwoBeadStructure::radiusKey = "AVLICMPFYWDNQHESTRKG";

/* Default constructor, doesn't do much.
 */
TwoBeadStructure::TwoBeadStructure():
	Structure()
{
}

/* Constructor from a collection of residues.
 */
TwoBeadStructure::TwoBeadStructure(vector<Residue> residues):
	Structure(residues)
{
}

/* Destructor, for gcc's benefit.
 */
TwoBeadStructure::~TwoBeadStructure(){
	/* DEBUG
	 *
	cout << "Destroying TwoBeadStructure class properties..." << endl;
	*/
}

// WONDER IF THIS IS ACTUALLY A GOOD IDEA...?
//
// ADD SOME COMMENTS ON THIS BEING LEVITT'S IDEA.
TwoBeadStructure::TwoBeadStructure(shared_ptr<AllAtomStructure> detailedStructure){

	vector<Residue> residues;
	vector<coord> beadCoords;
	vector<string> beadTypes;
	coord nowhere(3, 9999.0);
	int i, structLength;

	/* Go through all the Residues in the
	 * structure.
	 */
	structLength = detailedStructure->getNumResidues();
	for(i = 0; i < structLength; i++){

		Residue allAtomResidue = detailedStructure->getSingleResidue(i);

		/* Extract the Ca coordinates and use those.
		*/
		coord cAlphaBead = allAtomResidue.getAtomCoordsByType("CA");
		beadCoords.push_back(cAlphaBead);
		beadTypes.push_back("CA");

		/* Make sure they're not missing, though.
		 */
		if(cAlphaBead == nowhere){
			string message = "Can't convert to TwoBeadStructure: residue "+itos(i+1)+" has no C-alpha atom!";
			throw runtime_error(message);
		}

		/* Create a Cb bead for all residue except
		 * glycines.
		 */
		if(allAtomResidue.getType().compare("GLY") != 0){

			coord cBetaAtom = allAtomResidue.getAtomCoordsByType("CB");
			
			/* Make sure it's not missing before doing anything else.
			 */
			if(cBetaAtom == nowhere){
				string message = "Can't convert to TwoBeadStructure: residue "+itos(i+1)+" has no C-beta atom!";
				throw runtime_error(message);
			}

			/* Make a Cb bead coordinate by taking the centroid of all the side-chain atoms
			 * (Ca counts as side-chain), not weighted by mass.
			 */
			coord cBetaBead = NDCentroid( allAtomResidue.getAllSidechainCoords() );

			beadCoords.push_back(cBetaBead);
			beadTypes.push_back("CB");
		}

		/* Make a Residue from those beads and add it to
		 * the list of residues.
		 */
		Residue tempResidue(allAtomResidue.getType(), beadTypes, beadCoords);
		residues.push_back(tempResidue);
		beadCoords.clear();
		beadTypes.clear();

		/* DEBUG
		 *
		cout << "Residue " << i+1 << " Ca: " << tempResidue.getAtomCoordsByType("CA")[0] << "," << tempResidue.getAtomCoordsByType("CA")[1] << "," << tempResidue.getAtomCoordsByType("CA")[2] << ". Cb: " << tempResidue.getAtomCoordsByType("CB")[1] << "," << tempResidue.getAtomCoordsByType("CB")[1] << "," << tempResidue.getAtomCoordsByType("CB")[2] << endl;
		*/
	}

	residueList = residues;
}

/* Function to extract the backbone trace from a structure.
 * See Structure::getBackboneResidues() for more detail.
 */
shared_ptr<Structure> TwoBeadStructure::getBackbone() const{
	vector<Residue> residues = this->getBackboneResidues();

	return shared_ptr<Structure>(new TwoBeadStructure(residues));
}

/* Compute the contact map for this structure, with some
 * floating point distance threshold for contact between
 * C-beta beads in different residues. 
 *
 * Note that distances are counted from the edge of the bead, 
 * _not_ the center.
 *
 * Self-contacts (the diagonal) are always zero, and only
 * the upper triangular matrix is computed (the matrix is
 * symmetrical).
 */
vector< vector<int> > TwoBeadStructure::contactMap(float contactDistance) const{
	/* Did we do this already?
	 */
	if(this->contacts.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did the contact map." << endl;
		*/

		return this->contacts;
	}

	/* Create a suitably sized vector filled with 0's
	 */
	vector< vector<int> > map(residueList.size(), vector<int>(residueList.size(),0));
	unsigned int i,j;
	int resNum1, resNum2;
	float surfaceDistance;

	/* Compute the distance between all residue pairs,
	 * except self-contacts (always 0). I.e. do comparisons
	 * (i,j) but not (i,i) or (j,i) (matrix is symmetric).
	 */
	for(i = 0; i < residueList.size(); i++){
		resNum1 = radiusKey.find(residueList[i].getOneLetterType());
		for(j = i+1; j < residueList.size(); j++){
			resNum2 = radiusKey.find(residueList[j].getOneLetterType());

			/* Reduce distance by sum of radii (i.e. measure from
			 * surface of beads).
			 */
			surfaceDistance = NDDist(residueList[i].getAtomCoordsByType("CB"), residueList[j].getAtomCoordsByType("CB")) - (beadRadii[resNum1]+beadRadii[resNum2]);
			if(surfaceDistance < contactDistance){
				map[i][j] = 1;
			}
		}
	}

	return map;
}

// CALCULATING CONTACTS WITHIN A BINARY COMPLEX.
vector< vector<int> > TwoBeadStructure::interactionContactMap(float contactDistance, shared_ptr<Structure> interactor) const{
	
	/* Check type correctness of the interacting molecule.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*interactor)){
		throw runtime_error("Can only make contact map between two TwoBeadStructures.");
	}

	/* Create a suitably sized vector filled with 0's
	 */
	int interactorSize = interactor->getNumResidues();
	vector<Residue> interactorResidues = interactor->getAllResidues();
	vector< vector<int> > map(residueList.size(), vector<int>(interactorSize,0));
	unsigned int i,j;
	int resNum1, resNum2;
	float surfaceDistance;

	/* Compute the distance between all residue pairs in
	 * the complex.
	 */
	for(i = 0; i < residueList.size(); i++){
		resNum1 = radiusKey.find(residueList[i].getOneLetterType());
		for(j = 0; j < interactorSize; j++){
			resNum2 = radiusKey.find(interactorResidues[j].getOneLetterType());

			/* Reduce distance by sum of radii (i.e. measure from
			 * surface of beads).
			 */
			surfaceDistance = NDDist(residueList[i].getAtomCoordsByType("CB"), interactorResidues[j].getAtomCoordsByType("CB")) - (beadRadii[resNum1]+beadRadii[resNum2]);
			if(surfaceDistance < contactDistance){
				map[i][j] = 1;
			}
		}
	}

	return map;

}

/* Function to produce a distance map (Cbeta-Cbeta
 * distances between residues) for this structure. Note 
 * that distances are counted from the edge of each bead, 
 * _not_ the centers.
 * 
 * Self-distances (the diagonal) are always zero, and only
 * the lower triangular matrix is computed (the matrix is
 * symmetrical).
 */
vector< vector<float> > TwoBeadStructure::distanceMap() const{

	/* Create a suitably sized vector filled with 0's
	 */
	vector< vector<float> > map(residueList.size(), vector<float>(residueList.size(),0.0));
	unsigned int i,j;
	int resNum1, resNum2;

	/* Compute the distance between all residue pairs,
	 * except self-distances (always 0). I.e. do comparisons
	 * (i,j) but not (i,i) or (j,i) (matrix is symmetric).
	 */
	for(i = 0; i < residueList.size(); i++){
		resNum1 = radiusKey.find(residueList[i].getOneLetterType());
		for(j = i+1; j < residueList.size(); j++){
			resNum2 = radiusKey.find(residueList[j].getOneLetterType());

			/* Reduce distance by some of radii (i.e. measure
			 * between bead surfaces). Glycines need checking,
			 * as previously.
			 */
			coord atomOne = residueList[i].getOneLetterType().compare("G") == 0 ? residueList[i].getAtomCoordsByType("CA") : residueList[i].getAtomCoordsByType("CB");
			coord atomTwo = residueList[j].getOneLetterType().compare("G") == 0 ? residueList[j].getAtomCoordsByType("CA") : residueList[j].getAtomCoordsByType("CB");
			map[i][j] = NDDist(atomOne, atomTwo) - (beadRadii[resNum1]+beadRadii[resNum2]);
		}
	}

	return map;
}

/* Finds the number of contacts (based on bead-bead distance
 * measured from the bead surface) that are separated by
 * more than a certain number of residues in the sequence.
 */
unsigned int TwoBeadStructure::getNumContacts(float contactDistance, unsigned int sequenceSeparation){
	int i, j;
	unsigned int numContacts = 0;
	vector< vector<int> > contactMap = this->contactMap(contactDistance);

	/* Only accept contacts that are sufficiently separated in sequence.
	 */
	for(i = 0; i < contactMap.size(); i++){
		for(j = i+1; j < contactMap[i].size(); j++){
			if(abs(i-j) > sequenceSeparation && contactMap[i][j] == 1){
				numContacts++;
			}
		}
	}

	return numContacts;
}

vector< pair<float, float> > TwoBeadStructure::bondAngles() const{

	/* Have we already done this?
	 */
	if(this->angles.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did bond angles." << endl;
		*/

		return this->angles;
	}

	int seqLength =	this->getNumResidues();
	string sequence = this->getSequence();
	string resType;
	coord cAlpha1, cAlpha2, cBeta;
	vector< pair<float, float> > angles(seqLength, pair<float,float>(0.0,0.0));
	int i;

	for(i = 1; i < seqLength; i++){
		resType = sequence[i];

		if(resType.compare("G") != 0){
			cAlpha1 = this->getSingleResidue(i-1).getAtomCoordsByType("CA");
			cAlpha2 = this->getSingleResidue(i).getAtomCoordsByType("CA");
			cBeta = this->getSingleResidue(i).getAtomCoordsByType("CB");
			angles[i].first = NDAngle(cAlpha1, cAlpha2, cBeta);
		}
		else{
			angles[i].first = numeric_limits<float>::quiet_NaN();
		}
	}
	for(i = 0; i < seqLength-1; i++){
		resType = sequence[i];

		if(resType.compare("G") != 0){
			cAlpha1 = this->getSingleResidue(i).getAtomCoordsByType("CA");
			cAlpha2 = this->getSingleResidue(i+1).getAtomCoordsByType("CA");
			cBeta = this->getSingleResidue(i).getAtomCoordsByType("CB");
			angles[i].second = NDAngle(cBeta, cAlpha1, cAlpha2);
		}
		else{
			angles[i].second = numeric_limits<float>::quiet_NaN();
		}
	}

	return angles;
}

vector<float> TwoBeadStructure::caCbDistances() const{
	/* Did we already do this?
	 */
	if(this->cacbdists.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did Ca-Cb distances." << endl;
		*/

		return this->cacbdists;
	}

	int seqLength = this->getNumResidues();
	string sequence = this->getSequence();
	coord cAlpha1, cBeta2;
	string resType2;
	int i, j;
	vector<float> dists;

	/* Measure all Ca-Cb distances.
	 * 
	 * Bug fix 2013-05-21: Now measures all the distances, as opposed to
	 * previously (Nadia Bykova).
	 */
	for(i = 0; i < seqLength; i++){
		cAlpha1 = this->getSingleResidue(i).getAtomCoordsByType("CA");
		for(j = 0; j < seqLength; j++){

			if (abs(i - j) < SEQ_SEPARATION){
				continue;
			}
			/* End of edit
			 */

			resType2 = sequence[j];	

			if(resType2.compare("G") == 0){
				dists.push_back(numeric_limits<float>::quiet_NaN());
				continue;
			}
			else{
				cBeta2 = this->getSingleResidue(j).getAtomCoordsByType("CB");
			}
	
			dists.push_back(NDDist(cAlpha1, cBeta2));
		}
	}

	return dists;
}

vector<float> TwoBeadStructure::cbCbDistances() const{
	/* Did we already do this?
	 */
	if(this->cbcbdists.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did Cb-Cb distances." << endl;
		*/

		return this->cbcbdists;
	}

	int seqLength = this->getNumResidues();
	string sequence = this->getSequence();
	int i, j;
	string resType1, resType2;
	coord cBeta1, cBeta2;
	vector<float> dists;

	for(i = 0; i < seqLength; i++){
		resType1 = sequence.substr(i,1);

		if(resType1.compare("G") == 0){
			cBeta1 = this->getSingleResidue(i).getAtomCoordsByType("CA");
		}
		else{
			cBeta1 = this->getSingleResidue(i).getAtomCoordsByType("CB");
		}

		for(j = i + SEQ_SEPARATION; j < seqLength; j++){
			resType2 = sequence[j];
			if(resType2.compare("G") == 0){
				cBeta2 = this->getSingleResidue(j).getAtomCoordsByType("CA");
			}
			else{
				cBeta2 = this->getSingleResidue(j).getAtomCoordsByType("CB");
			}
			dists.push_back(NDDist(cBeta1, cBeta2));
		}
	}

	return dists;
}

vector<string> TwoBeadStructure::threeMers() const{
	/* Did we already do this?
	 */
	if(this->threemers.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did 3-mers." << endl;
		*/

		return this->threemers;
	}

	int seqLength = this->getNumResidues();
	int i;
	string resType1, resType2, resType3, resType4;
	vector<string> mers;

	for(i = 2; i < seqLength-3; i++){
		resType1 = this->getSingleResidue(i).getOneLetterType();
		resType2 = this->getSingleResidue(i+1).getOneLetterType();
		resType3 = this->getSingleResidue(i+2).getOneLetterType();
		mers.push_back(string(resType1+resType2+resType3));
	}

	return mers;
}

vector<string> TwoBeadStructure::helixFourMers() const{
	/* Did we already do this?
	 */
	if(this->hfourmers.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did helix 4-mers." << endl;
		*/

		return this->hfourmers;
	}

	int seqLength = this->getNumResidues();
	int i;
	string resType1, resType2, resType3, resType4;
	vector<string> mers;

	for(i = 2; i < seqLength-3; i++){
		resType1 = this->getSingleResidue(i).getOneLetterType();
		resType2 = this->getSingleResidue(i+1).getOneLetterType();
		resType3 = this->getSingleResidue(i+2).getOneLetterType();
		resType4 = this->getSingleResidue(i+3).getOneLetterType();
		mers.push_back(string(resType1+resType2+resType3+resType4));
	}

	return mers;
}

vector< pair<float,float> > TwoBeadStructure::backboneDistances() const{
	/* Did we already do this?
	 */
	if(this->bbdists.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did backbone distances." << endl;
		*/

		return this->bbdists;
	}

	int seqLength = this->getNumResidues();
	int i;
	coord cAlpha1, cAlpha2;
	vector< pair<float, float> > dists;

	for(i = 2; i < seqLength-3; i++){
		cAlpha1 = this->getSingleResidue(i).getAtomCoordsByType("CA");
		cAlpha2 = this->getSingleResidue(i+2).getAtomCoordsByType("CA");
		dists.push_back(pair<float, float>(NDDist(cAlpha1, cAlpha2), 0.0));

		cAlpha2 = this->getSingleResidue(i+3).getAtomCoordsByType("CA");
		dists.back().second = NDDist(cAlpha1, cAlpha2);
	}

	return dists;
}

vector<string> TwoBeadStructure::betaFourMers() const{
	/* Did we already do this?
	 */
	if(this->bfourmers.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did beta 4-mers." << endl;
		*/

		return this->bfourmers;
	}

	int seqLength = this->getNumResidues();
	int i;
	string resType1, resType2, resType3, resType4;
	vector<string> mers;
	for(i = 1; i < seqLength-3; i++){
		resType1 = this->getSingleResidue(i-1).getOneLetterType();
		resType2 = this->getSingleResidue(i).getOneLetterType();
		resType3 = this->getSingleResidue(i+1).getOneLetterType();
		resType4 = this->getSingleResidue(i+2).getOneLetterType();
		mers.push_back(string(resType1+resType2+resType3+resType4));
	}

	return mers;
}

vector<float> TwoBeadStructure::backboneTorsions() const{
	/* Did we already do this?
	 */
	if(this->bbtors.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did backbone torsions." << endl;
		*/

		return this->bbtors;
	}
	int seqLength = this->getNumResidues();
	int i;
	coord cAlpha1, cAlpha2, cAlpha3, cAlpha4;
	vector<float> torsions;

	for(i = 1; i < seqLength-3; i++){
		coord cAlpha1 = this->getSingleResidue(i-1).getAtomCoordsByType("CA");
		coord cAlpha2 = this->getSingleResidue(i).getAtomCoordsByType("CA");
		coord cAlpha3 = this->getSingleResidue(i+1).getAtomCoordsByType("CA");
		coord cAlpha4 = this->getSingleResidue(i+2).getAtomCoordsByType("CA");
		torsions.push_back(dihedralAngle(cAlpha1, cAlpha2, cAlpha3, cAlpha4));
	}

	return torsions;
}

vector<float> TwoBeadStructure::residueExposures() const{
	/* Did we already do this?
	 */
	if(this->exposures.size() > 0){
		/* DEBUG
		 *
		cerr << "Already did exposures." << endl;
		*/

		return this->exposures;
	}

	int seqLength = this->getNumResidues();
	int i;
	vector<float> exposures;

	for(i = 0; i < seqLength; i++){
		exposures.push_back(this->getSingleResidue(i).getApoExposure());
	}

	return exposures;
}

/* Function to read a 2-bead structure directly from file (as opposed
 * to creating one from an all-atom record). Will silently "fail" by
 * reading an all-atom record as composed of only Ca and Cb atoms,
 * if that's what you supply.
 *
 * Throws runtime_error on reading failure.
 */
TwoBeadStructure TwoBeadStructure::parseBeadFile(ifstream &inputFile){

	if(!inputFile.good()){
			string message = "Can't read from structure input file!";
			throw runtime_error(message);
	}

	int resNum;
	float x, y, z;
	string line, atomType, resType, lastResType;

	coord position;
	vector<coord> coords;
	vector<string> types;

	vector<Residue> parsedResidues;

	int lastResNum = -9999;

	/* Parse out the relevant records and construct Residue objects
	 * from them.
	 */
	while( getline(inputFile, line) ){

		/* Only ATOM records are interesting here.
		*/
		if(line.substr(0,4).compare("ATOM") == 0){

			/* Extract atom type.
			*/
			atomType = trimWhitespace(line.substr(12,4));

			/* And we're only interested in Ca and Cb atoms.
			*/
			if(atomType.compare("CA") == 0 || atomType.compare("CB") == 0){

				/* Extract and convert residue number.
				*/
				resNum = atoi(line.substr(22,4).c_str());

				/* End of a residue? Then construct the object. Note that
				 * the first residue is a bit of a special case.
				*/
				if(resNum != lastResNum && lastResNum != -9999){

					Residue tempResidue(lastResType, types, coords);
					parsedResidues.push_back(tempResidue);
					types.clear();
					coords.clear();
				}

				/* Parse residue type.
				*/
				resType = line.substr(17,3);

				/* Extract and convert the coordinates.
				*/
				x = atof(line.substr(30,8).c_str());
				y = atof(line.substr(38,8).c_str());
				z = atof(line.substr(46,8).c_str());

				/* Make a single 3D coordinate of them.
				*/
				position.push_back(x);
				position.push_back(y);
				position.push_back(z);

				/* Add the atom to the lists.
				*/
				coords.push_back(position);
				types.push_back(atomType);
				position.clear();

				lastResNum = resNum;
				lastResType = resType;
			}
		}
	}

	/* Save the last residue here (since there's no change in numbering
	 * after it).
	 */
	Residue tempResidue(lastResType, types, coords);
	parsedResidues.push_back(tempResidue);

	return TwoBeadStructure(parsedResidues);
}

// MCMC TAKE ON THREADING.
shared_ptr<Structure> TwoBeadStructure::threadSequence(string sequence){

	float temperature = 1.0;
	float meanStepSize = 1.0;
	float stepSizeVariance = 1.0;
	float neighborThreshold = 0.0;
	shared_ptr<StructureProposal> pBestState;
	shared_ptr<Proposal> pCurrent;
	shared_ptr<TwoBeadStructure> pFinalStructure;
	int i;

	/* DEBUG
	 *
	cout << "Setting up MCMC threading..." << endl;
	*/

	/* Setup the MCMC run.
	 */
	shared_ptr<TwoBeadStructure> startingStructure(new TwoBeadStructure(*this));
	shared_ptr<Proposal> pStartState(new StructureProposal(startingStructure, sequence, neighborThreshold));
	shared_ptr<MHMCMC> pOneChain(new MHMCMC(pStartState, temperature, MCMC_MAX_ITERATIONS));

	/* DEBUG
	 *
	cout << "MCMC threading running..." << endl;
	*/

	/* Run to convergence or the maximum number of iterations.
	 */
	for(i = 0; i < MCMC_MAX_ITERATIONS; i++){
		pOneChain->takeStep(meanStepSize, stepSizeVariance);

		/* Check current state.
		*/
		pCurrent = pOneChain->getCurrentState();

		/* DEBUG
		 *
		cout << pCurrent->toString() << endl;
		*/

		/* Finish if you've converged already.
		*/
		if(pCurrent->getProbDensity() == 0){
			break;
		}
	}

	/* DEBUG
	 *
	cout << "Ran " << i << " steps of the MCMC chain to adjust the structure." << endl;
	*/

	/* Get the best found structure.
	 */
	pBestState = dynamic_pointer_cast<StructureProposal>(pOneChain->getBestState()); // THIS MIGHT BE A REALLY BAD IDEA. dynamic_cast<> INSTEAD? DEREFERENCE TO OBJECT?
	pFinalStructure = pBestState->getStructure();

	return pFinalStructure;
}

// MCMC TAKE ON SIDECHAIN ADJUSTMENT. NOT PARTICULARLY
// EFFICIENT COMPARED TO E.G. SCWRL.
shared_ptr<Structure> TwoBeadStructure::adjustSidechains(){

	int i;
	float temperature = 1.0;
	float meanStepSize = 1.0;
	float stepSizeVariance = 1.0;

	shared_ptr<Proposal> pCurrent;
	shared_ptr<StructureProposal> pBestState;
	shared_ptr<TwoBeadStructure> pFinalStructure;
	
	/* DEBUG
	 *
	cout << "Setting up MCMC bead adjustment..." << endl;
	*/

	/* Setup the MCMC run based on all residues moving.
	 */
	vector<bool> moving(this->getNumResidues(), true);
	shared_ptr<TwoBeadStructure> startingStructure(new TwoBeadStructure(*this));
	shared_ptr<Proposal> pStartState(new StructureProposal(startingStructure, moving));
	shared_ptr<MHMCMC> pOneChain(new MHMCMC(pStartState, temperature, MCMC_MAX_ITERATIONS));

	/* DEBUG
	 *
	cout << "MCMC adjustment running..." << endl;
	*/

	/* Run to convergence or the maximum number of iterations.
	 */
	for(i = 0; i < MCMC_MAX_ITERATIONS; i++){
		pOneChain->takeStep(meanStepSize, stepSizeVariance);

		/* Check current state.
		*/
		pCurrent = pOneChain->getCurrentState();

		/* DEBUG
		 *
		 *
		cout << pCurrent->toString() << endl;
		*/

		/* Finish if you've converged already.
		*/
		if(pCurrent->getProbDensity() == 0){
			break;
		}
	}

	/* DEBUG
	 *
	cout << "Ran " << i << " steps of the MCMC chain to adjust the structure." << endl;
	*/

	/* Get the best found structure.
	 */
	pBestState = dynamic_pointer_cast<StructureProposal>(pOneChain->getBestState()); // THIS MIGHT BE A REALLY BAD IDEA. dynamic_cast<> INSTEAD? DEREFERENCE TO OBJECT?
	pFinalStructure = pBestState->getStructure();

	return pFinalStructure;
}

// MCMC TAKE ON COMPLEX ADJUSTMENT.
//
// SOMEWHAT FIDDLIER THAN ALLATOM CASE.
//
// DO A TYPE CHECK HERE BEFORE PROCEEDING?
vector< shared_ptr<Structure> > TwoBeadStructure::adjustComplex(vector< shared_ptr<Structure> > components){
	unsigned int i, j, k, l, numResidues;
	
	vector< shared_ptr<Structure> > finalComplex;
	vector<Residue> combinedResidues;
	vector<Residue> tempResidues;
	vector<Residue> splitResidues;
	vector<Residue> myResidues;
	vector<Residue> otherResidues;

	vector<int> startingPoint;
	
	float temperature = 1.0;
	float meanStepSize = 1.0;
	float stepSizeVariance = 1.0;
	double interactionDistance = 6.0; // ARBITRARY, COME UP WITH SOMETHING MORE REASONABLE

	shared_ptr<Proposal> pCurrent;
	shared_ptr<StructureProposal> pBestState;

	/* Complex should start out as null pointers
	 * in case something fails.
	 */
	for(i = 0; i < components.size(); i++){
		shared_ptr<TwoBeadStructure> pNothing;
		finalComplex.push_back(pNothing);
	}

	/* Construct the complex as a single structure,
	 * but keep track of where each piece starts.
	 */
	combinedResidues = components[0]->getAllResidues();
	startingPoint.push_back(0);
	for(i = 1; i < components.size(); i++){
		tempResidues = components[i]->getAllResidues();
		startingPoint.push_back(combinedResidues.size());
		combinedResidues.insert(combinedResidues.end(), tempResidues.begin(), tempResidues.end());
	}
	shared_ptr<TwoBeadStructure> combinedStruct(new TwoBeadStructure(combinedResidues));
	
	/* Figure out which residues should move around.
	 *
	 * In this case, any residue that is X Angstroms from
	 * any other component in the complex.
	 */

	/* DEBUG
	 *
	cout << "Determining interface residues..." << endl;
	*/

	vector<bool> moving(combinedStruct->getNumResidues(), false);
	for(i = 0; i < components.size(); i++){
		
		myResidues = components[i]->getAllResidues();

		for(j = i+1; j < components.size(); j++){

			otherResidues = components[j]->getAllResidues();

			/* Make pairwise distance comparisons between the residue
			 * sets to determine movement status.
			 */
			for(k = 0; k < myResidues.size(); k++){
				for(l = 0; l < otherResidues.size(); l++){

					if(myResidues[k].cBetaDistance(otherResidues[l]) < interactionDistance){
						moving[k + startingPoint[i]] = true;
						moving[l + startingPoint[j]] = true;

						/* DEBUG
						 *
						cout << k+1 << myResidues[k].getOneLetterType() << " in structure " << i+1 << " is close to " << l+1 << otherResidues[l].getOneLetterType() << " in structure " << j+1 << "." << endl;
						cout << "Pair " << k + startingPoint[i] << " and " << l + startingPoint[j] << " in the overall count (" << combinedStruct->getNumResidues()-1 << " total)." << endl << endl;
						*/
					}
				}
			}
		}
	}

	/* Adjust it via MCMC (mimics induced fit).
	 */
	shared_ptr<Proposal> pStartState(new StructureProposal(combinedStruct, moving));
	shared_ptr<MHMCMC> pOneChain(new MHMCMC(pStartState, temperature, MCMC_MAX_ITERATIONS));

	/* DEBUG
	 *
	cout << "MCMC adjustment running..." << endl;
	*/

	/* Run to convergence or the maximum number of iterations.
	 */
	for(i = 0; i < MCMC_MAX_ITERATIONS; i++){
		pOneChain->takeStep(meanStepSize, stepSizeVariance);

		/* Check current state.
		*/
		pCurrent = pOneChain->getCurrentState();

		/* DEBUG
		 *
		cout << pCurrent->toString() << endl;
		*/

		/* Finish if you've converged already.
		*/
		if(pCurrent->getProbDensity() == 0){
			break;
		}
	}

	/* DEBUG
	 *
	cout << "Ran " << i << " steps of the MCMC chain to adjust the structure." << endl;
	*/

	/* Get the best found structure.
	 */
	pBestState = dynamic_pointer_cast<StructureProposal>(pOneChain->getBestState()); // THIS MIGHT BE A REALLY BAD IDEA. dynamic_cast<> INSTEAD? DEREFERENCE TO OBJECT?
	shared_ptr<TwoBeadStructure> combinedAndAdjusted = pBestState->getStructure();

	/* Re-construct the components.
	 */
	combinedResidues = combinedAndAdjusted->getAllResidues();
	numResidues = 0;
	for(i = 0; i < components.size(); i++){
		splitResidues.clear();
		for(j = 0; j < components[i]->getNumResidues(); j++){
			splitResidues.push_back(combinedResidues[numResidues]);
			numResidues++;
		}
		shared_ptr<TwoBeadStructure> tempStructure(new TwoBeadStructure(splitResidues));
		finalComplex[i] = tempStructure;
	}
	
	/* Send them out again.
	 */
	return finalComplex;
}

// EXPERIMENTAL, SETS PERMANENT EXPOSURE VALUES FOR THE CONSTITUENT
// RESIDUES.
void TwoBeadStructure::calcExposure(bool inComplex){

	/* DEBUG
	 *
	cout << "Calculating the exposure (is in complex: " << inComplex << ").. " << endl;
	*/

	int i, j;
	string resType1, resType2;
	float cbDistance, neighborCount, neighborWeight, exposure;

	coord cBeta1, cBeta2, tempVector, exposureVector;
	coord origin(3,0.0);
	vector<coord> neighborVectors;
	Residue firstResidue;
	Residue secondResidue;

	// MODIFY AS NECESSARY, THESE ARE DURHAM'S DEFAULTS
	float lowerDistBound = 3.3;
	float upperDistBound = 11.1;

	/* Go over all beads and calculate the neighbor vector for
	 * each
	 *
	 * OPTIMIZATION: VERLET/CELL LISTS OF NEIGHBORS PRECOMPUTED.
	 */
	for(i = 0; i < this->getNumResidues(); i++){

		neighborVectors.clear();
		neighborCount = 0;
		firstResidue = this->getSingleResidue(i);
		resType1 = firstResidue.getOneLetterType();

		/* Special case for glycines: Cb bead is located at Ca position.
		*/
		if(resType1.compare("G") == 0){
			cBeta1 = firstResidue.getAtomCoordsByType("CA");
		}
		else{
			cBeta1 = firstResidue.getAtomCoordsByType("CB");
		}

		/* DEBUG
		 *
		 cout << (i+1) << resType1 << " has the following neighbors (and vectors): " << endl;
		 */

		/* Check if other residues are neighbors, and calculate the vector
		 * between them in that case.
		 */
		for(j = 0; j < this->getNumResidues(); j++){

			secondResidue = this->getSingleResidue(j);
			resType2 = secondResidue.getOneLetterType();

			/* Glycines are a special case.
			*/
			if(resType2.compare("G") == 0){
				cBeta2 = secondResidue.getAtomCoordsByType("CA");
			}
			else{
				cBeta2 = secondResidue.getAtomCoordsByType("CB");
			}

			/* Don't compare with self (you're not your own neighbor).
			*/
			if(i != j){

				/* Calculate Durham's NeighborWeight (add to NeighborCount) and vector for 
				 * this residue (if it happens to be a neighbor).
				 */
				cbDistance = NDDist(cBeta1, cBeta2);
				if(cbDistance < upperDistBound){

					if(cbDistance < lowerDistBound){
						neighborWeight = 1;
					}
					else{
						neighborWeight = 0.5*(cos( ( (cbDistance - lowerDistBound) / (upperDistBound-lowerDistBound) )*PI ) + 1);
					}
					neighborCount += neighborWeight;

					/* Assume that residue #1 is the center of the
					 * coordinate system and represent vectors as
					 * coordinates.
					 */
					coord tempVector = makeRelativeVector(cBeta1, cBeta2);

					/* Scale the vector to unit length and by
					 * "neighborness".
					 */
					tempVector = scaleVector(tempVector, (neighborWeight/NDDist(origin, tempVector)));
					neighborVectors.push_back(tempVector);

					/* DEBUG
					 *
					 cout << "\t" << (j+1) << resType2 << " ( ";
					 for(int k = 0; k < tempVector.size(); k++){
					 cout << tempVector[k] << " ";
					 }
					 cout << "), neighbor-ness " << neighborWeight << "." << endl;
					 */
				}
			}
		}

		/* Check to see that there were some neighbors. If not,
		 * just set the exposure to 1.
		 */
		if(neighborVectors.size() == 0){
			exposure = 1.0;

			/* DEBUG
			 *
			 cout << "No neighbors, assume full exposure (vector norm = 1)." << endl;
			 */
		}
		/* There are neighbors, calculate the norm of the
		 * exposure vector.
		 */
		else{

			/* Calculate the combined exposure vector
			*/
			exposureVector = vectorSum(neighborVectors);

			/* Scale it to [0,1] (e.g. by the "fractional count" of neighbors)
			*/
			exposureVector = scaleVector(exposureVector, (1.0/neighborCount));
			exposure = NDDist(origin, exposureVector);

			/* DEBUG
			 *
			 cout << "Total exposure vector: ( ";
			 for(int k = 0; k < exposureVector.size(); k++){
			 cout << exposureVector[k] << " ";
			 }
			 cout << "), norm " << exposure << endl;
			 */
		}

		if(inComplex){
			residueList[i].setHoloExposure(exposure);

			/* DEBUG
			 *
			 cout << "Set holo-exposure to " << residueList[i].getHoloExposure() << endl;
			 */
		}
		else{
			residueList[i].setApoExposure(exposure);

			/* DEBUG
			 *
			 cout << "Set apo-exposure to " << residueList[i].getApoExposure() << endl;
			 */
		}
	}
}

// EXPERIMENTAL, SETS EXPOSURE FOR RESIDUES WHEN IN A COMPLEX
vector< shared_ptr<TwoBeadStructure> > TwoBeadStructure::calcComplexExposure(vector< shared_ptr<TwoBeadStructure> > components){

	/* DEBUG
	 *
	cout << "Calculating holo-exposures for a complex..." << endl;
	*/

	unsigned int i;
	int j, numResidues;
	vector<Residue> combinedResidues;
	vector<Residue> tempResidues;
	vector<Residue> splitResidues;
	vector< shared_ptr<TwoBeadStructure> > finalComplex;

	/* Construct a single structure from the complex parts.
	 */
	combinedResidues = components[0]->getAllResidues();
	for(i = 1; i < components.size(); i++){
		tempResidues = components[i]->getAllResidues();
		combinedResidues.insert(combinedResidues.end(), tempResidues.begin(), tempResidues.end());
	}
	shared_ptr<TwoBeadStructure> combinedStruct(new TwoBeadStructure(combinedResidues));

	/* Calculate exposures in the entire complex.
	 */
	combinedStruct->calcExposure(true);

	/* Re-construct the components.
	 */
	combinedResidues = combinedStruct->getAllResidues();
	numResidues = 0;
	for(i = 0; i < components.size(); i++){
		splitResidues.clear();
		for(j = 0; j < components[i]->getNumResidues(); j++){
			splitResidues.push_back(combinedResidues[numResidues]);
			numResidues++;
		}
		shared_ptr<TwoBeadStructure> tempStructure(new TwoBeadStructure(splitResidues));
		finalComplex.push_back(tempStructure);
	}
	
	/* Send them out again.
	 */
	return finalComplex;

}

// CONVENIENCE FUNCTION TO TRIGGER ALL GEOMETRY CALCULATIONS
void TwoBeadStructure::calcAndFreezeGeometries(float contactDistance){
	contacts = this->contactMap(contactDistance);
	
	angles = this->bondAngles();
	cacbdists = this->caCbDistances();
	cbcbdists = this->cbCbDistances();
	threemers = this->threeMers();
	hfourmers = this->helixFourMers();
	bfourmers = this->betaFourMers();
	bbdists = this->backboneDistances();
	bbtors = this->backboneTorsions();
	exposures = this->residueExposures();
}


