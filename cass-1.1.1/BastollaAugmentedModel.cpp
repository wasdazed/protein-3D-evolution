#include "BastollaAugmentedModel.h"
#include "Structure.h"
#include "TwoBeadStructure.h"
#include <typeinfo>
#include <stdexcept>
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


/* Default constructor, doesn't do much.
 */
BastollaAugmentedModel::BastollaAugmentedModel():
	BastollaModel()
{
}

/* Destructor, required due to gcc being
 * finicky.
 */
BastollaAugmentedModel::~BastollaAugmentedModel(){
}

/* Calculating folding score/energy.
 */
float BastollaAugmentedModel::calcFoldScore(shared_ptr<Structure> protein){

	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only TwoBeadStructure is valid for BastollaAugmentedModel.");
	}

	/* Start the clock.
	*
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);
	*/

	float eng = 0.0;
	string res1, res2;
	unsigned int i, j, resNum1, resNum2;

	/* Construct a contact map at 4.5 A between heavy
	 * atoms.
	 *
	 * NOT VERY EFFICIENT. CAN THIS BE OPTIMIZED AWAY?
	 */
	vector< vector<int> > contactMap = protein->contactMap(4.5);
	string sequence = protein->getSequence();

	/* For residues separated by 3 or more positions in the
	 * structure and in contact, add their interaction energy
	 * to the total.
	 */
	for(i = 0; i < sequence.length(); i++){
		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){ 

			/* DEBUG
			 *
			 cout << sequence.substr(i,1) << i+1 << "-" << sequence.substr(j,1) << j+1 << endl;
			 */

			if(contactMap[i][j] == 1){

				/* Find the correct position in the Bastolla matrix for
				 * each residue.
				 */
				res1 = sequence.substr(i,1);
				res2 = sequence.substr(j,1);
				resNum1 = bastollaKey.find(res1);
				resNum2 = bastollaKey.find(res2);

				/* DEBUG
				 *
				 cout << res1 << " has position " << resNum1 << " in the matrix (start at 0)." << endl;
				 cout << res2 << " has position " << resNum2 << " in the matrix (start at 0)." << endl;
				 */

				/* Retrive the contact energy and add to
				 * total.
				 */
				eng += bastollaMatrix[resNum1][resNum2];

				/* DEBUG
				 *
				 cout << "Eng = " << bastollaMatrix[resNum1][resNum2] << endl;
				 */

			}
		}
	}

	/* Stop the clock and print.
	*
	gettimeofday(&tim, NULL);
	double timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
	cout << (timeStop - timeStart) << endl;
	*/

	return eng;
}

/* Calculating interaction score/energy.
 *
 * Based on assumption that Bastolla's contact energies
 * for residues work for inter- as well as intra-protein
 * contacts. Note that the model was NOT originially
 * parameterized for this purpose.
 */
float BastollaAugmentedModel::calcInteractionScore(shared_ptr<Structure> interactorOne, shared_ptr<Structure> interactorTwo){
	
	/* Crash if you're being passed the wrong type
	 * of pointer(s).
	 */
	TwoBeadStructure test;
	if(typeid(*interactorOne) != typeid(*interactorTwo)){
		throw runtime_error("Interacting structures must be of the same type.");
	}
	else{
		if(typeid(test) != typeid(*interactorOne) || typeid(test) != typeid(*interactorTwo)){
			throw runtime_error("Only TwoBeadStructure is valid for BastollaAugmentedModel.");
		}
	}

	/* Start the clock.
	 *
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);
	*/

	float bindEnergy = 0.0;
	unsigned int i, j;
	string res1, res2;
	int resNum1, resNum2;

	/* Examine contact status (at 4.5 A) of all the residues in one protein
	 * to all the residues in the other one, and calculate the 
	 * binding energy based on that.
	 */
	vector< vector<int> > contactMap = interactorOne->interactionContactMap(4.5, interactorTwo);
	vector<Residue> residuesOne = interactorOne->getAllResidues();
	vector<Residue> residuesTwo = interactorTwo->getAllResidues();

	for(i = 0; i < contactMap.size(); i++){
		Residue resOne = residuesOne[i];
		for(j = 0; j < contactMap[0].size(); j++){
			Residue resTwo = residuesTwo[j];

			/* Calculate the energy between all residues that are
			 * in contact.
			 *
			 * THIS IS GOING TO BE A SLIGHT PROBLEM... MIGHT WANT TO
			 * WRITE UP A INTER-STRUCTURE CONTACT MAP FUNCTION IN
			 * TWOBEADSTRUCTURE.
			 */
			if( contactMap[i][j] == 1){

				/* Find the correct position in the Bastolla matrix for
				 * each residue.
				 */
				res1 = resOne.getOneLetterType();
				res2 = resTwo.getOneLetterType();
				resNum1 = bastollaKey.find(res1);
				resNum2 = bastollaKey.find(res2);

				/* Add this to the binding energy.
				*/
				bindEnergy += bastollaMatrix[resNum1][resNum2];
			}
		}
	}

	/* Stop the clock and print.
	 *
	gettimeofday(&tim, NULL);
	double timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
	cout << (timeStop - timeStart) << endl;
	*/

	return bindEnergy;
}

// COPY COMMENTS FROM GrahnenModel VERSION OF THE SAME...
float BastollaAugmentedModel::calcFoldEnergyGap(shared_ptr<Structure> protein){
	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only TwoBeadStructure is valid for BastollaAugmentedModel.");
	}

	float natScore = 0.0;
	unsigned int i, j, k, resNum1, resNum2;
	float gapScore, randScore;
	vector<float> randEngDist;
	vector<float> quarts;

	/* Construct a contact map at 4.5 A between heavy
	 * atoms.
	 */
	vector< vector<int> > contactMap = protein->contactMap(4.5);
	string sequence = protein->getSequence();

	/* For residues separated by SEQ_SEPARATION or more positions in the
	 * structure and in contact, add their interaction energy
	 * to the total.
	 */
	for(i = 0; i < sequence.length(); i++){
		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){ 

			if(contactMap[i][j] == 1){

				/* Find the correct position in the Bastolla matrix for
				 * each residue.
				 */
				resNum1 = bastollaKey.find(sequence[i]);
				resNum2 = bastollaKey.find(sequence[j]);

				/* Retrive the contact energy and add to
				 * total.
				 */
				natScore += bastollaMatrix[resNum1][resNum2];
			}
		}
	}

	/* DEBUG
	 *
	cerr << "Native contact map:" << endl;
	for(i = 0; i < contactMap.size(); i++){
		for(j = 0; j < contactMap.size(); j++){
			cerr << contactMap[i][j];
		}
		cerr << endl;
	}
	cerr << endl;
	*/

	/* Sample a number of random structural contexts by
	 * randomizing the contact map.
	 *
	 * However, be sure to always do so in the _same_
	 * fashion for the _same_ sequence (ensuring a
	 * a stable score).
	 */
	std::hash<string> h;
	int seed = h(sequence);
	CRandomMersenne newRNG(seed);
	myRng random(newRNG);
	for(i = 0; i < NUM_RANDOM_STRUCTURES; i++){
		randScore = 0.0;
		/* Randomize the contact map.
		*/
		for(j = 0; j < contactMap.size()-SEQ_SEPARATION; j++){
			random_shuffle(contactMap[j].begin()+(j+SEQ_SEPARATION), contactMap[j].end(), random);
		}

		/* DEBUG
		*
		cerr << "Randomized contact map:" << endl;
		for(j = 0; j < contactMap.size(); j++){
			for(k = 0; k < contactMap.size(); k++){
				cerr << contactMap[j][k];
			}
			cerr << endl;
		}
		cerr << endl;
		*/

		/* Calculate the score of such a map.
		 */
		for(j = 0; j < contactMap.size(); j++){
			for(k = j + SEQ_SEPARATION; k < contactMap[j].size(); k++){
				if(contactMap[j][k] == 1){
					randScore += singleContact(sequence[i], sequence[j]);
				}
			}
		}
		randEngDist.push_back(randScore);
	}

	/* Calculate the gap score between native and the distribution
	 * of random structure scores.
	 */
	quarts = quartiles(randEngDist);
	//gapScore = (natScore - quarts[2]) / (quarts[2]-quarts[1]);
	gapScore = natScore - quarts[2];

	/* DEBUG
	 *
	cerr << "dG = " << natScore << " - " << quarts[2] << " = " << gapScore << endl;
	*/

	return gapScore;
}

// SUPPLIES POSSIBILITY OF KEEPING GEOMETRIC CONSTRAINTS (E.G.
// CONTACT MAP) CONSTANT INSTEAD OF RE-CALCULATING FROM SUPPLIED
// STRUCTURE.
float BastollaAugmentedModel::calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation){
	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only TwoBeadStructure is valid for BastollaAugmentedModel.");
	}

	float natScore = 0.0;
	unsigned int i, j, k, resNum1, resNum2;
	float gapScore, randScore;
	vector<float> randEngDist;
	vector<float> quarts;

	string sequence = protein->getSequence();
	vector< vector<int> > contactMap = protein->contactMap(4.5); // 4.5 A contact distance
	vector< vector<int> > otherContactMap = compactConformation->contactMap(4.5);

	/* For residues separated by SEQ_SEPARATION or more positions in the
	 * structure and in contact, add their interaction energy
	 * to the total.
	 */
	for(i = 0; i < sequence.length(); i++){
		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){ 

			if(contactMap[i][j] == 1){

				/* Find the correct position in the Bastolla matrix for
				 * each residue.
				 */
				resNum1 = bastollaKey.find(sequence[i]);
				resNum2 = bastollaKey.find(sequence[j]);

				/* Retrive the contact energy and add to
				 * total.
				 */
				natScore += bastollaMatrix[resNum1][resNum2];
			}
		}
	}

	/* DEBUG
	 *
	cerr << "Native contact map:" << endl;
	for(i = 0; i < contactMap.size(); i++){
		for(j = 0; j < contactMap.size(); j++){
			cerr << contactMap[i][j];
		}
		cerr << endl;
	}
	cerr << endl;
	*/

	/* Sample a number of random structural contexts by
	 * randomizing the supplied contact map.
	 *
	 * However, be sure to always do so in the _same_
	 * fashion for the _same_ sequence (ensuring a
	 * a stable score).
	 */
	std::hash<string> h;
	int seed = h(sequence);
	CRandomMersenne newRNG(seed);
	myRng random(newRNG);
	for(i = 0; i < NUM_RANDOM_STRUCTURES; i++){
		randScore = 0.0;

		/* Randomize the supplied contact map.
		*/
		for(j = 0; j < otherContactMap.size()-SEQ_SEPARATION; j++){
			random_shuffle(otherContactMap[j].begin()+(j+SEQ_SEPARATION), otherContactMap[j].end(), random);
		}

		/* Calculate the score of such a map.
		 */
		for(j = 0; j < otherContactMap.size(); j++){
			for(k = j + SEQ_SEPARATION; k < otherContactMap[j].size(); k++){
				if(otherContactMap[j][k] == 1){
					randScore += singleContact(sequence[i], sequence[j]);
				}
			}
		}
		randEngDist.push_back(randScore);
	}

	/* Calculate the gap score between native and the distribution
	 * of random structure scores.
	 */
	quarts = quartiles(randEngDist);
	//gapScore = (natScore - quarts[2]) / (quarts[2]-quarts[1]);
	gapScore = natScore - quarts[2];

	/* DEBUG
	 *
	cerr << "dG = " << natScore << " - " << quarts[2] << " = " << gapScore << endl;
	*/

	return gapScore;

}


// EXPERIMENTAL, DECOMPOSES ENERGY INTO CONTRIBUTION OF
// INDIVIDUAL RESIDUES
vector<float> BastollaAugmentedModel::residueFoldContributions(shared_ptr<Structure> protein){
	unsigned int i,j;
	float resEng;
	string res1, res2;
	int resNum1, resNum2;
	vector<float> energies;
	
	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	TwoBeadStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only TwoBeadStructure is valid for BastollaAugmentedModel.");
	}

	vector< vector<int> > contactMap = protein->contactMap(4.5);
	string sequence = protein->getSequence();
	
	for(i = 0; i < sequence.length(); i++){
		resEng = 0.0;
		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){ 

			/* DEBUG
			 *
			 cout << sequence.substr(i,1) << i+1 << "-" << sequence.substr(j,1) << j+1 << endl;
			 */

			if(contactMap[i][j] == 1){

				/* Find the correct position in the Bastolla matrix for
				 * each residue.
				 */
				res1 = sequence.substr(i,1);
				res2 = sequence.substr(j,1);
				resNum1 = bastollaKey.find(res1);
				resNum2 = bastollaKey.find(res2);

				/* DEBUG
				 *
				 cout << res1 << " has position " << resNum1 << " in the matrix (start at 0)." << endl;
				 cout << res2 << " has position " << resNum2 << " in the matrix (start at 0)." << endl;
				 */

				/* Retrive the contact energy and 
				 * total.
				 */
				resEng += bastollaMatrix[resNum1][resNum2];

				/* DEBUG
				 *
				 cout << "Eng = " << bastollaMatrix[resNum1][resNum2] << endl;
				 */

			}
		}
		energies.push_back(resEng);
	}

	return energies;
}

// EXPERIMENTAL DECOMPOSITION OF BINDING ENERGY: ONLY WORKS
// W R T TO PROTEIN, ONLY DELIVERS ENERGIES FOR BINDING
// SITE RESIDUES (E.G. THOSE WITHIN 4.5 A OF THE INTERACTOR).
vector<float> BastollaAugmentedModel::residueBindContributions(shared_ptr<Structure> protein, shared_ptr<Structure> interactor){
	
	/* Crash if you're being passed the wrong type
	 * of pointer(s).
	 */
	TwoBeadStructure test;
	if(typeid(*protein) != typeid(*interactor)){
		throw runtime_error("Interacting structures must be of the same type.");
	}
	else{
		if(typeid(test) != typeid(*protein) || typeid(test) != typeid(*interactor)){
			throw runtime_error("Only TwoBeadStructure is valid for BastollaAugmentedModel.");
		}
	}

	float bindEnergy = 0.0;
	vector<float> engs;
	unsigned int i, j;
	string res1, res2;
	int resNum1, resNum2;

	/* Examine contact status (at 4.5 A) of all the residues in one protein
	 * to all the residues in the other one, and calculate the 
	 * binding energy based on that.
	 */
	vector< vector<int> > contactMap = protein->interactionContactMap(4.5, interactor);
	vector<Residue> residuesOne = protein->getAllResidues();
	vector<Residue> residuesTwo = interactor->getAllResidues();

	for(i = 0; i < contactMap.size(); i++){
		Residue resOne = residuesOne[i];
		bindEnergy = 0.0;
		for(j = 0; j < contactMap[0].size(); j++){
			Residue resTwo = residuesTwo[j];

			/* Calculate the energy between all residues that are
			 * in contact.
			 */
			if( contactMap[i][j] == 1){

				/* Find the correct position in the Bastolla matrix for
				 * each residue.
				 */
				res1 = resOne.getOneLetterType();
				res2 = resTwo.getOneLetterType();
				resNum1 = bastollaKey.find(res1);
				resNum2 = bastollaKey.find(res2);

				/* Add this to the binding energy.
				*/
				bindEnergy += bastollaMatrix[resNum1][resNum2];
			}
		}
		engs.push_back(bindEnergy);
	}

	return engs;
}

/* Function to evaluate the energetic contribution of a single
 * contact of two specified residue types.
 *
 * Returns -9999.0 on error (e.g. improper residue types, etc).
 */
float BastollaAugmentedModel::singleContact(char residueOne, char residueTwo){
	float eng = -9999.0;

	size_t resNumOne = bastollaKey.find(residueOne);
	size_t resNumTwo = bastollaKey.find(residueTwo);

	if(resNumOne != string::npos && resNumTwo != string::npos){
		eng = bastollaMatrix[resNumOne][resNumTwo];
	}

	return eng;
}

