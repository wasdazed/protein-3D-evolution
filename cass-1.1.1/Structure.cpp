#include "Structure.h"
#include <stdexcept>

/* Default constructor.
 */
Structure::Structure(){
}

/* Constructor from parameters.
 */
Structure::Structure(vector<Residue> residues){
	residueList = residues;
}

/* Virtual destructor, required by gcc.
 */
Structure::~Structure(){
	/* DEBUG
	 *
	cout << "Destroying Structure class properties..." << endl;
	*/
}

/* Accessor methods.
 */
vector<Residue> Structure::getAllResidues() const{
	return residueList;
}

int Structure::getNumResidues() const{
	return residueList.size();
}
		
string Structure::getSequence() const{
	unsigned int i;
	string seq = "";

	for(i = 0; i < residueList.size(); i++){
		seq += residueList[i].getOneLetterType();
	}

	return seq;
}

Residue Structure::getSingleResidue(int residueNum) const{
	return residueList.at(residueNum); 
}

/* Function to extract only the backbone trace from a
 * structure, defined here as the "C", "N", "O" and "CA"
 * atoms.
 *
 * Throws a runtime Exception if there aren't any backbone
 * atoms in a particular residue.
 */
vector<Residue> Structure::getBackboneResidues() const{
	unsigned int i;
	int j;
	vector<Residue> residues;
	vector<coord> beadCoords;
	vector<string> beadTypes;

	/* Go through all the Residues in the
	 * structure, and re-construct them with just
	 * the backbone atoms.
	 */
	for(i = 0; i < residueList.size(); i++){

		Residue allAtomResidue = residueList[i];
		beadCoords.clear();
		beadTypes.clear();
		
		for(j = 0; j < allAtomResidue.getNumAtoms(); j++){

			if(allAtomResidue.getAtomTypeByNum(j).compare("CA") == 0
			|| allAtomResidue.getAtomTypeByNum(j).compare("C") == 0
			|| allAtomResidue.getAtomTypeByNum(j).compare("N") == 0
			|| allAtomResidue.getAtomTypeByNum(j).compare("O") == 0
			  ){
				beadCoords.push_back(allAtomResidue.getAtomCoordsByNum(j));
				beadTypes.push_back(allAtomResidue.getAtomTypeByNum(j));
			}
		}

		/* Make sure there was at least one atom available, otherwise
		 * throw an exception.
		 */
		if(beadCoords.size() == 0){
			string message = "Can't extract backbone trace: residue "+itos(i+1)+" has no backbone atoms!";
			throw runtime_error(message);
		}

		/* Make a Residue from those atoms only and add it to
		 * the list of residues.
		 */
		Residue tempResidue(allAtomResidue.getType(), beadTypes, beadCoords);
		residues.push_back(tempResidue);
	}

	return residues;
}

/* Mutator methods.
 */
void Structure::setSingleResidue(int residueNum, Residue newResidue){
	residueList[residueNum] = newResidue;
}

/* Compute the contact map for this structure, with some
 * floating point distance threshold for contact between
 * Residues. 
 * Self-contacts (the diagonal) are always zero, and only
 * the upper triangular matrix is computed (the matrix is
 * symmetrical).
 */
vector< vector<int> > Structure::contactMap(float contactDistance) const{

	/* Create a suitably sized vector filled with 0's
	 */
	vector< vector<int> > map(residueList.size(), vector<int>(residueList.size(),0));
	unsigned int i,j;

	/* Compute the distance between all residue pairs,
	 * except self-contacts (always 0). I.e. do comparisons
	 * (i,j) but not (i,i) or (j,i) (matrix is symmetric).
	 */
	for(i = 0; i < residueList.size(); i++){
		for(j = i+1; j < residueList.size(); j++){

			/* DEBUG
			 *
			cout << residueList[i].getOneLetterType() << "-" << residueList[j].getOneLetterType() << ":" << endl;
			*/

			if(residueList[i].isCloserThan(contactDistance, residueList[j])){
				/* DEBUG
				 *
				cout << i+1 << residueList[i].getOneLetterType() << "-" << j+1 << residueList[j].getOneLetterType() << endl;
				*/

				map[i][j] = 1;
			}
		}
	}

	return map;
}

vector< vector<int> > Structure::interactionContactMap(float contactDistance, shared_ptr<Structure> interactor) const{

	/* Create a suitably sized vector filled with 0's
	 */
	vector< vector<int> > map(residueList.size(), vector<int>(interactor->getNumResidues(),0));
	vector<Residue> interactorResidues = interactor->getAllResidues();
	unsigned int i,j;

	/* Compute the distance between all residue pairs
	 * between the structures.
	 */
	for(i = 0; i < residueList.size(); i++){
		for(j = 0; j < interactorResidues.size(); j++){

			/* DEBUG
			 *
			cout << residueList[i].getOneLetterType() << "-" << residueList[j].getOneLetterType() << ":" << endl;
			*/

			if(residueList[i].isCloserThan(contactDistance, interactorResidues[j])){
				/* DEBUG
				 *
				cout << i+1 << residueList[i].getOneLetterType() << "-" << j+1 << residueList[j].getOneLetterType() << endl;
				*/

				map[i][j] = 1;
			}
		}
	}

	return map;
}


/* Function to produce a distance map (residue-residue
 * distances) for this structure.
 *
 * Self-distances (the diagonal) are always zero, and only
 * the lower triangular matrix is computed (the matrix is
 * symmetrical).
 */
vector< vector<float> > Structure::distanceMap() const{
	
	/* Create a suitably sized vector filled with 0's
	 */
	vector< vector<float> > map(residueList.size(), vector<float>(residueList.size(),0.0));
	unsigned int i,j;

	/* Compute the distance between all residue pairs,
	 * except self-distances (always 0). I.e. do comparisons
	 * (i,j) but not (i,i) or (j,i) (matrix is symmetric).
	 */
	for(i = 0; i < residueList.size(); i++){
		for(j = i+1; j < residueList.size(); j++){

			/* DEBUG
			 *
			cout << residueList[i].getOneLetterType() << "-" << residueList[j].getOneLetterType() << ":" << endl;
			*/

			map[i][j] = residueList[i].distanceTo(residueList[j]);

			/* DEBUG
			 *
			cout << map[i][j] << endl;
			*/
		}
	}

	return map;
}


/* Function write the Structure in PDB format to a
 * string. Treats the Structure as a single chain
 * A and considers elements to only have a single
 * character as identifier (reasonable for most
 * structures).
 *
 * Just like the PDB format, it has issues with
 * Structures with more than 9,999 atoms. Deviates
 * slightly from the standard format for compatability
 * with SCWRL 3.0 (atom type field is 2 spaces and
 * 3 characters, instead of 1 space and 4 characters).
 *
 * BORKEN WITH SCWRL 4.0, FIX.
 */
string Structure::toPDBString() const{
	int i, j, seqLength, resSize, atomNum;
	vector<string> allAtomTypes;
	vector<coord> allAtomCoords;
	string elementSymbol;
	
	char tempString[20];  // Make room for 20 characters (no 2 fields are longer)
	string output = "";
	
	/* DEBUG
	 *
	cout << "Starting string creation..." << endl;
	*/

	/* Go through each Residue.
	 */
	seqLength = residueList.size();
	atomNum = 1;
	for(i = 0; i < seqLength; i++){
		allAtomTypes = residueList[i].getAllAtomTypes();
		allAtomCoords = residueList[i].getAllAtomCoords();
		
		/* Go through each atom in the Residue.
		 */
		resSize = allAtomTypes.size();
		for(j = 0; j < resSize; j++){

			/* ATOM identifier.
			 */
			sprintf(tempString, "ATOM  ");
			output += tempString;

			/* Atom number (5 characters), purely
			 * sequential.
			 */
			sprintf(tempString, "%5d", atomNum); 
			output += tempString;

			/* 2 spaces and atom name (4 characters), left
			 * adjusted.
			 */
			sprintf(tempString, "  %-4.4s", allAtomTypes[j].c_str());
			output += tempString;

			/* Residue name (3 characters)
			 */
			sprintf(tempString, "%3.3s", residueList[i].getType().c_str());
			output += tempString;

			/* Space and chain identifier (1 character,
			 * always A).
			 */
			sprintf(tempString, " A");
			output += tempString;

			/* Residue number (4 characters).
			 */
			sprintf(tempString, "%4d", i+1);
			output += tempString;

			/* Insertion code (not used, empty character), 3 spaces and
			 * X coordinate (8 characters, 3 decimal points,
			 * signed if negative).
			 */
			sprintf(tempString, "    %8.3f", allAtomCoords[j][0]);
			output += tempString;

			/* Y coordinate (8 characters, 3 decimal points,
			 * signed if negative).
			 */
			sprintf(tempString, "%8.3f", allAtomCoords[j][1]);
			output += tempString;

			/* Z coordinate (8 characters, 3 decimal points,
			 * signed if negative).
			 */
			sprintf(tempString, "%8.3f", allAtomCoords[j][2]);
			output += tempString;

			/* Occupancy (6 characters, 2 decimal points),
			 * always "1.00"
			 */
			sprintf(tempString, "%6.2f", 1.00);
			output += tempString;

			/* Temperature factor (6 characters, 2 decimal
			 * points), always "0.00".
			 */
			sprintf(tempString, "%6.2f", 0.00);
			output += tempString;

			/* 10 spaces and element symbol (2 characters) -- 
			 * stripped down to one character for readability.
			 */
			elementSymbol = allAtomTypes[j].substr(0,1);
			sprintf(tempString, "          %2.2s", elementSymbol.c_str());
			output += tempString;

			/* Atom charge (2 characters, signed if negative),
			 * always blank, and end of line.
			 */
			sprintf(tempString, "  \n");
			output += tempString;

			atomNum++;
		}

		/* DEBUG
		 *
		cout << "Finished residue " << i+1 << "." << endl;
		*/
	}

	/*
	cout << "Made a string with " << output.length() << " characters." << endl;
	*/

	return output;
}
