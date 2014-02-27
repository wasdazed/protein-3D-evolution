#include "Residue.h"

/* Default constructor.
 */
Residue::Residue(){
}

/* Constructor from values.
 */
Residue::Residue(string resType, vector<string> atoms, vector<coord> coords){
	type = resType;
	atomTypes = atoms;
	atomCoords = coords;
	apoExposure = 0.0;
	holoExposure = 0.0;
}

/* Copy constructor that switches out the atom type and coordinates.
 *
 * Doesn't play well with AllAtomStructure since it takes several
 * special cases from TwoBeadStructure into account.
 */
Residue::Residue(const Residue& org, string newResType, vector<coord> newCoords){
	atomTypes = org.atomTypes;
	apoExposure = org.apoExposure;
	holoExposure = org.holoExposure;

	/* If we're making a non-glycine residue from a glycine,
	 * or making a sidechain de novo, we'll need one more atom type 
	 * (and vice versa).
	 */
	if(org.type.compare("GLY") == 0 && newResType.compare("GLY") != 0){
		atomTypes.push_back("CB");
	}
	else if(org.type.compare("GLY") != 0 && atomTypes.size() < newCoords.size()){
		atomTypes.push_back("CB");
	}
	else if(org.type.compare("GLY") != 0 && newResType.compare("GLY") == 0){
		atomTypes.pop_back();
	}

	type = newResType;
	atomCoords = newCoords; 
}

/* Copy constructor that switches out the type but not coordinates.
 *
 * Quite dangerous to use on anything but OneBeadStructure
 * objects, since the others have varying coordinates for
 * different types.
 */
Residue::Residue(const Residue& org, string newResType){
	type = newResType;
	atomCoords = org.atomCoords;
	atomTypes = org.atomTypes;
}


/* Copy constructor that switches out the coordinates of a
 * particular atom by type. Operates in the same fashion as
 * getAtomCoordsByType() (e.g. only treats the first atom
 * of each type).
 */
Residue::Residue(const Residue& org, string atomType, coord newCoords){
	int i;
	type = org.type;
	atomTypes = org.atomTypes;
	apoExposure = org.apoExposure;
	holoExposure = org.holoExposure;
	
	/* Modify the coordinate list as appropriate.
	 */
	vector<coord> newCoordSet = org.atomCoords;
	for(i = 0; newCoordSet.size(); i++){
		if(atomTypes[i].compare(atomType) == 0){
			newCoordSet[i] = newCoords;
			break;
		}
	}
	atomCoords = newCoordSet;
}

/* Accessor methods.
 */
string Residue::getType() const{
	return type;
}

/* Stolen from Arne Elofsson and slightly
 * modified.
 */
string Residue::getOneLetterType() const{
	string onename = "X";
	if (type.compare("ALA") == 0) {onename = "A";}
	if (type.compare("ARG") == 0) {onename = "R";}
	if (type.compare("ASN") == 0) {onename = "N";}
	if (type.compare("ASP") == 0) {onename = "D";}
	if (type.compare("CYS") == 0) {onename = "C";}
	if (type.compare("GLN") == 0) {onename = "Q";}
	if (type.compare("GLU") == 0) {onename = "E";}
	if (type.compare("GLY") == 0) {onename = "G";}
	if (type.compare("HIS") == 0) {onename = "H";}
	if (type.compare("ILE") == 0) {onename = "I";}
	if (type.compare("LEU") == 0) {onename = "L";}
	if (type.compare("LYS") == 0) {onename = "K";}
	if (type.compare("MET") == 0) {onename = "M";}
	if (type.compare("PHE") == 0) {onename = "F";}
	if (type.compare("PRO") == 0) {onename = "P";}
	if (type.compare("SER") == 0) {onename = "S";}
	if (type.compare("THR") == 0) {onename = "T";}
	if (type.compare("TRP") == 0) {onename = "W";}
	if (type.compare("TYR") == 0) {onename = "Y";}
	if (type.compare("VAL") == 0) {onename = "V";}
	if (type.compare("MSE") == 0) {onename = "M";}
	
	return onename;
}

unsigned int Residue::getNumAtoms() const{
	return atomTypes.size();
}

vector<string> Residue::getAllAtomTypes() const{
	return atomTypes;
}

vector<coord> Residue::getAllAtomCoords() const{
	return atomCoords;
}

/* Function to retrieve all the atom coordinates for
 * only the sidechain portion of the residue, i.e.
 * excluding C, N and O, but including CA, CB and
 * so on.
 */
vector<coord> Residue::getAllSidechainCoords() const{
	vector<coord> sidechain;
	unsigned int i;

	/* Go through all the atoms and save
	 * only those that aren't in the backbone.
	 */
	for(i = 0; i < atomTypes.size(); i++){
		if( atomTypes[i].compare("C") != 0 && atomTypes[i].compare("N") != 0 && atomTypes[i].compare("O") != 0 ){
			sidechain.push_back(atomCoords[i]);
		}
	}

	return sidechain;
}

string Residue::getAtomTypeByNum(int atomNum) const{
	return atomTypes.at(atomNum);
}

coord Residue::getAtomCoordsByNum(int atomNum) const{
	return atomCoords.at(atomNum);
}

/* Returns the coordinates of an atom with a certain
 * type. If no atom of that type is found, returns
 * coordinates (9999,9999,9999).
 *
 * Always returns the first occurence of a particular
 * atom type, so use more round-about methods to
 * access later ones (e.g. getAllAtomCoords()).
 */
coord Residue::getAtomCoordsByType(string typeOfAtom) const{
	unsigned int i;
	coord foundAtom;

	/* Prep the default/error output.
	 */
	for(i = 0; i < 3; i++){
		foundAtom.push_back(9999.0);
	}

	/* DEBUG
	 */
	if(atomTypes.size() != atomCoords.size()){
		cout << "WARNING: This " << type << " has " << atomTypes.size() << " atom types and " << atomCoords.size() << " atom coordinates!" << endl;
		cout << "Types: ";
		for(i = 0; i < atomTypes.size(); i++){
			cout << atomTypes[i] << ",";
		}
		cout << endl;
	}

	/* Try to find the type of atom we're
	 * looking for.
	 */
	for(i = 0; i < atomTypes.size(); i++){
		if(atomTypes.at(i).compare(typeOfAtom) == 0){
			foundAtom = atomCoords.at(i);
			break;
		}
	}

	return foundAtom;
}

/* Function to compute the Euclidian distance between
 * this Residue and another Residue. Calculates the
 * distance between the two closest atoms in the two
 * Residues.
 *
 * If the Residues are more than 9999.0 A apart,
 * returns 9999.0;
 */
float Residue::distanceTo(const Residue &otherResidue) const{
	unsigned int i, j;
	float currentShortest = 9999.0;
	double distance;

	/* Compute distance between all atoms in the
	 * two residues, and keep track of the shortest
	 * one.
	 */
	for(i = 0; i < atomCoords.size(); i++){
		for(j = 0; j < otherResidue.getNumAtoms(); j++){

			/* Take the Euclidian 3D distance between
			 * the atoms.
			 */
			distance = NDDist( atomCoords.at(i), otherResidue.getAtomCoordsByNum(j) );

			/* DEBUG
			 *
			cout << atomTypes[i] << "-" << otherResidue.getAtomType(j) << " = " << distance << " A" << endl;
			*/

			/* Compare it to the current shortest distance,
			 * and record a new one if it's shorter.
			 */
			if(distance < currentShortest){
				currentShortest = distance;
			}
		}
	}

	/* Return the minimum distance between the
	 * residues.
	 */
	return currentShortest;
}

/* Function to calculate the distance between Cb atoms in
 * two residues.
 *
 * If one or both residues lack a Cb atom, returns "9999.0".
 */
float Residue::cBetaDistance(const Residue &otherResidue) const{
	
	float distance = 9999.0;

	coord cbOne = this->getAtomCoordsByType("CB");
	coord cbTwo = otherResidue.getAtomCoordsByType("CB");
	coord error(3, 9999.0);

	/* Only compute if both residues have a Cb atom.
	 */
	if(cbOne != error && cbTwo != error){
		distance = NDDist(cbOne, cbTwo);
	}

	return distance;
}


/* Function to determine if the Euclidean distance between
 * this residue and another is smaller than some threshold.
 * Is fail-fast, unlike distanceTo().
 */
bool Residue::isCloserThan(float distance, const Residue &otherResidue) const{
	unsigned int i, j;
	double distBetween;

	/* Compute distance between all atoms in the
	 * two residues, but break if one is closer than
	 * the threshold.
	 */
	for(i = 0; i < atomCoords.size(); i++){
		for(j = 0; j < otherResidue.getNumAtoms(); j++){

			/* Take the Euclidian 3D distance between
			 * the atoms.
			 */
			distBetween = NDDist( atomCoords.at(i), otherResidue.getAtomCoordsByNum(j) );

			/* If it's smaller than the threshold, finish now.
			 */
			if(distBetween < distance){
				return true;
			}
		}
	}

	/* No distance was smaller than the threshold.
	 */
	return false;

}

// EXPERIMENTAL, SETS EXPOSURE VALUE ONCE AND FOR ALL.
void Residue::setApoExposure(float exposure){
	apoExposure = exposure;
}

float Residue::getApoExposure() const{
	return apoExposure;
}

// EXPERIMENTAL, SETS EXPOSURE VALUE ONCE AND FOR ALL.
void Residue::setHoloExposure(float exposure){
	holoExposure = exposure;
}

float Residue::getHoloExposure() const{
	return holoExposure;
}
