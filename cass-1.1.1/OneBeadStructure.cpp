#include "OneBeadStructure.h"
#include <stdexcept>

OneBeadStructure::OneBeadStructure():
	Structure()
{
}

OneBeadStructure::OneBeadStructure(vector<Residue> residues):
	Structure(residues)
{
}

OneBeadStructure::OneBeadStructure(shared_ptr<Structure> detailedStructure){
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
		
		/* Make sure there are some Ca coordinates, though.
		 */
		if(cAlphaBead == nowhere){
			string message = "Can't convert to OneBeadStructure: residue "+itos(i+1)+" has no C-alpha atom!";
			throw runtime_error(message);
		}
		
		beadCoords.push_back(cAlphaBead);
		beadTypes.push_back("CA");

		/* Make a Residue from that bead only and add it to
		 * the list of residues.
		 */
		Residue tempResidue(allAtomResidue.getType(), beadTypes, beadCoords);
		residues.push_back(tempResidue);
		beadCoords.clear();
		beadTypes.clear();
	}

	/* Be sure that you actually found some Ca atoms. Otherwise, throw
	 * an exception.
	 */
	if(residues.size() == 0){
	}

	/* And save them.
	 */
	residueList = residues;
}

OneBeadStructure::~OneBeadStructure(){
	/* DEBUG
	 *
	cout << "Destroying OneBeadStructure class properties..." << endl;
	*/
}

shared_ptr<Structure> OneBeadStructure::getBackbone() const{
	return shared_ptr<Structure>(new OneBeadStructure(*this));
}

shared_ptr<Structure> OneBeadStructure::threadSequence(string sequence){

	vector<Residue> newResidues;
	unsigned int i;

	/* Switch out residue type for all the residues, but nothing
	 * else.
	 */
	for(i = 0; i < residueList.size(); i++){
		Residue tempResidue(residueList[i], threeLetterName(sequence.substr(i,1)));
		newResidues.push_back(tempResidue);
	}

	return shared_ptr<Structure>(new OneBeadStructure(newResidues));
}

shared_ptr<Structure> OneBeadStructure::adjustSidechains(){
	return shared_ptr<Structure>(this);
}

vector< shared_ptr<Structure> > OneBeadStructure::adjustComplex(vector< shared_ptr<Structure> > components){
	return components;
}
