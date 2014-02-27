/* Class to represent a single residue in a protein structure.
 *
 * A residue is considered to be a named collection of atoms with
 * 3D coordinates and name (typically PDB ATOM codes or similar).
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef RESIDUE_H
#define RESIDUE_H

#include "myUtils.h"

using namespace std;

class Residue{
	public:
		Residue();
		Residue(string resType, vector<string> atoms, vector<coord> atomCoords);

		Residue(const Residue& org, string newResType, vector<coord> newCoords); 		
		Residue(const Residue& org, string newResType);
		Residue(const Residue& org, string atomType, coord newCoords);

		string getType() const;
		string getOneLetterType() const;
		unsigned int getNumAtoms() const;

		vector<string> getAllAtomTypes() const;
		vector<coord> getAllAtomCoords() const;
		vector<coord> getAllSidechainCoords() const;

		string getAtomTypeByNum(int atomNum) const;
		coord getAtomCoordsByNum(int atomNum) const;
		coord getAtomCoordsByType(string type) const;

		float distanceTo(const Residue &otherResidue) const;
		float cBetaDistance(const Residue &otherResidue) const;
		bool isCloserThan(float distance, const Residue &otherResidue) const;

		void setApoExposure(float exposure);
		void setHoloExposure(float exposure);
		float getApoExposure() const;
		float getHoloExposure() const;

	private:
		string type;
		float apoExposure;
		float holoExposure;
		vector<string> atomTypes;
		vector<coord> atomCoords;
};

/* Idempotency end.
 */
#endif
