/* Base class for representing protein structures.
 *
 * A Structure is just considered to be an ordered collection of
 * residues, where each residue describes a set of atom coordinates
 * that have some amino acid name (e.g. ALA, VAL, LEU, etc.).
 */

/* Ensure idempotency (no multiple includes of this)
 */
#ifndef STRUCTURE_H
#define STRUCTURE_H

/* Forward declarations.
 */

/* Includes.
 */
#include "Residue.h"
#include <tr1/memory>

class Structure{
	public:
		Structure();
		Structure(vector<Residue> residues);
		virtual ~Structure();

		vector<Residue> getAllResidues() const;
		int getNumResidues() const;
		string getSequence() const;
		Residue getSingleResidue(int residueNum) const;
		virtual shared_ptr<Structure> getBackbone() const =0;

		void setSingleResidue(int residueNum, Residue newResidue);

		virtual vector< vector<int> > contactMap(float contactDistance) const;
		virtual vector< vector<float> > distanceMap() const;
		virtual vector< vector<int> > interactionContactMap(float contactDistance, shared_ptr<Structure> interactor) const;

		string toPDBString() const;

		virtual shared_ptr<Structure> threadSequence(string sequence) =0;
		virtual shared_ptr<Structure> adjustSidechains() =0;
		virtual vector< shared_ptr<Structure> > adjustComplex(vector< shared_ptr<Structure> > components) =0;
		
	protected:
		vector<Residue> getBackboneResidues() const;
		vector<Residue> residueList;
};

/* Idempotency end.
 */
#endif
