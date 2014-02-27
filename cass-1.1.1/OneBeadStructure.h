/* Class to represent a protein structure reduced to only
 * one bead per residue, as in (SOME REFERENCE).
 */

/* Idempotency.
 */
#ifndef ONEBEADSTRUCTURE_H
#define ONEBEADSTRUCTURE_H

/* Foward declarations.
 */
class AllAtomStructure;
class TwoBeadStructure;

/* Includes.
 */
#include "Structure.h"

/* Interface.
 */
class OneBeadStructure : public Structure{
	public:
		/* Constructors and destructors.
		 */
		OneBeadStructure();
		OneBeadStructure(shared_ptr<Structure> detailedStructure);
		virtual ~OneBeadStructure();

		virtual shared_ptr<Structure> getBackbone() const;

		// AT LEAST threadSequence() SHOULD HAVE SOME SORT OF EFFECT
		// (REPLACING RESIDUE TYPES)...
		shared_ptr<Structure> threadSequence(string sequence);
		shared_ptr<Structure> adjustSidechains();
		vector< shared_ptr<Structure> > adjustComplex(vector< shared_ptr<Structure> > components);

	protected:

		OneBeadStructure(vector<Residue> residues);
};

/* Idempotency end.
 */
#endif
