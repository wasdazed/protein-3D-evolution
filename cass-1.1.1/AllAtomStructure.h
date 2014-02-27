/* Class to represent a protein structure with all its
 * atoms, as derived from e.g. a PDB file.
 */

/* Idempotency.
 */
#ifndef ALLATOMSTRUCTURE_H
#define ALLATOMSTRUCTURE_H

/* Forward declarations.
 */

/* Includes.
 */
#include "Structure.h"

/* Definitions.
 */
#define SCWRL_MAX_TIME 10 // Seconds
#define SCWRL_CHECK_DELAY 100000 // Microseconds
//#define SCWRL_EXECUTABLE "scwrl4"
//#define SCWRL_OPTIONS "-h" // Turn off output of hydrogens
#define SCWRL_EXECUTABLE "scwrl3"
#define SCWRL_OPTIONS "" // SCWRL3 doesn't need options

/* Interface.
 */
class AllAtomStructure : public Structure{
	public:
		/* Constructors and destructors.
		 */
		AllAtomStructure();
		AllAtomStructure(vector<Residue> residues);
		virtual ~AllAtomStructure();
		
		virtual shared_ptr<Structure> getBackbone() const;

		/* Read a PDB structure file from disk or as a string.
		 */
		static AllAtomStructure parseScwrlOutput(ifstream &inputFile);
		static AllAtomStructure newFromPDBString(string pdbLines); // RETURN A POINTER?

		// SCWRL-BASED IMPLEMENTATIONS HERE.		
		shared_ptr<Structure> threadSequence(string sequence);
		shared_ptr<Structure> adjustSidechains();
		vector< shared_ptr<Structure> > adjustComplex(vector< shared_ptr<Structure> > components);

	protected:
};

/* Idempotency end.
 */
#endif
