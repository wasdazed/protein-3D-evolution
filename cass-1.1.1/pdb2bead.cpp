/* Utitility program to convert an all-atom PDB representation
 * into a 2-bead coarse-grained representation as in
 * Mukherjee and Bagchi, 2003, J Chem Phys 118(10) and/or
 *
 *
 * Run SCWRL 3.0 (or 4.0 with the "-h" switch) on your structure first to 
 * optimize and simplify parsing.
 *
 * Usage: ./pdb2bead [input PDB file] [output file]
 */

#include <iostream>
#include "AllAtomStructure.h"
#include "TwoBeadStructure.h"

using namespace std;

int main(int argc, char* argv[]){

	/* Check arguments.
	 */
	if(argc < 3){
		cout << "Usage: ./pdb2bead [input PDB file] [output file]" << endl;
		exit(1);
	}

	/* Open files.
	 */
	ifstream pdbFile(argv[1]);
	ofstream outBead(argv[2]);

	/* Convert
	 */
	shared_ptr<AllAtomStructure> original(new AllAtomStructure(AllAtomStructure::parseScwrlOutput(pdbFile)));
	shared_ptr<TwoBeadStructure> bead(new TwoBeadStructure(original));

	/* Output the bead representation.
	 */
	outBead << bead->toPDBString();

	pdbFile.close();
	outBead.close();

	return 0;
}
