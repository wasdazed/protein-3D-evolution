/* Program to adjust side chain angles near the interface of a protein-protein 
 * complex such that the L-J interaction energy is minimized. Operates
 * on 2-bead structure representations.
 */

#include "TwoBeadStructure.h"
#include "MHMCMC.h"
#include "ParameterRotation.h"
#include "myUtils.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]){

	string usage = " [protein 1 bead file] [protein 2 bead file] [output file 1] [output file 2]";

	/* Usage message on bad input.
	 */
	if(argc < 5){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Read arguments.
	 */
	ifstream protOneFile(argv[1]);
	ifstream protTwoFile(argv[2]);
	ofstream outOneFile(argv[3]);
	ofstream outTwoFile(argv[4]);

	/* Define variables.
	 */

	/* Seed the random number generators.
	 */
	CRandomMersenne mersenneRandGen(time(NULL));
	StochasticLib1 stochastRandGen(time(NULL));
	MHMCMC::setPRNG(mersenneRandGen);
	ParameterRotation::setPRNGs(mersenneRandGen, stochastRandGen);
	srand(time(NULL));

	/* Read structures into memory.
	 */
	shared_ptr<TwoBeadStructure> partOne(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(protOneFile)));
	shared_ptr<TwoBeadStructure> partTwo(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(protTwoFile)));

	/* Close input files.
	*/
	protOneFile.close();
	protTwoFile.close();

	/* Adjust them as a complex.
	 */
	vector< shared_ptr<Structure> > starting;
	starting.push_back(partOne);
	starting.push_back(partTwo);
	vector< shared_ptr<Structure> > adjusted = partOne->adjustComplex(starting);

	/* Output to disk.
	 */
	outOneFile << adjusted[0]->toPDBString();
	outTwoFile << adjusted[1]->toPDBString();

	/* Finish.
	*/
	outOneFile.close();
	outTwoFile.close();
	exit(0);
}
