/* Program to test the speed and workings of a number of energy
 * functions. See Rastogi and Liberles, 2006, Biophys Chem 126 and
 * [some paper of my own] for a description of some of the methods.
 *
 * Usage: [TBA]
 */
#include "TwoBeadStructure.h"
#include "BastollaAugmentedModel.h"
#include "GrahnenModel.h"
#include "MHMCMC.h"
#include "ParameterRotation.h"
#include "randomgen/randomc.h"
#include "randomgen/stocc.h"
#include <fstream>
#include <sys/time.h>
#include <dirent.h>
#include <typeinfo>

using namespace std;

int main(int argc, char *argv[]){

	string usage = " [protein bead file] [ligand bead file] [# test sequences] [energy function (Bastolla++/Grahnen] [energy output file] [opt: term weightings file]";

	/* Usage message on bad input.
	 */
	if(argc < 6){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Read arguments.
	 */
	ifstream protBeadFile(argv[1]);
	ifstream ligBeadFile(argv[2]);
	int decoyNum = atoi(argv[3]);
	string energyFunction = argv[4];
	ofstream decoyEngFile(argv[5]);
	ifstream weightsFile;
	if(argc == 7){
		weightsFile.open(argv[6]);
	}

	/* Prep a pointer for the energy model.
	 */
	shared_ptr<EnergyModel> pEngModel;

	/* Start the random number generators.
	 */
	srand(time(NULL));
	CRandomMersenne mersenneRandGen(time(NULL));
	StochasticLib1 stochastGen(time(NULL));
	ParameterRotation::setPRNGs(mersenneRandGen, stochastGen);
	MHMCMC::setPRNG(mersenneRandGen);
	
	/* Decide which energy model to use. 	 
	 */
	if(energyFunction.compare("Bastolla++") == 0){
		pEngModel = shared_ptr<EnergyModel>(new BastollaAugmentedModel());
	}
	/* The Grahnen model assigns weights to each individual term. Employ
	 * user-supplied weights if they exist.
	 */
	else if(energyFunction.compare("Grahnen") == 0){
		if(weightsFile.is_open()){
			float bendWeight, ljWeight, helixWeight, betaWeight, chargeWeight, solvationWeight, cysWeight;
			float bindLjWeight, bindChargeWeight, bindSolvWeight;
			string dummy;
			while(getline(weightsFile, dummy)){
				if(dummy.compare(0,1,"#") != 0){
					break;
				}
			}
			
			bindLjWeight = atof(dummy.c_str());
			weightsFile >> bindChargeWeight;
			weightsFile >> bindSolvWeight;

			bendWeight = 1.0;
			ljWeight = 1.0;
			helixWeight = 1.0;
			betaWeight = 1.0;
			chargeWeight = 1.0;
			solvationWeight = 1.0;
			cysWeight = 1.0;

			/* DEBUG
			 */
			cout << "#Parameters: ";
			cout << bendWeight << ", ";
			cout << ljWeight << ", ";
			cout << helixWeight << ", ";
			cout << betaWeight << ", ";
			cout << chargeWeight << ", ";
			cout << solvationWeight << ", ";
			cout << cysWeight << ", ";
			cout << bindLjWeight << ", ";
			cout << bindChargeWeight << ", ";
			cout << bindSolvWeight << endl;

			pEngModel = shared_ptr<EnergyModel>(new GrahnenModel(bendWeight, ljWeight, helixWeight, betaWeight, chargeWeight, solvationWeight, cysWeight, bindLjWeight, bindChargeWeight, bindSolvWeight));
		}
		else{
			pEngModel = shared_ptr<EnergyModel>(new GrahnenModel()); // No user input, all weights equal
		}
	}
	else{
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	if(argc == 7){
		weightsFile.close();
	}
	  

	/* Define variables.
	 */
	shared_ptr<Structure> nativeStructure;
	shared_ptr<Structure> nativeLigand;
	shared_ptr<Structure> decoyLigand;

	float nativeBindEng, decoyBindEng;
	int decoysBetterThanNative = 0;

	string aminoAcids = "AEQDNLKSVRTPIMFYCWHG";

	int i;
	
	bool goodSeq;

	string decoyLigSeq;

	/* Read the native structures into memory.
	 */
	shared_ptr<TwoBeadStructure> nativeProteinStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(protBeadFile)));
	shared_ptr<TwoBeadStructure> nativeLigandStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(ligBeadFile)));
	protBeadFile.close();
	ligBeadFile.close();

	/* Pre-calculate exposures.
	*/
	nativeProteinStructure->calcExposure(false);
	nativeLigandStructure->calcExposure(false);

	vector< shared_ptr<TwoBeadStructure> > wholeComplex;
	wholeComplex.push_back(nativeProteinStructure);
	wholeComplex.push_back(nativeLigandStructure);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	nativeStructure = wholeComplex[0];
	nativeLigand = wholeComplex[1];

	/* Get the binding energy of the original.
	 */
	nativeBindEng = pEngModel->calcInteractionScore(nativeStructure, nativeLigand);

	/* Print native structure scores to STDOUT.
	 */
	cout << "#Native structure sequence: " << nativeStructure->getSequence() << endl;
	cout << "#Native ligand sequence: " << nativeLigand->getSequence() << endl;
	cout << "#Native structure binding score: " << nativeBindEng << endl;

	/* Calculate folding energies of a number of decoy sequences
	 * threaded through the native structure.
	 */
	for(i = 0; i < decoyNum; i++){

		/* DEBUG
		 *
		cout << "Generating decoy sequence..." << endl;
		*/

		/* Generate a random protein sequence for the ligand and
		 * thread it.
		 */
		goodSeq = false;
		while(!goodSeq){

			/* Generate a random protein sequence of the correct length.
			*/
			decoyLigSeq = generateRandomSequence(nativeLigand->getNumResidues(), aminoAcids, mersenneRandGen);

			/* Use a threading method to make necessary side chain adjustments
			 * to the structure.
			 */
			decoyLigand = nativeLigand->threadSequence(decoyLigSeq);
			
			/* It didn't have any problems, accept it.
			 */
			if(decoyLigand){
				goodSeq = true;
			}
		}

		/* Get the calculated binding energy.
		 */
		decoyBindEng = pEngModel->calcInteractionScore(nativeStructure, decoyLigand);

		/* Save results to file and variable.
		 */
		decoyEngFile << decoyBindEng;
		decoyEngFile << "\t#" << decoyLigSeq << endl;

		/* Is the binding score in the native structure lower than that of
		 * the native sequence? Then count it.
		 */
		if(decoyBindEng < nativeBindEng){
			decoysBetterThanNative++;
		}
	}

	/* Print results to STDOUT.
	 */
	cout << "#% decoy sequences with better score than native: " << ((float) decoysBetterThanNative / (float) decoyNum)*100 << endl;
	cout << "#See " << (string) argv[5] << " for full distribution." << endl;

	/* Finish.
	 */
	decoyEngFile.close();
	exit(0);
}
