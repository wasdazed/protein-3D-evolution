/* Program to test the speed and workings of a number of energy
 * functions. See Rastogi and Liberles, 2006, Biophys Chem 126 and
 * Grahnen et al, 2011, BMC Evol Biol 11:361 for a description of the methods.
 *
 * Usage: [TBA]
 */
#include "TwoBeadStructure.h"
#include "BastollaAugmentedModel.h"
#include "MukherjeeModel.h"
#include "GrahnenModel.h"
#include "MHMCMC.h"
#include "ParameterRotation.h"
#include "randomgen/randomc.h"
#include "randomgen/stocc.h"
#include <fstream>
#include <sys/time.h>

using namespace std;

int main(int argc, char *argv[]){

	string usage = " [protein bead file] [# test sequences] [energy function (Mukherjee/Grahnen/Bastolla++] [energy output file] [opt: term weightings file]";

	/* Usage message on bad input.
	 */
	if(argc < 5){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Read arguments.
	 */
	ifstream protBeadFile(argv[1]);
	int decoyNum = atoi(argv[2]);
	string energyFunction = argv[3];
	ofstream decoyEngFile(argv[4]);
	ifstream weightsFile;
	if(argc == 6){
		weightsFile.open(argv[5]);
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
	else if(energyFunction.compare("Mukherjee") == 0){
		pEngModel = shared_ptr<EnergyModel>(new MukherjeeModel());
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
			bendWeight = atof(dummy.c_str());
			weightsFile >> ljWeight;
			weightsFile >> helixWeight;
			weightsFile >> betaWeight;
			weightsFile >> chargeWeight;
			weightsFile >> solvationWeight;
			weightsFile >> cysWeight;

			bindLjWeight = 1.0;
			bindChargeWeight = 1.0;
			bindSolvWeight = 1.0;

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

	if(weightsFile.is_open()){
		weightsFile.close();
	}
	  

	/* Define variables.
	 */
	int i;
	int decoysBetterThanNative = 0;
	
	float nativeFoldEng, decoyFoldEng;
	
	bool goodSeq;

	string nativeProteinSeq, decoyProteinSeq;
	string aminoAcids = "AEQDNLKSVRTPIMFYCWHG";

	vector<float> sampledComponents;

	shared_ptr<Structure> proteinStructure;
	shared_ptr<Structure> decoyStructure;
	shared_ptr<Structure> invariant;

	/* Read the native structure into memory.
	 */
	shared_ptr<TwoBeadStructure> nativeStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(protBeadFile)));
	nativeProteinSeq = nativeStructure->getSequence();
	protBeadFile.close();

	/* Pre-calculate exposures.
	 */
	nativeStructure->calcExposure(false);

	/* Have a reference geometry on hand to generate random
	 * structure samples from.
	 */
	shared_ptr<TwoBeadStructure> invariantStructure(new TwoBeadStructure(*nativeStructure));
	invariantStructure->calcAndFreezeGeometries(4.5);
	proteinStructure = nativeStructure;
	invariant = invariantStructure;

	/* Get the folding energy of the original.
	 */
	nativeFoldEng = pEngModel->calcFoldEnergyGap(proteinStructure, invariant);

	/* Print native structure scores to STDOUT.
	 */
	cout << "#Native sequence: " << proteinStructure->getSequence() << endl;
	cout << "#Native structure folding score: " << nativeFoldEng << endl;

	/* Calculate folding energies of a number of decoy sequences
	 * threaded through the native structure.
	 */
	for(i = 0; i < decoyNum; i++){

		/* DEBUG
		 *
		cout << "Generating decoy sequence..." << endl;
		*/

		/* Generate a random protein sequence of suitable length
		 * and thread it.
		 */
		goodSeq = false;
		while(!goodSeq){

			/* Generate a random protein sequence of the correct length.
			*/
			decoyProteinSeq = generateRandomSequence(proteinStructure->getNumResidues(), aminoAcids, mersenneRandGen);

			/* Thread it.
			 */
			decoyStructure = proteinStructure->threadSequence(decoyProteinSeq);
			
			/* It didn't have any problems, accept it.
			 */
			//if(decoyStructure != 0){
			if(decoyStructure){
				goodSeq = true;
			}
		}

		/* Get the calculated folding energy (is the
		 * fitness measure under this model).
		 */
		decoyFoldEng = pEngModel->calcFoldEnergyGap(decoyStructure, invariant);

		/* Save results to file and variable.
		 */
		decoyEngFile << decoyFoldEng;
		decoyEngFile << "\t#" << decoyProteinSeq << endl;

		/* Is the fold score in the native structure lower than that of
		 * the native sequence? Then count it.
		 */
		if(decoyFoldEng < nativeFoldEng){
			decoysBetterThanNative++;
		}
	}

	/* Print results to STDOUT.
	 */
	cout << "#% decoy sequences with better score than native: " << ((float) decoysBetterThanNative / (float) decoyNum)*100 << endl;
	cout << "#See " << (string) argv[4] << " for full distribution." << endl;

	/* Finish.
	 */
	decoyEngFile.close();
	exit(0);
}
