/* Program to test calculate values of individual terms of
 * energy functions for a native sequence, a collection of
 * random sequences, and a collection of randomly generated
 * "structures" (samples from geometry distributions, really).
 */

#include "GrahnenModel.h"
#include "AllAtomStructure.h"
#include "TwoBeadStructure.h"
#include "ParameterRotation.h"
#include "MHMCMC.h"
#include <fstream>
#include <limits>

using namespace std;

int main(int argc, char *argv[]){

	string usage = " [protein bead file] [ligand bead file] [# random samples] [energy function (folding/binding)]";

	/* Usage message on bad input.
	 */
	if(argc < 4){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Read arguments.
	 */
	ifstream proteinBeadFile(argv[1]);
	ifstream ligandBeadFile(argv[2]);
	unsigned int sampleNum = atoi(argv[3]);
	string energyFunction = argv[4];

	/* Error checking on input.
	 */
	if(energyFunction.compare("folding") != 0
	   && energyFunction.compare("binding") != 0
	  ){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Define variables.
	 */
	string newSequence;
	string aminoAcids = "AEQDNLGKSVRTPIMFYCWH";

	shared_ptr<Structure> randomStructure;

	vector<float> engTerms;
	vector<float> nativeEngTerms;
	vector<float> useTerms;

	unsigned int i, j; 
	unsigned int goodSeqs;
	
	/* Read the protein structure into memory.
	 */
	shared_ptr<TwoBeadStructure> nativeStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(proteinBeadFile)));
	proteinBeadFile.close();

	/* Do the same for the ligand.
	 */
	shared_ptr<TwoBeadStructure> ligandStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(ligandBeadFile)));
	ligandBeadFile.close();

	/* Pre-calculate exposure values.
	 */
	nativeStructure->calcExposure(false);
	ligandStructure->calcExposure(false);
	
	vector< shared_ptr<TwoBeadStructure> > wholeComplex;
	wholeComplex.push_back(nativeStructure);
	wholeComplex.push_back(ligandStructure);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	
	nativeStructure = wholeComplex[0];
	ligandStructure = wholeComplex[1];

	/* Seed and set the random number generators.
	 */
	CRandomMersenne mersenneRandGen(time(NULL));
	StochasticLib1 stochastRandGen(time(NULL));
	srand(time(NULL));
	ParameterRotation::setPRNGs(mersenneRandGen, stochastRandGen);
	MHMCMC::setPRNG(mersenneRandGen);

	/* Instantiate the energy model with default parameters.
	 */
	shared_ptr<GrahnenModel> pEngModel(new GrahnenModel());

	/* Get the unscaled terms from the energy function at
	 * hand for the native sequence.
	 *
	 * Also set which terms to actually use (1.0 for using
	 * a term, 0.0 for removing it).
	 */
	if(energyFunction.compare("folding") == 0){
		nativeEngTerms = pEngModel->getFoldScoreTerms(nativeStructure);
		useTerms = nativeEngTerms;
		useTerms[0] = 1.0; // Bond angle (bend)
		useTerms[1] = 1.0; // LJ
		useTerms[2] = 1.0; // Helix
		useTerms[3] = 1.0; // Beta
		useTerms[4] = 1.0; // Charge
		useTerms[5] = 1.0; // Solvation
		useTerms[6] = 1.0; // Disulfide
	}
	else if(energyFunction.compare("binding") == 0){
		nativeEngTerms = pEngModel->getBindScoreTerms(nativeStructure, ligandStructure);
		useTerms = nativeEngTerms;
		useTerms[0] = 1.0; // LJ
		useTerms[1] = 1.0; // Charge
		useTerms[2] = 1.0; // Solvation change
	}

	/* Zero-out terms as necessary.
	 */
	for(i = 0; i < nativeEngTerms.size(); i++){
		nativeEngTerms[i] *= useTerms[i];
	}

	/* Prep output.
	 */
	cout << "#Native protein sequence: " << nativeStructure->getSequence() << endl;
	cout << "#Original ligand sequence: " << ligandStructure->getSequence() << endl;
	cout << "#Native energy term values:" << endl;
	for(i = 0; i < nativeEngTerms.size()-1; i++){
		cout << nativeEngTerms[i] << "\t";
	}
	cout << nativeEngTerms[i] << endl;
	cout << "#Energy term values for " << sampleNum << " random sequences:" << endl;

	/* Calculate binding energies of a number of random sequences
	 * threaded through the native protein/ligand structure.
	 */
	goodSeqs = 0;
	while(goodSeqs < sampleNum){

		/* DEBUG
		 *
		cout << "Generating decoy sequence..." << endl;
		*/
		
		/* Generate a shuffled sequence of the appropriate length
		 * and get the term values for it.
		 */
		if(energyFunction.compare("folding") == 0){
			//newSequence = shuffleSequence(nativeStructure->getSequence());
			newSequence = generateRandomSequence(nativeStructure->getNumResidues(), aminoAcids, mersenneRandGen);
			randomStructure = nativeStructure->threadSequence(newSequence);

			/* If there's a threading failure, start over with a new sequence.
			 */
			if(randomStructure.get() == 0){
				continue;
			}
			else{
				goodSeqs++;
				engTerms = pEngModel->getFoldScoreTerms(randomStructure);
			}
		}
		else if(energyFunction.compare("binding") == 0){
			//newSequence = shuffleSequence(ligandStructure->getSequence());
			newSequence = generateRandomSequence(ligandStructure->getNumResidues(), aminoAcids, mersenneRandGen);
			randomStructure = ligandStructure->threadSequence(newSequence);
	
			/* If there's a threading failure, start over with a new sequence.
			*/
			if(randomStructure.get() == 0){
				continue;
			}
			else{
				goodSeqs++;
				engTerms = pEngModel->getBindScoreTerms(nativeStructure, randomStructure);
			}
		}

		/* Zero-out terms as necessary.
		 */
		for(j = 0; j < engTerms.size(); j++){
			engTerms[j] *= useTerms[j];
		}

		/* Print results to STDOUT.
		 */
		for(j = 0; j < engTerms.size()-1; j++){
			cout << engTerms[j] << "\t";
		}
		cout << engTerms[j] << "\t" << newSequence << endl;
	}

	/* That's it if we're parameterizing binding.
	 */
	if(energyFunction.compare("binding") == 0){
		exit(0);
	}

	/* Sample random conformations via the Random Energy Model.
	 */
	cout << "#Energy term values for " << sampleNum << " random structures:" << endl;

	/* First, extract the relevant geometric measurements.
	 */
	string sequence = nativeStructure->getSequence();
	vector< pair<float, float> > bondAngles = nativeStructure->bondAngles();
	vector<float> resBBDists = nativeStructure->caCbDistances();
	vector<float> resDists = nativeStructure->cbCbDistances();
	vector< pair<float, float> > backboneDists = nativeStructure->backboneDistances();
	vector<float> backboneTorsions = nativeStructure->backboneTorsions();
	vector<float> exposures = nativeStructure->residueExposures();
	vector<string> threeMers = nativeStructure->threeMers();
	vector<string> helixFourMers = nativeStructure->helixFourMers();
	vector<string> betaFourMers = nativeStructure->betaFourMers();
	srand(time(NULL));
	/* Sample random permutations of the geometry,
	 * representing random compact structures.
	 */
	for(i = 0; i < sampleNum; i++){
		engTerms = pEngModel->sampleFoldGapTerms(sequence, bondAngles, resBBDists, resDists, threeMers, helixFourMers, backboneDists, betaFourMers, backboneTorsions, exposures, mersenneRandGen);

		/* Zero-out terms as necessary.
		 */
		for(j = 0; j < engTerms.size(); j++){
			engTerms[j] *= useTerms[j];
		}

		/* Print results to STDOUT.
		 */
		for(j = 0; j < engTerms.size()-1; j++){
			cout << engTerms[j] << "\t";
		}
		cout << engTerms[j] << endl;
	}

	/* Finish.
	 */
	exit(0);
}
