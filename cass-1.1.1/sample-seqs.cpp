/* Program to get the sample the distribution of folding and binding 
 * energies in the portion of sequence space around a starting sequence.
 * See Rastogi and Liberles, 2006, Biophys Chem 126 and
 * [some paper of my own] for a description of the methods.
 *
 * Usage: [TBA]
 */

#include "GrahnenModel.h"
#include "BastollaAugmentedModel.h"
#include "TwoBeadStructure.h"
#include "MHMCMC.h"
#include "SequenceProposal.h"
#include "GoSequenceProposal.h"
#include "ParameterSequence.h"
#include "ParameterRotation.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]){

	string usage = " [protein bead file] [ligand bead file] [starting sequence] [energy function (folding/binding] [# of samples] [temperature] [max % divergence] [energy model (Grahnen/Bastolla++/Go)] [opt: parameter set file]";

	/* Usage message on bad input.
	 */
	if(argc < 9){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Read arguments.
	 */
	ifstream protFile(argv[1]);
	ifstream ligFile(argv[2]);
	ifstream startSeqFile(argv[3]);
	string energyFunction = argv[4];
	unsigned int numSamples = atoi(argv[5]);
	float temp = atof(argv[6]);
	float divergence = atof(argv[7]);
	string energyModel = argv[8];
	ifstream paramFile;

	if(argc == 10){
		paramFile.open(argv[9]);
	}

	/* Error checking on input.
	 */
	if(     (energyFunction.compare("folding") != 0 && energyFunction.compare("binding") != 0)
		||
		(energyModel.compare("Grahnen") != 0 && energyModel.compare("Bastolla++") != 0 && energyModel.compare("Go") != 0)
	  ){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Define variables.
	 */
	string inputSeq;

	unsigned int i;

	/* Dummy variables: step is always a single
	 * mutation.
	 */
	float meanStep = 1;
	float stepVariance = 1;
	
	shared_ptr<Structure> proteinStructure;
	shared_ptr<Structure> ligandStructure;
	shared_ptr<Structure> invariantStructure;

	shared_ptr<Proposal> pCurrent;
	shared_ptr<Proposal> pBest;

	shared_ptr<EnergyModel> engModel;

	vector< shared_ptr<Parameter> > startParams;
	
	/* Seed the random number generators.
	 */
	CRandomMersenne mersenneRandGen(time(NULL));
	StochasticLib1 stochastRandGen(time(NULL));
	MHMCMC::setPRNG(mersenneRandGen);
	ParameterRotation::setPRNGs(mersenneRandGen, stochastRandGen);
	srand(time(NULL));

	/* Initialize the energy model. The Go model uses BastollaAugmented.
	 * If no parameter values were supplied for Grahnen, use the
	 * hardcoded ones.
	 */
	if(energyModel.compare("Grahnen") == 0){
		if(!paramFile.is_open()){
			cout << "#No parameters supplied, using hardcoded ones..." << endl;
			engModel = shared_ptr<EnergyModel>(new GrahnenModel(0.213649, 0.255238, 0.0351669, 0.193426, 0.00636471, 0.296156, 0.0, 0.970, 0.011, 0.019)); //20110706
		}
		else{
			float bend, lj, helix, beta, charge, solv, ssbond, bindLj, bindCharge, bindSolv;
			if(energyFunction.compare("folding") == 0){
				paramFile >> bend;
				paramFile >> lj;
				paramFile >> helix;
				paramFile >> beta;
				paramFile >> charge;
				paramFile >> solv;
				paramFile >> ssbond;
				bindLj = 0.0;
				bindCharge = 0.0;
				bindSolv = 0.0;
			}
			else{
				bend = 0.0;
				lj = 0.0;
				helix = 0.0;
				beta = 0.0;
				charge = 0.0;
				solv = 0.0;
				ssbond = 0.0;
				paramFile >> bindLj;
				paramFile >> bindCharge;
				paramFile >> bindSolv;
			}
			cout << "#Parameters used: " << bend << "\t" << lj << "\t" << helix << "\t" << beta << "\t" << charge << "\t" << solv << "\t" << ssbond << "\t" << bindLj << "\t" << bindCharge << "\t" << bindSolv << endl;
			engModel = shared_ptr<EnergyModel>(new GrahnenModel(bend, lj, helix, beta, charge, solv, ssbond, bindLj, bindCharge, bindSolv));
			paramFile.close();
		}
	}
	else if(energyModel.compare("Bastolla++") == 0){
		engModel = shared_ptr<EnergyModel>(new BastollaAugmentedModel());
	}

	/* Read protein, ligands, sequence into memory.
	 */
	shared_ptr<TwoBeadStructure> proteinStruct = shared_ptr<TwoBeadStructure>(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(protFile)));
	shared_ptr<TwoBeadStructure> ligandStruct = shared_ptr<TwoBeadStructure>(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(ligFile)));
	startSeqFile >> inputSeq;
	protFile.close();
	ligFile.close();
	startSeqFile.close();

	/* Pre-calculate exposures as necessary.
	 */
	proteinStruct->calcExposure(false);
	ligandStruct->calcExposure(false);

	vector< shared_ptr<TwoBeadStructure> > wholeComplex;
	wholeComplex.push_back(proteinStruct);
	wholeComplex.push_back(ligandStruct);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	proteinStructure = wholeComplex[0];
	ligandStructure = wholeComplex[1];

	/* Make a copy of the initial conformation and pre-calculate the
	 * geometric measurements.
	 */
	shared_ptr<TwoBeadStructure> sourceOfGeometry(new TwoBeadStructure(*dynamic_pointer_cast<TwoBeadStructure>(proteinStructure)));
	sourceOfGeometry->calcAndFreezeGeometries(4.5);
	invariantStructure = sourceOfGeometry;

	cout << "#Calculating initial probability density..." << endl;

	/* Initialize the Markov chain.
	 */
	shared_ptr<Parameter> startingSeq(new ParameterSequence(inputSeq, inputSeq, divergence));
	startParams.push_back(startingSeq);
	shared_ptr<Proposal> pStartState;
	if(energyFunction.compare("folding") == 0 && (energyModel.compare("Grahnen") == 0 || energyModel.compare("Bastolla++") == 0)){
		pStartState = shared_ptr<Proposal>(new SequenceProposal(startParams, proteinStructure, engModel, invariantStructure));
	}
	else if(energyFunction.compare("binding") == 0 && (energyModel.compare("Grahnen") == 0 || energyModel.compare("Bastolla++") == 0)){
		pStartState = shared_ptr<Proposal>(new SequenceProposal(startParams, proteinStructure, ligandStructure, engModel));
	}
	else if(energyFunction.compare("folding") == 0 && energyModel.compare("Go") == 0){
		/* Calculate residue-wise contributions and make sure they're always
		 * attractive in the native state.
		 */
		BastollaAugmentedModel tempModel;
		vector<float> nativeEngs = tempModel.residueFoldContributions(proteinStructure);
		float largest = *max_element(nativeEngs.begin(), nativeEngs.end());
		for(i = 0; i < nativeEngs.size(); i++){
			nativeEngs[i] -= largest;
		}
		pStartState = shared_ptr<Proposal>(new GoSequenceProposal(startParams, proteinStructure->getSequence(), nativeEngs));
	}
	else if(energyFunction.compare("binding") == 0 && energyModel.compare("Go") == 0){
		cerr << "FATAL ERROR: Go model does not support sampling binding energies!" << endl;
		exit(1);
	}
	MHMCMC oneChain(pStartState, temp, numSamples);

	/* Report starting conditions to user.
	 */
	cout << "#MCMC intialized." << endl;
	cout << "#Steps: " << numSamples << endl;
	cout << "#Starting temperature: " << temp << endl;
	cout << "#Initial state (probability density + parameters): " << pStartState->toString() << endl;
	cout << "Step\tEnergy\tSeq" << endl;

	/* Run the chain for as long as necessary.
	 */
	for(i = 0; i < numSamples; i++){

		oneChain.takeStep(meanStep, stepVariance);

		/* Report current state.
		 */
		pCurrent = oneChain.getCurrentState();
		cout << (i+1) << "\t" << pCurrent->toString() << endl;
	}

	/* Output most stable state to STDOUT as well.
	 */
	pBest = oneChain.getBestState();
	cout << "#Best found state:" << endl;
	cout << "#" << pBest->toString() << endl;

	/* Finish.
	*/
	exit(0);
}
