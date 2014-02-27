/* Program to run an MCMC-based parameter optimization
 * of a set of floating-point parameters that govern the
 * behavior of an energy function.
 *
 * Relies heavily on pre-computed distributions of parts of
 * the function to do its job.
 */

/* Includes.
 */
#include "myUtils.h"
#include "MHMCMC.h"
#include "ParameterFloat.h"
#include "EnergyFuncProposal.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, const char *argv[]){
	string usage = " [parameter start file] [term distribution file] [energy function (folding/binding)] [# steps] [temp] [step length] [step variance]";


	/* Check arguments.
	 */
	if(argc < 8){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}
	
	/* Read arguments.
	 */
	ifstream paramFile(argv[1]);
	ifstream distribFile(argv[2]);
	string engFunc = argv[3];
	unsigned int numSteps = atoi(argv[4]);
	float startTemp = atof(argv[5]);
	float meanStep = atof(argv[6]);
	float stepVariance = atof(argv[7]);

	/* Define some variables.
	 */
	unsigned int i,j;
	float firstValue, paramSum;
	string line;

	vector< shared_ptr<Parameter> > parameters;
	vector< vector<float> > termDistributions;
	vector<float> nativeTerms;
	vector< vector<float> > dummy;

	vector<string> tokens;
	vector<string> sequences;
	
	shared_ptr<Parameter> param;
	shared_ptr<Proposal> pCurrent;
	shared_ptr<Proposal> pBest;

	vector<float> dummyFloats;
	vector<float> scaleFactors;
	vector<float> bestParams;

	/* Start random number generators.
	 */
	CRandomMersenne mersenneRandGen(time(NULL));
	StochasticLib1 stochastRandGen(time(NULL));
	srand(time(NULL));
	MHMCMC::setPRNG(mersenneRandGen);
	ParameterFloat::setPRNGs(mersenneRandGen, stochastRandGen);

	/* Read in starting parameter set.
	 */
	while( getline(paramFile, line) ){
		param = shared_ptr<Parameter>(new ParameterFloat(atof(line.c_str()))); 
		parameters.push_back(param);
	}
	paramFile.close();

	/* DEBUG
	 */
	cout << "#Reading energy distributions from file..." << endl;

	/* Read in the distributions of terms in the energy function.
	*/
	int lineCount = 0;
	int numDims = 0;
	while( getline(distribFile, line) ){

		/* Don't read lines starting with "#" or empty ones.
		*/
		if(line.substr(0,1).compare("#") != 0 && line.length() > 1){

			lineCount++;

			/* Make tokens.
			*/
			tokens = tokenize(line, "\t");

			/* Initialize vector and assign native sequence values.
			*/
			if(lineCount == 1){
				for(i = 0; i < tokens.size(); i++){
					termDistributions.push_back(dummyFloats);
				}
				
				/* Convert to floats and save.
				*/
				for(i = 0; i < tokens.size(); i++){
					nativeTerms.push_back(atof(tokens[i].c_str()));
					numDims++;
				}
			}
			/* The remaining data belongs to the random samples.
			 */
			else{
				/* Convert to floats and save.
				*/
				for(i = 0; i < numDims; i++){
					termDistributions[i].push_back(atof(tokens[i].c_str()));
				}
			}
		}

	}
	distribFile.close();

	/* Check for valid function to parameterize and
	 * make sure there are a sufficient number of
	 * parameters and previously known distributions.
	 */
	if(
	  (
		  engFunc.compare("folding") != 0
		  && engFunc.compare("binding") != 0 
	  )
	  ||
	  (
		  (engFunc.compare("folding") == 0 && parameters.size() != 7)
		  ||
		  (engFunc.compare("binding") == 0 && parameters.size() != 3)
	  )
	  ||
	  (
		  (engFunc.compare("folding") == 0 && termDistributions.size() != 7)
		  ||
		  (engFunc.compare("binding") == 0 && termDistributions.size() != 3)
	  )
	  ){
		cout << "Usage: " << argv[0] << usage << endl;
		cout << "Folding: 7 parameters, binding: 3 parameters." << endl;
		exit(1);
	}
	
	/* DEBUG
	 */
	cout << "#Range of the distributions from disk:" << endl;
	for(i = 0; i < termDistributions.size(); i++){
		cout << "#" << vectorRange(termDistributions[i]) << endl;
	}

	cout << "#Removing non-informative dimensions..." << endl;

	/* Special case: one or more parameters are information-less in the
	 * sense that they have no entropy (e.g. contain the same value for 
	 * every sample). We therefore ignore them during optimization.
	 */
	vector<bool> nonInformative(numDims, true);
	for(i = 0; i < termDistributions.size(); i++){
		for(j = 0; j < termDistributions[i].size(); j++){
			if(j == 0){
				firstValue = termDistributions[i][j];
			}
			else if(termDistributions[i][j] != firstValue){
				nonInformative[i] = false;
				break;
			}
		}
	}
	for(i = nonInformative.size(); i > 0; i--){
		if(nonInformative[i-1]){
			cout << "#Dimension " << i << " contains no information, removing it." << endl;
			nativeTerms.erase(nativeTerms.begin()+(i-1));
			termDistributions.erase(termDistributions.begin()+(i-1));
			parameters.erase(parameters.begin()+(i-1));
		}
	}

	cout << "#Scaling remaining dimensions to equal range..." << endl;

	/* Find the term distribution with the largest range.
	 */
	float currentRange;
	float largestRange = 0.0;
	for(i = 0; i < termDistributions.size(); i++){
		currentRange = vectorRange(termDistributions[i]);
		if(currentRange > largestRange){
			largestRange = currentRange;
		}
	}

	/* Scale the remaining distributions, and the native values,
	 * to that range.
	 */
	float scaleFactor;
	for(i = 0; i < termDistributions.size(); i++){
		scaleFactor = largestRange / vectorRange(termDistributions[i]);

		/* Scale the native value.
		 */
		nativeTerms[i] *= scaleFactor;

		/* And the random sequences.
		 */
		for(j = 0; j < termDistributions[i].size(); j++){
			termDistributions[i][j] *= scaleFactor;
		}

		/* And save it.
		 */
		scaleFactors.push_back(scaleFactor);
	}

	/* DEBUG
	 */
	cout << "#Range of the distributions after scaling:" << endl;
	for(i = 0; i < termDistributions.size(); i++){
		cout << "#" << vectorRange(termDistributions[i]) << endl;
	}

	/* Report progress.
	 */
	cout << "#Calculating initial probability density..." << endl;

	/* Formulate a starting state based on input parameters. No decoys
	 * involved here, so just feed in an empty vector for those terms.
	 */
	shared_ptr<Proposal> pStartState(new EnergyFuncProposal(parameters, nativeTerms, dummy, termDistributions));

	/* Initialize MCMC (standard Metropolis-Hastings algorithm).
	 */
	MHMCMC oneChain(pStartState, startTemp, numSteps);

	/* Report starting conditions to user.
	 */
	cout << "#MCMC intialized." << endl;
	cout << "#Steps: " << numSteps << endl;
	cout << "#Step size (variance): " << meanStep << " (" << stepVariance << ")" << endl;
	cout << "#Starting temperature: " << startTemp << endl;
	cout << "#Initial state (probability density + parameters): " << pStartState->toString() << endl;
	cout << "#Step\tProb density\tParameters" << endl;

	/* Run MCMC for required number of steps.
	 *
	 */
	for(i = 0; i < numSteps; i++){
		oneChain.takeStep(meanStep, stepVariance);

		/* Report current state.
		 */
		pCurrent = oneChain.getCurrentState();
		cout << "#" << (i+1) << "\t" << pCurrent->toString() << endl;
	}
	
	/* Output final and best states to STDOUT.
	 */
	pCurrent = oneChain.getCurrentState();
	pBest = oneChain.getBestState();
	cout << "#Current state of the chain:" << endl;
	cout << "#" << pCurrent->toString() << endl;
	cout << "#Best sampled state:" << endl;
	cout << "#" << pBest->toString() << endl;

	/* Calculate and output final parameter values,
	 * re-scaled such that the sum of parameters
	 * is always 1.
	 */
	bestParams = pBest->getFloatParams();
	paramSum = 0.0;
	for(i = 0; i < bestParams.size(); i++){
		paramSum += bestParams[i];
	}
	for(i = 0; i < bestParams.size(); i++){
		bestParams[i] /= paramSum;
	}
	
	/* Re-insert removed dimensions by simply setting them
	 * to 0.
	 */
	for(i = 0; i < nonInformative.size(); i++){
		if(nonInformative[i]){
			bestParams.insert(bestParams.begin()+i, 0.0);
			cout << "#Re-inserted non-informative dimension " << i+1 << "." << endl;
		}
	}

	cout << "#Final scaling factors:" << endl;
	for(i = 0; i < bestParams.size(); i++){
		cout << bestParams[i] << endl;
	}

	/* Finish up.
	 */
	return 0;
}
