#include "EnergyFuncProposal.h"
#include "Parameter.h"
#include "myUtils.h"
#include <cmath>
#include <iostream>
#include <algorithm>

/* Default constructor, does mostly nothing.
 */
EnergyFuncProposal::EnergyFuncProposal():
	Proposal()
{
}

/* Constructor from parameters, calls upward a bit.
 */
EnergyFuncProposal::EnergyFuncProposal(vector< shared_ptr<Parameter> > startParams, vector<float> nativeSeqTermValues, vector< vector<float> > decoyTermValues, vector< vector<float> > randomSeqTermValues):
	Proposal(startParams)
{
	/* Set some data members.
	 */
	randSeqTermVals = randomSeqTermValues;
	natSeqTermVals = nativeSeqTermValues;
	decoyTermVals = decoyTermValues;

	probDensity = calcProbDensity();
}

/* Virtual destructor, required by gcc.
 */
EnergyFuncProposal::~EnergyFuncProposal(){
}

/* Propose a new state from the current one
 * based on some size of the steps.
 */
shared_ptr<Proposal> EnergyFuncProposal::nextStep(float stepSize){
	unsigned int i;
	vector< shared_ptr<Parameter> > newParams;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->nextValue(stepSize));
	}

	return shared_ptr<Proposal>(new EnergyFuncProposal(newParams, natSeqTermVals, decoyTermVals, randSeqTermVals));
}

/* Propose a new state from the current one with
 * step sizes varying with some variance.
 */
shared_ptr<Proposal> EnergyFuncProposal::nextStep(float stepSize, float stepVar){
	unsigned int i;
	vector< shared_ptr<Parameter> > newParams;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->nextValue(stepSize, stepVar));
	}

	return shared_ptr<Proposal>(new EnergyFuncProposal(newParams, natSeqTermVals, decoyTermVals, randSeqTermVals));

}

/* Propose a completely random state based on a scaling
 * of the parameter space.
 */
shared_ptr<Proposal> EnergyFuncProposal::randomStep(float scaleFactor){
	unsigned int i;
	vector< shared_ptr<Parameter> > newParams;
	
	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->randomValue(scaleFactor));
	}

	return shared_ptr<Proposal>(new EnergyFuncProposal(newParams, natSeqTermVals, decoyTermVals, randSeqTermVals));
}

/* Function for calculating the probability density, BY SOME
 * METRIC THAT'S STILL UNDER DEVELOPMENT.
 *
 * Approximated by drawing samples from pre-computed random
 * sequences at each point and comparing to the native sequence.
 */
float EnergyFuncProposal::calcProbDensity(){
	unsigned int i, j, middle;
	float nativeEng = 0.0; 
	float density, dispersal, randomEng;
	vector<float> parameters;
	vector<float> randomEngs;
	vector<float> decoyEngs( decoyTermVals.size(), 0.0);
	vector<float> zeroes( decoyTermVals.size(), 0.0);
	vector<float> quarts;

	/* Convert current parameter set to floating point representation.
	 */
	for(i = 0; i < params.size(); i++){
		parameters.push_back( params[i]->toFloat() );
	}

	/* Calculate the native sequence energy under these
	 * parameters.
	 */
	for(i = 0; i < parameters.size(); i++){
		nativeEng += parameters[i]*natSeqTermVals[i];
	}

	/* And the respective decoy energies.
	 */
	for(i = 0; i < decoyTermVals.size(); i++){
		for(j = 0; j < parameters.size(); j++){
			decoyEngs[i] += parameters[i]*decoyTermVals[i][j];
		}
	}

	/* Calculate the scaled version of the energy function for
	 * all random sequences. Note that random energy values are
	 * organized by column (one for each term in the function),
	 * not by row (one for each sequence).
	 */
	for(i = 0; i < randSeqTermVals[0].size(); i++){
		
		randomEng = 0.0;

		/* Get the scaled energy value for this random sequence.
		 */
		for(j = 0; j < parameters.size(); j++){
			randomEng += parameters[j]*randSeqTermVals[j][i];
		}

		randomEngs.push_back(randomEng);
	}

	/* Calculate the quartiles of the distribution of random
	 * energies.
	 */
	quarts = quartiles(randomEngs);
	dispersal = quarts[2] - quarts[1];

	/* Calculate the enegy gap to the nearest explicit decoy
	 * conformation, or the energy gap to the middle of the distribution
	 * of random conformation.
	 */
	if(decoyEngs != zeroes){

		sort(decoyEngs.begin(), decoyEngs.end());

		/* Case 1: Best decoy is between native and random states,
		 * measure against it.
		 */
		if(randomEngs[0] > nativeEng && decoyEngs[0] > nativeEng && decoyEngs[0] < randomEngs[0]){
			/* DEBUG
			 *
			cerr << "A decoy is the closest misfold." << endl;
			*/

			density = (decoyEngs[0] - nativeEng) / dispersal;
		}
		/* Case 2: Best decoy is worse than best random state: measure 
		 * against random.
		 */
		else if(randomEngs[0] > nativeEng && decoyEngs[0] > randomEngs[0]){

			/* DEBUG
			 *
			cerr << "A random sequence is the closest misfold." << endl;
			*/

			density = (randomEngs[0] - nativeEng) / dispersal;
		}
		/* Either the best decoy or best random is better than native:
		 * not cool, sign should switch to negative to get out of this
		 * situation.
		 */
		else if(decoyEngs[0] < nativeEng || randomEngs[0] < nativeEng){

			/* DEBUG
			 *
			cerr << "DANGER: native not the most stable state, ";
			*/
			
			if(decoyEngs[0] < randomEngs[0]){

				/* DEBUG
				 *
				cout << "a decoy is!" << endl;
				*/

				density = (decoyEngs[0] - nativeEng) / dispersal;
			}
			else{
				/* DEBUG
				 *
				cout << "a random sequence is!" << endl;
				*/

				density = (randomEngs[0] - nativeEng) / dispersal;
			}
		}
	}
	else{
		density = -1.0*( (nativeEng - quarts[2]) / (quarts[2]-quarts[1]) );
	}

	return density;
}
