#include "MHMCMC.h"
#include "Proposal.h"
#include "randomgen/randomc.h"
#include <math.h>
#include <iostream>
#include <stdexcept>

/* Static shared PRNG.
 */
shared_ptr<CRandomMersenne> MHMCMC::pRandGen;

/* Default constructor, does nothing much.
 */
MHMCMC::MHMCMC():
	MCMC()
{
}

/* Constructor from parameters. Calls upward for pretty
 * much everything.
 */
MHMCMC::MHMCMC(shared_ptr<Proposal> startState, float startTemp, int numSteps):
	MCMC(startState, startTemp, numSteps)
{
	/* Make sure the class has a PRNG handy.
	 */
	if(MHMCMC::pRandGen.get() == 0){
		throw runtime_error("You must supply a random number generator via MHMCMC::setPRNG()!");
	}

}

/* Virtual destructor, required due to gcc finickyness.
 */
MHMCMC::~MHMCMC(){
	/* DEBUG
	 *
	cout << "Destroying MHMCMC class properties." << endl;
	*/
}

/* Function to evaluate moving to a proposed new state
 * (acceptance function). Based on the common acceptance
 * function
 * 
 * P(Pd(s), Pd(s'), T) = |---- 1, if Pd(s') > Pd(s)
 *                       |
 *		         |---- exp(-(Pd(s) - Pd(s'))/T), if Pd(s') < Pd(s)
 *
 * and a random number generator from the uniform
 * distribution.
 */
bool MHMCMC::moveToState(shared_ptr<Proposal> newState){

	/* DEBUG
	 *
	cout << "Calculating difference between proposed new state and current state." << endl;
	*/

	/* Get the differance in probabilty densities.
	 */
	float probDiff = currentState->getProbDensity() - newState->getProbDensity();
//cout << "ProbDiff: " << probDiff << endl;
	/* Definitely move if the new state is better.
	 */
	if(probDiff < 0){

		/* DEBUG
		 *
		cout << "Difference is " << probDiff << " (i.e. < 0), making the move." << endl;
		*/

		return true;
	}
	/* Maybe move if the new state is worse.
	 */
	else{
		/* DEBUG
		 *
		cout << "Difference is " << probDiff << ", might make the move." << endl;
		*/

		/* Get a random number between 0 and 1.
		 */
		float random = MHMCMC::pRandGen->Random();

		/* Accept the move with probability described above:
		 * moving probability decays quickly with increasing
		 * energy difference but is modulated by temperature.
		 */
		if( random < exp(-probDiff/temp) ){

			/* DEBUG
			 *
			cout << "Yes, moving to new state." << endl;
			*/

			return true;
		}
		else{
			/* DEBUG
			 *
			cout << "No, staying put." << endl;
			*/

			return false;
		}
	}
}

/* Function to take a single MCMC step with tunable
 * variance in step sizes.
 */
void MHMCMC::takeStep(float stepSize, float stepVar){

	/* Propose a new state.
	 */
	shared_ptr<Proposal> pNext = currentState->nextStep(stepSize, stepVar);

	/* Take a step if the acceptance function approves.
	 */
	if( moveToState(pNext) ){

		currentState = pNext;

		/* Is is a better state than any previously seen?
		 */
		if(currentState->getProbDensity() > bestState->getProbDensity() ){
			bestState = currentState;
		}
	}

	/* One more step taken.
	 */
	currentStep++;
}

void MHMCMC::setPRNG(CRandomMersenne &rangeGenerator){
	MHMCMC::pRandGen = shared_ptr<CRandomMersenne>(new CRandomMersenne(rangeGenerator));
}
