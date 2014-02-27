#include "MCMC.h"
#include "Proposal.h"
#include <iostream>

/* Default constructor, doesn't do much.
 */
MCMC::MCMC(){
}

/* Constructor from parameters.
 */
MCMC::MCMC(shared_ptr<Proposal> startState, float startTemp, int numSteps){
	currentState = startState;
	bestState = startState;

	temp = startTemp;
	
	maxSteps = numSteps;
	currentStep = 0;
}

/* Virtual destructor, required by gcc.
 */
MCMC::~MCMC(){
	/* DEBUG
	 *
	cout << "Destroying MCMC class properties." << endl;
	*/
}

/* Making a single step. Can be overridden in subclasses
 * that e.g. decrease the step length as they go (such as
 * in simulated annealing).
 *
 * Note that the algorithm is aiming toward the _lowest_
 * probability density here, _not_ the maximum.
 *
 * MEMORY MANAGEMENT IS ALL FUCKED HERE...
 */
void MCMC::takeStep(){

	shared_ptr<Proposal> pNext = currentState->nextStep(temp);

	/* Take a step if the acceptance function approves.
	 */
	if( moveToState(pNext) ){
		
		currentState = pNext;

		/* Is is a better state than any previously seen?
		 */
		if(currentState->getProbDensity() < bestState->getProbDensity() ){
			bestState = currentState;
		}

	}
	/* Not going anywhere -- deleting the proposed state (avoids
	 * memory leaks).
	 *
	 * KIND OF. WHAT HAPPENS IF YOU ACTUALLY _TAKE_ THE STEP?
	 *
	 * SHOULDN'T THIS HAPPEN AUTOMATICALLY AS THE POINTER GOES OUT
	 * OF SCOPE?
	 */
	else{
		//delete(pNext);
	}

	/* One more step taken.
	 */
	currentStep++;
}

/* Convenience function to do a full simulation at once.
 * Simply does takeStep() the required number of times.
 */
void MCMC::run(){
	int i;
	for(i = 0; i < maxSteps; i++){
		takeStep();
	}
}

/* Function to evaluate moving to a proposed new state
 * (acceptance function). Dummy implementation, should
 * be overridden in subclasses.
 */
bool MCMC::moveToState(shared_ptr<Proposal> newState){
	return true;
}

/* Accessor methods for the current state of the chain.
 */
shared_ptr<Proposal> MCMC::getCurrentState() const{
	return currentState;
}

/* Accessor method for the best known state of the
 * chain.
 */
shared_ptr<Proposal> MCMC::getBestState() const{
	return bestState;
}
