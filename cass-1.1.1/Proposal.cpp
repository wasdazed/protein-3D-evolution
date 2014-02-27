#include "Proposal.h"
#include "Parameter.h"
#include <sstream>
#include <iostream>

/* Default constructor, does mostly nothing.
 */
Proposal::Proposal(){
}

/* Constructor from parameters.
 */
Proposal::Proposal(vector< shared_ptr<Parameter> > startParams){
	params = startParams;
	probDensity = calcProbDensity();
}

/* Virtual destructor, must exist due to finicky
 * gcc.
 */
Proposal::~Proposal(){

}

/* Propose a new state from the current one.
 *
 * Dummy implementation, always returns the current
 * state (pointer to this object).
 */
shared_ptr<Proposal> Proposal::nextStep(float stepSize){
	return shared_ptr<Proposal>(this);
}

/* Propose a new state from the current one with
 * some variance.
 *
 * Dummy implementation, always returns the current
 * state (pointer to this object).
 */
shared_ptr<Proposal> Proposal::nextStep(float stepSize, float stepVar){
	return shared_ptr<Proposal>(this);
}

/* Propose a random state.
 *
 * Dummy implementation, always returns the current
 * state (pointer to this object).
 */
shared_ptr<Proposal> Proposal::randomStep(float scaleFactor){
	return shared_ptr<Proposal>(this);
}

/* Return relevant properties of the proposal,
 * namely the current probability density (always
 * a floating point number) and a string representation 
 * of the parameters (varies with parameter type),
 * all separated by "\t".
 */
string Proposal::toString(){
	unsigned int i;

	/* Somewhat tricky to string-ify a float.
	 * This seems to be the easiest option.
	 */
	stringstream os;
	os << probDensity;

	string object = os.str();

	for(i = 0; i < params.size(); i++){
		object += "\t";
		object += params[i]->toString();
	}

	return object;
}

/* Accessor method for the probability density.
 */
float Proposal::getProbDensity() const{
	return probDensity;
}
		
/* Acessor method for floating-point representation
 * of parameter values.
 *
 * Relies on the Parameter-derived class to define
 * some sensible value with the toFloat() function.
 */
vector<float> Proposal::getFloatParams() const{
	unsigned int i;
	vector<float> numbers;

	for(i = 0; i < params.size(); i++){
		numbers.push_back(params[i]->toFloat());
	}

	return numbers;
}

/* Function for calculating the probability density of
 * the current parameter set.
 *
 * Dummy implementation, should be overridden by 
 * subclasses. Depends on what distribution you're 
 * trying to sample.
 *
 * Always returns "-1.0".
 */
float Proposal::calcProbDensity(){

	return -1.0;
}
