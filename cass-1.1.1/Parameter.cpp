#include "Parameter.h"
#include "randomgen/randomc.h"
#include <iostream>

shared_ptr<CRandomMersenne> Parameter::pRandGen; 

/* Default constructor, doesn't do much.
 */
Parameter::Parameter(){
};

/* Virtual destructor, for the benefit of gcc.
 */
Parameter::~Parameter(){
	/* DEBUG
	 *
	cout << "Destroying Parameter class properties..." << endl;
	*/
};

/* Create a new value based on the current one
 * and a step size.
 *
 * Dummy implementation, always returns a pointer
 * to the current state (this object).
 */
shared_ptr<Parameter> Parameter::nextValue(float stepSize){
	return shared_ptr<Parameter>(this);
}

/* Create a new value based on the current one,
 * a step size and the variance in the step size.
 *
 * Dummy implementation, always returns a pointer
 * to the current state (this object).
 */
shared_ptr<Parameter> Parameter::nextValue(float stepSize, float stepVar){
	return shared_ptr<Parameter>(this);
}

/* Create a random new value.
 *
 * Dummy implementation, always returns a pointer
 * to the current state (this object).
 */
shared_ptr<Parameter> Parameter::randomValue(float scaleFactor){
	return shared_ptr<Parameter>(this);
}

/* Human-readable form of the parameter.
 *
 * Dummy implementation, should be overriden
 * in subclasses. 
 *
 * Always returns "DUMMY".
 */
string Parameter::toString() const{
	string dummy = "DUMMY";

	return dummy;
}

/* Floating point version of the parameter.
 *
 * Dummy implementation, should be overridden
 * in subclasses.
 *
 * Always returns "-1.0".
 */
float Parameter::toFloat() const{
	return -1.0;
}

void Parameter::setPRNG(CRandomMersenne &floatGenerator){
	Parameter::pRandGen = shared_ptr<CRandomMersenne>(new CRandomMersenne(floatGenerator));
}
