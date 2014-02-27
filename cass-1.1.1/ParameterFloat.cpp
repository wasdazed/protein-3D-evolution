#include "ParameterFloat.h"
#include "randomgen/randomc.h"
#include "randomgen/stocc.h"
#include <sstream>
#include <stdexcept>

shared_ptr<StochasticLib1> ParameterFloat::pStochastGen;

/* Default constructor, sets parameter to "-1.0"  
 */
ParameterFloat::ParameterFloat():
	Parameter()
{
	paramValue = -1.0;
}

/* Constructor from values.
 */
ParameterFloat::ParameterFloat(float value):
	Parameter()
{

	/* Make sure the class has a PRNG handy.
	 */
	if(Parameter::pRandGen.get() == 0 || ParameterFloat::pStochastGen.get() == 0){
		throw runtime_error("You must supply random number generators via ParameterFloat::setPRNGs()!");
	}

	paramValue = value;
}

/* Virtual destructor for gcc's benefit.
 */
ParameterFloat::~ParameterFloat(){
}

/* Create a new value based on the current one and a step size.
 *
 * In this case, step a positive value drawn from a distribution
 * Norm(stepSize, 1) up or down, as long as the new value falls
 * within (0,1].
 */
shared_ptr<Parameter> ParameterFloat::nextValue(float stepSize){
	return nextValue(stepSize, 1);
}

/* Create a new value based on the current one, a step size and 
 * a variance term.
 *
 * Generates positive-valued steps from a Gaussian distribution, 
 * centered around stepSize and with a variance stepVar,
 * within (0,1].
 */
shared_ptr<Parameter> ParameterFloat::nextValue(float stepSize, float stepVar){
	float step;
	float newValue = -1;

	/* Decide whether to step up or down on the real line,
	 * but never go outside (0,1].
	 */
	while(newValue <= 0 || newValue > 1){

		/* Draw a number from the relevant normal distribution.
		*/
		step = ParameterFloat::pStochastGen->Normal(stepSize, sqrt(stepVar));

		if(Parameter::pRandGen->IRandom(0,1) == 1){
			newValue = paramValue + step;
		}
		else{
			newValue = paramValue - step;
		}
	}

	return shared_ptr<Parameter>(new ParameterFloat(newValue));
}

/* Create a random new value based on some scaling
 * of the parameter space.
 *
 * In this case, return a pointer to a ParameterFloat
 * on [0, scaleFactor).
 */
shared_ptr<Parameter> ParameterFloat::randomValue(float scaleFactor){
	float newValue = Parameter::pRandGen->Random()*scaleFactor;

	return shared_ptr<Parameter>(new ParameterFloat(newValue));
}

/* Create a random new value within the parameter's natural
 * range of (0,1].
 */
shared_ptr<Parameter> ParameterFloat::randomValue(){
	float newValue = Parameter::pRandGen->Random();

	return shared_ptr<Parameter>(new ParameterFloat(newValue));
}

/* Human-readable form of the parameter.
 * Here, simply the parameter value in
 * string form.
 */
string ParameterFloat::toString() const{

	/* Somewhat tricky to string-ify a float.
	 * This seems to be the easiest option.
	 */
	ostringstream os;
	os << paramValue;
	string readable = os.str();

	return readable;
}

/* Floating-point form of the parameter.
 *
 * Here, simply the parameter value.
 */
float ParameterFloat::toFloat() const{

	return paramValue;
}

void ParameterFloat::setPRNGs(CRandomMersenne &rangeGenerator, StochasticLib1 &stochasticGenerator){
	ParameterFloat::pStochastGen = shared_ptr<StochasticLib1>(new StochasticLib1(stochasticGenerator));
	Parameter::setPRNG(rangeGenerator);
}
