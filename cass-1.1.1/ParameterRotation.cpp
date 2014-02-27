#include "ParameterRotation.h"
#include "randomgen/randomc.h"
#include "randomgen/stocc.h"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

/* Static shared PRNG.
 */
shared_ptr<StochasticLib1> ParameterRotation::pStochastGen;

/* Default constructor, doesn't do much.
 */
ParameterRotation::ParameterRotation():
	Parameter()
{
}

/* Constructor from values.
 *
 */
ParameterRotation::ParameterRotation(vector<double> rotations):
	Parameter()
{
	unsigned int i;

	/* Make sure the class has a PRNG handy.
	 */
	if(Parameter::pRandGen.get() == 0 || ParameterRotation::pStochastGen.get() == 0){
		throw runtime_error("You must supply random number generators via ParameterRotation::setPRNGs()!");
	}

	/* Verify that these are valid angles.
	 */
	for(i = 0; i < rotations.size(); i++){
		if(rotations[i] >= -90.0 && rotations[i] <= 90.0){
			rotationAngles.push_back(rotations[i]);
		}
		else{
			rotationAngles.push_back(-1.0);
			throw runtime_error("Angle out of range (-90 through 90 degrees)."); 
		}
	}
}

/* Virtual destructor for gcc's benefit.
 *
 */
ParameterRotation::~ParameterRotation(){
	/* DEBUG
	 *
	cout << "Destroying ParameterRotation class properties..." << endl;
	*/
}

/* Create a new value based on the current one
 * and a step size. Value of the step actually
 * drawn from a Gaussian distribution
 * specified by N(stepSize,1). Only does steps
 * in integer degrees.
 */
shared_ptr<Parameter> ParameterRotation::nextValue(float stepSize){
	return nextValue(stepSize,1);
}

/* Create a new value based on the current
 * one, a step size and a variance term.
 * Value of the step actually drawn from a Gaussian 
 * distribution specified by N(stepSize,sqrt(stepVar)). Only does 
 * steps in integer degrees.
 */
shared_ptr<Parameter> ParameterRotation::nextValue(float stepSize, float stepVar){
	double step;
	vector<double> newValues;
	unsigned int i;

	/* Generate a new value by stepping stepSize
	 * up or down within the allowed range.
	 */
	for(i = 0; i < rotationAngles.size(); i++){
		double newValue = -100.0;
		while(newValue < -90.0 || newValue > 90.0){

			/* Draw a number from the default Gaussian
			 * and make sure it's a whole integer.
			 */
			step = floor(ParameterRotation::pStochastGen->Normal(stepSize,sqrt(stepVar)));

			/* Decide whether to go up or down.
			*/
			if(Parameter::pRandGen->IRandom(0,1) == 1){
				newValue = rotationAngles[i] + step;
			}
			else{
				newValue = rotationAngles[i] - step;
			}
		}
		newValues.push_back(newValue);
	}

	return shared_ptr<Parameter>(new ParameterRotation(newValues));
}

/* Create a random new value on [-90, 90.0] for
 * each axis, always an integer.
 *
 * The scaleFactor input is meaningless here, and
 * only perserved for upwards compatability.
 */
shared_ptr<Parameter> ParameterRotation::randomValue(float scaleFactor){
	unsigned int i;
	vector<double> newValues;

	for(i = 0; i < rotationAngles.size(); i++){
		double newValue = -100.0;
		while(newValue < -90.0 || newValue > 90.0){
			newValue = floor(Parameter::pRandGen->Random()*90);
			if(Parameter::pRandGen->Random() > 0.5){
				newValue *= -1.0;
			}
		}
		newValues.push_back(newValue);
	}

	return shared_ptr<Parameter>(new ParameterRotation(newValues));
}

/* Human-readable form of the parameter.
 * Here, "XXX.XXX, YYY.YYY, ZZZ.ZZZ".
 */
string ParameterRotation::toString() const{

	int i;
	char out[9];
	string readable = "";
	
	/* Somewhat tricky to string-ify a collection
	 * of floats. Since we do want a specified precision,
	 * sprintf() should do the trick.
	 */
	for(i = 0; i < 2; i++){
		sprintf(out, "%03.3f, ", rotationAngles[i]);
		readable += out;
	}
	sprintf(out, "%03.3f", rotationAngles[i]);
	readable += out;

	return readable;
}


/* Return the rotation angles as a floating point
 * vector.
 */
vector<double> ParameterRotation::getAngles() const{
	return rotationAngles;
}

/* Set the shared PRNGs.
 */
void ParameterRotation::setPRNGs(CRandomMersenne &rangeGenerator, StochasticLib1 &stochasticGenerator){
	ParameterRotation::pStochastGen = shared_ptr<StochasticLib1>(new StochasticLib1(stochasticGenerator));
	Parameter::setPRNG(rangeGenerator);
}
