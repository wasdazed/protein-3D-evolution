/* Object to hold a constrained 3D rotation (-90 to 90 degrees
 * clockwise around the X, Y and Z axes) for an MCMC simulation.
 */

/* Idempotency.
 */
#ifndef PARAMETERROTATION_H
#define PARAMETERROTATION_H

/* Forward declarations.
 */
class CRandomMersenne;
class StochasticLib1;

/* Includes.
 */
#include "Parameter.h"
#include <vector>
#include <string>
#include <tr1/memory>

using namespace std;

/* Interface.
 */
class ParameterRotation : public Parameter{

	public:
		ParameterRotation();
		ParameterRotation(vector<double> rotations);
		virtual ~ParameterRotation();

		/* Create a new value based on the current one
		 * and a step size.
		 *
		 * Generates... POSSIBLE ANGLES? IN INCREMENTS
		 * OF 1, 5, 10 DEGREES?
		 */
		shared_ptr<Parameter> nextValue(float stepSize);
		
		/* Create a new value based on the current
		 * one, a step size and a variance term.
		 *
		 * Generates... ANGLES IN INCREMENTS OF 1 DEGREE?
		 * 5 DEGREES? 
		 */
		shared_ptr<Parameter> nextValue(float stepSize, float stepVar);

		/* Create a random new value on... AN APPROPRIATE
		 * RANGE.
		 */
		shared_ptr<Parameter> randomValue(float scaleFactor);

		/* Human-readable form of the parameter.
		 * Here, simply the parameter value in
		 * string form.
		 */
		string toString() const;

		/* Returns the rotation angles as a floating point
		 * vector.
		 */
		vector<double> getAngles() const;

		// EXPERIMENTING WITH STATIC PRNGS...
		static void setPRNGs(CRandomMersenne &rangeGenerator, StochasticLib1 &stochasticGenerator);

	protected:

		// EXPERIMENTING WITH STATIC PRNGS...
		static shared_ptr<StochasticLib1> pStochastGen;

		vector<double> rotationAngles;

};

/* Idempotency end.
 */
#endif
