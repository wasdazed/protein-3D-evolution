/* Object to hold a single floating point parameter for
 * an MCMC simulation. 
 *
 * Only generates parameter values on (0:1], so if you
 * want bigger values employ the scaled version of 
 * randomValue().
 */

/* Idempotency.
 */
#ifndef PARAMETERFLOAT_H
#define PARAMETERFLOAT_H

/* Forward declarations.
 */
class StochasticLib1;

/* Includes.
 */
#include "Parameter.h"

class ParameterFloat : public Parameter{
	public:
		ParameterFloat();
		ParameterFloat(float value);
		virtual ~ParameterFloat();

		/* Create a new value based on the current one
		 * and a step size.
		 *
		 * Generates positive-valued steps from a
		 * Gaussian distribution, centered around stepSize
		 * and with a variance of 1.
		 */
		shared_ptr<Parameter> nextValue(float stepSize);
		
		/* Create a new value based on the current
		 * one, a step size and a variance term.
		 *
		 * Generates positive-valued steps from a
		 * Gaussian distribution, centered around
		 * stepSize and with a variance stepVar,
		 * within the bounds of (0,1].
		 */
		shared_ptr<Parameter> nextValue(float stepSize, float stepVar);

		/* Create a random new value on (0, scaleFactor].
		 */
		shared_ptr<Parameter> randomValue(float scaleFactor);

		/* Create a random new value on (0,1].
		 */
		shared_ptr<Parameter> randomValue();

		/* Human-readable form of the parameter.
		 * Here, simply the parameter value in
		 * string form.
		 */
		string toString() const;

		/* Floating-point form of the parameter.
		 * Here, simply the parameter value.
		 */
		float toFloat() const;

		static void setPRNGs(CRandomMersenne &rangeGenerator, StochasticLib1 &stochasticGenerator);

	protected:
		float paramValue;

		/* Pointer to stochastic random number generator
		 * (used to generate e.g. normally distributed
		 * samples).
		 */
		static shared_ptr<StochasticLib1> pStochastGen;
};

/* Idempotency end.
 */
#endif
