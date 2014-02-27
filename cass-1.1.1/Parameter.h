/* Object that holds a single parameter (floating point number,
 * string, distribution, etc.) for e.g. the
 * current state of a Markov Chain Monte Carlo process.
 *
 * Mainly defines the interface for subclasses since
 * there are many things that could be a "parameter".
 */

/* Make source file idempotent (avoid multiple includes).
 */
#ifndef PARAMETER_H
#define PARAMETER_H

/* Forward declarations.
 */
class CRandomMersenne;

/* Includes.
 */
#include <string>
#include <tr1/memory>

using namespace std;

/* Interface.
 */
class Parameter{
	public:
		Parameter();
		virtual ~Parameter();

		/* Create a new value based on the current one
		 * and a step size (with or without a specified
		 * variance to the size).
		 */
		virtual shared_ptr<Parameter> nextValue(float stepSize);
		virtual shared_ptr<Parameter> nextValue(float stepSize, float stepVar);

		/* Create a random new value.
		 */
		virtual shared_ptr<Parameter> randomValue(float scaleFactor);

		/* Human-readable form of the parameter.
		 */
		virtual string toString() const;

		/* Floating point version of the parameter.
		 */
		virtual float toFloat() const;

		static void setPRNG(CRandomMersenne &floatGenerator);

	protected:

		/* Pointer to a random number generator.
		 */
		static shared_ptr<CRandomMersenne> pRandGen;

	private:
};

/* End idempotency statement.
 */
#endif
