/* Object to hold a particular protein sequence for an MCMC simulation.
 */

/* Idempotency.
 */
#ifndef PARAMETERSEQUENCE_H
#define PARAMETERSEQUENCE_H

/* Foward declarations.
 */

/* Includes.
 */
#include "Parameter.h"
#include <tr1/memory>

using namespace std;

/* Interface.
 */
class ParameterSequence : public Parameter{

	public:
		ParameterSequence();
		ParameterSequence(string inputSequence, string startingSequence, float maxPercentDivergence);
		virtual ~ParameterSequence();

		/* Create a sequence based on the current one.
		 */
		shared_ptr<Parameter> nextValue();
		shared_ptr<Parameter> nextValue(float stepSize);
		shared_ptr<Parameter> nextValue(float stepSize, float stepVar);

		/* Create a random new sequence.
		 */
		shared_ptr<Parameter> randomValue(float scaleFactor);

		/* Human-readable form of the parameter.
		 * Here, simply the sequence in
		 * string form.
		 */
		string toString() const;

	protected:

		string seq;
		string origSeq;
		float divThreshold;
};

/* Idempotency end.
 */
#endif
