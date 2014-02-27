#include "ParameterSequence.h"
#include "randomgen/randomc.h"
#include "myUtils.h"

using namespace std;

/* Default constructor, does lots of nothing.
 */
ParameterSequence::ParameterSequence():
	Parameter()
{
}

/* Constructor from values.
 */
ParameterSequence::ParameterSequence(string inputSequence, string startingSequence, float maxPercentDivergence){
	seq = inputSequence;
	origSeq = startingSequence;
	divThreshold = maxPercentDivergence;
}

/* Virtual destructor, required by for inheritance.
 */
ParameterSequence::~ParameterSequence(){
}

/* Create a sequence based on the current one by changing
 * a single position to some other residue. But don't
 * create one that diverges too much from the original.
 */
shared_ptr<Parameter> ParameterSequence::nextValue(){

	int pos, otherPos;
	string newSeq = seq;
	float divergence = 100*((float)hammingDist(origSeq, newSeq) / (float)origSeq.length());

	/* DEBUG
	 *
	cerr << "Current divergence: " << divergence << endl;
	*/

	string alphabet = "AEQDNLGKSVRTPIMFYCWH";
	string randLetter;

	/* Pick a position to mutate.
	 */
	pos = pRandGen->IRandom(0,seq.length()-1);

	/* Pick a random residue that isn't the previous
	 * one, and isn't too far from the original.
	 */
	while(newSeq.compare(seq) == 0){
		randLetter = alphabet.substr(pRandGen->IRandom(0,alphabet.length()-1), 1);
		newSeq.replace(pos,1,randLetter);
		divergence = 100*((float)hammingDist(origSeq, newSeq) / (float)origSeq.length());

		/* DEBUG
		 *
		cerr << "Mutation: " << (pos+1) << seq.substr(pos,1) << "->" << (pos+1) << newSeq.substr(pos,1) << endl;
		cerr << "Previous sequence: " << seq << endl;
		cerr << "Proposed sequence: " << newSeq << endl;
		cerr << "Proposed divergence: " << divergence << " (threshold " << divThreshold << ")" << endl;
		*/

		/* When the divergence gets too large, randomly revert a
		 * previously changed residue to the original state.
		 */
		if(divergence > divThreshold){
			otherPos = pRandGen->IRandom(0,seq.length()-1);
			while(newSeq.substr(otherPos,1).compare(origSeq.substr(otherPos,1)) == 0 || otherPos == pos){
				otherPos = pRandGen->IRandom(0,seq.length()-1);
			}
			newSeq.replace(otherPos,1,origSeq.substr(otherPos,1));
			divergence = 100*((float)hammingDist(origSeq, newSeq) / (float)origSeq.length());

			/* DEBUG
			 *
			cerr << "Reverted position " << (otherPos+1) << " to the original state (" << newSeq.substr(otherPos,1) << ", " << origSeq.substr(otherPos,1) << ")" << endl;
			cerr << "Previous sequence: " << seq << endl;
			cerr << "Proposed sequence: " << newSeq << endl;
			cerr << "Divergence now " << divergence << endl;
			*/
		}
	}

	return shared_ptr<Parameter>(new ParameterSequence(newSeq, origSeq, divThreshold));
}

/* Make more than one substitution at a time.
 *
 * CURRENTLY JUST MAKES ONE NO MATTER WHAT THE INPUT.
 */
shared_ptr<Parameter> ParameterSequence::nextValue(float stepSize){
	return nextValue();
}

/* Make some variable number of substitutions at
 * once.
 *
 * CURRENTLY JUST MAKE ONE NO MATTER WHAT THE INPUT.
 */
shared_ptr<Parameter> ParameterSequence::nextValue(float stepSize, float stepVar){
	return nextValue();
}

/* Create a random new sequence.
 *
 * SCALE FACTOR MEANINGLESS HERE, JUST GENERATES A RANDOM
 * SEQUENCE OF APPROPRIATE LENGTH.
 */
shared_ptr<Parameter> ParameterSequence::randomValue(float scaleFactor){
	return shared_ptr<Parameter>(new ParameterSequence(generateRandomSequence(seq.length(), "AEQDNLGKSVRTPIMFYCWH", *pRandGen), origSeq, divThreshold));
}

/* Human-readable form of the parameter.
 * Here, simply the sequence in string form.
 */
string ParameterSequence::toString() const{
	return seq;
}
