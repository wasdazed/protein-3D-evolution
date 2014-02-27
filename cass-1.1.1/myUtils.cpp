#include "myUtils.h"
#include <climits>

/* The scores for the BLOSUM62 matrix, and its key.
 * Retrieved from
 * http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
 */
string blosum62Key = "ARNDCQEGHILKMFPSTWYVBZX";
int blosum62Matrix[23][23] = {
{4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0},
{-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1},
{-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1},
{-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1},
{0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2},
{-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1},
{-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1},
{0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1},
{-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1},
{-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1},
{-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1},
{-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1},
{-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1},
{-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1},
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2},
{1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0},
{0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0},
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2},
{-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1},
{0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1},
{-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1},
{-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1},
{0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1}
};

/* Function to convert an integer to a string.
 *
 * Courtesy of Bjarne Soustrup. Or rather shamelessly stolen from 
 * http://www.research.att.com/~bs/bs_faq2.html#int-to-string
 */
string itos(int i){
	stringstream s;
	s << i;
	return s.str();
}

/* Function to tokenize a string into a vector of
 * strings based on some delimiting string (defaults
 * to " ", i.e. a space character).
 *
 * Stolen from http://www.oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
 * and slightly modified.
 */
vector<string> tokenize(const string &str, const string &delimiters){

	vector<string> tokens;

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));

		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);

		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;
}

/* REALLY SHITTY RANDOM INTEGER GENERATOR,
 * NEEDS IMPROVEMENT.
 *
 * STOLE FROM SHRUTI, WHO PROBABLY STOLE
 * IT ELSEWHERE.
 *
 * SHOULD PROBABLY BE REMOVED.
 */
int rand_int(int n)	// random number generator
{
	return rand()%n;
}

/* Function to place a floating point value in a bin
 * according to supplied limits of the bins. The
 * ranges are expressed as a vector of (min, max)
 * vectors, and the ranges are exclusive (no sense
 * in trying to find equal floating point numbers in 
 * binary representation).
 *
 * Returns a bin number starting from 0, or -1 if
 * the supplied value doesn't fit in a bin.
 */
int putInBin(float value, vector< vector<float> > binRanges){
	unsigned int i;
	int bin = -1;
	
	for(i = 0; i < binRanges.size(); i++){
		
		/* DEBUG
		 *
		cout << "Testing bin " << i+1 << ": " << binRanges[i][0] << "-" << binRanges[i][1] << "..." << endl;
		*/

		if(value > binRanges[i][0] && value < binRanges[i][1]){
			
			/* DEBUG
			 *
			cout << "Match!" << endl;
			*/

			bin = i;
			break;
		}
	}

	return bin;
}

/* Function to compute the geometric centroid of a collection of
 * points in N-dimensional Euclidian space.
 */
coord NDCentroid(vector<coord> points){
	coord centroid;
	unsigned int i, j;
	double avgDim;

	/* Take an average for each dimension...
	 */
	for(i = 0; i < points[0].size(); i++){

		/* ... by averaging the values for each point.
		 */
		avgDim = 0.0;
		for(j = 0; j < points.size(); j++){
			avgDim += points[j][i];
		}
		
		avgDim /= points.size();

		centroid.push_back(avgDim);
	}

	return centroid;
}

/* Function to compute the distance between two points in 
 * N-dimensional Euclidian space. The two points are
 * assumed to have the same dimensionality.
 */
double NDDist(coord point1, coord point2){
	double sumOfSquares = 0.0;
	unsigned int i;

	/* DEBUG
	 *
		cout << "Computing distance between points (" << point1[0] << "," << point1[1] << "," << point1[2] << ") and (" << point2[0] << "," << point2[1] << "," << point2[2] << ")" << endl;
	 */

	/* Compute the square of the distance in each
	 * dimension, and add them up.
	 */
	for(i = 0; i < point1.size(); i++){
		sumOfSquares += (point2.at(i) - point1.at(i))*(point2.at(i) - point1.at(i));

		/* DEBUG
		 *
		cout << "Dim: " << i << endl;
		cout << "SS: " << sumOfSquares << endl;
		*/
	}

	/* Return the square root of the sum of squares
	 * (i.e. the Euclidian metric).
	 */
	return sqrt(sumOfSquares);
}

/* Function to compute a N-dimensional angle in Euclidean 
 * space, i.e. between two vectors u and v. The vectors are 
 * represented as 3 points: vector u goes from point 1 to 2, 
 * and vector v goes from point 2 to 3, resulting in an angle 
 * 1-2-3.
 *
 * Returns the angle in degrees.
 */
double NDAngle(coord point1, coord point2, coord point3){
	/* In Euclidean space, for vectors u and v,
	 *
	 * u (dot) v = cos(theta) * norm(u) * norm(v)
	 *
	 * therefore
	 *
	 * theta = arccos( u (dot) v / [ norm(u) * norm(v) ] )
	 */

	double theta, normU, normV, uDotv;
	vector<double> u,v;
	unsigned int i;

	/* Construct u and v from points 1, 2
	 * and 3.
	 */
	for(i = 0; i < point1.size(); i++){
		u.push_back(point1[i] - point2[i]);
		v.push_back(point3[i] - point2[i]); // Scale u and v to start at origin
	}

	/* Calculate dot product and norms.
	 */
	uDotv = inner_product(u.begin(), u.end(), v.begin(), 0.0); 
	normU = NDDist(point1, point2);
	normV = NDDist(point2, point3);

	/* Calculate the angle.
	 */
	theta = acos(uDotv/(normU*normV)) * ((double)180.0/PI); // Convert to degrees from radians

	return theta;
}

/* Function to calculate the cross product of two vectors. 
 * Assumes a right-handed, orthogonal coordinate system, and
 * ignores dimensionality above 3.
 */
vector<double> crossProduct(coord vector1, coord vector2){

	/* The cross product of two vectors a and b:
	 *
	 * a x b = (a2*b3 − a3*b2, a3*b1 − a1*b3, a1*b2 − a2*b1)
	 */
	vector<double> crossProd;

	crossProd.push_back(vector1[1]*vector2[2] - vector1[2]*vector2[1]);
	crossProd.push_back(vector1[2]*vector2[0] - vector1[0]*vector2[2]);
	crossProd.push_back(vector1[0]*vector2[1] - vector1[1]*vector2[0]);

	return crossProd;
}

/* Function to compute the dihedral angle in 3
 * dimensions, defined here as the torsional angle
 * of the bond vector between point2 and point3 (point1
 * and point4 are needed for the calculation). Produces
 * a properly signed value due to method of calculation
 * (e.g. angles range from 0-180, -180-0 instead of
 * 0-360).
 *
 * Returns the angle in degrees.
 */
double dihedralAngle(coord point1, coord point2, coord point3, coord point4){

	/* For a dihedral angle of 4 points/3 vectors, 
	 *
	 * phi = atan2( (norm(b2)*b1) (dot) (b2 x b3), (b1 x b2) (dot) (b2 x b3) )
	 *
	 * where 
	 *
	 * b1 = point2 - point1
	 * b2 = point3 - point2
	 * b3 = point4 - point3
	 */
	double phi;
	vector<double> b1, b2, b3, b1_scaled;
	double norm_b2;
	vector<double> b2xb3, b1xb2;
	int i;

	/* Construct vectors b1, b2, and b3 from points.
	 */
	for(i = 0; i < 3; i++){
		b1.push_back(point2[i] - point1[i]);
		b2.push_back(point3[i] - point2[i]);
		b3.push_back(point4[i] - point3[i]);
	}

	/* Calculate norm and cross products.
	 */
	norm_b2 = NDDist(point2, point3);
	for(i = 0; i < 3; i++){
		b1_scaled.push_back(norm_b2*b1[i]);
	}

	b2xb3 = crossProduct(b2, b3);
	b1xb2 = crossProduct(b1, b2);

	/* Calculate phi as above.
	 */
	phi = atan2( inner_product(b1_scaled.begin(), b1_scaled.end(), b2xb3.begin(), 0) , inner_product(b1xb2.begin(), b1xb2.end(), b2xb3.begin(), 0) );

	/* Convert to degrees and return.
	 */
	return phi*((double)180.0/PI);
}

/* Function to find the range of a floating point vector,
 * here defined the same as the range of a distribution of
 * numbers (e.g. max(vector)-min(vector)).
 */
float vectorRange(vector<float> numbers){
	float maxVal = *max_element(numbers.begin(), numbers.end());
	float minVal = *min_element(numbers.begin(), numbers.end());
	return maxVal-minVal;
}

/* Function to construct a vector from two points in N-dimensional
 * space assuming that point #1 is the center of the coordinate
 * system (i.e. we can represent the vector a set of coordinates).
 */
coord makeRelativeVector(coord point1, coord point2){
	coord v;
	unsigned int i;

	/* Construct vector v from points 1 and 2.
	 */
	for(i = 0; i < point1.size(); i++){
		v.push_back(point2[i] - point1[i]); // i.e. put origin at point1
	}

	return v;
}

/* Function to compute the sum of a set of vectors going out
 * from the origin. Since they all start at the origin we
 * can represent them as coordinates in N-space.
 *
 * We assume that all vectors have the same dimensionality
 * as the first one (additional dimensions will not be
 * computed).
 */
coord vectorSum(vector<coord> vectors){
	unsigned int i,j;
	double dimensionSum;
	coord vectorSum;

	/* Construct the sum of vectors for all
	 * dimensions and combine.
	 */
	for(i = 0; i < vectors[0].size(); i++){
		
		dimensionSum = 0;

		for(j = 0; j < vectors.size(); j++){
			dimensionSum += vectors[j][i];
		}

		vectorSum.push_back(dimensionSum);
	}

	return vectorSum;
}

/* Scale the length of an N-dimensional vector by some factor.
 */
coord scaleVector(coord theVector, double scaleFactor){
	coord newVector;
	unsigned int i;

	/* Multiply all dimensions by the required factor.
	 */
	for(i = 0; i < theVector.size(); i++){
		newVector.push_back(theVector[i] * scaleFactor);
	}

	return newVector;
}

/* Function to rotate a vector (presumed to start at the origin
 * and extend to the supplied coordinates) some number of
 * degrees around the x, y and z axes. 
 *
 * Only works in 3D. Rotations are counter-clockwise for
 * positive angles and clockwise for negative angles. See
 * e.g. http://en.wikipedia.org/wiki/Rotation_matrix for
 * further information.
 *
 * On invalid degree input (must lie on -360.0 <= x <= 360.0), or
 * non-3D coordinates/angles, returns (0,0,0).
 *
 * Returns the "end" of the rotated vector (again assuming
 * it starts at the origin).
 */
coord rotateVector(coord theVector, vector<double> xyzRotationAngles){
	coord rotated(3, 0.0);
	double xRads, yRads, zRads;

	/* Check the input.
	 */
	if(xyzRotationAngles.size() != 3 || theVector.size() != 3 || xyzRotationAngles[0] < -360.0 || xyzRotationAngles[0] > 360.0 || xyzRotationAngles[1] < -360.0 || xyzRotationAngles[1] > 360.0 || xyzRotationAngles[2] < -360.0 || xyzRotationAngles[2] > 360.0){
		return rotated;
	}

	/* Convert degrees to radians.
	 */
	xRads = xyzRotationAngles[0]*(PI/(double)180.0);
	yRads = xyzRotationAngles[1]*(PI/(double)180.0);
	zRads = xyzRotationAngles[2]*(PI/(double)180.0);

	/* Rotate around the x axis.
	 */
	rotated[0] = theVector[0];
	rotated[1] = theVector[1]*cos(xRads) - theVector[2]*sin(xRads);
	rotated[2] = theVector[2]*cos(xRads) + theVector[1]*sin(xRads);

	/* And around the y axis.
	 */
	rotated[0] = rotated[0]*cos(yRads) + rotated[2]*sin(yRads);
	rotated[1] = rotated[1];
	rotated[2] = -1.0*rotated[0]*sin(yRads) + rotated[2]*cos(yRads);

	/* Finally around the z axis.
	 */
	rotated[0] = rotated[0]*cos(zRads) - rotated[1]*sin(zRads);
	rotated[1] = rotated[1]*cos(zRads) + rotated[0]*sin(zRads);
	rotated[2] = rotated[2];

	return rotated;
}

/* Function to reflect a vector (represented as coordinates,
 * as if its start lies at the origin) around a point in 3D space, 
 * according to
 *
 * Ref(vector) = 2*point - vector
 *
 * See e.g. http://en.wikipedia.org/wiki/Point_reflection
 * for more information.
 *
 * Returns (0,0,0) on error.
 */
coord reflectAroundPoint(coord theVector, coord point){
	int i;
	coord newVector(3,0.0);

	if(theVector.size() != 3 || point.size() != 3){
		return newVector;
	}

	/* Perform operations in 3-space.
	 */
	for(i = 0; i < 3; i++){
		newVector[i] = 2*point[i] - theVector[i];
	}

	return newVector;
}

/* Function to calculate the mean of a vector of
 * floating-point numbers, using Knuth's one-pass
 * algorithm.
 * 
 * Returns -9999.0 on error.
 */
float mean(vector<float> distribution){
	float mean = 0.0;
	float delta;
	unsigned int i;

	if(distribution.size() < 2){
		return -9999.0;
	}

	/* Use Knuth's one-pass algorithm for mean.
	 */
	for(i = 0; i < distribution.size(); i++){
		delta = distribution[i] - mean;
		mean += delta/(i+1);
	}

	return mean;
}


/* Function to calculation the variance of a vector of
 * floating-point numbers, using Knuth's one-pass
 * algorithm.
 * 
 * Returns -1.0 on error.
 */
float variance(vector<float> distribution){
	float mean = 0.0;
	float Mtwo = 0.0;
	float delta;
	float var = -1.0;
	unsigned int i;

	if(distribution.size() < 2){
		return var;
	}

	/* Use Knuth's one-pass algorithm for variance.
	 * This also requires calculating the mean,
	 * unfortunately.
	 */
	for(i = 0; i < distribution.size(); i++){

		delta = distribution[i] - mean;
		mean += delta/(i+1);
		Mtwo += delta*(distribution[i] - mean);
	}

	var = Mtwo/(i-1);

	return var;
}

/* Function to calculate the quartiles of a distribution
 * of floating point values.  Uses naive sorting 
 * implementation, doesn't deal well with NaN values
 * (exactly as badly as algorithm::sort(), in fact).
 *
 * Returns vector of (min, 25%, 50%, 75%, 100%) of
 * the distribution (i.e. minimum, 1Q, median, 3Q, maximum).
 */
vector<float> quartiles(vector<float> distribution){
	vector<float> values(5, 0.0);

	stable_sort(distribution.begin(), distribution.end());
	int middle = (int)floor(distribution.size()/2);
	int quarter = (int)floor(middle/2);
	values[0] = distribution.front(); //0%
	values[1] = distribution.at(quarter); //25%
	values[2] = distribution.at(middle); //50%
	values[3] = distribution.at(middle+quarter); //75%
	values[4] = distribution.back(); //100%
	
	return values;
}

/* Function to trim the whitespace from the start and
 * end of strings. Useful for e.g. parsing tasks.
 * Returns the trimmed string. DOES NOT trim tabs,
 * only space characters.
 */
string trimWhitespace(string input){
	string output = "";
	int i;
	int spacesFront = 0;
	int spacesBack = 0;

	/* Count whitespace characters up front.
	 */
	for(i = 0; i < input.length(); i++){
		if(input[i] != ' '){
			break;
		}
		spacesFront++;
	}

	/* Count whitespace characters in the
	 * back.
	 */
	for(i = input.length()-1; i > -1; i--){
		if(input[i] != ' '){
			break;
		}
		spacesBack++;
	}
	
	/* Construct a string without all those spaces.
	 */
	output = input.substr( (0+spacesFront), (input.length()-1-spacesBack) );

	return output;
}

/* Function to make a random transversion (A -> [C/T], C -> [A/G],
 * G -> [C/T], T -> [A/G]) for a single nucleotide.
 */
string makeTransversion(string nucleotide, CRandomMersenne &randGenerator){
	int randInt;

	/* Change can be to either of two other nucleotides,
	 * so it's a binary choice.
	 */
	randInt = randGenerator.IRandom(0,1);
	
	/* A and G transverse to C or T.
	 */
	if(nucleotide.compare("A") == 0 || nucleotide.compare("G") == 0){
		if(randInt == 0){
			return "C";
		}
		else{
			return "T";
		}
	}
	/* C and T transverse to A or G.
	 */
	else{
		if(randInt == 0){
			return "A";
		}
		else{
			return "G";
		}
	}
}

/* Function to make a transition (A <-> G, C <-> T) for a single nucleotide.
 *
 * On any input other than AGCT, returns "X".
 */
string makeTransition(string nucleotide){

	if(nucleotide.compare("A") == 0){
		return "G";
	}
	else if(nucleotide.compare("G") == 0){
		return "A";
	}
	else if(nucleotide.compare("C") == 0){
		return "T";
	}
	else if(nucleotide.compare("T") == 0){
		return "C";
	}
	else{
		return "X";
	}
}

/* Method for mutating DNA according to the Kimura 2-parameter
 * model (K80, with transitions being twice as probable as
 * transversions). Each nucleotide is mutated k times, with
 * k drawn from a Poisson distribution with a given lambda
 * (the average mutation rate per site per generation).
 */
string mutateDna(string dnaSequence, float lambda, StochasticLib1 &stochasticGenerator, CRandomMersenne &randGenerator){

	/* DEBUG
	 *
	cout << "Mutating sequence..." << endl;
	*/

	string newSequence = dnaSequence;
	
	unsigned int i, j, numMutations, randInt;

	/* For each nucleotide in the sequence, give it a possiblity of
	 * mutation.
	 */
	for(i = 0; i < newSequence.length(); i++){

		/* Randomly draw an integer number of mutations from a Poisson 
		 * distribution with lambda = mutation rate.
		 */
		numMutations = stochasticGenerator.Poisson(lambda); 

		/* If we're mutating something...
		*/
		if(numMutations != 0)
		{
			/* For each mutation, make a transition/transversion
			 * with probability ratio 2:1.
			 */
			for(j = 0; j < numMutations; j++){

				/* Draw a random number between 0 and 2.
				 */
				randInt = randGenerator.IRandom(0,2);

				/* 0 and 1 are transitions, 2 is a transversion.
				 */
				if(randInt < 2){
					newSequence.replace( i, 1, makeTransition(newSequence.substr(i,1)) );
				}
				else{
					newSequence.replace( i, 1, makeTransversion(newSequence.substr(i,1), randGenerator) );
				 }
			}
		}
	}

	return newSequence;
}


/* Function for translating a DNA sequence into protein,
 * according to the standard (nuclear) code. Will only translate an
 * integer number of codons (extra nucleotides on the "end"
 * will not have any effect). Only takes upper-case input.
 * Stop codons are translated to "_", unknown codons to "?".
 *
 * Translation table from Tisdall, "Beginning Perl for Bioinformatics",
 * 2001, O'Reilly & Associates, Sebastopol, CA.
 */
string translateDna(string dnaSeq){
	string proteinSeq = "";
	string codon = "";
	unsigned int i;

	for(i = 0; i < dnaSeq.length(); i = i + 3){
		codon = dnaSeq.substr(i,3);

		/* Translation table adapted from Tisdall's Perl solution.	
		 */
	    	if(codon.compare("TCA") == 0){
			proteinSeq += "S";	// Serine
		}
		else if(codon.compare("TCC") == 0){
			proteinSeq += "S";	// Serine
		}
		else if(codon.compare("TCG") == 0){
			proteinSeq += "S";	// Serine
		}
		else if(codon.compare("TCT") == 0){
			proteinSeq += "S";	// Serine
		}
		else if(codon.compare("TTC") == 0){
			proteinSeq += "F";	// Phenylalanine
		}
		else if(codon.compare("TTT") == 0){
			proteinSeq += "F";	// Phenylalanine
		}
		else if(codon.compare("TTA") == 0){
			proteinSeq += "L";	// Leucine
		}
		else if(codon.compare("TTG") == 0){
			proteinSeq += "L";	// Leucine
		}
		else if(codon.compare("TAC") == 0){
			proteinSeq += "Y";	// Tyrosine
		}
		else if(codon.compare("TAT") == 0){
			proteinSeq += "Y";	// Tyrosine
		}
		else if(codon.compare("TAA") == 0){
			proteinSeq += "_";	// Stop
		}
		else if(codon.compare("TAG") == 0){
			proteinSeq += "_";	// Stop
		}
		else if(codon.compare("TGC") == 0){
			proteinSeq += "C";	// Cysteine
		}
		else if(codon.compare("TGT") == 0){
			proteinSeq += "C";	// Cysteine
		}
		else if(codon.compare("TGA") == 0){
			proteinSeq += "_";	// Stop
		}
		else if(codon.compare("TGG") == 0){
			proteinSeq += "W";	// Tryptophan
		}
		else if(codon.compare("CTA") == 0){
			proteinSeq += "L";	// Leucine
		}
		else if(codon.compare("CTC") == 0){
			proteinSeq += "L";	// Leucine
		}
		else if(codon.compare("CTG") == 0){
			proteinSeq += "L";	// Leucine
		}
		else if(codon.compare("CTT") == 0){
			proteinSeq += "L";	// Leucine
		}
		else if(codon.compare("CCA") == 0){
			proteinSeq += "P";	// Proline
		}
		else if(codon.compare("CCC") == 0){
			proteinSeq += "P";	// Proline
		}
		else if(codon.compare("CCG") == 0){
			proteinSeq += "P";	// Proline
		}
		else if(codon.compare("CCT") == 0){
			proteinSeq += "P";	// Proline
		}
		else if(codon.compare("CAC") == 0){
			proteinSeq += "H";	// Histidine
		}
		else if(codon.compare("CAT") == 0){
			proteinSeq += "H";	// Histidine
		}
		else if(codon.compare("CAA") == 0){
			proteinSeq += "Q";	// Glutamine
		}
		else if(codon.compare("CAG") == 0){
			proteinSeq += "Q";	// Glutamine
		}
		else if(codon.compare("CGA") == 0){
			proteinSeq += "R";	// Arginine
		}
		else if(codon.compare("CGC") == 0){
			proteinSeq += "R";	// Arginine
		}
		else if(codon.compare("CGG") == 0){
			proteinSeq += "R";	// Arginine
		}
		else if(codon.compare("CGT") == 0){
			proteinSeq += "R";	// Arginine
		}
		else if(codon.compare("ATA") == 0){
			proteinSeq += "I";	// Isoleucine
		}
		else if(codon.compare("ATC") == 0){
			proteinSeq += "I";	// Isoleucine
		}
		else if(codon.compare("ATT") == 0){
			proteinSeq += "I";	// Isoleucine
		}
		else if(codon.compare("ATG") == 0){
			proteinSeq += "M";	// Methionine
		}
		else if(codon.compare("ACA") == 0){
			proteinSeq += "T";	// Threonine
		}
		else if(codon.compare("ACC") == 0){
			proteinSeq += "T";	// Threonine
		}
		else if(codon.compare("ACG") == 0){
			proteinSeq += "T";	// Threonine
		}
		else if(codon.compare("ACT") == 0){
			proteinSeq += "T";	// Threonine
		}
		else if(codon.compare("AAC") == 0){
			proteinSeq += "N";	// Asparagine
		}
		else if(codon.compare("AAT") == 0){
			proteinSeq += "N";	// Asparagine
		}
		else if(codon.compare("AAA") == 0){
			proteinSeq += "K";	// Lysine
		}
		else if(codon.compare("AAG") == 0){
			proteinSeq += "K";	// Lysine
		}
		else if(codon.compare("AGC") == 0){
			proteinSeq += "S";	// Serine
		}
		else if(codon.compare("AGT") == 0){
			proteinSeq += "S";	// Serine
		}
		else if(codon.compare("AGA") == 0){
			proteinSeq += "R";	// Arginine
		}
		else if(codon.compare("AGG") == 0){
			proteinSeq += "R";	// Arginine
		}
		else if(codon.compare("GTA") == 0){
			proteinSeq += "V";	// Valine
		}
		else if(codon.compare("GTC") == 0){
			proteinSeq += "V";	// Valine
		}
		else if(codon.compare("GTG") == 0){
			proteinSeq += "V";	// Valine
		}
		else if(codon.compare("GTT") == 0){
			proteinSeq += "V";	// Valine
		}
		else if(codon.compare("GCA") == 0){
			proteinSeq += "A";	// Alanine
		}
		else if(codon.compare("GCC") == 0){
			proteinSeq += "A";	// Alanine
		}
		else if(codon.compare("GCG") == 0){
			proteinSeq += "A";	// Alanine
		}
		else if(codon.compare("GCT") == 0){
			proteinSeq += "A";	// Alanine
		}
		else if(codon.compare("GAC") == 0){
			proteinSeq += "D";	// Aspartic Acid
		}
		else if(codon.compare("GAT") == 0){
			proteinSeq += "D";	// Aspartic Acid
		}
		else if(codon.compare("GAA") == 0){
			proteinSeq += "E";	// Glutamic Acid
		}
		else if(codon.compare("GAG") == 0){
			proteinSeq += "E";	// Glutamic Acid
		}
		else if(codon.compare("GGA") == 0){
			proteinSeq += "G";	// Glycine
		}
		else if(codon.compare("GGC") == 0){
			proteinSeq += "G";	// Glycine
		}
		else if(codon.compare("GGG") == 0){
			proteinSeq += "G";	// Glycine
		}
		else if(codon.compare("GGT") == 0){
			proteinSeq += "G";	// Glycine
		}
		else{
			proteinSeq += "?";
		}
	}

	return proteinSeq;
}

/* Generates a random sequence of characters of specified length from a
 * specified alphabet.
 */
string generateRandomSequence(unsigned int length, string alphabet, CRandomMersenne &randGen){
	string randSeq = "";
	unsigned int randomChar;

	/* Create the specified number of random characters
	 */
	while(randSeq.length() < length){

		/* Pick a random position in the alphabet.
		 */
		randomChar = randGen.IRandom(0,alphabet.length()-1);

		/* Add it to the string
		 */
		randSeq += alphabet.substr(randomChar,1);
	}

	return randSeq;
}

/* Generates a random sequence of characters by permuting the ones in
 * the existing string.
 *
 * Really just a wrapper for STL random_shuffle.
 */
string shuffleSequence(string inputString){

	vector<string> characters;
	unsigned int i;
	string outString = "";

	/* Convert to vector of characters.
	 */
	for(i = 0; i < inputString.length(); i++){
		characters.push_back(inputString.substr(i,1));
	}

	/* Shuffle the vector.
	 */
	random_shuffle(characters.begin(), characters.end());

	/* Reconstruct the string.
	 */
	for(i = 0; i < characters.size(); i++){
		outString += characters[i];
	}

	return outString;
}


/* Function to convert three-letter amino acid names to
 * their one-letter abbreviation. All-upper-case name is
 * expected as input.
 * Returns "X" for an unrecognized amino acid.
 */
string oneLetterName(string residue){
	string oneLetter = "X";

	if(residue.compare("ALA") == 0){
		oneLetter = "A";
	}
	else if(residue.compare("ARG") == 0){
		oneLetter = "R";
	}
	else if(residue.compare("ASN") == 0){
		oneLetter = "N";
	}
	else if(residue.compare("ASP") == 0){
		oneLetter = "D";
	}
	else if(residue.compare("CYS") == 0){
		oneLetter = "C";
	}
	else if(residue.compare("GLU") == 0){
		oneLetter = "E";
	}
	else if(residue.compare("GLN") == 0){
		oneLetter = "Q";
	}
	else if(residue.compare("GLY") == 0){
		oneLetter = "G";
	}
	else if(residue.compare("HIS") == 0){
		oneLetter = "H";
	}
	else if(residue.compare("ILE") == 0){
		oneLetter = "I";
	}
	else if(residue.compare("LEU") == 0){
		oneLetter = "L";
	}
	else if(residue.compare("LYS") == 0){
		oneLetter = "K";
	}
	else if(residue.compare("MET") == 0){
		oneLetter = "M";
	}
	else if(residue.compare("PHE") == 0){
		oneLetter = "F";
	}
	else if(residue.compare("PRO") == 0){
		oneLetter = "P";
	}
	else if(residue.compare("SER") == 0){
		oneLetter = "S";
	}
	else if(residue.compare("THR") == 0){
		oneLetter = "T";
	}
	else if(residue.compare("TRP") == 0){
		oneLetter = "W";
	}
	else if(residue.compare("TYR") == 0){
		oneLetter = "Y";
	}
	else if(residue.compare("VAL") == 0){
		oneLetter = "V";
	}

	return oneLetter;
}

/* Function to convert one-letter amino acid codes
 * to their three-letter representation. Essentially
 * inverts oneLetterName().
 *
 * Returns "XXX" on failure to recognize the one-letter
 * code.
 */
string threeLetterName(string oneLetterCode){
	string threeLetter = "XXX";

	if(oneLetterCode.compare("A") == 0){
		threeLetter = "ALA";
	}
	else if(oneLetterCode.compare("R") == 0){
		threeLetter = "ARG";
	}
	else if(oneLetterCode.compare("N") == 0){
		threeLetter = "ASN";
	}
	else if(oneLetterCode.compare("D") == 0){
		threeLetter = "ASP";
	}
	else if(oneLetterCode.compare("C") == 0){
		threeLetter = "CYS";
	}
	else if(oneLetterCode.compare("E") == 0){
		threeLetter = "GLU";
	}
	else if(oneLetterCode.compare("Q") == 0){
		threeLetter = "GLN";
	}
	else if(oneLetterCode.compare("G") == 0){
		threeLetter = "GLY";
	}
	else if(oneLetterCode.compare("H") == 0){
		threeLetter = "HIS";
	}
	else if(oneLetterCode.compare("I") == 0){
		threeLetter = "ILE";
	}
	else if(oneLetterCode.compare("L") == 0){
		threeLetter = "LEU";
	}
	else if(oneLetterCode.compare("K") == 0){
		threeLetter = "LYS";
	}
	else if(oneLetterCode.compare("M") == 0){
		threeLetter = "MET";
	}
	else if(oneLetterCode.compare("F") == 0){
		threeLetter = "PHE";
	}
	else if(oneLetterCode.compare("P") == 0){
		threeLetter = "PRO";
	}
	else if(oneLetterCode.compare("S") == 0){
		threeLetter = "SER";
	}
	else if(oneLetterCode.compare("T") == 0){
		threeLetter = "THR";
	}
	else if(oneLetterCode.compare("W") == 0){
		threeLetter = "TRP";
	}
	else if(oneLetterCode.compare("Y") == 0){
		threeLetter = "TYR";
	}
	else if(oneLetterCode.compare("V") == 0){
		threeLetter = "VAL";
	}

	return threeLetter;
}

/* Function to perform a standard alignment (no tweaked
 * options) between two sequences using CLUSTALW.
 *
 * Returns a vector of two sequences, now aligned with any 
 * gaps that CLUSTAL generated.
 *
 * Be sure that CLUSTALW 2.0.9 or better is in
 * your PATH before use.
 *
 * If CLUSTAL fails for any reason, returns a vector
 * of two empty strings.
 */
vector<string> clustalAlignment(string seq1, string seq2){
	string line;
	string aligned1 = "",  aligned2 = "";
	vector<string> alignedSeqs;

	/* Write sequences to file in FASTA format.
	 */
	ofstream seqInput("seqs");
	seqInput << ">seq1" << endl << seq1 << endl << ">seq2" << endl << seq2 << endl;
	seqInput.close();

	/* Run CLUSTALW and write to output file.
	 * Send all STDOUT and STDERR to the bit bucket.
	 */
	system("clustalw seqs -outfile=outaln > /dev/null 2>&1");

	/* Read the CLUSTAL results.
	 */
	ifstream resultFile("outaln");
	while( getline(resultFile, line) ){

		/* Parse aligned sequence 1
		 */
		if(line.length() > 4 && line.substr(0,4).compare("seq1") == 0){
			aligned1 += line.substr(16,60); // Just the sequence part
		}
		/* Parse aligned sequence 2
		 */
		else if(line.length() > 4 && line.substr(0,4).compare("seq2") == 0){
			aligned2 += line.substr(16,60); // Just the sequence part
		}

	}
	resultFile.close();

	/* Store them in a vector
	 */
	alignedSeqs.push_back(aligned1);
	alignedSeqs.push_back(aligned2);

	/* DEBUG
	 *
	cout << "Aligned sequences:" << endl << aligned1 << endl << aligned2 << endl;
	*/

	/* Remove the files.
	 */
	remove("seqs");
	remove("outaln");

	/* Return the sequences.
	 */
	return alignedSeqs;
}

/* Function to perform a non-gapped alignment between
 * two protein sequences using BLOSUM62 scores. If the 
 * sequences are of the same length, the raw score is 
 * returned. In any other case, the shorter sequence is 
 * slid along the longer to achieve a maximum score. The 
 * starting and ending gaps are not scored.
 *
 * Stores the starting and ending point of the
 * aligned region on the longer sequence of the two in
 * the supplied positions array.
 *
 * NOT QUITE KOSHER, BETTER APPROACHES EXIST THROUGH
 * DYNAMIC PROGRAMMING.
 */
int noGapAlignment(string seq1, string seq2, int (&positions)[2]){
	unsigned int i, j, resNum1, resNum2;
	string gappedSeq = "";
	int startPos = 0, endPos = 0;
	int bestScore = -32000, currentScore = 0;
	
	/* Only bother with the sliding procedure if the sequences
	 * are of unequal length.
	 */
	if(seq1.length() != seq2.length()){

		/* Sequence 1 is longer, slide sequence 2.
		 */
		if(seq1.length() > seq2.length()){

			/* Start with sequence 2 on the far left, gaps to the
			 * right, and move right until hitting the end of
			 * sequence 1.
			 */
			for(i = 0; i < (seq1.length() - seq2.length() + 1); i++){ // off-by-one?
	
				gappedSeq = seq2;

				/* Pad on the left
				 */
				for(j = 0; j < i; j++){
					gappedSeq = (string) "-" + gappedSeq;
				}

				/* Pad on the right
				 */
				for(j = 0; j < (seq1.length()-seq2.length()-i) ; j++){
					gappedSeq = gappedSeq + (string) "-";
				}

				/* DEBUG
				 *
				cout << "Padding attempt " << i+1 << ":" << endl << seq1 << endl << gappedSeq << endl;
				*/

				/* Compute score. Gaps don't count.
				 */
				for(j = 0; j < seq1.length(); j++){
					if(gappedSeq.substr(j,1).compare("-") != 0){
						resNum1 = blosum62Key.find(seq1[j]);
						resNum2 = blosum62Key.find(gappedSeq[j]);

						currentScore += blosum62Matrix[resNum1][resNum2];
					}
				}

				/* DEBUG
				 *
				cout << "Score: " << currentScore << endl;
				*/

				/* Compare with previous best score, save if
				 * better.
				 */
				if(currentScore > bestScore){
					bestScore = currentScore;
					// Avoid picking a position outside the string if no right padding exists
					startPos = max( (int) gappedSeq.find_first_of("ARNDCQEGHILKMFPSTWYVBZX") , 0 );
					endPos = min(gappedSeq.find_last_of("ARNDCQEGHILKMFPSTWYVBZX"), gappedSeq.length()-1); 
				}
			}

		}
		/* Sequence 2 is longer, slide sequence 1.
		 */
		else{
			/* Start with sequence 1 on the far left, gaps to the
			 * right, and move right until hitting the end of
			 * sequence 2.
			 */
			for(i = 0; i < (seq2.length() - seq1.length() + 1); i++){ // off-by-one?
				
				gappedSeq = seq1;

				/* Pad on the left
				 */
				for(j = 0; j < i; j++){
					gappedSeq = (string) "-" + gappedSeq;
				}

				/* Pad on the right
				 */
				for(j = 0; j < (seq2.length()-seq1.length()-i) ; j++){
					gappedSeq = gappedSeq + (string) "-";
				}

				/* DEBUG
				 *
				cout << "Padding attempt " << i+1 << ":" << endl << seq2 << endl << gappedSeq << endl;
				*/

				/* Compute score. Gaps don't count.
				 */
				for(j = 0; j < seq2.length(); j++){
					if(gappedSeq.substr(j,1).compare("-") != 0){
						resNum1 = blosum62Key.find(seq2[j]);
						resNum2 = blosum62Key.find(gappedSeq[j]);

						currentScore += blosum62Matrix[resNum1][resNum2];
					}
				}

				/* DEBUG
				 *
				cout << "Score: " << currentScore << endl;
				*/

				/* Compare with previous best score, save if necessary.
				 *
				 * FUCKING UP HERE: END POSITION IS CONSISTENLY 1 TOO
				 * HIGH.
				 *
				 * find_last_of() IS NOT RETURNING WHAT'S EXPECTED?
				 */
				if(currentScore > bestScore){

					/* DEBUG
					 *
					cout << "Updating score and positions." << endl;
					*/

					bestScore = currentScore;
					// Avoid picking a position outside the string if no padding exists
					startPos = max( (int) gappedSeq.find_first_of("ARNDCQEGHILKMFPSTWYVBZX") , 0 );
					endPos = min( gappedSeq.find_last_of("ARNDCQEGHILKMFPSTWYVBZX"), gappedSeq.length()-1 ); 
				}
			}
		}
	}

	/* Just calculate the raw score once for equilong
	 * sequences. The ends of the sequence will do
	 * nicely as a start and end.
	 */
	else{
		for(i = 0; i < seq1.length(); i++){
			resNum1 = blosum62Key.find(seq1[i]);
			resNum2 = blosum62Key.find(seq2[i]);

			bestScore += blosum62Matrix[resNum1][resNum2];
		}
		startPos = 0;
		endPos = seq1.length()-1;
	}

	/* Put the values in the output array.
	 */
	positions[0] = startPos;
	positions[1] = endPos;

	/* Return the score.
	 */
	return bestScore;
}

/* Function to compute the Hamming distance between two strings
 * of equal length.
 *
 * Returns the Hamming distance or -1 on error.
 */
int hammingDist(string stringOne, string stringTwo){
	int hdist = 0;
	unsigned int i;

	/* Strings are not of equal length: error.
	 */
	if(stringOne.length() != stringTwo.length()){
		return -1;
	}
	else{
		/* Compare each character in the strings
		 * and calculate the total Hamming distance.
		 */
		for(i = 0; i < stringOne.length(); i++){
			if(stringOne[i] != stringTwo[i]){
				hdist++;
			}
		}
	}

	return hdist;
}

/* Function to return all possible codons that encode a particular
 * amino acid in the standard code (specified as its one-letter
 * abbreviation in upper case). Stop codons are not differentiated 
 * from each other and should be specified as "_". "X" or unknown 
 * one-letter residue abbreviations result in an empty vector.
 */
vector<string> residueToCodons(char residue){

	vector<string> codons;

	/* Reverse of translation table adapated from 
	 * Tisdall, "Beginning Perl for Bioinformatics", 2001, 
	 * O'Reilly & Associates, Sebastopol, CA
	 */
	if(residue == 'S'){ // Serine
		codons.push_back("TCA");
		codons.push_back("TCC");
		codons.push_back("TCG");
		codons.push_back("TCT");
		codons.push_back("AGC");
		codons.push_back("AGT");
	}
	else if(residue == 'F'){ // Phenylalanine
		codons.push_back("TTC");
		codons.push_back("TTT");
	}
	else if(residue == 'L'){ // Leucine
		codons.push_back("TTA");
		codons.push_back("TTG");
		codons.push_back("CTA");
		codons.push_back("CTC");
		codons.push_back("CTG");
		codons.push_back("CTT");
	}
	else if(residue == 'Y'){ // Tyrosine
		codons.push_back("TAC");
		codons.push_back("TAT");
	}
	else if(residue == '_'){ // Stop
		codons.push_back("TAA");
		codons.push_back("TAG");
		codons.push_back("TGA");
	}
	else if(residue == 'C'){ // Cysteine
		codons.push_back("TGC");
		codons.push_back("TGT");
	}
	else if(residue == 'W'){
		codons.push_back("TGG");
	}
	else if(residue == 'P'){ // Proline
		codons.push_back("CCA");
		codons.push_back("CCC");
		codons.push_back("CCG");
		codons.push_back("CCT");
	}
	else if(residue == 'H'){ // Histidine
		codons.push_back("CAC");
		codons.push_back("CAT");
	}
	else if(residue == 'Q'){ // Glutamine
		codons.push_back("CAA");
		codons.push_back("CAG");
	}
	else if(residue == 'R'){ // Arginine
		codons.push_back("CGA");
		codons.push_back("CGC");
		codons.push_back("CGG");
		codons.push_back("CGT");
		codons.push_back("AGA");
		codons.push_back("AGG");
	}
	else if(residue == 'I'){ // Isoleucine
		codons.push_back("ATA");
		codons.push_back("ATC");
		codons.push_back("ATT");
	}
	else if(residue == 'M'){ // Methionine
		codons.push_back("ATG");
	}
	else if(residue == 'T'){ // Threonine
		codons.push_back("ACA");
		codons.push_back("ACC");
		codons.push_back("ACG");
		codons.push_back("ACT");
	}
	else if(residue == 'N'){ // Asparagine
		codons.push_back("AAC");
		codons.push_back("AAT");
	}
	else if(residue == 'K'){ // Lysine
		codons.push_back("AAA");
		codons.push_back("AAG");
	}
	else if(residue == 'V'){ // Valine
		codons.push_back("GTA");
		codons.push_back("GTC");
		codons.push_back("GTG");
		codons.push_back("GTT");
	}
	else if(residue == 'A'){ // Alanine
		codons.push_back("GCA");
		codons.push_back("GCC");
		codons.push_back("GCG");
		codons.push_back("GCT");
	}
	else if(residue == 'D'){ // Aspartic acid
		codons.push_back("GAC");
		codons.push_back("GAT");
	}
	else if(residue == 'E'){ // Glutamic acid
		codons.push_back("GAA");
		codons.push_back("GAG");
	}
	else if(residue == 'G'){ // Glycine
		codons.push_back("GGA");
		codons.push_back("GGC");
		codons.push_back("GGG");
		codons.push_back("GGT");
	}

	return codons;
}

/* Function to calculate the minimum number of
 * of DNA-level changes between two protein-coding 
 * sequences.
 *
 * Input should be upper-case sequences containing 
 * characters AGCT arranged in an integer number of 
 * codons.
 * 
 * Returns an integer number of required mutations,
 * or -1 on error. Doesn't deal realistically
 * with stop codons.
 */
int minDnaDist(string startDna, string endDna){

	unsigned int i, j, minDistance, distance;
	int numChanges = 0;
	string startCodon, endCodon;
	vector<string> codonsEnd;

	/* Are the sequences of equal length and
	 * an integer number of codons?
	 */
	if(startDna.length() % 3 != 0 || endDna.length() % 3 != 0 || startDna.length() != endDna.length()){
		return -1;
	}

	/* Translate sequences.
	*/
	string startProt = translateDna(startDna);
	string endProt = translateDna(endDna);

	// DEAL WITH STOP CODONS HERE ?

	/* Are they already completely synonymous?
	*/
	if(startProt.compare(endProt) == 0){
		return numChanges;
	}

	/* How far apart are the non-synonymous
	 * codons (minimum distance)?
	 */
	for(i = 0; i <= startDna.length()-3; i=i+3){
		startCodon = startDna.substr(i,3);
		endCodon = endDna.substr(i,3);

		/* For non-synonymous codons, calculate
		 * the minimum number of required changes 
		 * (order of changes doesn't matter).
		 */
		if(startProt[i/3] != endProt[i/3]){

			/* DEBUG
			 *
			cout << startProt[i/3] << " != " << endProt[i/3] << " (i.e. " << startCodon << " and " << endCodon << " translate to " << translateDna(startCodon) << " and " << translateDna(endCodon) << ")." << endl;
			*/

			codonsEnd = residueToCodons(endProt[i/3]);

			/* Find the shortest possible Hamming distance
			 * between the starting codon and one that's
			 * synonymous with the ending codon.
			 */
			minDistance = 3;
			for(j = 0; j < codonsEnd.size(); j++){
				distance = hammingDist(startCodon, codonsEnd[j]);

				/* DEBUG
				*
				cout << startCodon << " is " << distance << " steps from " << codonsEnd[j] << endl;
				*/

				if(distance < minDistance){
					minDistance = distance;

					/* DEBUG
					*
					cout << "New shortest distance: " << minDistance << endl;
					*/
				}
			}

			numChanges += minDistance;

			/* DEBUG
			 *
			cout << "Minimum number of changes so far: " << numChanges << endl;
			*/
		}
	}

	return numChanges;
}

/* Function to calculate a pairwise sequence distance
 * according go the BLOSUM62 matrix. Ignores gaps
 * completely. 
 *
 * Returns BLOSUM62 score or INT_MIN on error.
 */
int blosumDist(string protOne, string protTwo){
	int score = 0;
	unsigned int i, resNum1, resNum2;

	/* Strings are not of equal length: error.
	 */
	if(protOne.length() != protTwo.length()){
		return INT_MIN;
	}
	else{
		/* Compare each character in the strings
		 * and calculate the total BLOSUM62 score.
		 * Gaps are ignored.
		 */
		for(i = 0; i < protOne.length(); i++){
			if(protOne[i] != '-' && protTwo[i] != '-'){
				resNum1 = blosum62Key.find(protOne[i]);
				resNum2 = blosum62Key.find(protTwo[i]);
				score += blosum62Matrix[resNum1][resNum2];
			}
		}
	}

	return score;
}

/* Function to apply some number of non-synonymous mutations to
 * a DNA sequence. Each codon has an equal probability of being
 * mutated, and codons are drawn without replacement for each
 * mutation. Consequently a codon may not be mutated more than once,
 * and the maximum number of mutations applied is the
 * same as the number of codons. 
 *
 * Input DNA should have an integer number of codons.
 *
 * Returns an empty string on error.
 */
string mutateNonsyn(string inDna, int numMutations, CRandomMersenne &randGen){
	int i, codonNum, position, type;
	string outDna = inDna;
	string oldCodon, newCodon;
	vector<int> codonsUsed;

	/* Check for integer number of codons.
	 */
	if(inDna.length() % 3 != 0){
		return "";
	}

	/* Apply the indicated number of nonsynonymous
	 * mutations to the input sequence, if there
	 * are that many codons available.
	 */
	for(i = 0; i < numMutations; i++){

		/* Pick a random codon without replacement,
		 * until you run out of codons.
		 */
		if(codonsUsed.size() == inDna.length() / 3){

			/* DEBUG
			 *
			cout << "Made the maximum number of changes (" << (inDna.length() / 3) << "), stopping." << endl;
			*/

			return outDna;
		}
		else{
			codonNum = randGen.IRandom(0, (inDna.length()/3)-1);
			while( find(codonsUsed.begin(), codonsUsed.end(), codonNum) != codonsUsed.end() ){
				/* DEBUG
				 *
				cout << "Codon #" << (codonNum+1) << " already used, trying another..." << endl;
				*/

				codonNum = randGen.IRandom(0, (inDna.length()/3)-1);
			}
			codonsUsed.push_back(codonNum);

			/* Make a single non-synonymous change to it,
			 * at a random position in the codon.
			 *
			 * Most changes are non-synonymous, so just make
			 * random ones and check for changes in the
			 * translated residue.
			 */
			oldCodon = inDna.substr((codonNum*3),3);
			newCodon = outDna.substr((codonNum*3),3);
			while(translateDna(oldCodon).compare(translateDna(newCodon)) == 0){
				newCodon = outDna.substr((codonNum*3),3);
				position = randGen.IRandom(0,2); // All positions are equally likely
				type = randGen.IRandom(0,1); // Transitions and transversions are equally likely
				if(type == 0){
					newCodon.replace( position, 1, makeTransition(newCodon.substr(position,1)) );

					/* DEBUG
					 *
					cout << "Made a transition at position " << position+1 << " in codon " << codonNum+1 << " (" << oldCodon << " -> " << newCodon << ")" << endl;
					*/
				}
				else{
					newCodon.replace( position, 1, makeTransversion(newCodon.substr(position,1), randGen) );
					
					/* DEBUG
					 *
					cout << "Made a transversion at position " << position+1 << " in codon " << codonNum+1 << " (" << oldCodon << " -> " << newCodon << ")" << endl;
					*/
				}
			}
			outDna.replace((codonNum*3), 3, newCodon);
		}
	}

	return outDna;
}

/* Function to determine if two residues belong to the same biochemical
 * class or no. Classes are:
 *
 * Non-polar: A, V, M, I, P, L, G
 * Polar uncharged: S, T, C, N, Q
 * Aromatic: F, W, Y
 * Positive: R, K, H
 * Negative: D, E
 *
 * Unknown residues always return false.
 *
 * VERY UNPRETTY SET HANDLING HERE, DO UNORDERED MAPS
 * OR CHARACTER ARRAYS INSTEAD?
 */
bool sameBiochemClass(char residueOne, char residueTwo){
	bool same = false;
	string nonPolar = "AVMIPLG";
	string polar = "STCNQ";
	string aromatic = "FWY";
	string positive = "RKH";
	string negative = "DE";

	if( (nonPolar.find(residueOne) != string::npos && nonPolar.find(residueTwo) != string::npos)
	    ||
	    (polar.find(residueOne) != string::npos && polar.find(residueTwo) != string::npos)
	    ||
	    (aromatic.find(residueOne) != string::npos && aromatic.find(residueTwo) != string::npos)
	    ||
	    (positive.find(residueOne) != string::npos && positive.find(residueTwo) != string::npos)
	    ||
	    (negative.find(residueOne) != string::npos && negative.find(residueTwo) != string::npos)
	  ){
		same = true;
	}

	return same;
}

/* Position of individual upper-case alphanumerics in the
 * Greek alphabet. A (alpha) < B (beta) < G (gamma), etc.
 * Use only upper-case letter. Follows "beta code"
 * standard (see e.g. http://www.tlg.uci.edu/encoding/) 
 * but omits the "*" (implicit with only uppercase).
 * 
 * Returns the numerical position, 0-indexed. Currently
 * 
 * Returns -1 on error.
 */
int greekOrder(char latinLetter){
	int num = -1;

	switch(latinLetter){
		case 'A':
			num = 0;
			break;
		case 'B':
			num = 1;
			break;
		case 'G':
			num = 2;
			break;
		case 'D':
			num = 3;
			break;
		case 'E':
			num = 4;
			break;
		case 'Z':
			num = 5;
			break;
		case 'H':
			num = 6;
			break;
		case 'Q':
			num = 7;
			break;
		case 'I':
			num = 8;
			break;
		case 'K':
			num = 9;
			break;
		case 'L':
			num = 10;
			break;
		case 'M':
			num = 11;
			break;
		case 'C':
			num = 13;
			break;
		case 'O':
			num = 14;
			break;
		case 'P':
			num = 15;
			break;
		case 'R':
			num = 16;
			break;
		case 'S':
			num = 17;
			break;
		case 'T':
			num = 18;
			break;
		case 'U':
			num = 19;
			break;
		case 'F':
			num = 20;
			break;
		case 'X':
			num = 21;
			break;
		case 'Y':
			num = 22;
			break;
		case 'W':
			num = 23;
			break;
		default:
			break;
	}

	return num;
}

/* Provide a stable index converting residue names to integers.
 * Proceeds alphabetically from ALA = 0 to VAL = 19. Residue
 * identifier must be UPPERCASE, 3 letters.
 * 
 *  Returns -1 on error.
 */
int aaNum(string name){
	int number = -1;

	if(name.compare("ALA") == 0){
		number = 0;
	}
	else if(name.compare("ARG") == 0){
		number = 1;
	}
	else if(name.compare("ASN") == 0){
		number = 2;
	}
	else if(name.compare("ASP") == 0){
		number = 3;
	}
	else if(name.compare("CYS") == 0){
		number = 4;
	}
	else if(name.compare("GLN") == 0){
		number = 5;
	}
	else if(name.compare("GLU") == 0){
		number = 6;
	}
	else if(name.compare("GLY") == 0){
		number = 7;
	}
	else if(name.compare("HIS") == 0){
		number = 8;
	}
	else if(name.compare("ILE") == 0){
		number = 9;
	}
	else if(name.compare("LEU") == 0){
		number = 10;
	}
	else if(name.compare("LYS") == 0){
		number = 11;
	}
	else if(name.compare("MET") == 0){
		number = 12;
	}
	else if(name.compare("PHE") == 0){
		number = 13;
	}
	else if(name.compare("PRO") == 0){
		number = 14;
	}
	else if(name.compare("SER") == 0){
		number = 15;
	}
	else if(name.compare("THR") == 0){
		number = 16;
	}
	else if(name.compare("TRP") == 0){
		number = 17;
	}
	else if(name.compare("TYR") == 0){
		number = 18;
	}
	else if(name.compare("VAL") == 0){
		number = 19;
	}

	return number;
}
