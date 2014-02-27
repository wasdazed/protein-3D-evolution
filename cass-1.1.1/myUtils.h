/* FIGURE OUT SOME WAY OF PUTTING THIS IN THE SEARCH PATH FOR HEADERS,
 * JUST LIKE I HAVE FOR MY PERL PACKAGES.
 */

/* Make source file idempotent (avoid multiple includes).
 */
#ifndef MY_UTILS_CPP
#define MY_UTILS_CPP

/* Various includes.
 *
 * MOST OF THESE BELONG IN THE IMPLEMENTATION BIT...
 */
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>
#include "randomgen/randomc.h"
#include "randomgen/stocc.h"

using namespace std;

/* Definitions.
 */
//#define PI 3.1415926535
const double PI = atan(1.0)*4;

/* A point in N-dimensional space.
 */
typedef vector<double> coord;

/* Function to convert an integer to a string.
 *
 * Courtesy of Bjarne Soustrup. Or rather shamelessly stolen from 
 * http://www.research.att.com/~bs/bs_faq2.html#int-to-string
 */
string itos(int i);

/* Function to tokenize a string into a vector of
 * strings based on some delimiting string (defaults
 * to " ", i.e. a space character).
 *
 * Stolen from http://www.oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
 * and slightly modified.
 */
vector<string> tokenize(const string &str, const string &delimiters = " ");

/* REALLY SHITTY RANDOM INTEGER GENERATOR,
 * NEEDS IMPROVEMENT.
 *
 * STOLE FROM SHRUTI, WHO PROBABLY STOLE
 * IT ELSEWHERE.
 *
 * SHOULD PROBABLY BE REMOVED.
 */
int rand_int(int n);	// random number generator

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
int putInBin(float value, vector< vector<float> > binRanges);

/* Function to compute the geometric centroid of a collection of
 * points in N-dimensional Euclidian space.
 */
coord NDCentroid(vector<coord> points);

/* Function to compute the distance between two points in 
 * N-dimensional Euclidian space. The two points are
 * assumed to have the same dimensionality.
 */
double NDDist(coord point1, coord point2);

/* Function to compute a N-dimensional angle in Euclidean 
 * space, i.e. between two vectors u and v. The vectors are 
 * represented as 3 points: vector u goes from point 1 to 2, 
 * and vector v goes from point 2 to 3, resulting in an angle 
 * 1-2-3.
 *
 * Returns the angle in degrees.
 */
double NDAngle(coord point1, coord point2, coord point3);

/* Function to calculate the cross product of two vectors. 
 * Assumes a right-handed, orthogonal coordinate system, and
 * ignores dimensionality above 3.
 */
vector<double> crossProduct(coord vector1, coord vector2);

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
double dihedralAngle(coord point1, coord point2, coord point3, coord point4);

/* Function to find the range of a floating point vector,
 * here defined the same as the range of a distribution of
 * numbers (e.g. max(vector)-min(vector)).
 */
float vectorRange(vector<float> numbers);

/* Function to construct a vector from two points in N-dimensional
 * space assuming that point #1 is the center of the coordinate
 * system (i.e. we can represent the vector a set of coordinates).
 */
coord makeRelativeVector(coord point1, coord point2);

/* Function to compute the sum of a set of vectors going out
 * from the origin. Since they all start at the origin we
 * can represent them as coordinates in N-space.
 *
 * We assume that all vectors have the same dimensionality
 * as the first one (additional dimensions will not be
 * computed).
 */
coord vectorSum(vector<coord> vectors);

/* Scale the length of an N-dimensional vector by some factor.
 */
coord scaleVector(coord theVector, double scaleFactor);

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
coord rotateVector(coord theVector, vector<double> xyzRotationAngles);

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
coord reflectAroundPoint(coord theVector, coord point);

/* Function to calculate the mean of a vector of
 * floating-point numbers, using Knuth's one-pass
 * algorithm.
 * 
 * Returns -9999.0 on error.
 */
float mean(vector<float> distribution);

/* Function to calculation the variance of a vector of
 * floating-point numbers, using Knuth's one-pass
 * algorithm.
 * 
 * Returns -1.0 on error.
 */
float variance(vector<float> distribution);

/* Function to calculate the quartiles of a distribution
 * of floating point values.  Uses naive sorting 
 * implementation, doesn't deal well with NaN values
 * (exactly as badly as algorithm::sort(), in fact).
 *
 * Returns vector of (min, 25%, 50%, 75%, 100%) of
 * the distribution (i.e. minimum, 1Q, median, 3Q, maximum).
 */
vector<float> quartiles(vector<float> distribution);

/* Function to trim the whitespace from the start and
 * end of strings. Useful for e.g. parsing tasks.
 * Returns the trimmed string. DOES NOT trim tabs,
 * only space characters.
 */
string trimWhitespace(string input);

/* Function to make a random transversion (A -> [C/T], C -> [A/G],
 * G -> [C/T], T -> [A/G]) for a single nucleotide.
 */
string makeTransversion(string nucleotide, CRandomMersenne &randGenerator);

/* Function to make a transition (A <-> G, C <-> T) for a single nucleotide.
 *
 * On any input other than AGCT, returns "X".
 */
string makeTransition(string nucleotide);

/* Method for mutating DNA according to the Kimura 2-parameter
 * model (K80, with transitions being twice as probable as
 * transversions). Each nucleotide is mutated k times, with
 * k drawn from a Poisson distribution with a given lambda
 * (the average mutation rate per site per generation).
 */
string mutateDna(string dnaSequence, float lambda, StochasticLib1 &stochasticGenerator, CRandomMersenne &randGenerator);

/* Function for translating a DNA sequence into protein,
 * according to the standard (nuclear) code. Will only translate an
 * integer number of codons (extra nucleotides on the "end"
 * will not have any effect). Only takes upper-case input.
 * Stop codons are translated to "_", unknown codons to "?".
 *
 * Translation table from Tisdall, "Beginning Perl for Bioinformatics",
 * 2001, O'Reilly & Associates, Sebastopol, CA.
 */
string translateDna(string dnaSeq);

/* Generates a random sequence of characters of specified length from a
 * specified alphabet.
 */
string generateRandomSequence(unsigned int length, string alphabet, CRandomMersenne &randGen);

/* Generates a random sequence of characters by permuting the ones in
 * the existing string.
 *
 * Really just a wrapper for STL random_shuffle.
 */
string shuffleSequence(string inputString);

/* Function to convert three-letter amino acid names to
 * their one-letter abbreviation. All-upper-case name is
 * expected as input.
 * Returns "X" for an unrecognized amino acid.
 */
string oneLetterName(string residue);

/* Function to convert one-letter amino acid codes
 * to their three-letter representation. Essentially
 * inverts oneLetterName().
 *
 * Returns "XXX" on failure to recognize the one-letter
 * code.
 */
string threeLetterName(string oneLetterCode);

/* Function to perform a standard alignment (no tweaked
 * options) between two sequences using CLUSTALW.
 *
 * Returns a vector of two sequences, now aligned with any 
 * gaps that CLUSTAL generated.
 *
 * Be sure that CLUSTALW 2.0.9 or better is in
 * your PATH before use.
 */
vector<string> clustalAlignment(string seq1, string seq2);

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
int noGapAlignment(string seq1, string seq2, int (&positions)[2]);

/* Function to compute the Hamming distance between two strings
 * of equal length.
 *
 * Returns the Hamming distance or -1 on error.
 */
int hammingDist(string stringOne, string stringTwo);

/* Function to return all possible codons that encode a particular
 * amino acid in the standard code (specified as its one-letter
 * abbreviation in upper case). Stop codons are not differentiated 
 * from each other and should be specified as "_". "X" or unknown 
 * one-letter residue abbreviations result in an empty vector.
 */
vector<string> residueToCodons(char residue);

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
int minDnaDist(string startDna, string endDna);

/* Function to calculate a pairwise sequence distance
 * according go the BLOSUM62 matrix. Ignores gaps
 * completely. 
 *
 * Returns BLOSUM62 score or INT_MIN on error.
 */
int blosumDist(string protOne, string protTwo);

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
string mutateNonsyn(string inDna, int numMutations, CRandomMersenne &randGen);

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
 */
bool sameBiochemClass(char residueOne, char residueTwo);

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
int greekOrder(char latinLetter);

/* Provide a stable index converting residue names to integers.
 * Proceeds alphabetically from ALA = 0 to VAL = 19. Residue
 * identifier must be UPPERCASE, 3 letters.
 * 
 *  Returns -1 on error.
 */
int aaNum(string name);

/* Idempotency end.
 */
#endif
