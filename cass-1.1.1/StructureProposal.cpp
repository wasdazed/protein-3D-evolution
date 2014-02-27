#include "StructureProposal.h"
#include "TwoBeadStructure.h"
#include "GrahnenModel.h"
#include "ParameterRotation.h"
#include <iostream>

/* Default constructor, doesn't do much.
 */
StructureProposal::StructureProposal():
	Proposal()
{
}

/* Constructor from protein structure and desired
 * replacement sequence.
 *
 * THERE'S A GOOD ARGUMENT FOR PUTTING ALL THIS WORK
 * IN AN Init() METHOD: SEE http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml#Doing_Work_in_Constructors.
 *  
 * SO YOU CAN THROW AN EXCEPTION AND GET A QUASI-INFORMATIVE CRASH, AT LEAST.
 *
 * NOTE THAT startStruct IS NOW WHAT YOU HAVE _AFTER_ THE BEAD REPLACEMENT,
 * BUT _BEFORE_ POSITIONAL OPTIMIZATION.
 */
StructureProposal::StructureProposal(shared_ptr<TwoBeadStructure> inputStruct, string newSequence, float neighborDistance){

	int i, j, proteinLength;
	string oldSeq;
	double scaling;
	vector<double> noRotation(3,0.0);
	vector<bool> activeNums(inputStruct->getNumResidues(), false);
	vector< shared_ptr<Parameter> > rotPointers;
	vector<coord> newAtomCoordinates;
	coord oldAlpha, oldPreviousAlpha, oldNextAlpha, oldBeta, newBeta, firstBeta, secondBeta;
	coord nowhere(3, 9999.0);
	Residue oldResidue, firstResidue;

	oldSeq = inputStruct->getSequence();
	proteinLength = oldSeq.length();

	pStartStruct = shared_ptr<TwoBeadStructure>(new TwoBeadStructure(*inputStruct));	
	pCurrentStruct = shared_ptr<TwoBeadStructure>();
	targetSeq = newSequence;

	/* DEBUG
	 *
	cout << "Deciding which residues can move..." << endl;
	*/

	/* Decide which residues are "active" (allowed to move).
	 *
	 * Glycines never get to move (they don't have a CB
	 * atom anyway).
	 */
	for(i = 0; i < proteinLength; i++){

		oldResidue = pStartStruct->getSingleResidue(i);
		oldAlpha = oldResidue.getAtomCoordsByType("CA");
		oldBeta = oldResidue.getAtomCoordsByType("CB");

		/* Residue should be changed over if it's different,
		 * or if it's lacking a sidechain altogether.
		*/
		if(oldSeq[i] != newSequence[i] || (oldBeta == nowhere && oldSeq[i] != 'G')){

			/* Everything can move except glycines.
			 */
			if(newSequence[i] != 'G'){

				/* This residue can move.
				*/
				activeNums[i] = true;
			}

			/* Convert from a non-glycine to a non-glycine
			 * (if there are some coordinates to work with).
			 */
			if(oldSeq[i] != 'G' && newSequence[i] != 'G' && oldBeta != nowhere){
				/* DEBUG
				*
				cout << "Switching non-Gly " << oldSeq[i] << i+1 << " to non-Gly " << targetSeq[i] << i+1 << "." << endl;
				*/

				/* Adjust the bond length (Cb coordinates).
				 */
				scaling = GrahnenModel::getDefaultBondLength( targetSeq.substr(i,1) ) / GrahnenModel::getDefaultBondLength( oldSeq.substr(i,1) );

				newBeta = makeRelativeVector(oldAlpha, oldBeta); // Origin to oldAlpha
				newBeta = scaleVector(newBeta, scaling); // Scale bond length as needed
				newAtomCoordinates.push_back(oldAlpha);
				newAtomCoordinates.push_back(newBeta);
				newBeta = vectorSum(newAtomCoordinates); // Origin back to (0,0,0)
				newAtomCoordinates[1] = newBeta;
			}

			/* Convert from a glycine to a non-glycine, or from just
			 * backbone to something with a sidechain.
			 *
			 * Requires picking an orientation for the new Cb bead
			 * without having something old to go on.
			 */
			else if( (oldSeq[i] == 'G' && newSequence[i] != 'G') || oldBeta == nowhere){

				/* DEBUG
				*
				cout << "Switching " << oldSeq[i] << i+1 << " to " << targetSeq[i] << i+1 << ".";
				cout << "G->X transition, must invent starting Cb coordinates from scratch..." << endl;
				*/

				coord origin(3,0.0);

				/* Special case 1: We're at the N terminal.
				 * Place the bead pointing directly opposite
				 * to the Ca_i->Ca_i+1 vector (e.g. straight
				 * "backwards").
				 */
				if(i == 0){
					oldNextAlpha = pStartStruct->getSingleResidue(i+1).getAtomCoordsByType("CA");

					/* Move system to origin, centered on Ca_i.
					 */
					oldNextAlpha = makeRelativeVector(oldAlpha, oldNextAlpha);

					/* Reflect the Ca_i->Ca_i+1 vector around the
					 * origin.
					 */
					newBeta = reflectAroundPoint(oldNextAlpha, origin);

					/* Scale it appropriately.
					 */
					scaling = GrahnenModel::getDefaultBondLength( targetSeq.substr(i,1) ) / NDDist(origin, newBeta);
					newBeta = scaleVector(newBeta, scaling);

					/* Move system back.
					 */
					newAtomCoordinates.push_back(oldAlpha);
					newAtomCoordinates.push_back(newBeta);
					newBeta = vectorSum(newAtomCoordinates);
					newAtomCoordinates[1] = newBeta;
				}

				/* Special case 2: We're at the C terminal.
				 * Place the bead pointing in the same
				 * direction as the Ca_i-1->Ca_i vector,
				 * e.g. straight "forward".
				 */
				else if(i == proteinLength-1){
					oldPreviousAlpha = pStartStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");

					/* Move system to origin, centered on Ca_i.
					 */
					oldPreviousAlpha = makeRelativeVector(oldAlpha, oldPreviousAlpha);

					/* Move the start of the Ca_i-1->Ca_i vector
					 * to where the tip was, e.g. reflect it across
					 * the origin.
					 */
					newBeta = reflectAroundPoint(oldPreviousAlpha, origin);

					/* Scale it appropriately.
					 */
					scaling = GrahnenModel::getDefaultBondLength( targetSeq.substr(i,1) ) / NDDist(origin, newBeta);
					newBeta = scaleVector(newBeta, scaling);

					/* Move system back.
					 */
					newAtomCoordinates.push_back(oldAlpha);
					newAtomCoordinates.push_back(newBeta);
					newBeta = vectorSum(newAtomCoordinates);
					newAtomCoordinates[1] = newBeta;
				}

				/* We're somewhere in the middle of the protein.
				 * Place the bead pointing "out" from the internal
				 * angle between Ca_i-1, Ca_i and Ca_i+1.
				 */
				else{
					oldPreviousAlpha = pStartStruct->getSingleResidue(i-1).getAtomCoordsByType("CA");
					oldNextAlpha = pStartStruct->getSingleResidue(i+1).getAtomCoordsByType("CA");

					/* DEBUG
					 *
					cout << "Starting position of Ca_i: " << oldAlpha[0] << "," << oldAlpha[1] << "," << oldAlpha[2] << "." << endl;
					*/

					/* Move the whole system, centered on Ca_i,
					 * to the origin.
					 */
					oldPreviousAlpha = makeRelativeVector(oldAlpha, oldPreviousAlpha);
					oldNextAlpha = makeRelativeVector(oldAlpha, oldNextAlpha);

					/* DEBUG
					 *
					cout << "Set origin to Ca_i." << endl;
					*/

					/* Move the Ca_i-1->Ca_i vector to start
					 * at Ca_i (still points the same way),
					 * i.e. reflect it around the origin.
					 */
					newBeta = reflectAroundPoint(oldPreviousAlpha, origin);

					/* DEBUG
					 *
					cout << "Translated Ca_i-1->Ca_i vector tip to " << newBeta[0] << "," << newBeta[1] << "," << newBeta[2] << "." << endl;
					*/

					/* Add the Ca_i->Ca_i+1 vector, but reflected
					 * around Ca_i (simply the origin). The Ca_i->Cb vector now 
					 * points the right way but may not have the correct
					 * norm.
					 */
					newAtomCoordinates.push_back(newBeta);
					newAtomCoordinates.push_back(reflectAroundPoint(oldNextAlpha, origin));
					newBeta = vectorSum(newAtomCoordinates);

					/* DEBUG
					 *
					cout << "Added the reflected version of Ca_i->Ca_i+1, vector tip now at " << newBeta[0] << "," << newBeta[1] << "," << newBeta[2] << "." << endl;
					*/

					/* Scale the length of Cb vector as needed.
					 */
					scaling = GrahnenModel::getDefaultBondLength( targetSeq.substr(i,1) ) / NDDist(origin, newBeta);
					newBeta = scaleVector(newBeta, scaling);

					/* DEBUG
					 *
					cout << "Scaled vector norm, tip now at " << newBeta[0] << "," << newBeta[1] << "," << newBeta[2] << "." << endl;
					*/

					/* Move system back to original spot.
					 */	
					newAtomCoordinates[0] = oldAlpha;
					newAtomCoordinates[1] = newBeta;
					newBeta = vectorSum(newAtomCoordinates);
					newAtomCoordinates[1] = newBeta;
				
					/* DEBUG
					 *
					cout << "Final Cb coordinates: " << newBeta[0] << "," << newBeta[1] << "," << newBeta[2] << "." << endl;
					*/
				}
			}
			/* Convert from non-glycine to glycine.
			 */
			else{
				/* DEBUG
				*
				cout << "Switching non-Gly " << oldSeq[i] << i+1 << " to Gly " << targetSeq[i] << i+1 << "." << endl;
				*/

				/* No need to fiddle with the Cb bead: it doesn't exist anymore.
				 */
				newAtomCoordinates.push_back(oldAlpha);
			}

			/* DEBUG
			 *
			cout << "Constructing new residue " << targetSeq.substr(i,1) << i+1 << "..." << endl;
			*/

			/* Switch the residue. If there's no Cb atom in the original, this is
			 * taken care of automagically in the copy constructor.
			 *
			 * AUTOMAGIC IS APPARENTLY FAILING, CHECK ON THAT COPY CONSTRUCTOR.
			*/
			Residue tempResidue(oldResidue, threeLetterName(targetSeq.substr(i,1)), newAtomCoordinates);
			pStartStruct->setSingleResidue(i, tempResidue);
			newAtomCoordinates.clear();

			/* DEBUG
			 *
			cout << "Successfully processed this mutation." << endl;
			*/
		}
	}

	/* DEBUG
	 *
	cout << "Checking for neighbors that can move..." << endl;
	*/
	
	/* All the neighboring Cb beads (within some number of A) of moving residues can move, too.
	 *
	 * Again, no moving glycines.
	 */
	for(i = 0; i < proteinLength; i++){
		if(activeNums[i]){
			firstBeta = pStartStruct->getSingleResidue(i).getAtomCoordsByType("CB");

			for(j = 0; j < proteinLength; j++){

				if(i != j && newSequence[j] != 'G'){
					secondBeta = pStartStruct->getSingleResidue(j).getAtomCoordsByType("CB");

					/* Only consider Cb-Cb distances: Ca beads can't
					 * move anyway.
					 */
					if(NDDist(firstBeta, secondBeta) < neighborDistance){

						/* DEBUG
						 *
						cout << "Residue " << j+1 << " can also move around." << endl;
						*/

						activeNums[j] = true;
					}
				}
			}
		}
	}
	active = activeNums;

	/* Set the starting rotation to 0 on all
	 * axes. 
	 */
	for(i = 0; i < proteinLength; i++){
		if(active[i]){
			rotPointers.push_back( shared_ptr<ParameterRotation>(new ParameterRotation( noRotation )) );
		}
	}
	params = rotPointers;

	probDensity = calcProbDensity();
}

/* Constructor from a structure and a list of moving residues 
 *
 * NEEDS TO THROW EXCEPTIONS WHEN STRUCTURE AND VECTOR ARE NOT
 * OF EQUAL LENGTH.
 */
StructureProposal::StructureProposal(shared_ptr<TwoBeadStructure> inputStruct, vector<bool> movingResidues){
	int proteinLength, i;
	vector<double> noRotation(3,0.0);
	vector< shared_ptr<Parameter> > rotPointers;

	pStartStruct = shared_ptr<TwoBeadStructure>(new TwoBeadStructure(*inputStruct));
	pCurrentStruct = shared_ptr<TwoBeadStructure>();
	proteinLength = pStartStruct->getNumResidues();

	/* All marked residues should be active but unchanged.
	 */
	active = movingResidues;
	targetSeq = pStartStruct->getSequence();

	/* Set the starting rotation to 0 on all
	 * axes as necessary.
	 */
	for(i = 0; i < proteinLength; i++){
		if(active[i]){
			rotPointers.push_back( shared_ptr<ParameterRotation>(new ParameterRotation( noRotation )) );
		}
	}
	params = rotPointers;

	/* DEBUG
	 *
	cout << "Starting rotations of 0 assigned, calculating initial energy..." << endl;
	*/
	
	probDensity = calcProbDensity();
}

/* Constructor from parameter list. Pretty dangerous to
 * expose to user, so we don't.
 */
StructureProposal::StructureProposal(vector< shared_ptr<Parameter> > startParams, shared_ptr<TwoBeadStructure> scaffold, vector<bool> activeResidues, string novelSequence){

	params = startParams;

	/* Keep propagating the starting structure
	 * and which residues can move.
	 */
	pStartStruct = scaffold;
	pCurrentStruct = shared_ptr<TwoBeadStructure>();
	active = activeResidues;
	targetSeq = novelSequence;

	probDensity = calcProbDensity();
}

/* Virtual destructor, required for gcc to behave.
 */
StructureProposal::~StructureProposal(){
	/* DEBUG
	 *
	cout << "Destroying StructureProposal class properties..." << endl;
	*/

//	for(int i = 0; i < params.size(); i++){
//		/* DEBUG
//		 */
//		cout << "Killing parameter #" << i+1 << "..." << endl;
//		delete(params[i]);
//	}
//	cout << "StructureProposal: Parameters are dead." << endl;

}

/* Propose a new state from the current one
 * with varying-size steps.
*/
shared_ptr<Proposal> StructureProposal::nextStep(float stepSize){
	unsigned int i;
	vector< shared_ptr<Parameter> > newParams;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->nextValue(stepSize));
	}

	return shared_ptr<Proposal>(new StructureProposal(newParams, pStartStruct, active, targetSeq));
}

/* Propose a new state from the current one with
 * step sizes varying with some variance.
 */
shared_ptr<Proposal> StructureProposal::nextStep(float stepSize, float stepVar){
	unsigned int i;
	vector< shared_ptr<Parameter> > newParams;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->nextValue(stepSize, stepVar));
	}

	return shared_ptr<Proposal>(new StructureProposal(newParams, pStartStruct, active, targetSeq));
}

/* Propose a random state within parameter space.
 */
shared_ptr<Proposal> StructureProposal::randomStep(){
	unsigned int i;
	vector< shared_ptr<Parameter> > newParams;

	for(i = 0; i < params.size(); i++){
		newParams.push_back(params[i]->randomValue(0.0)); // Dummy value for implementation's sake
	}

	return shared_ptr<Proposal>(new StructureProposal(newParams, pStartStruct, active, targetSeq));
}

/* Function to update the old structure with the current
 * set of angle parameters.
 */
void StructureProposal::updateProteinStructure(){
	int i, num;
	vector<coord> twoCoords(2, coord(3,0.0));
	coord origin(3,0.0);
	coord oldAlpha, oldBeta, newBeta;
	double scaling, oldBondLength;
	Residue residueOld;
	
	pCurrentStruct = shared_ptr<TwoBeadStructure>(new TwoBeadStructure(*pStartStruct)); 

	/* Update the starting structure with the
	 * current parameter configuration.
	 */
	num = 0;
	for(i = 0; i < pCurrentStruct->getNumResidues(); i++){

		/* Move active residues around.
		 */
		if(active[i]){

			/* DEBUG
			 *
			cout << "Residue " << i+1 << " is allowed to move, changing coordinates..." << endl;
			*/

			/* Update the CB coordinates of the residue
			 * according to the rotation angles.
			 */
			residueOld = pCurrentStruct->getSingleResidue(i);
			oldBeta = residueOld.getAtomCoordsByType("CB");
			oldAlpha = residueOld.getAtomCoordsByType("CA");
			oldBondLength = NDDist(oldAlpha, oldBeta);

			/* DEBUG
			 *
			cout << "Old bond length for residue " << targetSeq.substr(i,1) << i+1 << ": " << oldBondLength << endl;
			*/

			newBeta = makeRelativeVector(oldAlpha, oldBeta); // Move origin of vector to 0,0,0
			newBeta = rotateVector(newBeta, dynamic_pointer_cast<ParameterRotation>(params[num])->getAngles()); // CHANGE TO DYNAMIC CAST?

			/* DEBUG
			 *
			cout << "After rotation: " << NDDist(origin, newBeta) << endl;
			*/

			scaling = oldBondLength / NDDist(origin, newBeta);
			newBeta = scaleVector(newBeta, scaling); // Re-normalize bond length after rotation

			/* DEBUG
			 *
			cout << "After re-normalization: " << NDDist(origin, newBeta) << endl;
			*/

			twoCoords[0] = oldAlpha;
			twoCoords[1] = newBeta;
			newBeta = vectorSum(twoCoords); // Reverse the previous move of origin

			/* DEBUG
			 *
			cout << "Old coords\tNew coords" << endl;
			for(int k = 0; k < oldBeta.size(); k++){
				cout << oldBeta[k] << "\t" << newBeta[k] << endl;
			}
			*/

			/* DEBUG
			 *
			cout << "Final new bond length for residue " << targetSeq.substr(i,1) << i+1 << ": " << NDDist(oldAlpha, newBeta) << endl;
			*/

			/* Modify the residue based on the new coordinates.
			 */

			// ISSUE OCCURS RIGHT HERE IF ATOM TYPE NOT DEFINED
			// FOR CB.
			Residue tempResidue(residueOld, "CB", newBeta);
			pCurrentStruct->setSingleResidue(i, tempResidue);

			num++;
		}
	}
}

/* Returns the current structure (based on angle
 * parameters).
 */
shared_ptr<TwoBeadStructure> StructureProposal::getStructure() const{
	return pCurrentStruct;
}

float StructureProposal::calcProbDensity(){
	return calcProbDensity(USE_SIMPLE_ENERGY_FUNCTION);
}

/* Function for calculating the probability density, i.e.
 * SOME_METRIC, for the current parameter set.
 */
float StructureProposal::calcProbDensity(bool useSimpleMethod){

	float energy = 0;

	/* DEBUG
	 *
	cout << "Updating the structure Cartesian coordinates..." << endl;
	*/

	/* Calculate what the current structure would look
	 * like.
	 */
	updateProteinStructure();

	/* DEBUG
	 *
	cout << "Calculating energy function..." << endl;
	*/

	/* Now calculate the relevant energy function for this
	 * structure.
	 */
	if(useSimpleMethod){
		energy = GrahnenModel::scwrlRepulsionEnergy(pCurrentStruct, active);
	}
	else{
		energy = GrahnenModel::LJRepulsionEnergy(pCurrentStruct, active);
	}

	/* DEBUG
	 *
	cout << "Energy is " << energy << endl;
	*/

	/* Return this as something negative, so we
	 * can optimize for larger values (e.g. closer to 0).
	 */
	return -1.0*energy;
}
