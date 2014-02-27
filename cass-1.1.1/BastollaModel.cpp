#include "BastollaModel.h"
#include "Structure.h"
#include "AllAtomStructure.h"
#include "TwoBeadStructure.h"
#include <typeinfo>
#include <stdexcept>

/* Static member: Bastolla's residue-residue interaction energy
 * matrix. All energies are in units of kT.
 */
float BastollaModel::bastollaMatrix[20][20] = 
{
{-0.0479, 0.1346, -0.0457, 0.1018, 0.1049, -0.1711, 0.1844, 0.0691, 0.0464, -0.1431, 0.1049, 0.0310, 0.1462, -0.0737, -0.0847, -0.1119, -0.1408, -0.1085, -0.0880, 0.0266},
{0.1346 , 0.1259, -0.0413, 0.1581, 0.1146, 0.0802, 0.2311, -0.2403, 0.0823, 0.1010, -0.3511, 0.0675, 0.2241, 0.1103, 0.0637, 0.0885, -0.0522, 0.1550, -0.0967, -0.0827},
{-0.0457 -0.0413, -0.0550, -0.0728, -0.0050, -0.0172, 0.1710, -0.0735, 0.1169, 0.1061, 0.0059, -0.0243, 0.1127, -0.0480, -0.1038, -0.0171, -0.1431, 0.0715, -0.0540, -0.0125},
{0.1018, 0.1581, -0.0728, 0.0840, 0.0192, 0.2673, 0.1115, -0.1154, 0.0424, 0.2728, -0.1859, 0.1043, 0.2386, 0.1892, -0.0197, 0.0827, -0.1165, 0.1169, -0.0124, -0.0749},
{0.1049, 0.1146, -0.0050, 0.0192, -0.0917, 0.0890, 0.1196, -0.0381, 0.1452, 0.1180, -0.0150, 0.0155, 0.1560, 0.1485, 0.0124, 0.0018, -0.1149, -0.0844, -0.0250, 0.0386},
{-0.1711, 0.0802, -0.0172, 0.2673, 0.0890, -0.5067, 0.0782, 0.0543, 0.0959, -0.4593, -0.0651, -0.0316, 0.0745, -0.5112, -0.1822, -0.5450, -0.2614, -0.1305, -0.2639, -0.0169},
{0.1844, 0.2311, 0.1710, 0.1115, 0.1196, 0.0782, 0.2219, 0.1963, 0.1075, 0.1859, -0.0251, 0.1763, 0.2131, 0.1174, -0.0573, 0.0789, -0.0176, -0.0982, -0.1567, 0.0979},
{0.0691, -0.2403, -0.0735, -0.1154, -0.0381, 0.0543, 0.1963, 0.1216, 0.1690, 0.0609, 0.0839, 0.0467, 0.1099, 0.0682, 0.0866, -0.0416, -0.1120, -0.0330, -0.1152, 0.0390},
{0.0464, 0.0823, 0.1169, 0.0424, 0.1452, 0.0959, 0.1075, 0.1690, 0.0941, 0.1766, 0.0442, 0.0228, 0.1626, 0.0332, 0.0185, 0.0398, 0.0214, -0.0132, -0.0802, -0.0005},
{-0.1431, 0.1010, 0.1061, 0.2728, 0.1180, -0.4593, 0.1859, 0.0609, 0.1766, -0.5193, 0.0475, 0.0119, 0.0868, -0.4223, -0.2127, -0.4001, -0.2792, -0.2349, -0.2898, -0.0039},
{0.1049, -0.3511, 0.0059, -0.1859, -0.0150, -0.0651, -0.0251, 0.0839, 0.0442, 0.0475, 0.0306, -0.0210, -0.0614, -0.0266, -0.0163, -0.0904, -0.1369, 0.0544, -0.2070, -0.0184},
{0.0310, 0.0675, -0.0243, 0.1043, 0.0155, -0.0316, 0.1763, 0.0467, 0.0228, 0.0119, -0.0210, 0.0150, 0.1908, -0.0700, 0.0018, -0.1120, -0.1445, -0.0013, 0.0052, 0.0681},
{0.1462, 0.2241, 0.1127, 0.2386, 0.1560, 0.0745, 0.2131, 0.1099, 0.1626, 0.0868, -0.0614, 0.1908, 0.1077, 0.0882, -0.0069, -0.0604, -0.1326, 0.0545, -0.0910, 0.0295},
{-0.0737, 0.1103, -0.0480, 0.1892, 0.1485, -0.5112, 0.1174, 0.0682, 0.0332, -0.4223, -0.0266, -0.0700, 0.0882, -0.5852, -0.2137, -0.3791, -0.3164, -0.2235, -0.1961, -0.0326},
{-0.0847, 0.0637, -0.1038, -0.0197, 0.0124, -0.1822, -0.0573, 0.0866, 0.0185, -0.2127, -0.0163, 0.0018, -0.0069, -0.2137, -0.1059, -0.1785, -0.1621, -0.0557, -0.0775, -0.0345},
{-0.1119, 0.0885, -0.0171, 0.0827, 0.0018, -0.5450, 0.0789, -0.0416, 0.0398, -0.4001, -0.0904, -0.1120, -0.0604, -0.3791, -0.1785, -0.3088, -0.4212, -0.3262, -0.3405, -0.1250},
{-0.1408, -0.0522, -0.1431, -0.1165, -0.1149, -0.2614, -0.0176, -0.1120, 0.0214 , -0.2792, -0.1369, -0.1445, -0.1326, -0.3164, -0.1621, -0.4212, -0.2793, -0.2444, -0.3209, -0.1976},
{-0.1085, 0.1550 , 0.0715 , 0.1169 , -0.0844, -0.1305, -0.0982, -0.0330, -0.0132, -0.2349, 0.0544 , -0.0013, 0.0545 , -0.2235, -0.0557, -0.3262, -0.2444, -1.0442, -0.1176, -0.0701},
{-0.0880, -0.0967, -0.0540, -0.0124, -0.0250, -0.2639, -0.1567, -0.1152, -0.0802, -0.2898, -0.2070, 0.0052 , -0.0910, -0.1961, -0.0775, -0.3405, -0.3209, -0.1176, -0.1066, -0.0200},
{0.0266, -0.0827, -0.0125, -0.0749, 0.0386, -0.0169, 0.0979, 0.0390, -0.0005, -0.0039, -0.0184, 0.0681, 0.0295, -0.0326, -0.0345, -0.1250, -0.1976, -0.0701, -0.0200, 0.0005}
};
string BastollaModel::bastollaKey = "AEQDNLGKSVRTPIMFYCWH";

/* Default constructor, doesn't do much.
 */
BastollaModel::BastollaModel():
	EnergyModel()
{
}

/* Destructor, required due to gcc being
 * finicky.
 */
BastollaModel::~BastollaModel(){
}

/* Calculating folding score/energy.
 *
 * NUKE SPEED CALCULATIONS AS NECESSARY.
 *
 * REQUIRE AN ALL-ATOM MODEL HERE?
 */
float BastollaModel::calcFoldScore(shared_ptr<Structure> protein){

	/* Crash if you're being passed the wrong type
	 * of pointer.
	 */
	AllAtomStructure test;
	if(typeid(test) != typeid(*protein)){
		throw runtime_error("Only AllAtomStructure is valid for BastollaModel.");
	}

	/* Start the clock.
	*
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);
	*/


	float eng = 0.0;
	string res1, res2;
	unsigned int i, j, resNum1, resNum2;

	/* Construct a contact map at 4.5 A between heavy
	 * atoms.
	 *
	 * NOT VERY EFFICIENT. CAN THIS BE OPTIMIZED AWAY?
	 */
	vector< vector<int> > contactMap = protein->contactMap(4.5);
	string sequence = protein->getSequence();

	/* For residues separated by 3 or more positions in the
	 * structure and in contact, add their interaction energy
	 * to the total.
	 */
	for(i = 0; i < sequence.length(); i++){
		for(j = i + SEQ_SEPARATION; j < sequence.length(); j++){ 

			/* DEBUG
			 *
			 cout << sequence.substr(i,1) << i+1 << "-" << sequence.substr(j,1) << j+1 << endl;
			 */

			if(contactMap[i][j] == 1){

				/* Find the correct position in the Bastolla matrix for
				 * each residue.
				 */
				res1 = sequence.substr(i,1);
				res2 = sequence.substr(j,1);
				resNum1 = bastollaKey.find(res1);
				resNum2 = bastollaKey.find(res2);

				/* DEBUG
				 *
				 cout << res1 << " has position " << resNum1 << " in the matrix (start at 0)." << endl;
				 cout << res2 << " has position " << resNum2 << " in the matrix (start at 0)." << endl;
				 */

				/* Retrive the contact energy and add to
				 * total.
				 */
				eng += bastollaMatrix[resNum1][resNum2];

				/* DEBUG
				 *
				 cout << "Eng = " << bastollaMatrix[resNum1][resNum2] << endl;
				 */

			}
		}
	}

	/* Stop the clock and print.
	*
	gettimeofday(&tim, NULL);
	double timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
	cout << (timeStop - timeStart) << endl;
	*/

	return eng;
}

// DUMMY IMPLEMENTATION FOR NOW...
float BastollaModel::calcFoldEnergyGap(shared_ptr<Structure> protein){
	return calcFoldScore(protein);
}

// DUMMY IMPLEMENTATION FOR NOW...
float BastollaModel::calcFoldEnergyGap(shared_ptr<Structure> protein, shared_ptr<Structure> compactConformation){
	return calcFoldScore(protein);
}

/* Calculating interaction score/energy.
 *
 * Based on assumption that Bastolla's contact energies
 * for residues work for inter- as well as intra-protein
 * contacts. Note that the model was NOT originially
 * parameterized for this purpose.
 */
float BastollaModel::calcInteractionScore(shared_ptr<Structure> interactorOne, shared_ptr<Structure> interactorTwo){
	
	/* Crash if you're being passed the wrong type
	 * of pointer(s).
	 */
	AllAtomStructure test;
	if(typeid(*interactorOne) != typeid(*interactorTwo)){
		throw runtime_error("Interacting structures must be of the same type.");
	}
	else{
		if( 	typeid(test) != typeid(*interactorOne) 
			|| typeid(test) != typeid(*interactorTwo)
		  ){
			throw runtime_error("Only AllAtomStructure is valid for BastollaModel.");
		}
	}

	/* Start the clock.
	 *
	timeval tim;
	gettimeofday(&tim, NULL);
	double timeStart = tim.tv_sec + (tim.tv_usec/1000000.0);
	*/

	float bindEnergy = 0.0;
	int i, j;
	string res1, res2;
	int resNum1, resNum2;

	/* Examine contact status (at 4.5 A) of all the residues in one protein
	 * to all the residues in the other one, and calculate the 
	 * binding energy based on that.
	 */
	for(i = 0; i < interactorOne->getNumResidues(); i++){
		Residue resOne = interactorOne->getSingleResidue(i);

		for(j = 0; j < interactorTwo->getNumResidues(); j++){
			Residue resTwo = interactorTwo->getSingleResidue(j);

			/* Calculate the energy between all residues that are within
			 * 4.5 A of each other.
			 */
			if( resOne.isCloserThan(4.5, resTwo) ){

				/* Find the correct position in the Bastolla matrix for
				 * each residue.
				 */
				res1 = resOne.getOneLetterType();
				res2 = resTwo.getOneLetterType();
				resNum1 = bastollaKey.find(res1);
				resNum2 = bastollaKey.find(res2);

				/* Add this to the binding energy.
				*/
				bindEnergy += bastollaMatrix[resNum1][resNum2];
			}
		}
	}

	/* Stop the clock and print.
	 *
	gettimeofday(&tim, NULL);
	double timeStop = tim.tv_sec + (tim.tv_usec/1000000.0);
	cout << (timeStop - timeStart) << endl;
	*/

	return bindEnergy;
}
