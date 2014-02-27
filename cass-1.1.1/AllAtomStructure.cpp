#include "AllAtomStructure.h"
#include <sys/stat.h>
#include <sys/wait.h>
#include <signal.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdexcept>
#include <unistd.h>
		
AllAtomStructure::AllAtomStructure():
	Structure()
{
}


AllAtomStructure::AllAtomStructure(vector<Residue> residues):
	Structure(residues)
{
}

AllAtomStructure::~AllAtomStructure(){
	/* DEBUG
	 *
	cout << "Destroying AllAtomStructure class properties..." << endl;
	*/
}

/* Function to extract the backbone trace from a structure.
 * See Structure::getBackboneResidues() for more detail.
 */
shared_ptr<Structure> AllAtomStructure::getBackbone() const{
	vector<Residue> residues = this->getBackboneResidues();

	return shared_ptr<Structure>(new AllAtomStructure(residues));
}

/* Function to parse SCWRL3 output into a Structure object. Provide
 * an opened ifstream as an argument. Parses the first chain or
 * model present in the input file, as terminated by "TER" or
 * "ENDMDL".
 *
 * Throws runtime_error exception on reading failure.
 */
AllAtomStructure AllAtomStructure::parseScwrlOutput(ifstream &inputFile){

	if(!inputFile.good()){
			string message = "Can't read from structure input file!";
			throw runtime_error(message);
	}


	string line;
	string pdbString = "";

	/* Read all the ATOM lines, delimited by "\n", into a
	 * string, stopping at the end of the chain or model.
	 */
	while( getline(inputFile, line) ){

		if(line.substr(0,4).compare("ATOM") == 0){
			pdbString += line;
			pdbString += "\n";
		}
		else if(line.substr(0,3).compare("TER") == 0 || line.substr(0,6).compare("ENDMDL") == 0){
			break;
		}
	}

	/* Create object from the parsed string
	 * and return it.
	 */
	return AllAtomStructure::newFromPDBString(pdbString);
}

/* Function to create a Structure from a string containing
 * only ATOM records separated by "\n".
 *
 * POSSIBLY BORKEN WITH SCWRL 4.0, FIX.
 */
AllAtomStructure AllAtomStructure::newFromPDBString(string pdbLines){

	string atomType, resType, lastResType;
	
	unsigned int resNum, i;
	int lastResNum = -9999;
	
	float x, y, z;
	
	coord position;
	vector<coord> coords;
	vector<string> types;
	vector<Residue> residueCollection;
	vector<string> allLines;

	/* Tokenize the string as lines.
	 */
	allLines = tokenize(pdbLines, "\n");

	/* Read all the lines in input string.
	 */
	for(i = 0; i < allLines.size(); i++){

		/* Extract atom type.
		 */
		atomType = trimWhitespace(allLines[i].substr(12,4));

		/* DEBUG
		 *
		 cout << "Current atom type: '" << atomType << "'" << endl;
		 */

		/* Extract and convert residue number.
		*/
		resNum = atoi(allLines[i].substr(22,4).c_str());

		/* DEBUG
		 *
		 cout << "Residue number: " << resNum << endl;
		 */

		/* If this is the start of a new residue, make a Residue of 
		 * the old values and stuff it in the list. First residue is
		 * a bit of a special case.
		 *
		 * TEMPORARY FIX: DON'T LOOK FOR OXT TERMINATOR, SCWRL4 OCCASIONALLY
		 * PRODUCES THIS ELSEWHERE.
		 */
		//if( resNum != lastResNum || atomType.compare("OXT") == 0){
		if( resNum != lastResNum && lastResNum != -9999){

			/* DEBUG
			 *
			 cout << "Hit a new residue or the end of the chain, time to store the old one." << endl;
			 */

			Residue tempResidue(lastResType, types, coords);
			residueCollection.push_back(tempResidue);
			types.clear();
			coords.clear();
		}

		/* Get the residue type we're working on.
		*/
		resType = allLines[i].substr(17,3);

		/* DEBUG
		 *
		 cout << "Current residue type: '" << resType << "'" << endl;
		 */

		/* Extract and convert the coordinates.
		*/
		x = atof(allLines[i].substr(30,8).c_str());
		y = atof(allLines[i].substr(38,8).c_str());
		z = atof(allLines[i].substr(46,8).c_str());

		/* DEBUG
		 *
		 cout << "Parsed coordinates: " << x << ", " << y << ", " << z << endl;
		 */

		/* Make a single 3D coordinate of them.
		*/
		position.push_back(x);
		position.push_back(y);
		position.push_back(z);

		/* Add the atom to the lists.
		*/
		coords.push_back(position);
		types.push_back(atomType);
		position.clear();

		lastResNum = resNum;
		lastResType = resType;

		/* DEBUG
		 *
		cout << "Read ATOM record " << (i+1) << "." << endl;
		*/
	}

	/* Save the last residue here (since there's no change in numbering
	 * after it).
	 */
	Residue tempResidue(lastResType, types, coords);
	residueCollection.push_back(tempResidue);

	/* DEBUG
	 *
	cout << "Saved the last residue (# " << residueCollection.size() << ")." << endl;
	*/

	/* Create Structure object.
	 */
	AllAtomStructure structure(residueCollection);

	/* DEBUG
	 *
	cout << "Sequence of parsed structure: " << structure.getSequence() << endl;
	*/

	/* And return it.
	 */
	return structure;
}

/* Threads the supplied sequence through the current structure, 
 * tweaks the sidechain orientation using SCWRL 3.0 (Canutescu et al,
 * 2003, Protein Science 12) and returns the tweaked structure.
 *
 * SCWRL is sensitive to upper- and lower-case characters in your
 * sequence -- see that paper for more information.
 * 
 * Random number generator is used to produce temporary file names.
 *
 * The progress of the SCWRL run is checked every SCWRL_CHECK_DELAY
 * microseconds (by checking if the output file exists), and if it doesn't 
 * terminate within SCWRL_MAX_TIME seconds, kills the process and
 * returns a null pointer.
 *
 * SWITCH TO SCWRL 4.0.
 */
shared_ptr<Structure> AllAtomStructure::threadSequence(string sequence){
	
	bool scwrlFinished = false;
	int i, statResult, randInt, waitStatus, killStatus;
	struct stat statInfo;

	/* Generate a random number for this series
	 * of files.
	 */
	randInt = rand();

	/* Write the sequence to disk.
	*/
	string seqFile = "newseq"+itos(randInt);
	ofstream newSeq(seqFile.c_str());
	newSeq << sequence;
	newSeq.close();

	/* Write this structure to disk.
	 */
	string structFileOrig = "origstruct"+itos(randInt);
	ofstream origStruct(structFileOrig.c_str());
	origStruct << toPDBString();
	origStruct.close();

	string structFileNew = "newstruct"+itos(randInt);
	
	/* Run SCWRL with the sequence and the original structure,
	 * forked off as a child process.
	 *
	 * THIS THOROUGHLY BREAKS WINDOWS COMPATABILITY. BUT FUCK
	 * 'EM -- MACS SHOULD STILL WORK.
	 */
	pid_t forkPid = fork();

	/* The child runs SCWRL.
	*/
	if(forkPid == 0){
		/* DEBUG
		 *
		cout << "Child process, running SCWRL." << endl;
		*/

		/* Open a pointer to /dev/null (write-only).
		 */
		int filePtr = open("/dev/null", O_WRONLY);

		/* Duplicate the STDOUT and STDERR pointers into that 
		 * pointer (redirects them to /dev/null). STDOUT is
		 * always pointer "1". STDIN is "0", STDERR is
		 * "2".
		 */
		dup2(filePtr, 1);
		dup2(filePtr, 2);

		/* Close the file pointer.
		 */
		close(filePtr);

		/* Run SCWRL with arguments.
		 */
		execlp(SCWRL_EXECUTABLE, SCWRL_EXECUTABLE, "-i", structFileOrig.c_str(), "-o", structFileNew.c_str(), "-s", seqFile.c_str(), SCWRL_OPTIONS, NULL);
	}
	/* Fork failed, quit.
	 */
	else if(forkPid < 0){
		cout << "Could not fork off child process for threading!" << endl;
		cout << "Error: " << strerror(errno) << endl;
		exit(1);
	}
	/* The parent checks for SCWRL finishing every once in a while.
	*/
	else{
		/* DEBUG
		 *
		cout << "Parent process, checking on SCWRL." << endl;
		*/

		/* Check on the output file every SCWRL_CHECK_DELAY microseconds (or so), 
		 * for SCWRL_MAX_TIME seconds, to see if it's there or not.
		 */
		for(i = 0; i < SCWRL_MAX_TIME*(1000000/SCWRL_CHECK_DELAY); i++){

			/* Sleep for SCWRL_CHECK_DELAY microseconds.
			*/
			usleep(SCWRL_CHECK_DELAY);

			/* Check if the SCWRL output file exists.
			 * If it does, we're finished.
			 */
			statResult = stat(structFileNew.c_str(), &statInfo);
			if(statResult == 0){
				scwrlFinished = true;
				wait(&waitStatus); // Be sure that SCWRL finishes before moving on

				/* DEBUG
				 *
				cout << "SCWRL processed finished and was reaped, moving on." << endl;
				*/

				break;
			}

			/* DEBUG
			 *
			 cout << i*(SCWRL_CHECK_DELAY/1000000.0) << " seconds." << endl;
			 */
		}

		/* Kill the child process, if it's still running, and
		 * reap the resulting zombie.
		 */
		if(!scwrlFinished){
			/* DEBUG
			 *
			cout << "SCWRL process did not finish on time, trying to kill and reap..." << endl;
			*/

			killStatus = kill(forkPid, 15);

			/* DEBUG
			 *
			cout << "Issued kill()." << endl;
			*/

			wait(&waitStatus);

			/* DEBUG
			 *
			cout << "Wait()ed." << endl;
			*/
		}

		/* DEBUG
		 *
		 cout << "SCWRL process dead and gone." << endl;
		 */
	}

	/* Create a new empty structure.
	 */
	shared_ptr<Structure> newStructure;
	
	/* Read the novel structure into memory if SCWRL finished
	 * on time.
	 */
	if(scwrlFinished){

		/* DEBUG
		*
		cout << "SCWRL finished after " << i*(SCWRL_CHECK_DELAY/1000000.0) << " seconds, reading the new file." << endl;
		*/

		ifstream newStruct(structFileNew.c_str());
		newStructure = shared_ptr<Structure>(new AllAtomStructure(AllAtomStructure::parseScwrlOutput(newStruct)));
		newStruct.close();

		remove(structFileNew.c_str());
	}
	else{
		/* DEBUG
		 *
		cerr << "SCWRL timed out." << endl;
		*/
	}

	/* DEBUG
	 *
	cout << "Done with reading/not reading novel structure." << endl;
	*/

	/* Delete the temporary files.
	 */
	remove(seqFile.c_str());
	remove(structFileOrig.c_str());

	/* DEBUG
	 *
	cout << "Temporary files deleted." << endl;
	*/

	/* Return the tweaked (or empty) structure.
	 */
	return newStructure;
}

/* Function to run SCWRL on the current structure, tweaking the
 * position/angle of the sidechains.
 *
 * Returns an adjusted Structure or a null pointer on error.
 *
 * See threadSequence() for more details.
 *
 * SWITCH TO SCWRL 4.0.
 */
shared_ptr<Structure> AllAtomStructure::adjustSidechains(){
	bool scwrlFinished = false;
	int i, statResult, randInt, waitStatus, killStatus;
	struct stat statInfo;
	
	/* Generate a random number for this series
	 * of files.
	 */
	randInt = rand();
	
	/* Write this structure to disk.
	 */
	string structFileOrig = "origstruct"+itos(randInt);
	ofstream origStruct(structFileOrig.c_str());
	origStruct << this->toPDBString();
	origStruct.close();

	string structFileNew = "newstruct"+itos(randInt);

	/* Run SCWRL with the original structure, forked off as a child process.
	 *
	 * THIS THOROUGHLY BREAKS WINDOWS COMPATABILITY. BUT FUCK
	 * 'EM -- MACS SHOULD STILL WORK.
	 */
	pid_t forkPid = fork();

	/* The child runs SCWRL.
	*/
	if(forkPid == 0){
		/* DEBUG
		 *
		cout << "Child process, running SCWRL." << endl;
		*/

		/* Open a pointer to /dev/null (write-only).
		 */
		int filePtr = open("/dev/null", O_WRONLY);

		/* Duplicate the STDOUT and STDERR pointers into that 
		 * pointer (redirects them to /dev/null). STDOUT is
		 * always pointer "1". STDIN is "0", STDERR is
		 * "2".
		 */
		dup2(filePtr, 1);
		dup2(filePtr, 2);

		/* Close the file pointer.
		 */
		close(filePtr);

		/* Run SCWRL with arguments.
		 */
		execlp(SCWRL_EXECUTABLE, SCWRL_EXECUTABLE, "-i", structFileOrig.c_str(), "-o", structFileNew.c_str(), SCWRL_OPTIONS, NULL);
	}
	/* Fork failed, quit.
	 */
	else if(forkPid < 0){
		cout << "Could not spawn child process!" << endl;
		cout << "Error: " << strerror(errno) << endl;
		exit(1);
	}
	/* The parent checks for SCWRL finishing every once in a while.
	*/
	else{
		/* DEBUG
		 *
		cout << "Parent process, checking on SCWRL." << endl;
		*/

		/* Check on the output file every SCWRL_CHECK_DELAY microseconds (or so), 
		 * for SCWRL_MAX_TIME seconds, to see if it's there or not.
		 */
		for(i = 0; i < SCWRL_MAX_TIME*(1000000/SCWRL_CHECK_DELAY); i++){

			/* Sleep for SCWRL_CHECK_DELAY microseconds.
			*/
			usleep(SCWRL_CHECK_DELAY);

			/* Check if the SCWRL output file exists.
			 * If it does, we're finished.
			 */
			statResult = stat(structFileNew.c_str(), &statInfo);
			if(statResult == 0){
				scwrlFinished = true;
				wait(&waitStatus); // Make sure SCWRL really finishes
				break;
			}

			/* DEBUG
			 *
			 cout << i*(SCWRL_CHECK_DELAY/1000000.0) << " seconds." << endl;
			 */

		}

		/* Kill the child process, if it's still running, and
		 * reap the resulting zombie.
		 */
		if(!scwrlFinished){

			killStatus = kill(forkPid, 15);

			wait(&waitStatus);
		}
	}

	/* Create a new empty structure.
	 */
	shared_ptr<Structure> newStructure;
	
	/* Read the novel structure into memory if SCWRL finished
	 * on time.
	 */
	if(scwrlFinished){

		/* DEBUG
		*
		cout << "SCWRL finished after " << i*(SCWRL_CHECK_DELAY/1000000.0) << " seconds, reading the new file." << endl;
		*/

		ifstream newStruct(structFileNew.c_str());
		newStructure = shared_ptr<Structure>(new AllAtomStructure(AllAtomStructure::parseScwrlOutput(newStruct)));
		newStruct.close();

		remove(structFileNew.c_str());
	}

	/* Delete the temporary file.
	 */
	remove(structFileOrig.c_str());

	/* Return the tweaked (or empty) structure.
	 */
	return newStructure;
}

/* Function to optimially adjust side chain positions in
 * a complex of several Structures (mimicing a limited
 * degree of "induced fit" in the binding mechanism).
 *
 * Returns a vector of null pointers if SCWRL fails.
 */
vector< shared_ptr<Structure> > AllAtomStructure::adjustComplex(vector< shared_ptr<Structure> > components){
	vector< shared_ptr<Structure> > adjusted;
	vector<Residue> combinedResidues, splitResidues, tempResidues;
	shared_ptr<Structure> combinedAndAdjusted;
	int i, j, numResidues;

	/* Combine all components into one complex.
	*/
	combinedResidues = components[0]->getAllResidues();
	for(i = 1; i < components.size(); i++){
		tempResidues = components[i]->getAllResidues();
		combinedResidues.insert(combinedResidues.end(), tempResidues.begin(), tempResidues.end());
	}
	AllAtomStructure combinedStruct(combinedResidues);

	/* DEBUG
	 *
	 cout << "Structure before SCWRL: " << combinedStruct.getNumResidues() << " residues." << endl;
	 cout << "Sequence: " << combinedStruct.getSequence() << endl;
	 cout << "Structure:\n" << combinedStruct.toPDBString() << endl;
	 */

	/* Minimize the steric clashes using SCWRL on the whole
	 * complex.
	 */
	combinedAndAdjusted = combinedStruct.adjustSidechains();

	/* DEBUG
	 *
	 cout << "After SCWRL: " << combinedStruct.getNumResidues() << endl;
	 cout << "Sequence: " << combinedStruct.getSequence() << endl;
	 cout << "Structure: " << combinedStruct.toPDBString() << endl;
	 */
	
	/* Check for SCWRL failure.
	 */
	if(!combinedAndAdjusted){
		for(i = 0; i < components.size(); i++){
			shared_ptr<Structure> pNothing;
			adjusted.push_back( pNothing );
		}
		return adjusted;
	}

	/* Split up the pieces again.
	*/
	combinedResidues = combinedAndAdjusted->getAllResidues();

	/* Rebuild the protein structures.
	*/
	numResidues = 0;
	for(i = 0; i < components.size(); i++){
		splitResidues.clear();
		for(j = 0; j < components[i]->getNumResidues(); j++){
			splitResidues.push_back(combinedResidues[numResidues]);
			numResidues++;
		}
		shared_ptr<Structure> tempStructure(new AllAtomStructure(splitResidues));
		adjusted.push_back(tempStructure);
	}

	/* Send them back out.
	 */
	return adjusted;
}
