#include "Organism.h"
#include "myUtils.h"
#include "Structure.h"
#include <stdexcept>

/* Map of what sequences correspond to what folding energies
 * given a certain Structure, if an object of the class has
 * encountered that combination before.
 *
 * Is re-initialized to proper size (same as structs vector) by 
 * first Organism object that is constructed.
 */
vector<string_float_map> Organism::knownFoldEnergies;
vector<string_structure_map> Organism::knownStructures;

/* Some other shared properties.
 */
shared_ptr<CRandomMersenne> Organism::pMersenneGen;
shared_ptr<StochasticLib1> Organism::pStochastGen;
unsigned long int Organism::totalCount = 0;

/* Default constructor,  doesn't do much.
 */
Organism::Organism(){
	id = Organism::totalCount;
	parentId = Organism::totalCount;
	Organism::totalCount++;
}

/* Default copy constructor.
 */
Organism::Organism(const Organism& org){
	mutRate = org.mutRate;
	fitness = org.fitness;
	dnaSeqs = org.dnaSeqs;
	structs = org.structs;
	cumNumMutations = org.cumNumMutations;
	foldThreshold = org.foldThreshold;
	phenotype = org.phenotype;
	changeHappened = org.changeHappened;

	pEngModel = org.pEngModel;

	id = Organism::totalCount;
	parentId = org.id;
	Organism::totalCount++;
}

/* Destructor. Needs definition due to finicky
 * feature of gcc.
 */
Organism::~Organism(){
}

/* Copy constructor with possibility of mutation.
 */
Organism::Organism(const Organism& org, bool evolve){

	string oldProt, newProt;
	unsigned int i, j;

	/* DEBUG
	 *
	 cout << "Base class copy constructor." << endl;
	 */

	mutRate = org.mutRate;
	fitness = org.fitness;
	phenotype = org.phenotype;
	changeHappened = org.changeHappened;
	structs = org.structs;

	cumNumMutations = org.cumNumMutations;
	foldThreshold = org.foldThreshold;

	pEngModel = org.pEngModel;

	/* If there's a possiblity of mutation,
	 * go ahead and try to mutate the DNA of
	 * all the sequences.
	 */
	if(evolve){
		for(i = 0; i < org.dnaSeqs.size(); i++){

			/* DEBUG
			 *
			cout << "Mutating the DNA..." << endl;
			*/

			dnaSeqs.push_back( mutateDna(org.dnaSeqs[i], mutRate, *Organism::pStochastGen, *Organism::pMersenneGen ) );

			/* Count the mutations.
			 */
			for(j = 0; j < dnaSeqs[i].length(); j++){
				if(dnaSeqs[i][j] != org.dnaSeqs[i][j]){
					cumNumMutations++;
				}
			}

			oldProt = translateDna(org.dnaSeqs[i]);
			newProt = translateDna(dnaSeqs[i]);

			/* If there's been a non-synonymous substitution, thread the
			 * novel sequence (unless there's a stop codon).
			 */
			if(newProt.compare(oldProt) != 0){

				/* DEBUG
				 *
				cout << "Non-synonymous change occured, evaluating..." << endl;
				*/

				changeHappened = true;

				/* Check for stop codons in mutated sequence. If they're found,
				 * the resulting Structure* is NULL.
				 */
				if(newProt.find("_") != string::npos){
					/* DEBUG
					*
					cout << "Stop codon appeared: automatic folding failure." << endl;
					*/

					shared_ptr<Structure> pNothing;
					structs[i] = pNothing;
				}
				/* Thread the new sequence. Whatever method is used
				 * may fail: needs to be taken care of downstream.
				 */
				else{

					/* DEBUG
					 *
					 cout << "Re-threading structure " << (i+1) << "..." << endl;
					 */

					/* Have we done this thread before?
					 */
					if(Organism::knownStructures[i].find(newProt) == Organism::knownStructures[i].end()){
						/* DEBUG
						 *
						cout << "Haven't threaded this before, doing so..." << endl;
						*/

						structs[i] = structs[i]->threadSequence(newProt);

						/* DEBUG
						 *
						cout << "Finished, saving it for later use at address " << structs[i].get() << "." << endl;
						*/

						Organism::knownStructures[i].insert( stringAndStructure(newProt, structs[i]) );
					}
					else{
						/* DEBUG
						 *
						cout << "Seen this sequence in this conformation before, using the stored one..." << endl;
						*/

						structs[i] = Organism::knownStructures[i].find(newProt)->second;
					}

				}
			}
		}
	}
	/* Otherwise, just copy the DNA sequences.
	*/
	else{
		dnaSeqs = org.dnaSeqs;
	}
	
	id = Organism::totalCount;
	parentId = org.id;
	Organism::totalCount++;

	/* DEBUG
	 *
	cout << "Finished copying organism " << org.id << "." << endl;
	*/
}

/* Copy constructor with possibility of mutation.
 * Mainly for use with derived classes or the
 * heap.
 */
shared_ptr<Organism> Organism::Clone(bool evolve){

	/* DEBUG
	 *
	cout << "Base clone." << endl;
	*/

	return shared_ptr<Organism>(new Organism(*this, evolve));
}

/* Constructor from genome of coding sequences and associated
 * structures.
 */
Organism::Organism(float mutationRate, vector<string> codingSequences, vector< shared_ptr<Structure> > structures, shared_ptr<EnergyModel> engModel, float foldingThreshold){

	unsigned int i;
	string structProt, geneProt;

	mutRate = mutationRate;
	cumNumMutations = 0;
	pEngModel = engModel;
	foldThreshold = foldingThreshold;

	dnaSeqs = codingSequences;
	structs = structures;
	changeHappened = true; // Being created sure was a big change...

	/* Make sure the class has its PRNGs handy.
	 */
	if(Organism::pMersenneGen.get() == 0 || Organism::pStochastGen.get() == 0){
		throw runtime_error("You must supply random number generators via Organism::setPRNGs()!");
	}

	/* Initialize the list of known structure-
	 * sequence combinations if this is the first
	 * Organism constructed.
	 */
	if(Organism::knownFoldEnergies.empty()){
		for(i = 0; i < structs.size(); i++){
			Organism::knownFoldEnergies.push_back(string_float_map());
		}
	}	
	if(Organism::knownStructures.empty()){
		for(i = 0; i < structs.size(); i++){
			Organism::knownStructures.push_back(string_structure_map());
		}
	}

	pEngModel = engModel;

	/* If the translated input DNA is different from the supplied structures,
	 * thread the genomic sequences. 
	 *
	 * This is expected for e.g. restarting from sequences on disk.
	 */
	for(i = 0; i < dnaSeqs.size(); i++){
		structProt = structs[i]->getSequence();
		geneProt = translateDna(dnaSeqs[i]);
		if(structProt.compare(geneProt) != 0){
			if(geneProt.find("_") == string::npos){
				structs[i] = structs[i]->threadSequence(geneProt);
				Organism::knownStructures[i].insert( stringAndStructure(geneProt, structs[i]) );
			}
			else{
				shared_ptr<Structure> pNothing;
				structs[i] = pNothing;
			}
		}
	}

	/* Insert some dummy values for fitness
	 * and phenotype.
	 */
	calcFitnessAndPhenotype(); 

	id = Organism::totalCount;
	parentId = Organism::totalCount;
	Organism::totalCount++;
}

/* Accessor method for mutation rate
 */
float Organism::getMutationRate() const{
	return mutRate;
}

/* Accessor method for all the coding sequences.
 */
vector<string> Organism::getDnaSequences() const{
	vector<string> seqs;
	unsigned int i;

	for(i = 0; i < dnaSeqs.size(); i++){
		seqs.push_back(dnaSeqs[i]);
	}

	return seqs;
}

/* Accessor method for all the protein sequences.
 */
vector<string> Organism::getProteinSequences() const{
	vector<string> prots;
	unsigned int i;

	for(i = 0; i < dnaSeqs.size(); i++){
		prots.push_back(translateDna(dnaSeqs[i]));
	}

	return prots;
}

/* Accessor method for all the structures.
 */
vector< shared_ptr<Structure> > Organism::getStructures() const{
	return structs;
}

/* Accessor method for cumulative number of mutations
 */
int Organism::getNumberOfMutations() const{
	return cumNumMutations;
}

/* Accessor method for fitness. Calculates the
 * fitness if not already set.
 */
float Organism::getFitness(){

	/* DEBUG
	 *
	cout << "Retrieving fitness..." << endl;
	*/

	/* Only the dummy fitness value has been set,
	 * or there's been a recent mutation: calculate the real one.
	 */
	if(fitness == -1.0 || changeHappened){
		calcFitnessAndPhenotype();
		changeHappened = false;
	}

	return fitness;
}

/* Accessor method for the phenotype. Calculates
 * the phenotype if not already set.
 */
string Organism::getPhenotype(){

	/* DEBUG
	 *
	cout << "Retrieving phenotype..." << endl;
	*/

	/* Only the dummy phenotype value has been set,
	 * or there's been a recent mutation: calculate the real one.
	 */
	if(phenotype.compare("dummy") == 0 || changeHappened){
		calcFitnessAndPhenotype();
		changeHappened = false;
	}

	return phenotype;
}

/* Generates a string representation of some of the data members 
 * in this object.
 *
 * Currently: id, parent id, phenotype, dna seqs, protein seqs
 *
 * Useful to call when using subclass versions of
 * the same function.
 */
string Organism::toString() const{

	string stringRep = "";
	char tempString[500];
	unsigned int i;

	/* Each simple data member needed to represent the object
	 * gets its own line.
	 */
	sprintf(tempString, "%u\t", id);	
	stringRep += tempString; 
	
	sprintf(tempString, "%u\t", parentId);
	stringRep += tempString;

	stringRep += phenotype;

	/* Each DNA sequence gets its own column.
	 */
	for(i = 0; i < dnaSeqs.size(); i++){
		stringRep += "\t";
		stringRep += dnaSeqs[i];
	}

	/* Same for the proteins.
	 */
	for(i = 0; i < dnaSeqs.size(); i++){
		stringRep += "\t";
		stringRep += translateDna(dnaSeqs[i]);
	}

	return stringRep;
}

/* Default implementation for fitness function:
 * always sets fitness to "-1.0" and phenotype to
 * "dummy".
 */
void Organism::calcFitnessAndPhenotype(){

	/* DEBUG
	 *
	cout << "In Organism fitness calculation." << endl;
	*/

	fitness = -1.0;
	phenotype = "dummy";
}

/* Default folding function: always returns
 * a vector with a single 1e10.
 */
vector<float> Organism::calcFoldingEnergies(){
	vector<float> dummy;
	dummy.push_back(1000000000.0);

	return dummy;
}

/* Default binding function: always returns
 * 1e10.
 */
float Organism::calcBindingEnergy(shared_ptr<Structure> bindingStructure, shared_ptr<Structure> proteinStructure){
	return 1000000000.0;
}

void Organism::setPRNGs(CRandomMersenne &rangeGenerator, StochasticLib1 &stochasticGenerator){
	Organism::pMersenneGen = shared_ptr<CRandomMersenne>(new CRandomMersenne(rangeGenerator));
	Organism::pStochastGen = shared_ptr<StochasticLib1>(new StochasticLib1(stochasticGenerator));
}
