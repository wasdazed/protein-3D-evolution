/* Program to run evolution simulations of an SH2 fold binding a
 * ligand. See Rastogi et al., 2006, Biophys chem 126 for a
 * description of some of the methods: Grahnen et al, 2011, BMC Evol
 * Biol 11:361 for a more detailed description.
 *
 */

/* Includes.
 */
#include "TwoBeadStructure.h"
#include "SingleSiteSim.h"
#include "MultipleStructureOrganism.h"
#include "BastollaAugmentedModel.h"
#include "GrahnenModel.h"
#include "ParameterRotation.h"
#include "MHMCMC.h"
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>

using namespace std;

/* Main loop: run the simulation and collect data.
*/
int main(int argc, char *argv[]){
	string usage = " [phenotype output file] [# of generations] [sample size] [mutation rate] [fold threshold] [bind threshold] [protein bead file] [protein coding sequence] [ligand bead file] [decoy ligand bead file] [opt: sequence file to restore simulation from]";

	/* Check for sufficient number of arguments
	 */
	if(argc < 11){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}
	
	/* Read a bunch of arguments
	 */
	ofstream outFile(argv[1]);
	int generations = atoi(argv[2]);
	int samples = atoi(argv[3]);
	float mutRate = atof(argv[4]);
	float foldThreshold = atof(argv[5]);
	float bindThreshold = atof(argv[6]);
	ifstream proteinBeadFile(argv[7]);
	ifstream proteinSeqFile(argv[8]);
	ifstream ligandBeadFile(argv[9]);
	ifstream decoyLigBeadFile(argv[10]);

	/* Declare some variables.
	 */
	unsigned int i,j;
	int numDead, numNormal, numDecoy;
	int currentGeneration, lastGeneration;

	string proteinCodingSeq;
	string decoySeq;
	string line;

	vector<string> currentDna;
	vector<string> currentProtein;
	vector<string> currentPhenotypes;

	vector<string> startingDna;
	vector< shared_ptr<Structure> > startingStructures;
	vector< shared_ptr<Structure> > ligandStructures;
	vector< shared_ptr<Structure> > decoyLigandStructures;
	vector< shared_ptr<Organism> > allOrganisms;

	shared_ptr<SingleSiteSim> pSimulation;
	shared_ptr<MultipleStructureOrganism> pOriginalSample;

	/* Seed and set the random number generators.
	 */
	CRandomMersenne mersenneRandGen(time(NULL));
	StochasticLib1 stochastRandGen(time(NULL));
	srand(time(NULL));
	Organism::setPRNGs(mersenneRandGen, stochastRandGen);
	Simulation::setPRNG(mersenneRandGen);
	ParameterRotation::setPRNGs(mersenneRandGen, stochastRandGen);
	MHMCMC::setPRNG(mersenneRandGen);
		
	/* Read protein, ligand and sequences into memory.
	*/
	shared_ptr<TwoBeadStructure> proteinStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(proteinBeadFile)));
	shared_ptr<TwoBeadStructure> ligandStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(ligandBeadFile)));
	shared_ptr<TwoBeadStructure> decoyLigandStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(decoyLigBeadFile)));
	proteinBeadFile.close();
	ligandBeadFile.close();
	decoyLigBeadFile.close();

	/* Pre-calculate exposure values as necessary.
	*/
	proteinStructure->calcExposure(false);
	ligandStructure->calcExposure(false);
	decoyLigandStructure->calcExposure(false);

	vector< shared_ptr<TwoBeadStructure> > wholeComplex;
	wholeComplex.push_back(proteinStructure);
	wholeComplex.push_back(ligandStructure);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	proteinStructure = wholeComplex[0];
	ligandStructure = wholeComplex[1];

	wholeComplex.clear();
	wholeComplex.push_back(proteinStructure);
	wholeComplex.push_back(decoyLigandStructure);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	decoyLigandStructure = wholeComplex[1];

	/* And pre-calculate starting geometry values.
	 */
	shared_ptr<TwoBeadStructure> geometry(new TwoBeadStructure(*proteinStructure));
	geometry->calcAndFreezeGeometries(4.5);

	/* No sequence input: starting from a homogenous population.
	 */
	if(argc == 11){

		cout << "#Creating a new homogenous population of " << samples << " organisms from parameters..." << endl;

		proteinSeqFile >> proteinCodingSeq;
		proteinSeqFile.close();

		/* Create a starting organism from input data
		*/
		startingDna.push_back(proteinCodingSeq);

		startingStructures.push_back(proteinStructure);
		ligandStructures.push_back(ligandStructure);
		decoyLigandStructures.push_back(decoyLigandStructure);

		//shared_ptr<EnergyModel> engModel(new BastollaAugmentedModel());
		shared_ptr<EnergyModel> engModel(new GrahnenModel(0.213649, 0.255238, 0.0351669, 0.193426, 0.00636471, 0.296156, 0.0, 0.970, 0.011, 0.019)); // 20110725

		pOriginalSample = shared_ptr<MultipleStructureOrganism>(new MultipleStructureOrganism( mutRate, startingDna, startingStructures, ligandStructures, decoyLigandStructures, engModel, foldThreshold, bindThreshold, geometry ));

		/* Initiate the simulation.
		 */
		pSimulation = shared_ptr<SingleSiteSim>(new SingleSiteSim(generations, pOriginalSample, samples));
	}

	/* Sequence input is given: starting from a population with 
	 * non-identical sequences but otherwise identical parameters.
	 */
	else if(argc == 12){
		
		cout << "#Creating a population from parameters and population of sequences..." << endl;

		startingStructures.push_back(proteinStructure);
		ligandStructures.push_back(ligandStructure);
		decoyLigandStructures.push_back(decoyLigandStructure);

		//shared_ptr<EnergyModel> engModel(new BastollaAugmentedModel());
		shared_ptr<EnergyModel> engModel(new GrahnenModel(0.213649, 0.255238, 0.0351669, 0.193426, 0.00636471, 0.296156, 0.0, 0.970, 0.011, 0.019)); //20110725

		/* Figure out which generation was the last in the previous
		 * experiment.
		 */
		ifstream seqCollectionFile(argv[11]);
		while( getline(seqCollectionFile, line) ){
			/* Don't read comment lines or empty lines.
			 */
			if(line.substr(0,1).compare("#") != 0 && line.length() > 1){
				vector<string> tokens = tokenize(line, "\t");
				lastGeneration = atoi(tokens[0].c_str());
			}
		}
		seqCollectionFile.close();

		/* DEBUG
		 *
		cout << "Last generation in previous experiment: " << lastGeneration << endl;
		*/

		/* Read each sequence in the last generation and create an organism from
		 * it and the other parameters.
		 */
		seqCollectionFile.open(argv[11]);
		while( getline(seqCollectionFile, line) ){
			/* Don't read comment lines or empty lines.
			 */
			if(line.substr(0,1).compare("#") != 0 && line.length() > 1){
				/* Split line on tabs. First column holds the 
				 * generation number, fifth one the sequences.
				 */
				vector<string> tokens = tokenize(line, "\t");
				currentGeneration = atoi(tokens[0].c_str());

				if(currentGeneration == lastGeneration){
					/* DEBUG
					 *
					cout << "Input sequence: " << tokens[4] << endl;
					*/

					/* Create a starting organism from input data.
					*/
					startingDna.clear();
					startingDna.push_back(tokens[4]);

					allOrganisms.push_back( shared_ptr<MultipleStructureOrganism>(new MultipleStructureOrganism(mutRate, startingDna, startingStructures, ligandStructures, decoyLigandStructures, engModel, foldThreshold, bindThreshold, geometry)) );
				}
			}
		}

		/* Close input file.
		*/
		seqCollectionFile.close();

		/* Initiate the simulation.
		 */
		pSimulation = shared_ptr<SingleSiteSim>(new SingleSiteSim(generations, allOrganisms));

		/* DEBUG
		 *
		cout << "Starting sequence set: " << endl;
		currentProtein = pSimulation->getPopProtein();
		for(i = 0; i < currentProtein.size(); i++){
			cout << currentProtein[i] << endl;
		}
		*/
	}
	
	/* Prep the outfile.
	 */
	outFile << "#Gen\tDead\tNormal\tDecoy" << endl;

	/* DEBUG
	 */
	cout << "#Instantiated simulation, going to run " << pSimulation->getNumGenerations() << " steps with " << pSimulation->getNumSamples() << " organisms." << endl;
	cout << "Gen\tID\tParentID\tPhenotype\tDNA\tProtein\tCount" << endl;

	/* In each generation...
	*/
	for(i = 0; i < pSimulation->getNumGenerations(); i++){
		/* Do one whole step of the simulation algorithm
		 */
		pSimulation->doTimeStep();

		/* Get some phenotypic data from the simulation.
		 */
		numDead = 0;
		numNormal = 0;
		numDecoy = 0;
		currentPhenotypes = pSimulation->getPopPhenotypes();
		for(j = 0; j < currentPhenotypes.size(); j++){
			/* DEBUG
			 *
			cout << "Phenotype " << j+1 << ": '" << currentPhenotypes[j] << "'" << endl;
			*/

			if(currentPhenotypes[j].compare("Dead") == 0){
				numDead++;
			}
			else if (currentPhenotypes[j].compare("Normal") == 0){
				numNormal++;
			}
			else if (currentPhenotypes[j].compare("Decoy") == 0){
				numDecoy++;
			}
		}

		/* Print current status to file
		 */
		outFile << i+1 << "\t" << numDead << "\t" << numNormal << "\t" << numDecoy << endl;

		/* Write relevant properties to file every generation.
		 *
		 * In this case, not the whole population, but only the
		 * unique alleles and their counts.
		 */
		unordered_map<string, pair< shared_ptr<Organism>, unsigned int > > uniqueAlleles;
		vector< shared_ptr<Organism> > currentPop(pSimulation->getPop());
		for(j = 0; j < currentPop.size(); j++){
			vector<string> dna = currentPop[j]->getDnaSequences();
			if(uniqueAlleles.find(dna[0]) == uniqueAlleles.end()){ // Should only contain one gene
				uniqueAlleles[dna[0]] = pair< shared_ptr<Organism>, unsigned int >(currentPop[j], 1);
			}
			else{
				uniqueAlleles[dna[0]].second++;
			}
		}
		for(unordered_map< string, pair<shared_ptr<Organism>, unsigned int> >::iterator it = uniqueAlleles.begin(); it != uniqueAlleles.end(); it++){
			cout << i+1 << "\t" << it->second.first->toString() << "\t" << it->second.second << endl;
		}

		/* Print status message to STDOUT.
		 *
		cout << "Generation " << i+1 << " done." << endl;
		*/
	}

	/* Print final results to file
	*/
	outFile << "#Simulation finished." << endl;

	/* Finish up
	*/
	outFile.close();
	exit(0);
}
