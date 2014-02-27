/* Program to generate and calculate binding energy for a
 * number of possible variations on a particular motif
 * under an unparameterized Grahnen model and a particular
 * protein complex. Motifs are for the "ligand" part of
 * the complex, and are not scanned for folding energy.
 *
 * Motifs should be specified as e.g.
 *
 * XXXTIYXX(VI)X(KHR)
 *
 * where X is replaced by all possible residues and positions
 * within parentheses are replaced by any of the specified
 * options. Do NOT put X within parentheses, and do NOT
 * nest parentheses.
 */

/* Includes.
 */
#include "TwoBeadStructure.h"
#include "AllAtomStructure.h"
#include "GrahnenModel.h"
#include "BastollaAugmentedModel.h"
#include "ParameterRotation.h"
#include "MHMCMC.h"

int main(int argc, char* argv[]){

	string usage = " [protein bead file] [current ligand bead file] [motif] [# samples] [energy model (Grahnen/Bastolla++)] [protein DNA file] [binding residue position file] [# mutations to make]";

	/* Check number of arguments.
	 */
	if(argc < 9){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Read arguments.
	 */
	ifstream proteinFile(argv[1]);
	ifstream ligandFile(argv[2]);
	string motif = argv[3];
	unsigned int numSamples = atoi(argv[4]);
	string energyModel = argv[5];
	ifstream dnaFile(argv[6]);
	ifstream residueFile(argv[7]);
	unsigned int numMutsMade = atoi(argv[8]);

	/* Check input.
	 */
	if(energyModel.compare("Grahnen") != 0 && energyModel.compare("Bastolla++") != 0){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Define some variables.
	 */
	unsigned int i, j, k, numBetterThanNative;
	unsigned int numMutSamples = 100;
	bool inOption;
	float bindScore, mutNativeScore, mutDecoyScore;
	string motifSequence;
	string protDna, siteDna, mutSiteDna, mutProtDna, mutProtSeq;

	shared_ptr<Structure> mutStructure;

	vector<string> motifComposition;
	vector<float> engTerms;
	
	/* Read protein and ligand structures, and make a "fake" all-atom
	 * version to use for deterministic threading.
	 */
	shared_ptr<TwoBeadStructure> proteinStructure(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(proteinFile)));
	shared_ptr<TwoBeadStructure> beadLigand(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(ligandFile)));
	proteinFile.close();
	ligandFile.close();

	/* And the DNA.
	 */
	dnaFile >> protDna;
	dnaFile.close();

	/* And the binding site residue positions.
	 */
	vector<int> bindResidues;
	int site;
	while(!residueFile.eof()){
		residueFile >> site;
		bindResidues.push_back(site);
	}

	/* Also read in a "fake" all-atom version to use
	 * for deterministic threading.
	 */
	ligandFile.open(argv[2]);
	shared_ptr<AllAtomStructure> threadingTarget(new AllAtomStructure(AllAtomStructure::parseScwrlOutput(ligandFile)));
	ligandFile.close();

	/* Get the binding site DNA.
	 */
	siteDna = "";
	for(i = 0; i < bindResidues.size(); i++){
		siteDna += protDna.substr((bindResidues[i]-1)*3, 3);
	}

	/* Pre-calculate exposures.
	 */
	proteinStructure->calcExposure(false);
	beadLigand->calcExposure(false);
	vector< shared_ptr<TwoBeadStructure> > wholeComplex;
	wholeComplex.push_back(proteinStructure);
	wholeComplex.push_back(beadLigand);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	proteinStructure = wholeComplex[0];
	beadLigand = wholeComplex[1];

	/* Start random number generators.
	 */
	StochasticLib1 stochasticGenerator(time(NULL));
	CRandomMersenne randGenerator(time(NULL));
	ParameterRotation::setPRNGs(randGenerator, stochasticGenerator);
	MHMCMC::setPRNG(randGenerator);
	
	/* Instantiate the relevant energy model.
	 */
	shared_ptr<EnergyModel> engModel;
	if(energyModel.compare("Grahnen") == 0){
		engModel = shared_ptr<EnergyModel>(new GrahnenModel(0.213649, 0.255238, 0.0351669, 0.193426, 0.00636471, 0.296156, 0.0, 0.970, 0.011, 0.019)); // Replace parameters as necessary
	}
	else{
		engModel = shared_ptr<EnergyModel>(new BastollaAugmentedModel());
	}

	/* Prep output.
	 */
	cout << "#Making " << numSamples << " variations on motif " << motif << "..." << endl;
	cout << "#Original ligand sequence: " << beadLigand->getSequence() << endl;
	cout << "Num\tMotif\tEnergy\tNumBetter" << endl;

	/* The native motif energy.
	 */
	cout << 0 << "\t" << beadLigand->getSequence() << "\t" << engModel->calcInteractionScore(proteinStructure, beadLigand) << "\tNA" << endl;

	/* Parse the motif.
	 */
	inOption = false;
	for(i = 0; i < motif.length(); i++){
		
		/* Start of a multi-character optional position.
		 */
		if(!inOption && motif.substr(i,1).compare("(") == 0){
			inOption = true;
			motifComposition.push_back("");
		}
		/* End of a multi-character optional position.
		 */
		else if(inOption && motif.substr(i,1).compare(")") == 0){
			inOption = false;
		}
		/* A character option for a multi-character position.
		 */
		else if( inOption && motif.substr(i,1).compare("(") != 0 && motif.substr(i,1).compare(")") != 0 ){
			motifComposition.back() += motif.substr(i,1);
		}
		/* An optional character (i.e. any residue).
		 */
		else if(!inOption && motif.substr(i,1).compare("X") == 0){
			motifComposition.push_back("AEQDNLGKSVRTPIMFYCWH");
		}
		/* A fixed character.
		 */
		else if(!inOption && motif.substr(i,1).compare("X") != 0){
			motifComposition.push_back(motif.substr(i,1));
		}
	}

	/* DEBUG
	 *
	cout << "Motif composition:" << endl;
	for(i = 0; i < motifComposition.size(); i++){
		cout << motifComposition[i] << endl;
	}
	*/

	/* Generate a sufficient number of the possible variations.
	 * and calculate energies.
	 */
	for(i = 0; i < numSamples; i++){

		motifSequence = "";

		/* Generate a random sequence within the motif
		 * bounds.
		 */
		for(j = 0; j < motifComposition.size(); j++){

			/* Pick a random character for this position.
			 */
			if(motifComposition[j].length() > 1){
				motifSequence += motifComposition[j].substr(randGenerator.IRandom(0,motifComposition[j].length()-1),1);
			}
			else{
				motifSequence += motifComposition[j];
			}
		}

		/* Thread new sequence with SCWRL (deterministic).
		 */
		shared_ptr<Structure> newLigand = threadingTarget->threadSequence(motifSequence);
		
		/* Check for successs, go back to start if SCWRL
		 * fails.
		 */
		if(newLigand.get() == 0){
			i--;
			continue;
		}

		/* Convert to two-bead representation.
		 */
		shared_ptr<TwoBeadStructure> motifLigand(new TwoBeadStructure(dynamic_pointer_cast<AllAtomStructure>(newLigand)));

		/* Re-calc exposure.
		 */
		motifLigand->calcExposure(false);
		wholeComplex.clear();
		wholeComplex.push_back(proteinStructure);
		wholeComplex.push_back(motifLigand);
		wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
		motifLigand = wholeComplex[1];

		/* Get binding energy for this configuration.
		*/
		bindScore = engModel->calcInteractionScore(proteinStructure, motifLigand);

		/* Perform some number of mutations to the binding site, and
		 * re-measure binding scores.
		 */
		numBetterThanNative = 0;
		for(j = 0; j < numMutSamples; j++){
			/* Mutate the binding site.
			 */
			mutSiteDna = mutateNonsyn(siteDna, numMutsMade, randGenerator);
			mutProtDna = protDna;
			for(k = 0; k < bindResidues.size(); k++){
				mutProtDna.replace((bindResidues[k]-1)*3, 3, mutSiteDna.substr(k*3,3));
			}

			/* Ensure that it didn't produce a stop codon.
			 */
			mutProtSeq = translateDna(mutProtDna);
			if(mutProtSeq.find("_") != string::npos){
				j--;
				continue;
			}

			/* Thread the resulting sequence.
			 */
			mutStructure = proteinStructure->threadSequence(translateDna(mutProtDna));

			/* Make sure the threading worked.
			 */
			if(mutStructure.get() == 0){
				j--;
				continue;
			}

			/* Measure the new ligand binding scores.
			 */
			mutDecoyScore = engModel->calcInteractionScore(mutStructure, motifLigand);
			mutNativeScore = engModel->calcInteractionScore(mutStructure, beadLigand);

			/* Was it better than native (e.g. outcompetes it)?
			 */
			if(mutDecoyScore < mutNativeScore){
				numBetterThanNative++;
			}
		}

		/* Print motif and stats to STDOUT.
		 */
		cout << (i+1) << "\t" << motifSequence << "\t" << bindScore << "\t" << numBetterThanNative << endl;
	}

	/* Finish.
	 */
	exit(0);
}
