/* Program to get the distribution of folding and binding energies from
 * a single structure given a set of ligands. See 
 * Rastogi and Liberles, 2006, Biophys Chem 126 and
 * Grahnen et al, 2011, BMC Evol Biol 11:361 for a description of the methods.
 *
 * Usage: [TBA]
 */

#include "BastollaModel.h"
#include "BastollaAugmentedModel.h"
#include "GrahnenModel.h"
#include "AllAtomStructure.h"
#include "TwoBeadStructure.h"
#include "myUtils.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

using namespace std;

int main(int argc, char *argv[]){

	string usage = " [structure bead file] [ligand bead file] [neofunc ligand bead file] [decoy ligand bead file] [energy function (Grahnen/Bastolla++]";

	/* Usage message on bad input.
	 */
	if(argc < 6){
		cout << "Usage: " << argv[0] << usage << endl;
		exit(1);
	}

	/* Read arguments.
	 */
	ifstream protFile(argv[1]);
	ifstream ligFile(argv[2]);
	ifstream neofuncFile(argv[3]);
	ifstream decoyLigFile(argv[4]);
	string energyFunction = argv[5];

	/* Error checking on input.
	 */
	if(energyFunction.compare("Grahnen") != 0 
	   &&
	   energyFunction.compare("Bastolla++") != 0
	  ){
		cout << "Usage: ./" << argv[0] << usage << endl;
		exit(1);
	}

	/* Define variables.
	 */
	string inputProtein;

	unsigned int i;

	float foldO, bindO, bindD, bindN;
	
	vector<float> otherScores;
	vector<float> bindingScores;

	shared_ptr<TwoBeadStructure> proteinStructure;
	shared_ptr<Structure> ligandStructure;
	shared_ptr<Structure> neoLigandStructure;
	shared_ptr<Structure> decoyLigandStructure;

	shared_ptr<EnergyModel> engModel;
	
	/* Seed the random number generators.
	 */
	CRandomMersenne mersenneRandGen(time(NULL));
	StochasticLib1 stochastRandGen(time(NULL));
	srand(time(NULL));

	/* Initialize the energy model.
	 */
	if(energyFunction.compare("Grahnen") == 0){
		engModel = shared_ptr<EnergyModel>(new GrahnenModel(0.213649, 0.255238, 0.0351669, 0.193426, 0.00636471, 0.296156, 0.0, 0.970, 0.011, 0.019)); // Best set from paper
	}
	else{
		engModel = shared_ptr<EnergyModel>(new BastollaAugmentedModel());
	}

	/* Read protein, ligands into memory.	 
	 */
	shared_ptr<TwoBeadStructure> proteinStruct(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(protFile)));
	shared_ptr<TwoBeadStructure> ligandStruct(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(ligFile)));
	shared_ptr<TwoBeadStructure> neoLigandStruct(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(neofuncFile)));
	shared_ptr<TwoBeadStructure> decoyLigandStruct(new TwoBeadStructure(TwoBeadStructure::parseBeadFile(decoyLigFile)));
	inputProtein = proteinStruct->getSequence();
	protFile.close();
	ligFile.close();
	neofuncFile.close();
	decoyLigFile.close();

	/* Pre-calculate exposures.
	 */
	proteinStruct->calcExposure(false);
	ligandStruct->calcExposure(false);
	neoLigandStruct->calcExposure(false);
	decoyLigandStruct->calcExposure(false);

	vector< shared_ptr<TwoBeadStructure> > wholeComplex;
	wholeComplex.push_back(proteinStruct);
	wholeComplex.push_back(ligandStruct);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	proteinStructure = wholeComplex[0];
	ligandStructure = wholeComplex[1];

	wholeComplex.clear();
	wholeComplex.push_back(proteinStruct);
	wholeComplex.push_back(neoLigandStruct);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	neoLigandStructure = wholeComplex[1];

	wholeComplex.clear();
	wholeComplex.push_back(proteinStruct);
	wholeComplex.push_back(decoyLigandStruct);
	wholeComplex = TwoBeadStructure::calcComplexExposure(wholeComplex);
	decoyLigandStructure = wholeComplex[1];

	/* Make a "frozen" copy for geometry sampling.
	 */
	shared_ptr<TwoBeadStructure> invariant(new TwoBeadStructure(*proteinStructure));
	invariant->calcAndFreezeGeometries(4.5);

	/* Do the folding and three binding energies.
	*/
	foldO = engModel->calcFoldEnergyGap(proteinStructure, invariant);
	bindO = engModel->calcInteractionScore(proteinStructure, ligandStructure);
	bindD = engModel->calcInteractionScore(proteinStructure, decoyLigandStructure);
	bindN = engModel->calcInteractionScore(proteinStructure, neoLigandStructure);

	/* Common info to all models.
	 */
	cout << "Seq\tFoldO\tBindO\tBindD\tBindN";

	/* If using the Grahnen model, also output individual terms of the
	 * energy function for original folding and function.
	 */
	if(energyFunction.compare("Grahnen") == 0){

		/* Model-specific info will be outputted.
		 */
		cout << "\tV_bend\tV_LJ\tV_helix\tV_beta\tV_charge\tV_solv\tV_ssbond\tV_bindLJ\tV_bindCharge\tV_bindSolv";
		shared_ptr<GrahnenModel> actualModel = dynamic_pointer_cast<GrahnenModel>(engModel);
		otherScores = actualModel->getFoldScoreTerms(proteinStructure);
		bindingScores = actualModel->getBindScoreTerms(proteinStructure, ligandStructure);
		otherScores.insert(otherScores.end(), bindingScores.begin(), bindingScores.end());
	}
	cout << endl;

	/* DEBUG
	 *
	 cerr << "Calculated energies." << endl;
	 */

	/* Print results
	*/
	cout << inputProtein << "\t" << foldO << "\t" << bindO << "\t" << bindD << "\t" << bindN;
	for(i = 0; i < otherScores.size(); i++){
		cout << "\t" << otherScores[i];
	}
	cout << endl;

	cout << inputProtein << "\t" << foldO << "\t" << bindO << "\t" << bindD << "\t" << bindN;
  cout << foldO << "\t" << bindO;
	for(i = 0; i < otherScores.size(); i++){
		cout << "\t" << otherScores[i];
	}
  cout << endl;




	/* Finish.
	*/
	exit(0);
}
