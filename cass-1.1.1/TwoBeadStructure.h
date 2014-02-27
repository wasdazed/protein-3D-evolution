/* Class to represent a protein structure reduced to the
 * two-bead model used in Mukherjee and Bagchi, 2003, J Chem Phys
 * 118:10.
 */

/* Idempotency.
 */
#ifndef TWOBEADSTRUCTURE_H
#define TWOBEADSTRUCTURE_H

/* Forward declarations.
 */
class AllAtomStructure;
class OneBeadStructure;

/* Includes.
 */
#include "Structure.h"

/* Definitions.
 */
#define MCMC_MAX_ITERATIONS 100
#define SEQ_SEPARATION 4 // Keep this synched to EnergyModel.h

/* Interface.
 */
class TwoBeadStructure : public Structure{
	public:
		/* Constructors and destructors.
		 */
		TwoBeadStructure();
		TwoBeadStructure(shared_ptr<AllAtomStructure> detailedStructure);
		TwoBeadStructure(vector<Residue> residues); // SHOULD THIS IN FACT BE protected?
		virtual ~TwoBeadStructure();

		virtual shared_ptr<Structure> getBackbone() const;
		virtual vector< vector<int> > contactMap(float contactDistance) const;
		virtual vector< vector<float> > distanceMap() const;
		virtual vector< vector<int> > interactionContactMap(float contactDistance, shared_ptr<Structure> interactor) const;
		virtual unsigned int getNumContacts(float contactDistance, unsigned int sequenceSeparation);

		vector< pair<float, float> > bondAngles() const;
		vector<float> caCbDistances() const;
		vector<float> cbCbDistances() const;
		vector<string> threeMers() const;
		vector<string> helixFourMers() const;
		vector<string> betaFourMers() const;
		vector< pair<float,float> > backboneDistances() const;
		vector<float> backboneTorsions() const;
		vector<float> residueExposures() const;

		static TwoBeadStructure parseBeadFile(ifstream &inputFile);

		shared_ptr<Structure> threadSequence(string sequence);
		shared_ptr<Structure> adjustSidechains();
		vector< shared_ptr<Structure> > adjustComplex(vector< shared_ptr<Structure> > components);

		void calcExposure(bool inComplex);
		static vector< shared_ptr<TwoBeadStructure> > calcComplexExposure(vector< shared_ptr<TwoBeadStructure> > components);

		void calcAndFreezeGeometries(float contactDistance);

	protected:

		static float beadRadii[20];
		static string radiusKey;

		vector< vector<int> > contacts;
		
		vector< pair<float, float> > angles;
		vector<float> cacbdists;
		vector<float> cbcbdists;
		vector<string> threemers;
		vector<string> hfourmers;
		vector<string> bfourmers;
		vector< pair<float,float> > bbdists;
		vector<float> bbtors;
		vector<float> exposures;

};

/* Idempotency end.
 */
#endif
