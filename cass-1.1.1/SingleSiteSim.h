/* Class to run an evolutionary simulation under the single binding
 * site model (see e.g. Massey et al, 2008, Gene 418(1-2):22-6 for
 * a description of the basic model).
 */

/* Header guards.
 */
#ifndef SINGLESITESIM_H
#define SINGLESITESIM_H

/* Forward declarations.
 */

/* Includes.
 */
#include "Simulation.h"

class SingleSiteSim : public Simulation{
	public:
		SingleSiteSim();
		SingleSiteSim(int generations, shared_ptr<Organism> origSample, int sampleSize);
		SingleSiteSim(int generations, vector< shared_ptr<Organism> > inputPopulation);

		vector<string> getPopDna() const;
		vector<string> getPopProtein() const;
		vector<string> getPopPhenotypes() const;

	private:

		void calcPopulationStats();
};

/* Header guard end.
 */
#endif
