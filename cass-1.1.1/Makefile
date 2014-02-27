#
# Options
#
VPATH = ./randomgen
CC = g++
CFLAGS = -c -std=c++0x

#
# Compile and link everything
#
all: sh2-evo-decoys sh2-evo-neo engfunc-nodecoys bindfunc-nodecoys to-bead samp-seq motifs complex termdist mcmc-opt structeng pos

# 
# sh2-evolution-decoys target
#
sh2-evo-decoys-objects = sh2-evolution-decoys.o AllAtomStructure.o TwoBeadStructure.o Structure.o Residue.o myUtils.o mersenne.o stoc1.o userintf.o StructureProposal.o Proposal.o ParameterRotation.o Parameter.o MHMCMC.o MCMC.o SingleSiteSim.o Simulation.o MultipleStructureOrganism.o Organism.o GrahnenModel.o RastogiModel.o MukherjeeModel.o BastollaModel.o BastollaAugmentedModel.o EnergyModel.o

sh2-evo-decoys: $(sh2-evo-decoys-objects)
	$(CC) -o sh2-evolution-decoys $(sh2-evo-decoys-objects)

sh2-evolution-decoys.o: sh2-evolution-decoys.cpp AllAtomStructure.h TwoBeadStructure.h SingleSiteSim.h MultipleStructureOrganism.h BastollaModel.h GrahnenModel.h ParameterRotation.h MHMCMC.h
	$(CC) $(CFLAGS) $<

#
# sh2-evolution-neo target
#
sh2-evo-neo-objects = sh2-evolution-decoys+neo.o AllAtomStructure.o TwoBeadStructure.o Structure.o Residue.o myUtils.o mersenne.o stoc1.o userintf.o StructureProposal.o Proposal.o ParameterRotation.o Parameter.o MHMCMC.o MCMC.o SingleSiteSim.o Simulation.o MultipleStructureOrganism.o MultipleStructureOrganismNeo.o Organism.o GrahnenModel.o RastogiModel.o MukherjeeModel.o EnergyModel.o

sh2-evo-neo: $(sh2-evo-neo-objects)
	$(CC) -o sh2-evolution-decoys+neo $(sh2-evo-neo-objects)

sh2-evolution-decoys+neo.o: sh2-evolution-decoys+neo.cpp AllAtomStructure.h TwoBeadStructure.h SingleSiteSim.h MultipleStructureOrganismNeo.h GrahnenModel.h ParameterRotation.h MHMCMC.h
	$(CC) $(CFLAGS) $<



#
# test-engfunc-nodecoys 
#
engfunc-nodecoys-objects = test-engfunc-nodecoys.o Structure.o Residue.o EnergyModel.o AllAtomStructure.o TwoBeadStructure.o BastollaModel.o BastollaAugmentedModel.o MukherjeeModel.o RastogiModel.o GrahnenModel.o mersenne.o stoc1.o userintf.o myUtils.o MHMCMC.o MCMC.o StructureProposal.o Proposal.o ParameterRotation.o Parameter.o

engfunc-nodecoys: $(engfunc-nodecoys-objects)
	$(CC) -o test-engfunc-nodecoys $(engfunc-nodecoys-objects)

test-engfunc-nodecoys.o: test-engfunc-nodecoys.cpp TwoBeadStructure.h BastollaAugmentedModel.h GrahnenModel.h randomc.h MHMCMC.h ParameterRotation.h
	$(CC) $(CFLAGS) test-engfunc-nodecoys.cpp

#
# test-bindfunc-nodecoys
#
bindfunc-nodecoys-objects = test-bindfunc-nodecoys.o Structure.o Residue.o EnergyModel.o AllAtomStructure.o TwoBeadStructure.o BastollaModel.o BastollaAugmentedModel.o MukherjeeModel.o RastogiModel.o GrahnenModel.o mersenne.o stoc1.o userintf.o myUtils.o MHMCMC.o MCMC.o StructureProposal.o Proposal.o ParameterRotation.o Parameter.o

bindfunc-nodecoys: $(bindfunc-nodecoys-objects)
	$(CC) -o test-bindfunc-nodecoys $(bindfunc-nodecoys-objects)

test-bindfunc-nodecoys.o: test-bindfunc-nodecoys.cpp AllAtomStructure.h TwoBeadStructure.h BastollaModel.h GrahnenModel.h randomc.h stocc.h MHMCMC.h ParameterRotation.h
	$(CC) $(CFLAGS) test-bindfunc-nodecoys.cpp

#
# PDB-to-2bead converter
#
to-bead-objects = pdb2bead.o AllAtomStructure.o TwoBeadStructure.o Structure.o Residue.o myUtils.o mersenne.o stoc1.o userintf.o MHMCMC.o MCMC.o StructureProposal.o Proposal.o ParameterRotation.o Parameter.o GrahnenModel.o RastogiModel.o MukherjeeModel.o EnergyModel.o

to-bead: $(to-bead-objects)
	$(CC) -o pdb2bead $(to-bead-objects)

pdb2bead.o: pdb2bead.cpp AllAtomStructure.h TwoBeadStructure.h
	$(CC) $(CFLAGS) $<

#
# Scoring calculator for extant structure
#
structeng-objects = engs-from-struct.o TwoBeadStructure.o AllAtomStructure.o Structure.o Residue.o myUtils.o mersenne.o stoc1.o userintf.o MHMCMC.o MCMC.o StructureProposal.o Proposal.o ParameterRotation.o Parameter.o BastollaModel.o BastollaAugmentedModel.o GrahnenModel.o RastogiModel.o MukherjeeModel.o EnergyModel.o

structeng: $(structeng-objects)
	$(CC) -o engs-from-struct $(structeng-objects)

engs-from-struct.o: engs-from-struct.cpp BastollaModel.h BastollaAugmentedModel.h GrahnenModel.h AllAtomStructure.h TwoBeadStructure.h myUtils.h
	$(CC) $(CFLAGS) $<

#
# Sequence sampler based on energy
#
samp-seq-objects = sample-seqs.o TwoBeadStructure.o AllAtomStructure.o Structure.o Residue.o myUtils.o mersenne.o stoc1.o userintf.o MHMCMC.o MCMC.o SequenceProposal.o GoSequenceProposal.o StructureProposal.o Proposal.o ParameterRotation.o ParameterSequence.o Parameter.o GrahnenModel.o RastogiModel.o MukherjeeModel.o BastollaAugmentedModel.o BastollaModel.o EnergyModel.o

samp-seq: $(samp-seq-objects)
	$(CC) -o sample-seqs $(samp-seq-objects)

sample-seqs.o: sample-seqs.cpp GrahnenModel.h AllAtomStructure.h TwoBeadStructure.h MHMCMC.h SequenceProposal.h ParameterSequence.h ParameterRotation.h
	$(CC) $(CFLAGS) $<

#
# Ligand motif energy sampler
#
motifs-objects = test-motifs.o TwoBeadStructure.o AllAtomStructure.o Structure.o Residue.o myUtils.o mersenne.o stoc1.o userintf.o MHMCMC.o MCMC.o StructureProposal.o Proposal.o ParameterRotation.o Parameter.o GrahnenModel.o RastogiModel.o MukherjeeModel.o BastollaModel.o BastollaAugmentedModel.o EnergyModel.o

motifs: $(motifs-objects)
	$(CC) -o test-motifs $(motifs-objects)

test-motifs.o: test-motifs.cpp AllAtomStructure.h TwoBeadStructure.h GrahnenModel.h BastollaAugmentedModel.h
	$(CC) $(CFLAGS) $<

#
# PPI energy minimizer
#
complex-objects = minimize-complex.o TwoBeadStructure.o Structure.o Residue.o myUtils.o mersenne.o stoc1.o userintf.o MHMCMC.o MCMC.o StructureProposal.o Proposal.o ParameterRotation.o Parameter.o GrahnenModel.o RastogiModel.o MukherjeeModel.o EnergyModel.o

complex: $(complex-objects)
	$(CC) -o minimize-complex $(complex-objects)

minimize-complex.o: minimize-complex.cpp TwoBeadStructure.h MHMCMC.h ParameterRotation.h myUtils.h
	$(CC) $(CFLAGS) $<

#
# Term distribution calculator
#
termdist-objects = calc-termdists.o Structure.o Residue.o EnergyModel.o AllAtomStructure.o TwoBeadStructure.o GrahnenModel.o RastogiModel.o MukherjeeModel.o StructureProposal.o mersenne.o stoc1.o userintf.o myUtils.o Parameter.o ParameterRotation.o Proposal.o MCMC.o MHMCMC.o 

termdist: $(termdist-objects)
	$(CC) -o calc-termdists $(termdist-objects)

calc-termdists.o: calc-termdists.cpp GrahnenModel.h AllAtomStructure.h TwoBeadStructure.h ParameterRotation.h MHMCMC.h 
	$(CC) $(CFLAGS) $<

#
# pos.cpp by Nadia
#
pos-objects = pos.o Structure.o Residue.o EnergyModel.o AllAtomStructure.o TwoBeadStructure.o GrahnenModel.o RastogiModel.o MukherjeeModel.o StructureProposal.o mersenne.o stoc1.o userintf.o myUtils.o Parameter.o ParameterRotation.o Proposal.o MCMC.o MHMCMC.o 

pos: $(pos-objects)
	$(CC) -o pos $(pos-objects)

pos.o: pos.cpp GrahnenModel.h AllAtomStructure.h TwoBeadStructure.h ParameterRotation.h MHMCMC.h Structure.h
	$(CC) $(CFLAGS) $<

#
# Weight optimization, using MCMC
#
mcmc-opt-objects = run-mcmc-engparams.o MHMCMC.o MCMC.o myUtils.o mersenne.o stoc1.o userintf.o ParameterFloat.o Parameter.o EnergyFuncProposal.o Proposal.o

mcmc-opt: $(mcmc-opt-objects)
	$(CC) -o run-mcmc-engparams $(mcmc-opt-objects)

run-mcmc-engparams.o: run-mcmc-engparams.cpp myUtils.h MHMCMC.h ParameterFloat.h EnergyFuncProposal.h
	$(CC) $(CFLAGS) $<

#
# Phony targets
#
clean:
	rm -f *.o

# Targets that are part of the simulation system
AllAtomStructure.o: AllAtomStructure.cpp AllAtomStructure.h TwoBeadStructure.h OneBeadStructure.h
	$(CC) $(CFLAGS) $<

TwoBeadStructure.o: TwoBeadStructure.cpp TwoBeadStructure.h AllAtomStructure.h OneBeadStructure.h
	$(CC) $(CFLAGS) $<

OneBeadStructure.o: OneBeadStructure.cpp OneBeadStructure.h AllAtomStructure.h TwoBeadStructure.h
	$(CC) $(CFLAGS) $<

Structure.o: Structure.cpp Structure.h
	$(CC) $(CFLAGS) $<

Residue.o: Residue.cpp Residue.h
	$(CC) $(CFLAGS) $<

myUtils.o: myUtils.cpp myUtils.h
	$(CC) $(CFLAGS) $<

BastollaModel.o: BastollaModel.cpp Structure.h AllAtomStructure.h BastollaModel.h 
	$(CC) $(CFLAGS) $<

BastollaAugmentedModel.o: BastollaAugmentedModel.cpp Structure.h AllAtomStructure.h TwoBeadStructure.h BastollaAugmentedModel.h 
	$(CC) $(CFLAGS) $<

MukherjeeModel.o: MukherjeeModel.cpp Structure.h TwoBeadStructure.h MukherjeeModel.h
	$(CC) $(CFLAGS) $<

RastogiModel.o: RastogiModel.cpp Structure.h TwoBeadStructure.h RastogiModel.h
	$(CC) $(CFLAGS) $<

GrahnenModel.o: GrahnenModel.cpp Structure.h TwoBeadStructure.h GrahnenModel.h
	$(CC) $(CFLAGS) $<

EnergyModel.o: EnergyModel.cpp EnergyModel.h
	$(CC) $(CFLAGS) $<

Simulation.o: Simulation.cpp Simulation.h Organism.h
	$(CC) $(CFLAGS) $<

SingleSiteSim.o: SingleSiteSim.cpp SingleSiteSim.h Organism.h
	$(CC) $(CFLAGS) $<

Organism.o: Organism.cpp Organism.h myUtils.h Structure.h
	$(CC) $(CFLAGS) $<

MultipleStructureOrganism.o: MultipleStructureOrganism.cpp MultipleStructureOrganism.h Structure.h EnergyModel.h
		$(CC) $(CFLAGS) $<

MultipleStructureOrganismNeo.o: MultipleStructureOrganismNeo.cpp MultipleStructureOrganismNeo.h Structure.h EnergyModel.h
	$(CC) $(CFLAGS) $<

#
# Targets for random number generation
#
mersenne.o: mersenne.cpp randomc.h
	$(CC) $(CFLAGS) $<

stoc1.o: stoc1.cpp stocc.h
	$(CC) $(CFLAGS) $<

userintf.o: userintf.cpp
	$(CC) $(CFLAGS) $<

#
# Targets for parameter optimization
#
Proposal.o: Proposal.cpp Proposal.h Parameter.h
	$(CC) $(CFLAGS) $<

Parameter.o: Parameter.cpp Parameter.h 
	$(CC) $(CFLAGS) $<

ParameterRotation.o: ParameterRotation.cpp ParameterRotation.h randomc.h stocc.h
	$(CC) $(CFLAGS) $<

ParameterFloat.o: ParameterFloat.cpp ParameterFloat.h randomc.h stocc.h
	$(CC) $(CFLAGS) $<

ParameterSequence.o: ParameterSequence.cpp ParameterSequence.h randomc.h myUtils.h
	$(CC) $(CFLAGS) $<

MCMC.o: MCMC.cpp MCMC.h Proposal.h
	$(CC) $(CFLAGS) $<

MHMCMC.o: MHMCMC.cpp MHMCMC.h Proposal.h
	$(CC) $(CFLAGS) $<

StructureProposal.o: StructureProposal.cpp StructureProposal.h Residue.h TwoBeadStructure.h Parameter.h ParameterRotation.h GrahnenModel.h myUtils.h 
	$(CC) $(CFLAGS) $<

SequenceProposal.o: SequenceProposal.cpp SequenceProposal.h Parameter.h EnergyModel.h Structure.h
	$(CC) $(CFLAGS) $<

GoSequenceProposal.o: GoSequenceProposal.cpp GoSequenceProposal.h Parameter.h
	$(CC) $(CFLAGS) $<

EnergyFuncProposal.o: EnergyFuncProposal.cpp EnergyFuncProposal.h Parameter.h
	$(CC) $(CFLAGS) $<
