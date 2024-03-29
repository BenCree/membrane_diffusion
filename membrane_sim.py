from openmm.app import ForceField, Modeller, PME, Simulation, PDBFile, PDBReporter, StateDataReporter, PDBReporter
from openmm import unit, Platform
from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule, Topology
from rdkit import Chem
from rdkit.Chem import AllChem
import mdtraj as md
import openmm as mm
import numpy as np
import sys
 
mol = Chem.MolFromSmiles('OCCCCO')
hmol = Chem.AddHs(mol)
AllChem.EmbedMolecule(hmol)
 
ofmol = Molecule.from_rdkit(hmol)
ff = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
 
smirnoff = SMIRNOFFTemplateGenerator(ofmol)
ff.registerTemplateGenerator(smirnoff.generator)
 
md_ligand_top = ofmol.to_topology()
omm_top = md_ligand_top.to_openmm()
pos = ofmol.conformers[0]
print(pos)
print(type(pos))
pos = pos.to('nanometers')
print(pos)
print(type(pos))
#pos = unit.Quantity(np.zeros([len(ofmol.atoms), 3]), unit=unit.nanometers)
#md_ligand_top = md.Topology.from_openmm(omm_top)
#md_ligand_top = md_ligand_top.to_openmm()
md_ligand_pos = pos.to_openmm()
md_ligand_pos = md_ligand_pos + [3,3,5] * unit.nanometers
 
modeller = Modeller(omm_top, md_ligand_pos)
modeller.addMembrane(ff,
                    lipidType='POPC',
                    #minimumPadding=1.0*unit.nanometers,
                    membraneCenterZ=3.5*unit.nanometer,
                   )
#platform = Platform.getPlatformByName('CUDA') 
system = ff.createSystem(modeller.topology, nonbondedMethod=PME)
integrator = mm.LangevinIntegrator(
    300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
)
simulation = Simulation(modeller.topology, system, integrator,)
simulation.context.setPositions(modeller.positions)
 
simulation.minimizeEnergy()
with open( "mol_topology.pdb", "w") as pdb_file:
    PDBFile.writeFile(
        simulation.topology,
        simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
        file=pdb_file,
        keepIds=True,
    )
steps = 500000  # corresponds to 20 fs
write_interval = 100  # write every 2 fs
log_interval = 1  # log progress to stdout every 2 fs
simulation.reporters.append(
    md.reporters.XTCReporter(file=str( "trajectory.xtc"), reportInterval=write_interval)
)
simulation.reporters.append(
    StateDataReporter(
        sys.stdout,
        log_interval,
        step=True,
        potentialEnergy=True,
        temperature=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=steps,
        separator="\t",
    )
)
simulation.reporters.append(PDBReporter('ligand_md.pdb', 1000))
 
 
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
simulation.step(steps)  # perform the simulation
