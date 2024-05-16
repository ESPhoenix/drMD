## BASIC LIBS
import os
from os import path as p
## OPEN MM LIBS
import simtk.openmm.app as app
from simtk.openmm.app import metadynamics
import simtk.openmm as openmm
import  simtk.unit  as unit
## CUSTOM LIBS
import instruments.drSim as drSim
import instruments.drConstraints as drConstraints
import instruments.drCheckup as drCheckup
from pdbUtils import pdbUtils
########################################################################################################
def run_metadynamics(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile):
    print("Running MetaDynamics!")
    ## initilaise new system
    system = drSim.init_system(prmtop)
    ## deal with restraints (clear all lurking restraints and constants)
    system, clearRestraints = drConstraints.constraints_handler(system, prmtop, inpcrd, sim, saveFile)
    # Add a Monte Carlo Barostat to maintain constant pressure
    barostat = openmm.MonteCarloBarostat(1.0*unit.atmospheres, 300*unit.kelvin)  # Set pressure and temperature
    system.addForce(barostat)

    ## create bias variable
    biasVariables = []
    if sim["biasVar"].upper() == "RMSD":
        biasVariable = gen_rmsd_bias_variable(prmtop, inpcrd, sim)
        biasVariables.append(biasVariable)
    elif sim["biasVar"].upper() == "DIHEDRAL":
        biasVariable = gen_dihedral_bias_variable(prmtop, inpcrd, sim, pdbFile)
        biasVariables.append(biasVariable)

    ## create metadynamics object -- adds bias variable as a force to the system
    meta = metadynamics.Metadynamics(system = system,
                                     variables = biasVariables,
                                     temperature = sim["temp"],
                                     biasFactor = 5,
                                     height = 1,
                                     frequency = 50,
                                     saveFrequency = 50,
                                     biasDir = simDir)
    ## set up intergrator
    integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    ## create new simulation
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # load state from previous simulation (or continue from checkpoint)
    simulation = drSim.load_simulation_state(simulation, saveFile)
    ## set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    reporters = drSim.init_reporters(simDir = simDir,
                                nSteps =  totalSteps,
                                nLogSteps = sim["nLogSteps"])
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])
    ## run metadynamics simulation
    meta.step(simulation, sim["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"Meta_step.xml")
    simulation.saveState(saveXml)
    # save result as pdb
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir,"Meta_final_geom.pdb"), 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    ## run drCheckup to assess simulation health (will this work for METADYNAMICS???)
    drCheckup.check_vitals(simDir,reporters["vitals"][1], reporters["progress"][1])

    
    ## return checkpoint file for continuing simulation
    return saveXml
########################################################################################################
def gen_dihedral_bias_variable(prmtop, inpcrd, sim, pdbFile):
    topology = prmtop.topology
    positions = inpcrd.positions
    ## get atom indecies
    dihedralAtomsInput = sim["biasAtoms"]
    atomIndecies = []

    pdbDf = pdbUtils.pdb2df(pdbFile)
    for atomInfo in dihedralAtomsInput:
        atomDf = pdbDf[(pdbDf["CHAIN_ID"].astype(str) == atomInfo[0]) &
                        (pdbDf["RES_NAME"].astype(str) == atomInfo[1]) &
                        (pdbDf["RES_ID"].astype(str) == atomInfo[2]) &
                       (pdbDf["ATOM_NAME"].astype(str) == atomInfo[3]) ]
        atomIndex  = atomDf.index.values[0]
        atomIndecies.append(atomIndex)

    ## generate a dihedral bias force:
    dihedralForce = openmm.CustomTorsionForce("theta")

    dihedralForce.addTorsion(atomIndecies[0],
                              atomIndecies[1],
                                atomIndecies[2],
                                  atomIndecies[3])

    dihedralBiasVariable = metadynamics.BiasVariable(force = dihedralForce,
                                                     minValue = -180 * unit.degrees,
                                                     maxValue = 180 * unit.degrees,
                                                     biasWidth = 1 * unit.degrees,
                                                     periodic = True)
    return dihedralBiasVariable




########################################################################################################
def gen_rmsd_bias_variable(prmtop, inpcrd, sim):
    topology = prmtop.topology
    positions = inpcrd.positions

    ## we need to pass a list of atom indecies to create our RMSDForce object
    ##  we need to also pass a list of all original atom coords
    ## depending on our user inputs stored in sim, we will select different atoms

    ## get all atom coords
    atomCoords = []
    for i, atom in enumerate(topology.atoms()):
        atomCoords.append(positions[i])

    atomIndecies = []
    if "biasAtoms" in sim:
        biasAtoms = sim["biasAtoms"]
        for biasAtom in biasAtoms:
            for chain in topology.chains():
                if not chain == biasAtom[0]:
                    continue
                for residue in chain.residues():
                    if not residue.name == biasAtom[1]:
                        continue
                    if not residue.index == biasAtom[2]:
                        continue
                    for atom in residue.atoms():
                        if not atom.name == biasAtom[3]:
                            continue
                        atomIndecies.append(atom.index)

    if "biasResidues" in sim:
        biasResidues = sim["biasResidues"]
        for biasResidue in biasResidues:
            for chain in topology.chains():
                if not chain == biasResidue[0]:
                    continue
                for residue in chain.residues():
                    if not residue.name == biasResidue[1]:
                        continue
                    if not residue.index == biasResidue[2]:
                        continue
                    for atom in residue.atoms():
                        atomIndecies.append(atom.index)
    if "biasKeywords" in sim:
        biasKeywords = sim["biasKeywords"]
        for biasKeyword in biasKeywords:
            if biasKeyword == "backbone":
                backboneAtomNames = ["N","CA","C","O"]
                for atom in topology.atoms():
                    if atom.name in backboneAtomNames:
                        atomIndecies.append(i)
    ## deal with repeated indecies from selection
    atomIndecies = list(set(atomIndecies))

    
    rmsdForce = openmm.RMSDForce(atomCoords, atomIndecies)
    rmsdBiasVariable = metadynamics.BiasVariable(force = rmsdForce,
                                                 minValue = 10 * unit.angstrom,
                                                 maxValue = 100 * unit.angstrom,
                                                 biasWidth = 1 * unit.angstrom,
                                                 periodic = False)
    return rmsdBiasVariable



########################################################################################################
def add_funnel_potential(system, ):
    """https://github.com/jeff231li/funnel_potential"""
