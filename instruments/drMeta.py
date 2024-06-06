## BASIC LIBS
import os
from os import path as p
## OPEN MM LIBS
import simtk.openmm.app as app
from simtk.openmm.app import metadynamics
import simtk.openmm as openmm
import  simtk.unit  as unit
## CUSTOM drMD LIBS
from instruments import drSim
from instruments import drRestraints
from instruments import drCheckup 
from instruments import drClusterizer
from instruments import drSelector
## generic pdb <-> df utils
from pdbUtils import pdbUtils
########################################################################################################
def run_metadynamics(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile):
    print("Running MetaDynamics!")
    ## initilaise new system
    system = drSim.init_system(prmtop)
    ## deal with restraints (clear all lurking restraints and constants)
    system = drRestraints.constraints_handler(system, prmtop, inpcrd, sim, saveFile, pdbFile)
    # Add a Monte Carlo Barostat to maintain constant pressure
    barostat = openmm.MonteCarloBarostat(1.0*unit.atmospheres, 300*unit.kelvin)  # Set pressure and temperature
    system.addForce(barostat)

    ## read biases from sim config and create bias variables
    biases = sim["biases"]
    biasVariables = []
    for bias in biases:
        print(bias)
        print(bias["selection"])
        atomIndexes = drSelector.get_atom_indexes(bias["selection"], pdbFile)
        atomCoords = get_atom_coords_for_metadynamics(prmtop, inpcrd)
        ## create bias variable
        if bias["biasVar"].upper() == "RMSD":
            biasVariable = gen_rmsd_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)

        if bias["biasVar"].upper() == "DIHEDRAL":
            biasVariable = gen_dihedral_bias_variable(bias, atomCoords, atomIndexes)
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
    # set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    reportInterval = sim["logInterval"]
    reporters = drSim.init_reporters(simDir = simDir,
                                nSteps =  totalSteps,
                                reportInterval= reportInterval)
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

    ## run trajectory clustering 
    if "clusterTrajectory" in sim:
        if sim["clusterTrajectory"]["clusterBool"]:
            drClusterizer.rmsd_clustering_protocol(simDir, sim["clusterTrajectory"])

    ## return checkpoint file for continuing simulation
    return saveXml
########################################################################################################
def get_atom_coords_for_metadynamics(prmtop, inpcrd):
    topology = prmtop.topology
    positions = inpcrd.positions
    ## get all atom coords
    atomCoords = []
    for i, atom in enumerate(topology.atoms()):
        atomCoords.append(positions[i])
    return atomCoords


########################################################################################################
def gen_dihedral_bias_variable(bias, atomCoords, atomIndexes):


    ## generate a dihedral bias force:
    dihedralForce = openmm.CustomTorsionForce("theta")

    dihedralForce.addTorsion(atomIndexes[0],
                              atomIndexes[1],
                                atomIndexes[2],
                                  atomIndexes[3])

    dihedralBiasVariable = metadynamics.BiasVariable(force = dihedralForce,
                                                     minValue = -180 * unit.degrees,
                                                     maxValue = 180 * unit.degrees,
                                                     biasWidth = 1 * unit.degrees,
                                                     periodic = True)
    return dihedralBiasVariable




########################################################################################################
def gen_rmsd_bias_variable(bias, atomCoords, atomIndexes):
  
    ## generate a rmsd bias variable
    rmsdForce = openmm.RMSDForce(atomCoords, atomIndexes)
    rmsdBiasVariable = metadynamics.BiasVariable(force = rmsdForce,
                                                 minValue = 10 * unit.angstrom,
                                                 maxValue = 100 * unit.angstrom,
                                                 biasWidth = 1 * unit.angstrom,
                                                 periodic = False)
    return rmsdBiasVariable



########################################################################################################
def add_funnel_potential(system, ):
    """https://github.com/jeff231li/funnel_potential"""
