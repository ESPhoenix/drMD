## BASIC PYTHON LIBRARIES
import os
from os import path as p
import numpy as np

## OPENMM LIBRARIES
import openmm.app as app
from openmm.app import metadynamics
import openmm as openmm
import  simtk.unit  as unit

## drMD LIBRARIES
from Surgery import drSim, drRestraints, drFirstAid
from ExaminationRoom import drLogger, drCheckup
from UtilitiesCloset import drSelector

########################################################################################################
##TODO: allow for custom parameters to be passed from config file ---> bias force generation and metadynamics construction

########################################################################################################
@drLogger.monitor_progress_decorator()
@drFirstAid.firstAid_handler(drFirstAid.run_firstAid_energy_minimisation)
@drCheckup.check_up_handler()
def run_metadynamics(prmtop: app.Topology,
                      inpcrd: any,
                        sim: dict,
                          saveFile: str,
                            outDir: str,
                              platform: openmm.Platform,
                                refPdb: str,
                                config:dict) -> str:
    """
    Run a simulation at constant pressure (NpT) step with biases.

    Args:
        prmtop (str): The path to the topology file.
        inpcrd (str): The path to the coordinates file.
        sim (dict): The simulation parameters.
        saveFile (str): The path to the checkpoint or XML file.
        simDir (str): The path to the simulation directory.
        platform (openmm.Platform): The simulation platform.
        pdbFile (str): The path to the PDB file.

    Returns:
        str: The path to the XML file containing the final state of the simulation.

    This function runs a simulation at constant pressure (NpT) step with biases.
    It initializes the system, handles any constraints, adds a constant pressure force,
    sets up the integrator and system, loads the state from a checkpoint or XML file,
    sets up reporters, runs the simulation, saves the final geometry as a PDB file,
    resets the chain and residue Ids, runs drCheckup, performs trajectory clustering
    if specified in the simulation parameters, and saves the simulation state as an
    XML file.
    """
    stepName = sim["stepName"]
    drLogger.log_info(f"-->{' '*4}Running MetaDynamics Step: {stepName}",True)
    ## make a simulation directory
    simDir: str = p.join(outDir, stepName)
    os.makedirs(simDir, exist_ok=True)


    sim = drSim.process_sim_data(sim)
    # Define the nonbonded method and cutoff.
    nonbondedMethod: openmm.NonbondedForce = app.PME
    nonbondedCutoff: unit.Quantity = 1 * unit.nanometer

    # Define the restraints.
    hBondconstraints: openmm.Force = app.HBonds

    # Create the system.
    system: openmm.System = prmtop.createSystem(nonbondedMethod=nonbondedMethod,
                                                nonbondedCutoff=nonbondedCutoff,
                                                constraints=hBondconstraints)
    # Deal with restraints (clear all lurking restraints and constants)
    system: openmm.System = drRestraints.restraints_handler(system, prmtop, inpcrd, sim, saveFile, refPdb)
    # Add a Monte Carlo Barostat to maintain constant pressure
    barostat: openmm.MonteCarloBarostat = openmm.MonteCarloBarostat(1.0*unit.atmospheres, sim["temperature"])  # Set pressure and temperature
    system.addForce(barostat)
    # Read metaDynamicsInfo from sim config
    metaDynamicsInfo: dict = sim["metaDynamicsInfo"]

    # Read biases from sim config and create bias variables
    biases: list = metaDynamicsInfo["biases"]
    biasVariables: list = []
    for bias in biases:
        # Get atom indexes and coordinates for the biases
        atomIndexes: list = drSelector.get_atom_indexes(bias["selection"], refPdb)
        atomCoords: list = get_atom_coords_for_metadynamics(prmtop, inpcrd)
        # Create bias variable based on the type of bias
        if bias["biasVar"].upper() == "RMSD":
            biasVariable: metadynamics.BiasVariable = gen_rmsd_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)
        elif bias["biasVar"].upper() == "DIHEDRAL":
            biasVariable: metadynamics.BiasVariable = gen_dihedral_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)
        elif bias["biasVar"].upper() == "DISTANCE":
            biasVariable: metadynamics.BiasVariable = gen_distance_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)
        elif bias["biasVar"].upper() == "ANGLE":
            biasVariable: metadynamics.BiasVariable = gen_angle_bias_variable(bias, atomCoords, atomIndexes)
            biasVariables.append(biasVariable)
    # Create metadynamics object and add bias variables as forces to the system
    meta: metadynamics.Metadynamics = metadynamics.Metadynamics(system=system,
                                     variables=biasVariables,
                                     temperature=sim["temperature"],
                                     biasFactor=metaDynamicsInfo["biasFactor"],
                                     height=metaDynamicsInfo["height"],
                                     frequency=50,
                                     saveFrequency=50,
                                     biasDir=simDir)
    # Set up integrator
    integrator: openmm.LangevinMiddleIntegrator = openmm.LangevinMiddleIntegrator(sim["temperature"], 1/unit.picosecond, sim["timestep"])
    # Create new simulation
    simulation: app.Simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # Load state from previous simulation (or continue from checkpoint)
    simulation: app.Simulation = drSim.load_simulation_state(simulation, saveFile)
    # Set up reporters
    totalSteps: int = simulation.currentStep + sim["nSteps"]
    reportInterval: int = sim["logInterval"]
    simulation: app.Simulation = drSim.init_reporters(simDir=simDir,
                                     nSteps=totalSteps,
                                     reportInterval=reportInterval,
                                       simulation=simulation)
    # Run metadynamics simulation
    meta.step(simulation, sim["nSteps"])
    # Save simulation as XML
    saveXml: str = p.join(simDir, f"{stepName}.xml")
    simulation.saveState(saveXml)
    # Save result as pdb
    state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir, f"{stepName}.pdb"), 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)

    # Return checkpoint file for continuing simulation
    return saveXml
########################################################################################################
def get_atom_coords_for_metadynamics(prmtop: app.Topology, inpcrd: any) -> list:
    """
    Get the coordinates of all atoms in the system.

    Parameters
    ----------
    prmtop : file
        The topology file.
    inpcrd : file
        The coordinates file.

    Returns
    -------
    atomCoords : list
        The coordinates of all atoms in the system.
    """
    # Get the topology and positions from the input files
    topology: app.Topology = prmtop.topology
    positions: np.ndarray = inpcrd.positions

    # Initialize an empty list to store atom coordinates
    atomCoords: list = []

    # Loop over all atoms in the topology
    for i, atom in enumerate(topology.atoms()):
        # Append the coordinates of the current atom to the atomCoords list
        atomCoords.append(positions[i])

    # Return the list of atom coordinates
    return atomCoords
########################################################################################################
def gen_angle_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate an angle bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the angle parameters.
    atomCoords : np.ndarray
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the angle calculation.

    Returns
    -------
    angleBiasVariable : openmm.CustomTorsionForce
        The angle bias variable.
    """

    # Create an angle bias force
    angleForce: openmm.CustomAngleForce = openmm.CustomAngleForce("theta")

    # Add the atom indexes for the angle
    angleForce.addAngle(atomIndexes[0],
                          atomIndexes[1],
                          atomIndexes[2])
    
    angleBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(force = angleForce,
                                                     minValue = 0 * unit.degrees,
                                                     maxValue = 180 * unit.degrees,
                                                     biasWidth = 1 * unit.degrees,
                                                     periodic = False)
    return angleBiasVariable

########################################################################################################
def gen_dihedral_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate a dihedral bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the dihedral angle parameters.
    atomCoords : list
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the dihedral angle.

    Returns
    -------
    dihedralBiasVariable : openmm.CustomTorsionForce
        The dihedral bias variable.
    """
    ## express dihedral as a harmonic potential
    dihedralEnergyExpression: str = "theta"
    
    ## generate a dihedral bias force:
    # Create a custom torsion force object
    dihedralForce: openmm.CustomTorsionForce = openmm.CustomTorsionForce(dihedralEnergyExpression)

    
    # Add the atom indexes for the dihedral angle
    dihedralForce.addTorsion(atomIndexes[0],
                              atomIndexes[1],
                                atomIndexes[2],
                                  atomIndexes[3])

    # Create a dihedral bias variable
    dihedralBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(force = dihedralForce,
                                                     minValue = -179 * unit.degrees,
                                                     maxValue = 179 * unit.degrees,
                                                     biasWidth = 1 * unit.degrees,
                                                     periodic = True)
    
    return dihedralBiasVariable

########################################################################################################
def gen_distance_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate a distance bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the distance parameters.
    atomCoords : np.ndarray
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the distance calculation.

    Returns
    -------
    distanceBiasVariable : openmm.CustomTorsionForce
        The distance bias variable.
    """

    # Create a distance bias force
    distanceForce: openmm.CustomTorsionForce = openmm.CustomBondForce("0.5*k*(r-r0)^2")
    distanceForce.addPerBondParameter("k")
    distanceForce.addPerBondParameter("r0") 
    # Add the atom indexes for the distance
    distanceForce.addBond(atomIndexes[0],
                             atomIndexes[1])
    
    distanceBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(force = distanceForce,
                                                     minValue = 0 * unit.angstrom,
                                                     maxValue = 100 * unit.angstrom,
                                                     biasWidth = 1 * unit.angstrom,
                                                     periodic = False)
    
    return distanceBiasVariable

########################################################################################################
def gen_rmsd_bias_variable(bias: dict, atomCoords: np.ndarray, atomIndexes: list) -> metadynamics.BiasVariable:
    """
    Generate a RMSD (Root Mean Square Deviation) bias variable.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing the RMSD parameters.
    atomCoords : np.ndarray
        The coordinates of all atoms in the system.
    atomIndexes : list
        The indexes of the atoms involved in the RMSD calculation.

    Returns
    -------
    rmsdBiasVariable : openmm.RMSDForce
        The RMSD bias variable.
    """

    # Generate a RMSD bias force
    rmsdForce: openmm.RMSDForce = openmm.RMSDForce(atomCoords, atomIndexes)

    # Create a RMSD bias variable
    rmsdBiasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(
        force=rmsdForce,
        minValue=10 * unit.angstrom,
        maxValue=100 * unit.angstrom,
        biasWidth=1 * unit.angstrom,
        periodic=False)

    return rmsdBiasVariable



########################################################################################################
##TODO implement funnel potential
# def add_funnel_potential(system, ):
#     """https://github.com/jeff231li/funnel_potential"""
