## BASIC LIBS
import os
from os import path as p
## OPEN MM LIBS
import openmm.app as app
import openmm as openmm
import  simtk.unit  as unit
## CUSTOM LIBS
from  instruments import drRestraints
from  instruments import drCheckup 
from  instruments import drClusterizer
from  instruments import drFixer
from  instruments import drFirstAid
###########################################################################################
def init_system(prmtop: app.Topology) -> openmm.System:
    """
    Create an openmm.System from a prmtop.

    Parameters
    ----------
    prmtop : app.Topology
        The topology of the system.

    Returns
    -------
    system : openmm.System
        The initialized system.

    """

    # Define the nonbonded method and cutoff.
    nonbondedMethod: openmm.NonbondedForce = app.PME
    nonbondedCutoff: unit.Quantity = 1 * unit.nanometer

    # Define the restraints.
    hBondconstraints: openmm.Force = app.HBonds

    # Create the system.
    system: openmm.System = prmtop.createSystem(nonbondedMethod=nonbondedMethod,
                                                nonbondedCutoff=nonbondedCutoff,
                                                constraints=hBondconstraints)

    return system

###########################################################################################
def process_sim_data(sim: dict) -> dict:
    """
    Reads simulation variables from config file and processes them.

    Parameters
    ----------
    sim : dict
        Dictionary containing simulation variables.

    Returns
    -------
    sim : dict
        Updated dictionary with processed simulation variables.
    """

    if "processed" in sim:
        return sim

    # Set up unit translator
    timescale: dict = {"fs":unit.femtoseconds,
                 "ps":unit.picoseconds,
                 "ns":unit.nanoseconds}
    
    # Read timestep data and process it
    timestepData = sim["timestep"].split()
    timestep: unit.Quantity = float(timestepData[0]) * timescale[timestepData[1]]

    # Read duration data and process it
    durationData = sim["duration"].split()
    duration: int = int(durationData[0]) * timescale[durationData[1]]

    # Read log interval data and process it
    logIntervalData = sim["logInterval"].split()
    logInterval: unit.Quantity = float(logIntervalData[0]) * timescale[logIntervalData[1]]
    logIntervalInSteps: int = int(round(logInterval / timestep))


    # Calculate number of steps
    nSteps: int = int(duration / timestep)

    # Calculate number of log steps (every 500 steps)
    nLogSteps: int = round(nSteps / 500)

    # Process temperature
    temp: unit.Quantity = sim["temp"] * unit.kelvin

    # Update sim dictionary with processed variables
    sim.update({
        "timestep": timestep,
        "duration": duration,
        "nSteps": nSteps,
        "nLogSteps": nLogSteps,
        "temp": temp,
        "logInterval": logIntervalInSteps
    })

    sim["processed"] = True
    return sim
###########################################################################################
def init_reporters(simDir: str, nSteps: int, reportInterval) -> dict:
    """
    Initializes and returns a dictionary of reporters for a simulation.

    Args:
        simDir (str): The directory where the simulation files will be saved.
        nSteps (int): The total number of simulation steps.
        reportInterval (int): The interval at which the reporters will report.

    Returns:
        dict: A dictionary containing the initialized reporters. The keys are the names of the reporters and the values are lists containing the reporter object and the corresponding file.

    The function initializes and returns a dictionary of reporters for a simulation. It takes in the simulation directory, the total number of simulation steps, and the report interval. The reporters are initialized with the appropriate file paths and report intervals. The dictionary contains the initialized reporters, with the keys being the names of the reporters and the values being lists containing the reporter object and the corresponding file.

    Example usage:
    ```
    reporters = init_reporters("/path/to/simulation", 1000, 100)
    ```
    """
    ## as well as temperature, volume and density
    vitalsCsv = p.join(simDir, "vitals_report.csv")
    vitalsStateReporter: app.StateDataReporter = app.StateDataReporter(file = vitalsCsv, reportInterval = reportInterval, step = True,
                                            time = True, potentialEnergy = True, kineticEnergy = True,
                                            totalEnergy = True, temperature = True, volume = True,
                                            density = True)
    ## progresCsv gives you information on how long the simulation has been running and
    ## how long is left
    progressCsv: str = p.join(simDir, "progress_report.csv")
    progressStateReporter: app.StateDataReporter = app.StateDataReporter(file = progressCsv, progress = True, remainingTime = True,
                                            speed = True, elapsedTime = True, totalSteps = nSteps,
                                            reportInterval = reportInterval)
    ## dcdFile is the trajectory of the simulation
    dcdFile: str = p.join(simDir, "trajectory.dcd")
    dcdTrajectoryReporter: app.DCDReporter = app.DCDReporter(file = dcdFile, reportInterval = reportInterval, append = False)
    ## chkFile works as a checkpoint so that the simulation can be resumed if something goes wrong
    chkFile: str = p.join(simDir, "checkpoint.chk")
    chkReporter: app.CheckpointReporter = app.CheckpointReporter(file = chkFile, reportInterval = reportInterval, writeState = False)
    reporters: dict = {"vitals":[vitalsStateReporter,vitalsCsv],
                "progress": [progressStateReporter,progressCsv],
                 "trajectory": [dcdTrajectoryReporter, dcdFile],
                 "checkpoint": [chkReporter, chkFile]}

    return  reporters
###########################################################################################
def load_simulation_state(simulation: app.Simulation, saveFile: str) -> app.Simulation:
    """
    Load the simulation state from a checkpoint or XML file.

    Args:
        simulation (app.Simulation): The simulation object.
        saveFile (str): The path to the checkpoint or XML file.

    Returns:
        app.Simulation: The modified simulation object.

    This function loads the simulation state from a checkpoint or XML file.
    The function determines the file extension of the `saveFile` and calls
    the appropriate method to load the state.

    If the file extension is ".chk", the function calls the `loadCheckpoint`
    method of the simulation object.

    If the file extension is ".xml", the function calls the `loadState`
    method of the simulation object.
    """
    # Determine the file extension of the saveFile
    saveFileExt: str = p.splitext(saveFile)[1]

    # Load the simulation state from a checkpoint file
    if saveFileExt == ".chk":
        simulation.loadCheckpoint(saveFile)
    # Load the simulation state from an XML file
    elif saveFileExt == ".xml":
        simulation.loadState(saveFile)

    ## reset time and step count of simulation to zero
    simulation.context.setTime(0.0)
    simulation.context.setStepCount(0)

    # Return the modified simulation object
    return simulation
###########################################################################################
@drFirstAid.firstAid_handler(drFirstAid.run_firstAid_energy_minimisation, max_retries=10)
@drCheckup.check_up_handler()
def run_molecular_dynamics(prmtop: str, inpcrd: str, sim: dict, saveFile: str, outDir: str, platform: openmm.Platform, refPdb: str) -> str:
    """
    Run a simulation at constant pressure (NpT) step.

    Args:
        prmtop (str): The path to the topology file.
        inpcrd (str): The path to the coordinates file.
        sim (dict): The simulation parameters.
        saveFile (str): The path to the checkpoint or XML file.
        simDir (str): The path to the simulation directory.
        platform (openmm.Platform): The simulation platform.
        refPdb (str): The path to the reference PDB file.

    Returns:
        str: The path to the XML file containing the final state of the simulation.

    This function runs a simulation at constant pressure (NpT) step. The function
    initializes the system, handles any restraints, adds a constant pressure force,
    sets up the integrator and system, loads the state from a checkpoint or XML file,
    sets up reporters, runs the simulation, saves the final geometry as a PDB file,
    resets the chain and residue Ids, runs drCheckup, performs trajectory clustering
    if specified in the simulation parameters, and saves the simulation state as an
    XML file.
    """
    stepName = sim["stepName"]
    simulationType = sim["simulationType"]
    print(f"-->\tRunning {simulationType} Step:\t{stepName}")
    sim = process_sim_data(sim)

    ## create simluation directory
    simDir: str = p.join(outDir, sim["stepName"])
    os.makedirs(simDir, exist_ok=True)

    ## initialise a new system from parameters
    system: openmm.System = init_system(prmtop)
    ## deal with any restraints
    system: openmm.System = drRestraints.restraints_handler(system, prmtop, inpcrd, sim, saveFile, refPdb)
    if sim["simulationType"].upper() == "NPT":
        # add constant pressure force to system (makes this an NpT simulation)
        system.addForce(openmm.MonteCarloBarostat(1*unit.bar, sim["temp"]))
    integrator: openmm.Integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timestep"])
    simulation: app.simulation.Simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # set up intergrator and system
    # load state from previous simulation (or continue from checkpoint)
    simulation: app.Simulation = load_simulation_state(simulation, saveFile)
    # set up reporters
    totalSteps: int = simulation.currentStep + sim["nSteps"]
    reportInterval: int = sim["logInterval"]
    reporters: dict = init_reporters(simDir=simDir,
                                nSteps=totalSteps,
                                reportInterval=reportInterval)
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])
    # run NVT / NPT simulation
    simulation.step(sim["nSteps"])

    # find name to call outFiles
    protName = p.basename(p.dirname(simDir))
    # save result as pdb - reset chain and residue Ids
    state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
    nptPdb: str = p.join(simDir, f"{protName}.pdb")
    with open(nptPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    drFixer.reset_chains_residues(refPdb, nptPdb)


    # run trajectory clustering
    if "clusterTrajectory" in sim and sim["clusterTrajectory"]["clusterBool"]:
        drClusterizer.rmsd_clustering_protocol(simDir, sim["clusterTrajectory"])

    # save simulation as XML
    saveXml: str = p.join(simDir, f"{stepName}.xml")
    simulation.saveState(saveXml)
    return saveXml

###########################################################################################
def run_energy_minimisation(prmtop: str, inpcrd: str, sim: dict, outDir: str, platform: openmm.Platform, refPdb: str, saveFile: str) -> str:

    """
    Run energy minimisation on a system.

    Args:
        prmtop (str): The path to the topology file.
        inpcrd (str): The path to the coordinates file.
        sim (dict): The simulation parameters.
        simDir (str): The path to the simulation directory.
        platform (openmm.Platform): The simulation platform.
        refPdb (str): The path to the reference PDB file.

    Returns:
        str: The path to the XML file containing the final state of the simulation.

    This function runs energy minimisation on a system. The function initializes the
    system, sets up the integrator and simulation, runs the energy minimisation, saves
    the result as a PDB file, resets the chain and residue Ids, saves the simulation
    state as an XML file, and returns the path to the XML file.
    """
    stepName = sim["stepName"]
    print(f"-->\tRunning Energy Minimisation Step:\t{stepName}")
    ## create simluation directory
    simDir: str = p.join(outDir, sim["stepName"])
    os.makedirs(simDir, exist_ok=True)


    # Initialize the system
    system: openmm.System = init_system(prmtop)


    if isinstance(sim["temp"], int):
        sim["temp"] = sim["temp"] * unit.kelvin

    # Set up integrator and simulation
    integrator: openmm.Integrator = openmm.LangevinMiddleIntegrator(sim["temp"],
                                                                  1/unit.picosecond,
                                                                  0.004*unit.picoseconds)
    simulation: app.simulation.Simulation = app.simulation.Simulation(prmtop.topology,
                                                                    system,
                                                                    integrator,
                                                                    platform)
    simulation.context.setPositions(inpcrd.positions)

    # Set box vectors if provided
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    # Run energy minimisation
    simulation.minimizeEnergy(maxIterations=sim['maxIterations'])

    # find name to call outFiles
    protName = p.basename(p.dirname(simDir))
    # save result as pdb - reset chain and residue Ids
    state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
    minimisedPdb: str = p.join(simDir, f"{protName}.pdb")
    with open(minimisedPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology,
                                     state.getPositions(),
                                     output)
    drFixer.reset_chains_residues(refPdb, minimisedPdb)

    # Save simulation as XML
    saveXml: str = p.join(simDir, "energy_minimisation.xml")
    simulation.saveState(saveXml)
    return saveXml


