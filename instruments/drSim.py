## BASIC LIBS
import os
from os import path as p
## OPEN MM LIBS
import simtk.openmm.app as app
import simtk.openmm as openmm
import  simtk.unit  as unit
## CUSTOM LIBS
from  instruments import drRestraints
from  instruments import drCheckup 
from  instruments import drClusterizer
from instruments import drFixer

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

    # Define the constraints.
    constraints: openmm.Force = app.HBonds

    # Create the system.
    system: openmm.System = prmtop.createSystem(nonbondedMethod=nonbondedMethod,
                                                nonbondedCutoff=nonbondedCutoff,
                                                constraints=constraints)

    return system

###########################################################################################
def process_sim_data(sim: dict, timescale: dict) -> dict:
    """
    Reads simulation variables from config file and processes them.

    Parameters
    ----------
    sim : dict
        Dictionary containing simulation variables.
    timescale : dict
        Dictionary mapping time units to their corresponding values.

    Returns
    -------
    sim : dict
        Updated dictionary with processed simulation variables.
    """
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

    # Print log interval in steps
    print(logIntervalInSteps)

    # Calculate number of steps
    nSteps: int = int(duration / timestep)

    # Calculate number of log steps (every 500 steps)
    nLogSteps: int = round(nSteps / 500)

    # Process temperature
    temp: unit.Quantity = sim["temp"] * unit.kelvin

    # Update sim dictionary with processed variables
    sim.update({
        "timeStep": timestep,
        "duration": duration,
        "nSteps": nSteps,
        "nLogSteps": nLogSteps,
        "temp": temp,
        "logInterval": logIntervalInSteps
    })

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

    # Return the modified simulation object
    return simulation
###########################################################################################
def run_npt(prmtop: str, inpcrd: str, sim: dict, saveFile: str, simDir: str, platform: openmm.Platform, refPdb: str) -> str:
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
    initializes the system, handles any constraints, adds a constant pressure force,
    sets up the integrator and system, loads the state from a checkpoint or XML file,
    sets up reporters, runs the simulation, saves the final geometry as a PDB file,
    resets the chain and residue Ids, runs drCheckup, performs trajectory clustering
    if specified in the simulation parameters, and saves the simulation state as an
    XML file.
    """
    print("Running NpT Step!")
    ## initialise a new system from parameters
    system: openmm.System = init_system(prmtop)
    ## deal with any constraints
    system: openmm.System = drRestraints.constraints_handler(system, prmtop, inpcrd, sim, saveFile, refPdb)
    # add constant pressure force to system (makes this an NpT simulation)
    system.addForce(openmm.MonteCarloBarostat(1*unit.bar, sim["temp"]))
    integrator: openmm.Integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
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
    # run NVT simulation
    simulation.step(sim["nSteps"])

    # save result as pdb - reset chain and residue Ids
    state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
    nptPdb: str = p.join(simDir, "NpT_final_geom.pdb")
    with open(nptPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    drFixer.reset_chains_residues(refPdb, nptPdb)

    # run drCheckup
    drCheckup.check_vitals(simDir, reporters["vitals"][1], reporters["progress"][1])

    # run trajectory clustering
    if "clusterTrajectory" in sim and sim["clusterTrajectory"]["clusterBool"]:
        drClusterizer.rmsd_clustering_protocol(simDir, sim["clusterTrajectory"])

    # save simulation as XML
    saveXml: str = p.join(simDir, "NpT_step.xml")
    simulation.saveState(saveXml)
    return saveXml
###########################################################################################
def run_nvt(
        prmtop: str,
        inpcrd: str,
        sim: dict,
        saveFile: str,
        simDir: str,
        platform: openmm.Platform,
        refPdb: str
) -> str:
    """
    Run a simulation at constant volume (NVT) step.

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

    This function runs a simulation at constant volume (NVT) step. The function
    initializes the system, handles any constraints, checks the forces applied to the
    system, sets up the integrator and system, loads the state from a checkpoint or XML
    file, sets up reporters, runs the simulation, saves the final geometry as a PDB file,
    resets the chain and residue Ids, runs drCheckup, performs trajectory clustering if
    specified in the simulation parameters, and saves the simulation state as an XML file.
    """
    print("Running NVT Step!")

    # Initialize a new system from parameters
    system: openmm.System = init_system(prmtop)

    # Deal with any constraints
    system: openmm.System = drRestraints.constraints_handler(system, prmtop, inpcrd, sim, saveFile, refPdb)

    # Check forces
    for i in range(system.getNumForces()):
        force: openmm.Force = system.getForce(i)
        print("Force", i, "type:", type(force).__name__)
        # If you want to explore properties of specific forces, you can add more detailed checks
        if isinstance(force, openmm.CustomExternalForce):
            print("  - Energy expression:", force.getEnergyFunction())
            print("  - Number of global parameters:", force.getNumGlobalParameters())
            for j in range(force.getNumGlobalParameters()):
                parameterName = force.getGlobalParameterName(j)
                print("    - Global parameter:", parameterName)

    # Set up integrator and system
    integrator: openmm.Integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    simulation: app.simulation.Simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)

    # Load state from previous simulation (or continue from checkpoint)
    simulation: app.Simulation = load_simulation_state(simulation, saveFile)

    # Set up reporters
    totalSteps: int = simulation.currentStep + sim["nSteps"]
    reportInterval: int = sim["logInterval"]
    reporters: dict = init_reporters(simDir=simDir,
                                     nSteps=totalSteps,
                                     reportInterval=reportInterval)
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])

    # Run NVT simulation
    simulation.step(sim["nSteps"])

    # Save result as pdb - reset chain and residue Ids
    state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
    nvtPdb: str = p.join(simDir, "NVT_final_geom.pdb")
    with open(nvtPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    drFixer.reset_chains_residues(refPdb, nvtPdb)

    # Run drCheckup
    drCheckup.check_vitals(simDir, reporters["vitals"][1], reporters["progress"][1])

    # Run trajectory clustering
    if "clusterTrajectory" in sim and sim["clusterTrajectory"]["clusterBool"]:
        drClusterizer.rmsd_clustering_protocol(simDir, sim["clusterTrajectory"])

    # Save simulation as XML
    saveXml: str = p.join(simDir, "NVT_step.xml")
    simulation.saveState(saveXml)

    return saveXml

###########################################################################################
def run_energy_minimisation(prmtop: str, inpcrd: str, sim: dict, simDir: str, platform: openmm.Platform, refPdb: str) -> str:
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
    print("Running Energy Minimisation!")
    # Initialize the system
    system: openmm.System = init_system(prmtop)

    # Set up integrator and simulation
    integrator: openmm.Integrator = openmm.LangevinMiddleIntegrator(sim["temp"]*unit.kelvin,
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

    # Save result as pdb - reset chain and residue Ids
    state: openmm.State = simulation.context.getState(getPositions=True,
                                                      getEnergy=True)
    minimisedPdb: str = p.join(simDir, "minimised_geom.pdb")
    with open(minimisedPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology,
                                     state.getPositions(),
                                     output)
    drFixer.reset_chains_residues(refPdb, minimisedPdb)

    # Save simulation as XML
    saveXml: str = p.join(simDir, "energy_minimisation.xml")
    simulation.saveState(saveXml)
    return saveXml


