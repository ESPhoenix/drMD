## BASIC LIBS
import os
from os import path as p
from shutil import rmtree
## OPEN MM LIBS
import simtk.openmm.app as app
import simtk.openmm as openmm
import simtk.unit  as unit
## drMD UTILS
from instruments import drPrep
from instruments import drSim
from instruments import drMeta
from instruments import drConfigInspector
## BASIC PDB <-> DF UTILS
from pdbUtils import pdbUtils
#####################################################################################
def drMD_protocol(configYaml: str) -> None:
    """
    Run the drMD protocol for a given configuration file.

    Args:
        configYaml (str): The path to the YAML configuration file.

    Returns:
        None

    This function reads the configuration file, prepares the protocol, and runs the simulation.
    """
    # Read the configuration file
    config: dict = drConfigInspector.read_config(configYaml)

    # Create the output directory if it doesn't exist
    outDir: str = config["pathInfo"]["outputDir"]
    os.makedirs(outDir, exist_ok=True)

    # Prepare the protocol
    mergedPdb, inputCoords, amberParams = drPrep.prep_protocol(config)

    # Run the simulation
    run_simulation(config, outDir, inputCoords, amberParams, mergedPdb)
###########################################################################################
def run_simulation(config: dict, outDir: str, inputCoords: str, amberParams: str, pdbFile: str) -> None:
    """
    Run the simulation according to the given configuration.

    Args:
        config (dict): The configuration dictionary.
        outDir (str): The output directory.
        inputCoords (str): The path to the input coordinates file.
        amberParams (str): The path to the Amber parameters file.
        pdbFile (str): The path to the PDB file.

    Returns:
        None
    """
    # Set up unit translator
    timescale: dict = {"fs":unit.femtoseconds,
                 "ps":unit.picoseconds,
                 "ns":unit.nanoseconds}

    # Set up platform
    usePlatform: str = config["generalInfo"]["platform"]
    if usePlatform == "CUDA":
        platform=openmm.Platform.getPlatformByName("CUDA")
    elif usePlatform == "OpenCL":
        platform=openmm.Platform.getPlatformByName("OpenCL")
    elif usePlatform == "CPU":
        platform=openmm.Platform.getPlatformByName("CPU")

    # Load Amber files and create system
    prmtop: app.Topology = app.AmberPrmtopFile(amberParams)
    inpcrd: app.InpcrdFile = app.AmberInpcrdFile(inputCoords)

    # Loop over simulations
    simulations = config["simulationInfo"]
    for i in range(len(simulations)):
        sim: dict = simulations[i]
        simDir: str = p.join(outDir,sim["stepName"])

        # Decide whether to skip, resume, or start a new simulation
        skipResumeSim, saveFile = skip_resume_or_simulate(simDir=simDir,
                                                           simulations = simulations,
                                                           i = i, 
                                                           outDir=outDir)

        # Skip or resume simulation
        if skipResumeSim == "skip":
            stepName: str = sim["stepName"]
            print(f"Skipping {stepName}")
            continue
        if skipResumeSim == "resume":
            print(f"Resuming {stepName} from checkpoint file")

        # Set up directories
        os.makedirs(simDir,exist_ok=True)
        os.chdir(simDir)

        # Run the simulation
        if sim["simulationType"].upper() == "EM":
            saveFile = drSim.run_energy_minimisation(prmtop, inpcrd, sim, simDir, platform, pdbFile)
        elif sim["simulationType"].upper() == "NVT":
            sim = drSim.process_sim_data(sim, timescale)
            saveFile = drSim.run_nvt(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile)
        elif sim["simulationType"].upper() == "NPT":
            sim = drSim.process_sim_data(sim, timescale)
            saveFile = drSim.run_npt(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile)
        elif sim["simulationType"].upper() == "META":
            sim = drSim.process_sim_data(sim, timescale)
            saveFile = drMeta.run_metadynamics(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile)
###########################################################################################
def skip_resume_or_simulate(simDir: str, simulations: list, i: int, outDir: str) -> tuple:
    """
    Check if the simulation directory exists and decide whether to skip, resume or start a new simulation.

    Args:
        simDir (str): The path to the simulation directory.
        simulations (list): List of simulation information dictionaries.
        i (int): Index of the current simulation in the list.
        outDir (str): The path to the output directory.

    Returns:
        tuple: A tuple containing the action to be taken (simulate, skip or resume) and the path to the save file.
    """
    # run simulation if simDir doesn't exist
    if not p.isdir(simDir):
        # if it's the first simulation, return "simulate"
        if i == 0:
            return "simulate", "foo"
        ## find previous saveXml file
        previousSim: dict = simulations[i-1]
        previousSimDir: str = p.join(outDir, previousSim["stepName"])
        # look for previous saveXml file
        for file in os.listdir(previousSimDir):
            if p.splitext(file)[1] == ".xml":
                saveXml: str = p.join(previousSimDir,file)
                return "simulate", saveXml
    ## look for save.xml file that is written once a sim finishes
    ## skip over sim if this exists
    for file in os.listdir(simDir):
        if p.splitext(file)[1] == ".xml":
            saveXml: str = p.join(simDir,file)
            return "skip", saveXml
    ## look for checkpoint.chk file that is written during simulation
    ## resume simulation from this checkpoint file
    for file in os.listdir(simDir):
        if p.splitext(file)[1] == ".chk":
            saveChk: str = p.join(simDir,file)
            return "resume", saveChk
    ## if it's got here, we have an empty simDir
    ## remove and simulate
    rmtree(simDir)
    previousSim: dict = simulations[i-1]
    previousSimDir: str = p.join(outDir, previousSim["stepName"])
    # look for previous saveXml file
    for file in os.listdir(previousSimDir):
        if p.splitext(file)[1] == ".xml":
            saveXml: str = p.join(previousSimDir,file)
            return "simulate", saveXml
