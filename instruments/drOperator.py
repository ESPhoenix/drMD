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
def drMD_protocol(configYaml):
    config = drConfigInspector.read_config(configYaml)
    outDir = config["pathInfo"]["outputDir"]
    os.makedirs(outDir,exist_ok=True)
        
    mergedPdb, inputCoords, amberParams = drPrep.prep_protocol(config)

    run_simulation(config, outDir, inputCoords, amberParams, mergedPdb)
###########################################################################################
def run_simulation(config, outDir, inputCoords, amberParams, pdbFile):
    # set up simple unit translator
    timescale = {"fs":unit.femtoseconds,
                 "ps":unit.picoseconds,
                 "ns":unit.nanoseconds}
    ## set up paltform

    usePlatform = config["generalInfo"]["platform"]

    if usePlatform == "CUDA":
        platform=openmm.Platform.getPlatformByName("CUDA")
    elif usePlatform == "OpenCL":
        platform=openmm.Platform.getPlatformByName("OpenCL")
    elif usePlatform == "CPU":
        platform=openmm.Platform.getPlatformByName("CPU")

    omp_num_threads = os.environ.get('OMP_NUM_THREADS')
    
    # load amber files, create system
    prmtop = app.AmberPrmtopFile(amberParams)
    inpcrd = app.AmberInpcrdFile(inputCoords)  
    simulations = config["simulationInfo"]
    for i in range(len(simulations)):
        sim = simulations[i]
        # sort out directories
        simDir = p.join(outDir,sim["stepName"])
        ###############################
        ## this bit deals with whether we continue from a saveXml or saveChk file
        skipResumeSim, saveFile = skip_resume_or_simulate(simDir=simDir,
                                                           simulations = simulations,
                                                             i = i, 
                                                             outDir=outDir)

        if skipResumeSim == "skip":
            stepName = sim["stepName"]
            print(f"skipping {stepName}")
            continue
        if skipResumeSim == "resume":
            print(f"resuming {stepName} from checkpoint file")

        ################################
        # sort out directories
        os.makedirs(simDir,exist_ok=True)
        os.chdir(simDir)
        # Check for simulation type, run as needed:
        if sim["type"].upper() == "EM":
            saveFile = drSim.run_energy_minimisation(prmtop, inpcrd, sim, simDir, platform, pdbFile)
        elif sim["type"].upper() == "NVT":
            sim = drSim.process_sim_data(sim,timescale)
            saveFile = drSim.run_nvt(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile)
        elif sim["type"].upper() == "NPT":
            sim = drSim.process_sim_data(sim,timescale)
            saveFile = drSim.run_npt(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile)
        elif sim["type"].upper() == "META":
            sim = drSim.process_sim_data(sim, timescale)
            saveFile = drMeta.run_metadynamics(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile)
###########################################################################################
def skip_resume_or_simulate(simDir, simulations, i, outDir):
    # run simulation if simDir doesn't exist
    if not p.isdir(simDir):
        if i == 0:
            return "simulate", "foo"
        ## find previous saveXml file
        previousSim = simulations[i-1]
        previousSimDir = p.join(outDir, previousSim["stepName"])
        for file in os.listdir(previousSimDir):
            if p.splitext(file)[1] == ".xml":
                saveXml = p.join(previousSimDir,file)
                return "simulate", saveXml
    ## look for save.xml file that is written once a sim finishes
    ## skip over sim if this exists
    for file in os.listdir(simDir):
        if p.splitext(file)[1] == ".xml":
            saveXml = p.join(simDir,file)
            return "skip", saveXml
    ## look for checkpoint.chk file that is written during simulation
    ## resume simulation from this checkpoint file
    for file in os.listdir(simDir):
        if p.splitext(file)[1] == ".chk":
            saveChk = p.join(simDir,file)
            return "resume", saveChk
    ## if it's got here, we have an empty simDir
    ## remove and simulate
    rmtree(simDir)
    previousSim = simulations[i-1]
    previousSimDir = p.join(outDir, previousSim["stepName"])
    for file in os.listdir(previousSimDir):
        if p.splitext(file)[1] == ".xml":
            saveXml = p.join(previousSimDir,file)
            return "simulate", saveXml