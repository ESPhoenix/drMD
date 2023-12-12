## BASIC LIBS
import os
from os import path as p
## OPEN MM LIBS
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout



###########################################################################################
def run_simulation(config, outDir, inputCoords, amberParams):

    timescale = {"fs":femtoseconds,
                 "ps":picoseconds,
                 "ns":nanoseconds}
    ## set up gpu
    platform=Platform.getPlatformByName("CUDA")

    # load amber files, create system
    prmtop = AmberPrmtopFile(amberParams)
    inpcrd = AmberInpcrdFile(inputCoords)  
    system = prmtop.createSystem(nonbondedMethod=PME,
                                    nonbondedCutoff=1*nanometer,
                                    constraints = HBonds)
    # set up simple boolean to check for whether we need to continue
    isInitialSim = True
    simulations = config["simulationInfo"]
    for sim in simulations:
        # sort out directories
        simDir = p.join(outDir,sim["stepName"])
        os.makedirs(simDir,exist_ok=True)
        os.chdir(simDir)
        # Check for simulation type, run as needed:
        if sim["type"] == "EM":
            saveXml = run_energy_minimisation(prmtop, inpcrd, system, simDir)

        elif sim["type"] == "NVT":
            simData = read_sim_data(sim,timescale)
            saveXml = run_nvt(system, prmtop, simData, saveXml, simDir)

        elif sim["type"] in ["NpT","NPT"]:
            simData = read_sim_data(sim,timescale)
            saveXml = run_npt(system, prmtop, simData, saveXml, simDir)
###########################################################################################
def read_sim_data(sim,timescale):
    # Reads simulation variables from config file and processes them
    timestepData = sim["timestep"].split()
    timestep = int(timestepData[0]) * timescale[timestepData[1]]
    durationData = sim["duration"].split()
    duration = int(durationData[0]) * timescale[durationData[1]]
    nSteps = int(duration / timestep)
    nLogSteps = round(nSteps / 500)
    temp = sim["temp"] * kelvin


    simData = {"timeStep":timestep,
                "duration": duration,
                "nSteps":nSteps,
                "nLogSteps":nLogSteps,
                "temp":temp}
    return simData
###########################################################################################
def run_npt(system, prmtop, simData, saveXml, simDir):
    # add constant pressure force to system (makes this an NpT simulation)
    system.addForce(MonteCarloBarostat(1*bar, simData["temp"]))
    # set up intergrator and system
    integrator = LangevinMiddleIntegrator(simData["temp"], 1/picosecond, simData["timeStep"])
    simulation = Simulation(prmtop.topology, system, integrator)
    # load state from previous simulation
    simulation.loadState(saveXml)
    # set up reporters
    simulation.reporters.append(PDBReporter(p.join(simDir,'NpT_step.pdb'), 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                                potentialEnergy=True, temperature=True))
    # run NVT simulation
    simulation.step(simData["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"NpT_step.xml")
    simulation.saveState(saveXml)
    return saveXml
###########################################################################################
def run_nvt(system, prmtop, simData, saveXml, simDir):
    # set up intergrator and system
    integrator = LangevinMiddleIntegrator(simData["temp"], 1/picosecond, simData["timeStep"])
    simulation = Simulation(prmtop.topology, system, integrator)
    # load state from previous simulation
    simulation.loadState(saveXml)
    # set up reporters
    simulation.reporters.append(DCDReporter(p.join(simDir,'NVT_step.dcd'), 1000))
    simulation.reporters.append(StateDataReporter(p.join(simDir,'NVT_step.csv'), 1000, time=True,
        kineticEnergy=True, potentialEnergy=True))
    # run NVT simulation
    simulation.step(simData["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"NVT_step.xml")
    simulation.saveState(saveXml)
    return saveXml

###########################################################################################
def run_energy_minimisation(prmtop, inpcrd, system, simDir):
    # set up intergrator and simulation
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    # run energy minimisation
    simulation.minimizeEnergy()
    # save result as pdb
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir,"minimised_geom.pdb"), 'w') as output:
        PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    # save simulation as XML
    saveXml = p.join(simDir,"energy_minimisation.xml")
    simulation.saveState(saveXml)
    return saveXml

