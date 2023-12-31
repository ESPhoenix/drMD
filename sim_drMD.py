## BASIC LIBS
import os
from os import path as p
## OPEN MM LIBS
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
## CUSTOM LIBS
from constraints_drMD import heavy_atom_position_restraints
###########################################################################################
def run_simulation(config, outDir, inputCoords, amberParams):
    # set up simple unit translator
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
            sim = process_sim_data(sim,timescale)
            saveXml = run_nvt(system, prmtop, inpcrd, sim, saveXml, simDir)

        elif sim["type"] in ["NpT","NPT"]:
            sim = process_sim_data(sim,timescale)
            saveXml = run_npt(system, prmtop, inpcrd, sim, saveXml, simDir)
###########################################################################################
def process_sim_data(sim,timescale):
    # Reads simulation variables from config file and processes them
    timestepData = sim["timestep"].split()
    timestep = int(timestepData[0]) * timescale[timestepData[1]]
    durationData = sim["duration"].split()
    duration = int(durationData[0]) * timescale[durationData[1]]
    nSteps = int(duration / timestep)
    nLogSteps = round(nSteps / 500)
    temp = sim["temp"] * kelvin


    sim.update({"timeStep":timestep,
                "duration": duration,
                "nSteps":nSteps,
                "nLogSteps":nLogSteps,
                "temp":temp})
    return sim
###########################################################################################
def run_npt(system, prmtop, inpcrd, sim, saveXml, simDir):
    # add constant pressure force to system (makes this an NpT simulation)
    system.addForce(MonteCarloBarostat(1*bar, sim["temp"]))
    if sim["freezeHeavy"]:
        print("Adding position restraints to heavy atoms!")
        system = heavy_atom_position_restraints(system,prmtop,inpcrd)
    # set up intergrator and system
    integrator = LangevinMiddleIntegrator(sim["temp"], 1/picosecond, sim["timeStep"])
    simulation = Simulation(prmtop.topology, system, integrator)
    # load state from previous simulation
    simulation.loadState(saveXml)
    # set up reporters
    simulation.reporters.append(DCDReporter(p.join(simDir,'NpT_step.dcd'), 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                                potentialEnergy=True, temperature=True))
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    # run NVT simulation
    simulation.step(sim["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"NpT_step.xml")
    simulation.saveState(saveXml)
    return saveXml
###########################################################################################
def run_nvt(system, prmtop, inpcrd, sim, saveXml, simDir):
    ## Add position restraints to heavy atoms if instructed
    if sim["freezeHeavy"]:
        print("Adding position restraints to heavy atoms!")
        system = heavy_atom_position_restraints(system,prmtop,inpcrd)
      

    # set up intergrator and system
    integrator = LangevinMiddleIntegrator(sim["temp"], 1/picosecond, sim["timeStep"])
    simulation = Simulation(prmtop.topology, system, integrator)
    # load state from previous simulation
    simulation.loadState(saveXml)
    # set up reporters
    simulation.reporters.append(DCDReporter(p.join(simDir,'NVT_step.dcd'), 1000))
    simulation.reporters.append(StateDataReporter(p.join(simDir,'NVT_step.csv'), 1000, time=True,
        kineticEnergy=True, potentialEnergy=True))
    # box vectors
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    # run NVT simulation
    simulation.step(sim["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"NVT_step.xml")
    simulation.saveState(saveXml)
    return saveXml

###########################################################################################
def run_energy_minimisation(prmtop, inpcrd, system, simDir):
    print("Running Energy Minimisation!")

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

