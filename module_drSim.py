## BASIC LIBS
import os
from os import path as p
## OPEN MM LIBS
import simtk.openmm.app as app
import simtk.openmm as openmm
import  simtk.unit  as unit
from sys import stdout
## CUSTOM LIBS
import module_drConstraints as drConstraints
import module_drCheckup as drCheckup
###########################################################################################
def run_simulation(config, outDir, inputCoords, amberParams):
    # set up simple unit translator
    timescale = {"fs":unit.femtoseconds,
                 "ps":unit.picoseconds,
                 "ns":unit.nanoseconds}
    ## set up gpu
    platform=openmm.Platform.getPlatformByName("CUDA")
    # load amber files, create system

    prmtop = app.AmberPrmtopFile(amberParams)
    inpcrd = app.AmberInpcrdFile(inputCoords)  
    simulations = config["simulationInfo"]
    for sim in simulations:
        # sort out directories
        simDir = p.join(outDir,sim["stepName"])
        os.makedirs(simDir,exist_ok=True)
        os.chdir(simDir)
        # Check for simulation type, run as needed:
        if sim["type"] == "EM":
            saveXml = run_energy_minimisation(prmtop, inpcrd, sim, simDir, platform)
        elif sim["type"] == "NVT":
            sim = process_sim_data(sim,timescale)
            saveXml = run_nvt(prmtop, inpcrd, sim, saveXml, simDir, platform)

        elif sim["type"] in ["NpT","NPT"]:
            sim = process_sim_data(sim,timescale)
            saveXml = run_npt(prmtop, inpcrd, sim, saveXml, simDir, platform)
###########################################################################################
def init_system(prmtop):
    system = prmtop.createSystem(nonbondedMethod=app.PME,
                                    nonbondedCutoff=1*unit.nanometer,
                                    constraints = app.HBonds)
    return system
###########################################################################################
def init_reporters(simDir, nSteps, nLogSteps):
    vitalsCsv = p.join(simDir, "vitals_report.csv")
    vitalsStateReporter = app.StateDataReporter(file = vitalsCsv, reportInterval = nLogSteps, step = True,
                                            time = True, potentialEnergy = True, kineticEnergy = True,
                                            totalEnergy = True, temperature = True, volume = True,
                                            density = True)
    progressCsv = p.join(simDir, "progress_report.csv")
    progressStateReporter = app.StateDataReporter(file = progressCsv, progress = True, remainingTime = True,
                                            speed = True, elapsedTime = True, totalSteps = nSteps,
                                            reportInterval = nLogSteps)

    dcdFile = p.join(simDir, "trajectory.dcd")
    dcdTrajectoryReporter = app.DCDReporter(file = dcdFile, reportInterval = nLogSteps, append = False)


    return vitalsStateReporter, progressStateReporter, dcdTrajectoryReporter, vitalsCsv, progressCsv
###########################################################################################
def process_sim_data(sim,timescale):
    # Reads simulation variables from config file and processes them
    timestepData = sim["timestep"].split()
    timestep = int(timestepData[0]) * timescale[timestepData[1]]
    durationData = sim["duration"].split()
    duration = int(durationData[0]) * timescale[durationData[1]]
    nSteps = int(duration / timestep)
    nLogSteps = round(nSteps / 500)
    temp = sim["temp"] *unit.kelvin

    sim.update({"timeStep":timestep,
                "duration": duration,
                "nSteps":nSteps,
                "nLogSteps":nLogSteps,
                "temp":temp})
    return sim
###########################################################################################
def run_npt(prmtop, inpcrd, sim, saveXml, simDir, platform):
    print("Running NpT Step!")
    ## initialise a new system from parameters
    system = init_system(prmtop)
    # add constant pressure force to system (makes this an NpT simulation)
    system.addForce(openmm.MonteCarloBarostat(1*unit.bar, sim["temp"]))

    ## deal with any constraints
    system = drConstraints.constraints_handler(system, prmtop, inpcrd, sim, saveXml)

    # set up intergrator and system
    integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # load state from previous simulation
    simulation.loadState(saveXml)
    # set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    vitalsStateReporter, progressStateReporter, dcdTrajectoryReporter, vitalsCsv, progressCsv = init_reporters(simDir = simDir,
                                                                                        nSteps =  totalSteps,
                                                                                        nLogSteps = sim["nLogSteps"])

    for rep in [vitalsStateReporter, progressStateReporter, dcdTrajectoryReporter]:
        simulation.reporters.append(rep)

    # set up periodic boundary conditions
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    # run NVT simulation
    simulation.step(sim["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"NpT_step.xml")
    simulation.saveState(saveXml)

    # save result as pdb
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir,"NpT_final_geom.pdb"), 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)

    drCheckup.check_vitals(simDir,vitalsCsv, progressCsv)
    return saveXml
###########################################################################################
def run_nvt(prmtop, inpcrd, sim, saveXml, simDir, platform):
    print("Running NVT Step!")
    ## initialise a new system from parameters
    system = init_system(prmtop)

    ## deal with any constraints
    system = drConstraints.constraints_handler(system, prmtop, inpcrd, sim, saveXml)
      
    # set up intergrator and system
    integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # load state from previous simulation
    simulation.loadState(saveXml)
    # set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    vitalsStateReporter, progressStateReporter, dcdTrajectoryReporter, vitalsCsv, progressCsv = init_reporters(simDir = simDir,
                                                                                        nSteps =  totalSteps,
                                                                                        nLogSteps = sim["nLogSteps"])

    for rep in [vitalsStateReporter, progressStateReporter, dcdTrajectoryReporter]:
        simulation.reporters.append(rep)

    # box vectors
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    # run NVT simulation
    simulation.step(sim["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"NVT_step.xml")
    simulation.saveState(saveXml)

    # save result as pdb
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir,"NVT_final_geom.pdb"), 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    drCheckup.check_vitals(simDir,vitalsCsv, progressCsv)
    return saveXml

###########################################################################################
def run_energy_minimisation(prmtop, inpcrd, sim, simDir,platform):
    print("Running Energy Minimisation!")
    system = init_system(prmtop)
    # set up intergrator and simulation
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    # run energy minimisation
    simulation.minimizeEnergy(maxIterations=sim['maxIterations']) 
    # save result as pdb
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir,"minimised_geom.pdb"), 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    # save simulation as XML
    saveXml = p.join(simDir,"energy_minimisation.xml")
    simulation.saveState(saveXml)
    return saveXml

