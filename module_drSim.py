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
    platform=openmm.Platform.getPlatformByName("OpenCL")
    # load amber files, create system
    prmtop = app.AmberPrmtopFile(amberParams)
    inpcrd = app.AmberInpcrdFile(inputCoords)  
    simulations = config["simulationInfo"]
    for sim in simulations:
        # sort out directories
        simDir = p.join(outDir,sim["stepName"])
        ###############################
        ## skip if complete:
        skipToNextSim = False
        if p.isdir(simDir):
            for file in os.listdir(simDir):
                if p.splitext(file)[1] == ".xml":
                    saveXml = p.join(simDir,file)
                    skipToNextSim = True
        if skipToNextSim:
            stepName = sim["stepName"]
            print(f"skipping {stepName}")
            continue
        ################################
        # sort out directories
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

    reporters = {"vitals":[vitalsStateReporter,vitalsCsv],
                "progress": [progressStateReporter,progressCsv],
                 "trajectory": [dcdTrajectoryReporter, dcdFile]}

    return  reporters
###########################################################################################
def process_sim_data(sim,timescale):
    # Reads simulation variables from config file and processes them
    timestepData = sim["timestep"].split()
    timestep = float(timestepData[0]) * timescale[timestepData[1]]
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
    ## deal with any constraints
    system, clearRestraints = drConstraints.constraints_handler(system, prmtop, inpcrd, sim, saveXml)
    # add constant pressure force to system (makes this an NpT simulation)
    system.addForce(openmm.MonteCarloBarostat(1*unit.bar, sim["temp"]))
    integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)

    # if clearRestraints:
    #     print("WONK")
    #     simulation.context.setParameter("k", (0 *unit.kilocalories_per_mole/unit.angstroms**2))

    # set up intergrator and system
    # load state from previous simulation
    simulation.loadState(saveXml)
    # set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    reporters = init_reporters(simDir = simDir,
                                nSteps =  totalSteps,
                                nLogSteps = sim["nLogSteps"])
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])
    # set restraints constant to 0 


    # # set up periodic boundary conditions
    # if inpcrd.boxVectors is not None:
    #     simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    # run NVT simulation
    simulation.step(sim["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"NpT_step.xml")
    simulation.saveState(saveXml)

    # save result as pdb
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir,"NpT_final_geom.pdb"), 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)

    drCheckup.check_vitals(simDir,reporters["vitals"][1], reporters["progress"][1])
    return saveXml
###########################################################################################
def run_nvt(prmtop, inpcrd, sim, saveXml, simDir, platform):
    print("Running NVT Step!")
    ## initialise a new system from parameters
    system = init_system(prmtop)

    ## deal with any constraints
    system, clearRestraints = drConstraints.constraints_handler(system, prmtop, inpcrd, sim, saveXml)
      
    # set up intergrator and system
    integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # load state from previous simulation
    simulation.loadState(saveXml)
    # set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    reporters = init_reporters(simDir = simDir,
                                nSteps =  totalSteps,
                                nLogSteps = sim["nLogSteps"])
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])

    # run NVT simulation
    simulation.step(sim["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"NVT_step.xml")
    simulation.saveState(saveXml)

    # save result as pdb
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir,"NVT_final_geom.pdb"), 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    drCheckup.check_vitals(simDir,reporters["vitals"][1], reporters["progress"][1])
    return saveXml

###########################################################################################
def run_energy_minimisation(prmtop, inpcrd, sim, simDir,platform):
    print("Running Energy Minimisation!")
    system = init_system(prmtop)
    # set up intergrator and simulation
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)

    # set up reporters
    totalSteps = simulation.currentStep + sim["maxIterations"]
    reporters = init_reporters(simDir = simDir,
                                nSteps =  totalSteps,
                                nLogSteps = 100)
    
    simulation.reporters.append(reporters["vitals"][0])
    # box vectors
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
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

