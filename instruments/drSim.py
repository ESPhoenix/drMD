## BASIC LIBS
import os
from os import path as p
from shutil import rmtree
## OPEN MM LIBS
import simtk.openmm.app as app
import simtk.openmm as openmm
import  simtk.unit  as unit
## CUSTOM LIBS
import instruments.drConstraints as drConstraints
import instruments.drCheckup as drCheckup
import instruments.drMeta as drMeta
from instruments import drFixer
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
            saveFile = run_energy_minimisation(prmtop, inpcrd, sim, simDir, platform, pdbFile)
        elif sim["type"].upper() == "NVT":
            sim = process_sim_data(sim,timescale)
            saveFile = run_nvt(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile)
        elif sim["type"].upper() == "NPT":
            sim = process_sim_data(sim,timescale)
            saveFile = run_npt(prmtop, inpcrd, sim, saveFile, simDir, platform, pdbFile)
        elif sim["type"].upper() == "META":
            sim = process_sim_data(sim, timescale)
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

###########################################################################################
def init_system(prmtop):
    system = prmtop.createSystem(nonbondedMethod=app.PME,
                                    nonbondedCutoff=1*unit.nanometer,
                                    constraints = app.HBonds)
    return system
###########################################################################################
def init_reporters(simDir, nSteps, reportInterval):
    ## vitalsCsv gives you information on potential, kinetic and total energy 
    ## as well as temperature, volume and density
    vitalsCsv = p.join(simDir, "vitals_report.csv")
    vitalsStateReporter = app.StateDataReporter(file = vitalsCsv, reportInterval = reportInterval, step = True,
                                            time = True, potentialEnergy = True, kineticEnergy = True,
                                            totalEnergy = True, temperature = True, volume = True,
                                            density = True)
    ## progresCsv gives you information on how long the simulation has been running and
    ## how long is left
    progressCsv = p.join(simDir, "progress_report.csv")
    progressStateReporter = app.StateDataReporter(file = progressCsv, progress = True, remainingTime = True,
                                            speed = True, elapsedTime = True, totalSteps = nSteps,
                                            reportInterval = reportInterval)
    ## dcdFile is the trajectory of the simulation
    dcdFile = p.join(simDir, "trajectory.dcd")
    dcdTrajectoryReporter = app.DCDReporter(file = dcdFile, reportInterval = reportInterval, append = False)
    ## chkFile works as a checkpoint so that the simulation can be resumed if something goes wrong
    chkFile = p.join(simDir, "checkpoint.chk")
    chkReporter = app.CheckpointReporter(file = chkFile, reportInterval = reportInterval, writeState = False)
    ## restraintsCsv contains all custom forces


    reporters = {"vitals":[vitalsStateReporter,vitalsCsv],
                "progress": [progressStateReporter,progressCsv],
                 "trajectory": [dcdTrajectoryReporter, dcdFile],
                 "checkpoint": [chkReporter, chkFile]}



    return  reporters
###########################################################################################
def process_sim_data(sim,timescale):
    # Reads simulation variables from config file and processes them
    timestepData = sim["timestep"].split()
    timestep = float(timestepData[0]) * timescale[timestepData[1]]
    durationData = sim["duration"].split()
    duration = int(durationData[0]) * timescale[durationData[1]]
    logIntervalData = sim["logInterval"].split()
    logInterval = float(logIntervalData[0]) * timescale[logIntervalData[1]]
    logIntervalInSteps = int(round(logInterval / timestep))
    print(logIntervalInSteps)
    nSteps = int(duration / timestep)
    nLogSteps = round(nSteps / 500)
    temp = sim["temp"] *unit.kelvin

    sim.update({"timeStep":timestep,
                "duration": duration,
                "nSteps":nSteps,
                "nLogSteps":nLogSteps,
                "temp":temp,
                "logInterval": logIntervalInSteps})
    return sim
###########################################################################################
def load_simulation_state(simulation, saveFile):
    saveFileExt = p.splitext(saveFile)[1]
    if saveFileExt == ".chk":
        simulation.loadCheckpoint(saveFile)
    elif saveFileExt == ".xml":
        simulation.loadState(saveFile)
    return simulation
###########################################################################################
def run_npt(prmtop, inpcrd, sim, saveFile, simDir, platform, refPdb):
    print("Running NpT Step!")
    ## initialise a new system from parameters
    system = init_system(prmtop)
    ## deal with any constraints
    system = drConstraints.constraints_handler(system, prmtop, inpcrd, sim, saveFile, refPdb)
    # add constant pressure force to system (makes this an NpT simulation)
    system.addForce(openmm.MonteCarloBarostat(1*unit.bar, sim["temp"]))
    integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # set up intergrator and system
    # load state from previous simulation (or continue from checkpoint)
    simulation = load_simulation_state(simulation, saveFile)
    # set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    reportInterval = sim["logInterval"]
    reporters = init_reporters(simDir = simDir,
                                nSteps =  totalSteps,
                                reportInterval= reportInterval)
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])
    # run NVT simulation
    simulation.step(sim["nSteps"])

    # save result as pdb - reset chain and residue Ids
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    nptPdb = p.join(simDir,"NpT_final_geom.pdb")
    with open(nptPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    drFixer.reset_chains_residues(refPdb,nptPdb)

    # run drCheckup
    drCheckup.check_vitals(simDir,reporters["vitals"][1], reporters["progress"][1])


    # save simulation as XML
    saveXml = p.join(simDir,"NpT_step.xml")
    simulation.saveState(saveXml)
    return saveXml
###########################################################################################
def run_nvt(prmtop, inpcrd, sim, saveFile, simDir, platform, refPdb):
    print("Running NVT Step!")
    ## initialise a new system from parameters
    system = init_system(prmtop)

    ## deal with any constraints
    system = drConstraints.constraints_handler(system, prmtop, inpcrd, sim, saveFile, refPdb)

    ## check forces
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        print("Force", i, "type:", type(force).__name__)
        # If you want to explore properties of specific forces, you can add more detailed checks
        if isinstance(force, openmm.CustomExternalForce):
            print("  - Energy expression:", force.getEnergyFunction())
            print("  - Number of global parameters:", force.getNumGlobalParameters())
            for j in range(force.getNumGlobalParameters()):
                parameterName = force.getGlobalParameterName(j)
                print("    - Global parameter:", parameterName)


    # set up intergrator and system
    integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # load state from previous simulation (or continue from checkpoint)
    simulation = load_simulation_state(simulation, saveFile)
    # set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    reportInterval = sim["logInterval"]
    reporters = init_reporters(simDir = simDir,
                                nSteps =  totalSteps,
                                reportInterval= reportInterval)
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])

    # run NVT simulation
    simulation.step(sim["nSteps"])

    # save result as pdb - reset chain and residue Ids
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    nvtPdb = p.join(simDir,"NVT_final_geom.pdb")
    with open(nvtPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    drFixer.reset_chains_residues(refPdb, nvtPdb)

    ## run drCheckup    
    drCheckup.check_vitals(simDir,reporters["vitals"][1], reporters["progress"][1])

    # save simulation as XML
    saveXml = p.join(simDir,"NVT_step.xml")
    simulation.saveState(saveXml)

    return saveXml

###########################################################################################
def run_energy_minimisation(prmtop, inpcrd, sim, simDir,platform, refPdb):
    print("Running Energy Minimisation!")
    system = init_system(prmtop)
    # set up intergrator and simulation
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    # box vectors
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    # run energy minimisation
    simulation.minimizeEnergy(maxIterations=sim['maxIterations']) 

    # save result as pdb - reset chain and residue Ids
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    minimisedPdb = p.join(simDir,"minimised_geom.pdb")
    with open(minimisedPdb, 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    drFixer.reset_chains_residues(refPdb, minimisedPdb)
    
    # save simulation as XML
    saveXml = p.join(simDir,"energy_minimisation.xml")
    simulation.saveState(saveXml)
    return saveXml

