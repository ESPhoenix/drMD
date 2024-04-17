## BASIC LIBS
import os
from os import path as p
## OPEN MM LIBS
import simtk.openmm.app as app
from simtk.openmm.app import metadynamics
import simtk.openmm as openmm
import  simtk.unit  as unit
## CUSTOM LIBS
import instruments.drSim as drSim
import instruments.drConstraints as drConstraints
import instruments.drCheckup as drCheckup
import instruments.pdbUtils as pdbUtils
########################################################################################################
def run_metadynamics(prmtop, inpcrd, sim, saveXml, simDir, platform, pdbFile):
    print("Running MetaDynamics!")
    ## initilaise new system
    system = drSim.init_system(prmtop)
    ## deal with restraints (clear all lurking restraints and constants)
    system, clearRestraints = drConstraints.constraints_handler(system, prmtop, inpcrd, sim, saveXml)
    ## create bias variable
    biasVariables = []
    if sim["biasVar"].upper() == "RMSF":
        biasVariable = gen_rmsd_bias_variable(prmtop, inpcrd, sim)
        biasVariables.append(biasVariable)
    elif sim["biasVar"].upper() == "DIHEDRAL":
        biasVariable = gen_dihedral_bias_variable(prmtop, inpcrd, sim, pdbFile)
    ## create metadynamics object -- adds bias variable as a force to the system
    meta = metadynamics.Metadynamics(system = system,
                                     variables = biasVariables,
                                     temperature = sim["temp"],
                                     biasFactor = 5,
                                     height = 1,
                                     width = 0.05, 
                                     frequency = 50,
                                     saveFrequency = 50,
                                     biasDir = simDir)
    ## set up intergrator
    integrator = openmm.LangevinMiddleIntegrator(sim["temp"], 1/unit.picosecond, sim["timeStep"])
    ## create new simulation
    simulation = app.simulation.Simulation(prmtop.topology, system, integrator, platform)
    # load state from previous simulation
    simulation.loadState(saveXml)    
    ## set up reporters
    totalSteps = simulation.currentStep + sim["nSteps"]
    reporters = drSim.init_reporters(simDir = simDir,
                                nSteps =  totalSteps,
                                nLogSteps = sim["nLogSteps"])
    for rep in reporters:
        simulation.reporters.append(reporters[rep][0])
    ## run metadynamics simulation
    meta.step(simulation, sim["nSteps"])
    # save simulation as XML
    saveXml = p.join(simDir,"Meta_step.xml")
    simulation.saveState(saveXml)
    # save result as pdb
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(p.join(simDir,"Meta_final_geom.pdb"), 'w') as output:
        app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    ## run drCheckup to assess simulation health (will this work for METADYNAMICS???)
    drCheckup.check_vitals(simDir,reporters["vitals"][1], reporters["progress"][1])
    ## return checkpoint file for continuing simulation


    return saveXml
########################################################################################################
def gen_dihedral_bias_variable(prmtop, inpcrd, sim, pdbFile):
    topology = prmtop.topology
    positions = inpcrd.positions
    ## get atom indecies
    dihedralAtomsInput = sim["biasAtoms"]
    atomIndecies = []

    pdbDf = pdbUtils.pdb2df(pdbFile)
    for atomInfo in dihedralAtomsInput:
        print(atomInfo)
        atomDf = pdbDf[(pdbDf["CHAIN_ID"].astype(str) == atomInfo[0]) &
                        (pdbDf["RES_NAME"].astype(str) == atomInfo[1]) &
                        (pdbDf["RES_ID"].astype(str) == atomInfo[2])]# &
        #                (pdbDf["ATOM_NAME"].astype(str) == atomInfo[3]) ]
        print(atomDf)
        continue

    #         if (atom.residue.chain.id == chain_id and 
    #             atom.residue.name == res_name and 
    #             atom.residue.id == res_id):# and 
    #             # atom.name == atom_name):
    #             atomIndecies.append(i)
    # print(atomIndecies)
    exit()
    ## generate a dihedral bias force:
    dihedralForce = openmm.CustomTorsionForce("k*(theta-theta0)^2")
    dihedralForce.addTorsion(atomIndecies[0],
                              atomIndecies[1],
                                atomIndecies[2],
                                  atomIndecies[3],
                                    [], 0, 0)

    dihedralBiasVariable = metadynamics.BiasVariable(force = dihedralForce,
                                                     minValue = -180 * unit.degrees,
                                                     maxValue = 180 * unit.degrees,
                                                     biasWidth = 1 * unit.degrees,
                                                     periodic = True)
    return dihedralBiasVariable




########################################################################################################
def gen_rmsd_bias_variable(prmtop, inpcrd, sim):
    topology = prmtop.topology
    positions = inpcrd.positions

    ## generate a backbone bias force:
    backboneAtomNames = ["N","CA","C","O"]
    atomIndecies = []
    atomCoords = []
    for i, atom in enumerate(topology.atoms()):
        if atom.name in backboneAtomNames:
            atomIndecies.append(i)
        atomCoords.append(positions[i])

    rmsdForce = openmm.RMSDForce(atomCoords, atomIndecies)
    rmsdBiasVariable = metadynamics.BiasVariable(force = rmsdForce,
                                                 minValue = 10 * unit.angstrom,
                                                 maxValue = 100 * unit.angstrom,
                                                 biasWidth = 1 * unit.angstrom,
                                                 periodic = False)
    return rmsdBiasVariable



########################################################################################################
def add_funnel_potential(system, ):
    """https://github.com/jeff231li/funnel_potential"""
