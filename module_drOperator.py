## BASIC LIBS
import os
from os import path as p
from shutil import copy
## INPUT LIBS
import yaml
import argpass
## drMD UTILS
import module_drPrep as drPrep
## drMD simulation
import module_drSim as drSim
from module_pdbUtils import mergePdbs
#####################################################################################
def read_inputs(configYaml):
    ## Read config.yaml into a dictionary
    with open(configYaml,"r") as yamlFile:
        config = yaml.safe_load(yamlFile) 
    return config
#####################################################################################
def check_inputs(config):
    if not config["proteinInfo"]["nProteins"] == len(config["proteinInfo"]["proteins"]):
        print("Number of proteins in config does not match nProteins")
        exit()
    if "ligandInfo" in config and not "ligandInfo" == None:
        if not config["ligandInfo"]["nLigands"] == len(config["ligandInfo"]["ligands"]):
            print("Number of ligands in config does not match nLigands")
            exit()
    if not p.isfile(config["pathInfo"]["inputPdb"]):
        print("Input PDB does not exist")
        exit()       

#####################################################################################
#####################################################################################
def drMD_protocol(configYaml):
    config = read_inputs(configYaml)
    #check_inputs(config=config)
    ## MAKE OUTPUT DIRECTORY
    outDir = config["pathInfo"]["outputDir"]
    prepDir = p.join(outDir,"00_prep")

    amberParams = False
    inputCoords = False

    ## skip prep if complete
    wholeDir = p.join(prepDir,"WHOLE")
    if p.isdir(p.join(wholeDir)):
        for file in os.listdir(wholeDir):
            if p.splitext(file)[1] == ".prmtop":
                amberParams = p.join(wholeDir,file)
            elif p.splitext(file)[1] == ".inpcrd":
                inputCoords = p.join(wholeDir,file)

    if amberParams and inputCoords:
        drSim.run_simulation(config = config,
                outDir = outDir,
                inputCoords=inputCoords,
                amberParams=amberParams)
        return

    ## MAIN PREP PROTOCOL
    os.makedirs(outDir,exist_ok=True)
    os.makedirs(prepDir,exist_ok=True)
    
    prepLog = p.join(prepDir,"prep.log")

    if "ligandInfo" in config:
        ## SPLIT INPUT PDB INTO PROT AND ONE FILE PER LIGAND
        inputPdb = config["pathInfo"]["inputPdb"]
        drPrep.split_input_pdb(inputPdb =inputPdb,
                        config = config,
                        outDir=prepDir)
        ## PREPARE LIGAND PARAMETERS, OUTPUT LIGAND PDBS
        ligandPdbs,ligandFileDict = drPrep.prepare_ligand_parameters(config = config, outDir = prepDir, prepLog = prepLog)
        ## PREPARE PROTEIN STRUCTURE
        proteinPdbs = drPrep.prepare_protein_structure(config=config, outDir = prepDir, prepLog=prepLog)
        ## RE-COMBINE PROTEIN AND LIGAND PDB FILES
        wholePrepDir = p.join(prepDir,"WHOLE")
        os.makedirs(wholePrepDir,exist_ok=True)
        allPdbs = proteinPdbs + ligandPdbs
        outName = config["pathInfo"]["outputName"]
        mergedPdb = p.join(wholePrepDir,f"{outName}.pdb")
        mergePdbs(pdbList=allPdbs, outFile = mergedPdb)
        ## MAKE AMBER PARAMETER FILES WITH TLEAP
        inputCoords, amberParams = drPrep.make_amber_params(outDir = wholePrepDir,
                            ligandFileDict=ligandFileDict,
                            pdbFile= mergedPdb,
                            outName= outName,
                            prepLog= prepLog)

    else:
        ## PREPARE PROTEIN STRUCTURE
        proteinPdbs = drPrep.prepare_protein_structure(config=config, outDir = prepDir, prepLog = prepLog)  
        ## MERGE PROTEIN PDBS
        outName = config["pathInfo"]["outputName"]
        mergedPdb = p.join(p.join(prepDir,"PROT",f"{outName}.pdb"))
        mergePdbs(pdbList=proteinPdbs, outFile = mergedPdb)
        ## MAKE AMBER PARAMETER FILES WITH TLEAP
        inputCoords, amberParams = drPrep.make_amber_params(outDir = p.join(prepDir,"PROT"),
                                                        pdbFile= mergedPdb,
                                                        outName= outName,
                                                        prepLog = prepLog)

    drSim.run_simulation(config = config,
                   outDir = outDir,
                   inputCoords=inputCoords,
                   amberParams=amberParams)
#####################################################################################
#####################################################################################
