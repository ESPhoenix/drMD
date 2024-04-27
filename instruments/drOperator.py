## BASIC LIBS
import os
from os import path as p
## drMD UTILS
from instruments import drPrep
## drMD simulation
from instruments import drSim
from pdbUtils import pdbUtils
import instruments.drConfigInspector as drConfigInspector
#####################################################################################
def drMD_protocol(configYaml):
    config = drConfigInspector.read_config(configYaml)
    # validate_file_config()
    ## MAKE OUTPUT DIRECTORY
    outDir = config["pathInfo"]["outputDir"]
    prepDir = p.join(outDir,"00_prep")

    amberParams = False
    inputCoords = False

    ###### skip prep if complete ######
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
    #####################################
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
        mergedPdb = p.join(wholePrepDir,"MERGED.pdb")
        pdbUtils.mergePdbs(pdbList=allPdbs, outFile = mergedPdb)
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
        mergedPdb = p.join(p.join(prepDir,"PROT","MERGED.pdb"))
        pdbUtils.mergePdbs(pdbList=proteinPdbs, outFile = mergedPdb)
        ## MAKE AMBER PARAMETER FILES WITH TLEAP
        inputCoords, amberParams = drPrep.make_amber_params(outDir = p.join(prepDir,"PROT"),
                                                        pdbFile= mergedPdb,
                                                        outName= outName,
                                                        prepLog = prepLog)

    drSim.run_simulation(config = config,
                   outDir = outDir,
                   inputCoords=inputCoords,
                   amberParams=amberParams,
                   pdbFile = mergedPdb)
#####################################################################################
