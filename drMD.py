## BASIC LIBS
import os
from os import path as p
## INPUT LIBS
import yaml
import argpass
## drMD UTILS
from drMDUtils import *
#####################################################################################
def read_inputs():
    ## create an argpass parser, read config file, snip off ".py" if on the end of file
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()
    configYaml=args.config

    ## Read config.yaml into a dictionary
    with open(configYaml,"r") as yamlFile:
        config = yaml.safe_load(yamlFile) 
    return config
#####################################################################################
def check_inputs(config):
    if not config["proteinInfo"]["nProteins"] == len(config["proteinInfo"]["proteins"]):
        print("Number of proteins in config does not match nProteins")
        exit()
    elif not config["ligandInfo"]["nLigands"] == len(config["ligandInfo"]["ligands"]):
        print("Number of ligands in config does not match nLigands")
        exit()
    elif not p.isfile(config["pathInfo"]["inputPdb"]):
        print("Input PDB does not exist")
        exit()       

#####################################################################################
#####################################################################################
def main():
    config = read_inputs()
    check_inputs(config=config)
    ## MAKE OUTPUT DIRECTORY
    outDir = config["pathInfo"]["outputDir"]
    prepDir = p.join(outDir,"01_prep")
    os.makedirs(outDir,exist_ok=True)
    os.makedirs(prepDir,exist_ok=True)
    ## SPLIT INPUT PDB INTO PROT AND ONE FILE PER LIGAND
    inputPdb = config["pathInfo"]["inputPdb"]
    split_input_pdb(inputPdb =inputPdb,
                    config = config,
                    outDir=prepDir)
    ## PREPARE LIGAND PARAMETERS, OUTPUT LIGAND PDBS
    ligandPdbs,ligandFileDict = prepare_ligand_parameters(config = config, outDir = prepDir)
    #ligandPdbs,ligandMol2s = ["/home/esp/scriptDevelopment/drMD/01_outputs/01_prep/FAD/FAD_amber.pdb",
    #              "/home/esp/scriptDevelopment/drMD/01_outputs/01_prep/PLM/PLM_amber.pdb"]
    ## PREPARE PROTEIN STRUCTURE
    proteinPdbs = prepare_protein_structure(config=config, outDir = prepDir)
    ## RE-COMBINE PROTEIN AND LIGAND PDB FILES
    wholePrepDir = p.join(prepDir,"WHOLE")
    os.makedirs(wholePrepDir,exist_ok=True)
    allPdbs = proteinPdbs + ligandPdbs
    outName = config["pathInfo"]["outputName"]
    mergedPdb = p.join(wholePrepDir,f"{outName}.pdb")
    mergePdbs(pdbList=allPdbs, outFile = mergedPdb)
    ## MAKE AMBER PARAMETER FILES WITH TLEAP
    make_amber_params(outDir = wholePrepDir,
                      ligandFileDict=ligandFileDict,
                       pdbFile= mergedPdb,
                         outName= outName)

#####################################################################################
#####################################################################################
main()