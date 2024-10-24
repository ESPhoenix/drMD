## BASIC PYTHON LIBRARIES
import os
from os import path as p
import pandas as pd
import yaml

## drMD LIBRARIES
from Surgery import drPrep
from UtilitiesCloset import drListInitiator

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## CLEAN CODE
from typing import Optional, Dict, List
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath
######################################################################################
def make_per_protein_config(
    pdbFile: FilePath,
    batchConfig: dict,
) -> tuple[dict, Optional[dict], dict]:
    """
    This function writes a config.yaml file for each MD simulation

    Args:
        pdbDf (pd.DataFrame): DataFrame containing PDB information
        protName (str): Name of the protein
        batchConfig (dict): Batch configuration dictionary

    Returns:
        tuple[dict, Optional[dict], dict]: A tuple containing proteinInfo, ligandInfo, and hardwareInfo
    """    
    # Skip if not a PDB file
    fileData: List[str] = p.splitext(pdbFile)
    if fileData[1] != ".pdb":
        return
    
    ## read info from batchConfig
    outDir = batchConfig["pathInfo"]["outputDir"]
    yamlDir = p.join(outDir, "00_configs")


    # Get filename of pdb file and use that to make a run directory
    protName: str = p.splitext(p.basename(pdbFile))[0]
    outDir = batchConfig["pathInfo"]["outputDir"]
    runDir: DirectoryPath = p.join(outDir, protName)
    os.makedirs(runDir, exist_ok=True)

    ## if config file has already been made, skip and return it
    configYaml = p.join(yamlDir, f"{protName}_config.yaml")
    if p.exists(configYaml):
        return configYaml


    ## load pdb file into DataFrame
    pdbDf = pdbUtils.pdb2df(pdbFile)
    ## generate infomation on the protein portion of the pdb file
    proteinInfo = make_proteinInfo(pdbDf, protName)
    ## generate information on the ligand in the pdb file
    inputDir = p.dirname(pdbFile)

    ligandInfo = make_ligandInfo(pdbDf, inputDir, yamlDir, batchConfig)
    # construct pathInfo
    pathInfo = {
        "inputDir": inputDir,
        "inputPdb": pdbFile,
        "outputDir": runDir,
        "outputName": protName
    }

    runConfig: Dict[Dict] = {
        "pathInfo": pathInfo,
        "hardwareInfo": batchConfig["hardwareInfo"],
        "proteinInfo": proteinInfo,
        "simulationInfo": batchConfig["simulationInfo"],
        "miscInfo": batchConfig["miscInfo"],
        "loggingInfo": batchConfig["loggingInfo"]

    }

    if bool(ligandInfo):
        runConfig["ligandInfo"] = ligandInfo
    # Write config to YAML
    configYaml = p.join(yamlDir, f"{protName}_config.yaml")
    with open(configYaml, "w") as f:
        yaml.dump(runConfig, f, default_flow_style=False)

    return configYaml
######################################################################################

def make_proteinInfo(
    pdbDf: pd.DataFrame,
    protName: str,
) -> Dict:
    """
    This function writes a config.yaml file for each MD simulation

    Args:
        pdbDf (pd.DataFrame): DataFrame containing PDB information
        pdbDir (str): Path to the directory containing PDB files
        protName (str): Name of the protein
        batchConfig (dict): Batch configuration dictionary

    Returns:
        tuple[dict, Optional[dict], dict]: A tuple containing proteinInfo, ligandInfo, and hardwareInfo
    """
    
    ## GET PROTEIN AND ION ATOMS IN INPUT GEOMETRY
    aminoAcidNames = drListInitiator.get_amino_acid_residue_names()
    ionAtomNames = drListInitiator.get_ion_names()

    protDf = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidNames) |
                   pdbDf["ATOM_NAME"].str.upper().isin(ionAtomNames)]
    
    ## CHECK TO SEE IF PROTEIN HAS HYDROGENS
    isProteinProtonated = False
    if (protDf["ELEMENT"] == "H").any():
        isProteinProtonated = True
    else:
        isProteinProtonated = False
    ## CREATE proteinInfo
    proteinInfo = {"proteinName":  protName,
                    "protons": isProteinProtonated}   
    
    return proteinInfo
######################################################################################
def make_ligandInfo(
    pdbDf: pd.DataFrame,
    pdbDir: str,
    yamlDir: str,
    batchConfig: dict,
    ) -> Optional[Dict]:
    ## USE ligandInfo IF SUPPLIED IN BATCH CONFIG
    if "ligandInfo" in batchConfig:    
        ligandInfo = batchConfig["ligandInfo"]
        return ligandInfo
    ## GET LIGAND ATOMS IN INPUT GEOMETRY 
    aminoAcidNames = drListInitiator.get_amino_acid_residue_names()
    ionAtomNames = drListInitiator.get_ion_names()

    ligandDf = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcidNames) &
                    ~pdbDf["ATOM_NAME"].str.upper().isin(ionAtomNames)]
    ## GET NAMES OF LIGANDS
    ligNames = ligandDf["RES_NAME"].unique().tolist()

    ## SKIP IF NOT LIGAND
    if len(ligNames) == 0:
        ligandInfo = None
        return ligandInfo
    

    
    ## CREATE ligandInfo AUTOMATICALLY (WORKS FOR SIMPLE LIGANDS)
    else:
        ligandInfo = []
        for ligName in ligNames:
            thisLigandDf = ligandDf[ligandDf["RES_NAME"] == ligName]
            # detect protons in ligand
            isLigandProtonated = False
            if (thisLigandDf["ELEMENT"] == "H").any():
                isLigandProtonated = True
            # check for mol2 file in input pdb directory
            ligMol2 = p.join(pdbDir, f"{ligName}.mol2")
            isLigandMol2 = False
            if p.isfile(ligMol2):
                isLigandMol2 = True
            ## checl for lib file in input pdb directory
            ligLib = p.join(pdbDir, f"{ligName}.lib")
            isLigandMol2 = False
            if p.isfile(ligLib):
                isLigandMol2 = True
            # check for frcmod file in input pdb directory
            ligFrcmod = p.join(pdbDir, f"{ligName}.frcmod")
            isLigandFrcmod = False
            if p.isfile(ligFrcmod):
                isLigandFrcmod = True  
            # deal with charge
            charge = drPrep.find_ligand_charge(thisLigandDf, ligName, yamlDir, pH=7.4)
            # write to temporary dict, then to ligandInfo for config
            tmpDict = {"ligandName": ligName,
                       "protons": isLigandProtonated,
                       "mol2": isLigandMol2,
                       "toppar": isLigandFrcmod,
                       "charge": charge}
            ligandInfo.append(tmpDict)

        return ligandInfo
    
######################################################################################
