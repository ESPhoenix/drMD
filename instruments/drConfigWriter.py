## BASIC LIBS
import os
from os import path as p
import pandas as pd
## INPUT LIBS
import yaml
import argpass
## CUSTOM DR MD MODULES
from pdbUtils import pdbUtils
from  instruments import drPrep 
from  instruments import drCleanup 
from  instruments import drOperator 
from  instruments import drConfigInspector 
from  instruments import drSplash
from  instruments import drPdbTriage

## Multiprocessing
import concurrent.futures as cf
from tqdm import tqdm
from subprocess import run
import multiprocessing as mp
## CLEAN CODE
from typing import Optional, Dict, List, Tuple, Union
from instruments.drCustomClasses import FilePath, DirectoryPath
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
        tuple[dict, Optional[dict], dict]: A tuple containing proteinInfo, ligandInfo, and generalInfo
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


    ## copy generalInfo from batchConfig to new per run config
    generalInfo = batchConfig["generalInfo"]
    ## copy simInfo from batchConfig to new per run config
    simInfo = batchConfig["simulationInfo"]
    ## load pdb file into DataFrame
    pdbDf = pdbUtils.pdb2df(pdbFile)
    ## generate infomation on the protein portion of the pdb file
    proteinInfo = make_proteinInfo(pdbDf, protName)
    ## generate information on the ligands in the pdb file
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
        "generalInfo": generalInfo,
        "proteinInfo": proteinInfo,
        "simulationInfo": simInfo
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
        tuple[dict, Optional[dict], dict]: A tuple containing proteinInfo, ligandInfo, and generalInfo
    """
    
    ## GET PROTEIN AND ION ATOMS IN INPUT GEOMETRY
    aminoAcids, monovalentIons, multivalentIons = init_residue_name_lists()
    protDf = pdbDf[pdbDf["RES_NAME"].isin(aminoAcids) |
                   pdbDf["ATOM_NAME"].str.upper().isin(monovalentIons) |
                   pdbDf["ATOM_NAME"].str.upper().isin(multivalentIons)]
    
    ## CHECK TO SEE IF PROTEIN HAS HYDROGENS
    isProteinProtonated = False
    if (protDf["ELEMENT"] == "H").any():
        isProteinProtonated = True

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
    
    ## GET LIGAND ATOMS IN INPUT GEOMETRY 
    aminoAcids, monovalentIons, multivalentIons = init_residue_name_lists()

    ligandsDf = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcids) &
                    ~pdbDf["ATOM_NAME"].str.upper().isin(monovalentIons) &
                    ~pdbDf["ATOM_NAME"].str.upper().isin(multivalentIons)]

    ## GET NAMES OF LIGANDS
    ligNames = ligandsDf["RES_NAME"].unique().tolist()

    ## SKIP IF NOT LIGAND
    if len(ligNames) == 0:
        ligandInfo = None
        return None
    
    ## USE ligandInfo IF SUPPLIED IN BATCH CONFIG
    if "ligandInfo" in batchConfig:
        ligandInfo = batchConfig["ligandInfo"]
        return ligandInfo
    
    ## CREATE ligandInfo AUTOMATICALLY (WORKS FOR SIMPLE LIGANDS)
    else:
        ligandInfo = []
        for ligName in ligNames:
            ligandDf = ligandsDf[ligandsDf["RES_NAME"] == ligName]
            # detect protons in ligand
            isLigandProtonated = False
            if (ligandDf["ELEMENT"] == "H").any():
                isLigandProtonated = True
            # check for mol2 file in input pdb directory
            ligMol2 = p.join(pdbDir, f"{ligName}.mol2")
            isLigandMol2 = False
            if p.isfile(ligMol2):
                isLigandMol2 = True
            # check for frcmod file in input pdb directory
            ligFrcmod = p.join(pdbDir, f"{ligName}.frcmod")
            isLigandFrcmod = False
            if p.isfile(ligFrcmod):
                isLigandFrcmod = True  
            # deal with charge
            charge = drPrep.find_ligand_charge(ligandDf, ligName, yamlDir, pH=7.4)
            # write to temporary dict, then to ligandInfo for config
            tmpDict = {"ligandName": ligName,
                       "protons": isLigandProtonated,
                       "mol2": isLigandMol2,
                       "toppar": isLigandFrcmod,
                       "charge": charge}
            ligandInfo.append(tmpDict)

        return ligandInfo
    
######################################################################################
def init_residue_name_lists():
    aminoAcids =   ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL']
    
    
    monovalentIons = ["LI", "NA", "K", "RB", "CS", "TL", "CU", "AG", "NH4", "H3O", "F", "CL", "BR", "I"]

    multivalentIons = [
        "BE2", "CU2", "NI2", "PT2", "ZN2", "CO2", "PD2", "AG2", "CR2", "FE2", 
        "MG2", "V2", "MN2", "HG2", "CD2", "YB2", "CA2", "SN2", "PB2", "EU2", 
        "SR2", "SM2", "BA2", "RA2", "AL3", "FE3", "CR3", "IN3", "TL3", "Y3", 
        "LA3", "CE3", "PR3", "ND3", "SM3", "EU3", "GD3", "TB3", "DY3", "ER3", 
        "TM3", "LU3", "HF4", "ZR4", "CE4", "U4", "PU4", "TH4"
    ]

    return aminoAcids, monovalentIons, multivalentIons