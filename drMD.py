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
from  instruments import drSpash
from  instruments import drPdbTriage
## Multiprocessing
import concurrent.futures as cf
from tqdm import tqdm
from subprocess import run
import multiprocessing as mp
## CLEAN CODE
from typing import Optional

######################################################################################################
def main():
    '''
    Main function for drMD
    Unpacks batchConfig dictionary 
    Based on desired CPU useage, decides to use multiprocessing or not
    Based on desired CPU useage, manages CPU useage per run

    Args:
        Nothing
    Returns:
        Nothing
    '''
    ## print drMD logo
    drSpash.print_drMD_logo()

    ## get batchConfig
    batchConfig = drConfigInspector.read_and_validate_config()
    ## unpack batchConfig into variables for this function
    outDir = batchConfig["pathInfo"]["outputDir"]
    yamlDir = p.join(outDir,"00_configs")
    pdbDir = batchConfig["pathInfo"]["inputDir"]
    simInfo = batchConfig["simulationInfo"]
    parallelCPU = batchConfig["generalInfo"]["parallelCPU"]
    subprocessCpus = batchConfig["generalInfo"]["subprocessCpus"]


    ## run pdbTriage to detect commmon problems with pdb files
    drPdbTriage.pdb_triage(pdbDir, batchConfig)
    exit()
    ## set envorment variables for OpenMP and OpenMM
    manage_cpu_usage_for_subprocesses("ON",subprocessCpus)

    ## create yamlDir if it doesn't exist, this will be used to store per-run yaml files
    os.makedirs(yamlDir,exist_ok=True)
    ## run simulations in serial or paralell
    if parallelCPU == 1:
        run_serial(batchConfig, pdbDir, outDir, yamlDir, simInfo)
    elif parallelCPU > 1:
        run_parallel(parallelCPU, batchConfig, pdbDir, outDir, yamlDir, simInfo)

    ## unset envorment variables for OpenMP and OpenMM
    manage_cpu_usage_for_subprocesses("OFF")

######################################################################################################
def manage_cpu_usage_for_subprocesses(mode, subprocessCpus=None):
    '''
    In ON mode, sets thread usage for OpenMP and OpenMM 
    In OFF mode, unsets thread usage

    Args:
        mode (string): "ON" or "OFF"
        subprocessCpus (int): will set thread usage if mode == "ON"
    Returns:
        Nothing
    '''
    if mode == "ON":
        if subprocessCpus is not None:
            os.environ['OMP_NUM_THREADS'] = str(subprocessCpus)
            os.environ['OPENMM_CPU_THREADS'] = str(subprocessCpus)
        else:
            raise ValueError("subprocessCpus must be provided when mode is 'ON'")
    elif mode == "OFF":
        os.environ.pop('OMP_NUM_THREADS', None)
        os.environ.pop('OPENMM_CPU_THREADS', None)
    else:
        raise ValueError("mode must be 'ON' or 'OFF'")

#####################################################################################################
def process_pdb_file(pdbFile, pdbDir, outDir, yamlDir, simInfo, batchConfig):
    """
    Process a PDB file and run a sequence of MD simulations.

    Args:
        pdbFile (str): name of the PDB file to process
        pdbDir (str): input directory with PDB files
        outDir (str): output directory for simulation results
        yamlDir (str): directory to write YAML configuration files
        simInfo (dict): simulation information
        batchConfig (dict): batch configuration information
    """
    # Skip if not a PDB file
    fileData = p.splitext(pdbFile)
    if fileData[1] != ".pdb":
        return
    
    # Extract basic info, make dirs
    protName = fileData[0]
    pdbPath = p.join(pdbDir, pdbFile)
    runDir = p.join(outDir, protName)
    os.makedirs(runDir, exist_ok=True)

    # Convert to DataFrame, extract rest of info
    pdbDf = pdbUtils.pdb2df(pdbPath)
    proteinInfo, ligandInfo, generalInfo = make_per_protein_config(pdbDf, pdbDir, protName, yamlDir, batchConfig)

    # Get path info
    pathInfo = {
        "inputDir": pdbDir,
        "inputPdb": pdbPath,
        "outputDir": runDir,
        "outputName": protName
    }

    # Combine infos into one dict (ignore ligInfo if empty)
    if ligandInfo is None:
        config = {
            "pathInfo": pathInfo,
            "generalInfo": generalInfo,
            "proteinInfo": proteinInfo,
            "simulationInfo": simInfo
        }
    else:
        config = {
            "pathInfo": pathInfo,
            "generalInfo": generalInfo,
            "proteinInfo": proteinInfo,
            "ligandInfo": ligandInfo,
            "simulationInfo": simInfo
        }

    # Write config to YAML
    configYaml = p.join(yamlDir, f"{protName}_config.yaml")
    with open(configYaml, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    # Pass config YAML to drOperator to run a sequence of MD simulations
    drOperator.drMD_protocol(configYaml)

###################################################################################################### 
def run_serial(
    batchConfig: dict,
    pdbDir: str,
    outDir: str,
    yamlDir: str,
    simInfo: dict
) -> None:
    """
    Process each PDB file in the given directory serially.

    Args:
        batchConfig (dict): Batch configuration dictionary.
        pdbDir (str): Path to the directory containing PDB files.
        outDir (str): Path to the output directory.
        yamlDir (str): Path to the directory to write YAML configuration files.
        simInfo (dict): Simulation information dictionary.

    Returns:
        None
    """
    # Iterate over each file in the PDB directory
    for pdbFile in os.listdir(pdbDir):
        # Skip if the file is not a PDB file
        fileData = p.splitext(pdbFile)
        if fileData[1] != ".pdb":
            continue  
        # Process the PDB file
        process_pdb_file(pdbFile, pdbDir, outDir, yamlDir, simInfo, batchConfig)
    ## CLEAN UP
    drCleanup.clean_up_handler(batchConfig)


######################################################################################################
def run_parallel(
    parallelCPU: int,
    batchConfig: dict,
    pdbDir: str,
    outDir: str,
    yamlDir: str,
    simInfo: dict
) -> None:
    """
    Process each PDB file in the given directory in parallel using multiple worker threads.

    Args:
        parallelCPU (int): Number of worker threads to use.
        batchConfig (dict): Batch configuration dictionary.
        pdbDir (str): Path to the directory containing PDB files.
        outDir (str): Path to the output directory.
        yamlDir (str): Path to the directory to write YAML configuration files.
        simInfo (dict): Simulation information dictionary.

    Returns:
        None
    """
    # Get list of PDB files in the directory
    pdbFiles: list[str] = [pdbFile for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]

    inputArgs: list[tuple] = [(pdbFile, pdbDir, outDir, yamlDir, simInfo, batchConfig) for pdbFile in pdbFiles]

# Function to update the progress bar
    def update_progress_bar(*args):
        with progress.get_lock():
            progress.value += 1
            pbar.update()


            

    # Create a Pool with the desired number of worker processes
    with mp.Pool(processes=parallelCPU) as pool:
        # Create a shared Value for progress tracking
        progress = mp.Value('i', 0)
        
        # Create a tqdm progress bar
        with tqdm(total=len(inputArgs)) as pbar:
            # Submit tasks to the pool
            results = [pool.apply_async(process_pdb_file, args=args, callback=update_progress_bar) for args in inputArgs]
            
            # Close the pool and wait for the work to finish
            pool.close()
            pool.join()


    # CLEAN UP
    drCleanup.clean_up_handler(batchConfig)
######################################################################################################
def make_per_protein_config(
    pdbDf: pd.DataFrame,
    pdbDir: str,
    protName: str,
    yamlDir: str,
    batchConfig: dict,
) -> tuple[dict, Optional[dict], dict]:
    """
    This function writes a config.yaml file for each MD simulation

    Args:
        pdbDf (pd.DataFrame): DataFrame containing PDB information
        pdbDir (str): Path to the directory containing PDB files
        protName (str): Name of the protein
        yamlDir (str): Path to the directory to write YAML configuration files
        batchConfig (dict): Batch configuration dictionary

    Returns:
        tuple[dict, Optional[dict], dict]: A tuple containing proteinInfo, ligandInfo, and generalInfo
    """

    generalInfo = batchConfig["generalInfo"]
    
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
    ## GET LIGAND ATOMS IN INPUT GEOMETRY 
    ligandsDf = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcids) &
                    ~pdbDf["ATOM_NAME"].str.upper().isin(monovalentIons) &
                    ~pdbDf["ATOM_NAME"].str.upper().isin(multivalentIons)]

    ## GET NAMES OF LIGANDS
    ligNames = ligandsDf["RES_NAME"].unique().tolist()

    ## SKIP IF NOT LIGAND
    if len(ligNames) == 0:
        ligandInfo = None
        return proteinInfo, ligandInfo, generalInfo 
    
    ## USE ligandInfo IF SUPPLIED IN BATCH CONFIG
    if "ligandInfo" in batchConfig:
        ligandInfo = batchConfig["ligandInfo"]
        return proteinInfo, ligandInfo, generalInfo
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

        
        return proteinInfo, ligandInfo, generalInfo
######################################################################################################

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

main()