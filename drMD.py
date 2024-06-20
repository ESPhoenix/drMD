## BASIC LIBS
import os
from os import path as p
import pandas as pd
## INPUT LIBS
import yaml
import argpass
## CUSTOM DR MD MODULES
from pdbUtils import pdbUtils
import instruments.drPrep as drPrep
import instruments.drCleanup as drCleanup
import instruments.drOperator as drOperator
import instruments.drConfigInspector as drConfigInspector

## Multiprocessing
import concurrent.futures as cf
from tqdm import tqdm
from subprocess import run

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
    
    ## Read config.yaml into a dictionary and run checks on it
    batchConfig = drConfigInspector.read_and_validate_config()

    ## unpack batchConfig into variables for this function
    outDir = batchConfig["pathInfo"]["outputDir"]
    yamlDir = p.join(outDir,"00_configs")
    pdbDir = batchConfig["pathInfo"]["inputDir"]
    simInfo = batchConfig["simulationInfo"]
    parallelCPU = batchConfig["generalInfo"]["parallelCPU"]
    subprocessCpus = batchConfig["generalInfo"]["subprocessCpus"]
    
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
def manage_cpu_usage_for_subprocesses(mode , subprocessCpus = None):
    '''
    In ON mode, sets thread useage for OpenMP and OpenMM 
    In OFF mode, unsets thread useage

    Args:
        mode (string): "ON" or "OFF"
        subprocessCpus (int): will set thread useage if mode == "ON"
    Returns:
        Nothing
    '''
    if mode == "ON":
        # Set environment variables directly in Python
        os.environ['OMP_NUM_THREADS'] = str(subprocessCpus)
        os.environ['OPENMM_CPU_THREADS'] = str(subprocessCpus)

        omp_num_threads = os.environ.get('OMP_NUM_THREADS')
        print("OMP_NUM_THREADS:", omp_num_threads)
    elif mode == "OFF":
        # remove cpu useage limits
        run(f"unset OMP_NUM_THREADS", shell=True)
        run(f"unset OPENMM_CPU_THREADS", shell=True)

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
    proteinInfo, ligandInfo, generalInfo = extract_info(pdbDf, pdbDir, protName, yamlDir, batchConfig)

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

    # Create a ThreadPoolExecutor with the desired number of worker threads
    with cf.ThreadPoolExecutor(max_workers=parallelCPU) as executor:
        # Prepare the arguments for process_pdb_file
        inputArgs: list[tuple] = [(pdbFile, pdbDir, outDir, yamlDir, simInfo, batchConfig) for pdbFile in pdbFiles]

        # Use starmap to apply process_pdb_file to each PDB file
        # Use tqdm to display progress bar
        list(tqdm(executor.map(lambda args: process_pdb_file(*args), inputArgs), total=len(inputArgs)))

    # CLEAN UP
    drCleanup.clean_up_handler(batchConfig)
######################################################################################################
def extract_info(
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
    
    ## GET PROTEIN INFORMATION
    aminoAcids =   ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL']
    protDf = pdbDf[pdbDf["RES_NAME"].isin(aminoAcids)]
    ## CHECK TO SEE IF PROTEIN HAS HYDROGENS
    protH = False
    if (protDf["ELEMENT"] == "H").any():
        protH = True

    ## CREATE proteinInfo
    proteinInfo = {"nProteins": 1,
                   "proteins": [{"proteinName": f"{protName}",
                                 "protons": protH}]}   
    ## GET LIGAND INFORMATION 
    ligsDf = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcids)]
    ligNames = ligsDf["RES_NAME"].unique().tolist()

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
        ligList = []
        for ligName in ligNames:
            ligDf = ligsDf[ligsDf["RES_NAME"] == ligName]
            # deal with protonation
            ligH = False
            if (ligDf["ELEMENT"] == "H").any():
                ligH = True

            # deal with mol2
            ligMol2 = p.join(pdbDir, f"{ligName}.mol2")
            mol2 = False
            if p.isfile(ligMol2):
                mol2 = True
            # deal with frcmod
            ligFrcmod = p.join(pdbDir, f"{ligName}.frcmod")
            frcmod = False
            if p.isfile(ligFrcmod):
                frcmod = True  
            # deal with charge
            charge = drPrep.find_ligand_charge(ligDf, ligName, yamlDir, pH=7.4)
            # write to temporary dict, then to ligandInfo for config
            tmpDict = {"ligandName": ligName,
                       "protons": ligH,
                       "mol2": mol2,
                       "toppar": frcmod,
                       "charge": charge}
            ligList.append(tmpDict)
        nLigands = len(ligList)
        ligandInfo = {"nLigands": nLigands,
                      "ligands": ligList}

        return proteinInfo, ligandInfo, generalInfo
######################################################################################################
main()