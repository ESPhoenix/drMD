## BASIC LIBS
import os
from os import path as p
import pandas as pd
import numpy as np
from shutil import move
## INPUT LIBS
import yaml
import argpass
## CUSTOM DR MD MODULES
from pdbUtils import pdbUtils
from  instruments import drPrep 
from  instruments import drCleanup 
from  instruments import drOperator 
from  instruments import drConfigTriage 
from  instruments import drSplash
from  instruments import drPdbTriage
from  instruments import drConfigWriter
## Multiprocessing
import concurrent.futures as cf
from tqdm import tqdm
from subprocess import run
import multiprocessing as mp
from tqdm.contrib.concurrent import process_map
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

## CLEAN CODE
from typing import Optional, Dict, List, Tuple, Union, Any
from instruments.drCustomClasses import FilePath, DirectoryPath

######################################################################################################
def main() -> None:
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
    drSplash.print_drMD_logo()

    ## get batchConfig
    batchConfig, configTriageLog  = drConfigTriage.read_and_validate_config()
    ## unpack batchConfig into variables for this function
    outDir: DirectoryPath = batchConfig["pathInfo"]["outputDir"]
    yamlDir: DirectoryPath = p.join(outDir,"00_configs")
    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    parallelCPU: int = batchConfig["generalInfo"]["parallelCPU"]
    subprocessCpus: int = batchConfig["generalInfo"]["subprocessCpus"]

    ## create logDir if it doesn't exist
    logDir: DirectoryPath = p.join(outDir, "00_drMD_logs")
    os.makedirs(logDir, exist_ok=True)
    os.replace(configTriageLog, p.join(logDir,"config_triage.log"))
    exit()

    ## run pdbTriage to detect commmon problems with pdb files
    drPdbTriage.pdb_triage(pdbDir, batchConfig)

    ## set envorment variables for OpenMP and OpenMM
    manage_cpu_usage_for_subprocesses("ON",subprocessCpus)

    ## create yamlDir if it doesn't exist, this will be used to store per-run yaml files
    os.makedirs(yamlDir,exist_ok=True)
    ## run simulations in serial or paralell
    if parallelCPU == 1:
        run_serial(batchConfig)
    elif parallelCPU > 1:
        run_parallel(batchConfig)

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
def process_pdb_file(pdbFile: FilePath, batchConfig: Dict):
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
    fileData: List[str] = p.splitext(pdbFile)
    if fileData[1] != ".pdb":
        return

    ## read batchConfig
    outDir: DirectoryPath = batchConfig["pathInfo"]["outputDir"]

    ## create a per-protein config
    runConfigYaml: FilePath = drConfigWriter.make_per_protein_config(pdbFile, batchConfig)

    # Pass config YAML to drOperator to run a sequence of MD simulations
    drOperator.drMD_protocol(runConfigYaml)

###################################################################################################### 
def run_serial(batchConfig: Dict) -> None:
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

    ## unpack batchConfig to get pdbDir
    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    ## create a list of PDB files
    pdbFiles = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
    # Iterate over each file in the PDB directory
    for pdbFile in pdbFiles:
        # Process the PDB file
        process_pdb_file(pdbFile, batchConfig)
    ## CLEAN UP
    drCleanup.clean_up_handler(batchConfig)


######################################################################################################
def run_parallel(batchConfig: Dict) -> None:
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

    ## read input directory from batchConfig
    pdbDir = batchConfig["pathInfo"]["inputDir"]
    parallelCpus: int = batchConfig["generalInfo"]["parallelCPU"]

    # Get list of PDB files in the directory
    pdbFiles: list[str] = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
    ## construct inputArgs for multiprocessing
    inputArgs: list[tuple] = [(pdbFile, batchConfig) for pdbFile in pdbFiles]
    ## create batched inputs
    batchedArgsWithPos = [(batch, pos) for pos, batch in enumerate(np.array_split(inputArgs, parallelCpus))]
    ## add a dummy batch to be used for printing logging
    batchedArgsWithPos = [(["dummy"], -1)] + batchedArgsWithPos

    ## run simulations in parallel
    process_map(per_core_worker, batchedArgsWithPos, 
                max_workers=parallelCpus)
    # CLEAN UP
    drCleanup.clean_up_handler(batchConfig)
######################################################################################################
def per_core_worker(batchedArgsWithPos: Tuple[Dict, int]) -> None:
    """
    Each core is passed a batch of arguments to process
    This function unpacks the arguments and handles the progress bar

    Args:
        batchedArgsWithPos List[(Tuple[Dict, int])]: 
            Alist of tuples containing input argunents for process_pdb_file and 
            the position of the core for the loading bar
    Returns:
        None
    """

    ## unpack batchedArgsWithPos into the batch of arguments for 
    batchedArgs, pos = batchedArgsWithPos

    cmap = plt.get_cmap('coolwarm', 32)
    colors = [mcolors.rgb2hex(cmap(i)) for i in range(32)]
    if pos == -1:
        with tqdm(total=1, position=0, bar_format='{desc}', 
                  colour="#000000", leave=True) as dummy_progress:
            dummy_progress.set_description_str("Logs:")
            dummy_progress.refresh()
    else:
        with tqdm(desc=f"Core {str(pos)}", total=len(batchedArgs), 
                position=pos+1, colour=colors[pos % len(colors)], 
                leave=False) as progress:
            for args in batchedArgs:
                pdbFile, batch_config = args
                process_pdb_file(pdbFile, batch_config)
                progress.update(1)
            progress.close()  

######################################################################################################

if __name__ == "__main__":
    main()
