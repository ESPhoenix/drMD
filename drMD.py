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
from  instruments import drConfigWriter
## Multiprocessing
import concurrent.futures as cf
from tqdm import tqdm
from subprocess import run
import multiprocessing as mp
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
    batchConfig: Dict = drConfigInspector.read_and_validate_config()
    ## unpack batchConfig into variables for this function
    outDir: DirectoryPath = batchConfig["pathInfo"]["outputDir"]
    yamlDir: DirectoryPath = p.join(outDir,"00_configs")
    pdbDir: DirectoryPath = batchConfig["pathInfo"]["inputDir"]
    parallelCPU: int = batchConfig["generalInfo"]["parallelCPU"]
    subprocessCpus: int = batchConfig["generalInfo"]["subprocessCpus"]


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
    # Get list of PDB files in the directory
    pdbFiles: list[str] = [p.join(pdbDir, pdbFile) for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]
    ## construct inputArgs for multiprocessing
    inputArgs: list[tuple] = [(pdbFile, batchConfig) for pdbFile in pdbFiles]

    # Function to update the progress bar
    def update_progress_bar(*args):
        with progress.get_lock():
            progress.value += 1
            pbar.update()

    parallelCpus: int = batchConfig["generalInfo"]["parallelCPU"]
    # Create a Pool with the desired number of worker processes
    with mp.Pool(processes=parallelCpus) as pool:
        # Create a shared Value for progress tracking
        progress: int = mp.Value('i', 0)
        
        # Create a tqdm progress bar
        with tqdm(total=len(inputArgs), colour="CYAN") as pbar:
            # Submit tasks to the pool
            results: Any = [pool.apply_async(process_pdb_file, args=args, callback=update_progress_bar) for args in inputArgs]
            
            # Close the pool and wait for the work to finish
            pool.close()
            pool.join()


    # CLEAN UP
    drCleanup.clean_up_handler(batchConfig)
######################################################################################################

main()