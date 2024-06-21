import os
from os import path as p 
from shutil import copy, rmtree
from pdbUtils import pdbUtils
from typing import List, Dict, Union, Any, Optional
import pandas as pd
######################################################################################################
def clean_up_handler(batchConfig: dict) -> None:
    """
    Handler for performing clean-up operations based on the batch configuration.

    Args:
        batchConfig (dict): The batch configuration dictionary.

    Returns:
        None
    """
    # READ FROM batchConfig
    simulation_info: List[Dict[str, Union[str, int]]] = batchConfig["simulationInfo"]
    path_info: Dict[str, str] = batchConfig["pathInfo"]

    # RUN THROUGH OPTIONS IN cleanUpInfo
    if "cleanUpInfo" not in batchConfig:
        return  # If cleanUpInfo is not present in batchConfig, return

    clean_up_info: Dict[str, Union[bool, List[str]]] = batchConfig["cleanUpInfo"]

    # If getEndpointPdbs is True in cleanUpInfo, call get_endpoint_pdbs
    if "getEndpointPdbs" in clean_up_info and clean_up_info["getEndpointPdbs"]:
        get_endpoint_pdbs(simulation_info, path_info)

        # If either removeWaters or removeIons is True in cleanUpInfo, call remove_atoms_from_pdb
        if any(key in clean_up_info for key in ["removeWaters", "removeIons"]):
            remove_atoms_from_pdb(simulation_info, clean_up_info, path_info)
######################################################################################################
def get_endpoint_pdbs(
    simulationInfo: List[Dict[str, Union[str, int]]],
    pathInfo: Dict[str, str]
) -> None:
    """
    Gets the pdb files at the end of each simulation step and collates them into a single directory.

    Args:
        simulationInfo (List[Dict[str, Union[str, int]]]): List of dictionaries containing simulation information.
        pathInfo (Dict[str, str]): Dictionary containing path information.

    Returns:
        None
    """
    # Get input and output directories from pathInfo
    inputDir: str = pathInfo["inputDir"]
    outDir: str = pathInfo["outputDir"]

    # Get the names of all pdb files in the input directory
    inputNames: List[str] = [p.splitext(pdbFile)[0] 
                             for pdbFile in os.listdir(inputDir) 
                             if p.splitext(pdbFile)[1] == ".pdb"]

    # Create a new directory to collate the pdbs into
    collatedPdbDir: str = p.join(outDir, "collatedPdbs")
    os.makedirs(collatedPdbDir,exist_ok=True)

    # List of directories to exclude from collation
    excudeDirNames: List[str] = ["00_configs","01_ligand_parameters"]

    # Iterate over each simulation and collate the pdbs
    for sim in simulationInfo:
        stepName: str = sim["stepName"]
        tag: str = sim["simulationType"]
        collateSubDir: str = p.join(collatedPdbDir,stepName)
        os.makedirs(collateSubDir,exist_ok=True)
        excudeDirNames.append(stepName)

        # Iterate over each input name and collate the pdbs
        for inputName in inputNames:
            if inputName in excudeDirNames:
                continue
            stepDir: str = p.join(outDir,inputName,stepName)
            for file in os.listdir(stepDir):
                if p.splitext(file)[1] == ".pdb":
                    # Copy the pdb file to the collated directory
                    copy(p.join(stepDir,file),
                         p.join(collateSubDir,f"{inputName}_{tag}.pdb"))
######################################################################################################
def remove_atoms_from_pdb(
    simulationInfo: List[Dict[str, Union[str, int]]],
    cleanUpInfo: Dict[str, Union[bool, List[str]]],
    pathInfo: Dict[str, str]
) -> None:
    """
    Remove atoms from PDB files based on the cleanUpInfo dictionary.

    Args:
        simulationInfo (List[Dict[str, Union[str, int]]]): List of dictionaries containing simulation information.
        cleanUpInfo (Dict[str, Union[bool, List[str]]]): Dictionary containing information about atoms to remove.
        pathInfo (Dict[str, str]): Dictionary containing path information.

    Returns:
        None
    """
    # Get output directory from pathInfo
    outDir: str = pathInfo["outputDir"]

    # Create path to collated PDB directory
    collatedPdbDir: str = p.join(outDir, "collatedPdbs")

    # Iterate over each simulation
    for sim in simulationInfo:
        # Get step name from simulation info
        stepName: str = sim["stepName"]

        # Create path to subdirectory in collated PDB directory
        collateSubDir: str = p.join(collatedPdbDir, stepName)

        # Iterate over each PDB file in subdirectory
        for file in os.listdir(collateSubDir):
            # Check if file is a PDB file
            if not p.splitext(file)[1] == ".pdb":
                continue

            # Create path to PDB file
            pdbFile: str = p.join(collateSubDir, file)

            # Read PDB file into a DataFrame
            pdbDf: pd.DataFrame = pdbUtils.pdb2df(pdbFile)

            # Check if waters should be removed
            if "removeWaters" in cleanUpInfo:
                if cleanUpInfo["removeWaters"]:
                    # Remove waters from DataFrame
                    pdbDf = pdbDf[~pdbDf["RES_NAME"].isin(["HOH"])].copy()

            # Check if ions should be removed
            if "removeIons" in cleanUpInfo:
                if cleanUpInfo["removeIons"]:
                    # Remove ions from DataFrame
                    pdbDf = pdbDf[~pdbDf["RES_NAME"].isin(["Na+", "Cl-", "Mg2+", "F-"])].copy()

            # Write modified DataFrame back to PDB file
            pdbUtils.df2pdb(pdbDf, pdbFile)
######################################################################################################

