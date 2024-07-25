import os
from os import path as p 
from shutil import copy, rmtree
import pandas as pd

## clean code
from typing import List, Dict, Union, Any, Optional
from os import PathLike
## WellsWood modules
from pdbUtils import pdbUtils

## drMD modules
from instruments import drClusterizer
from instruments import drSelector
######################################################################################################
def clean_up_handler(batchConfig: dict) -> None:
    """
    Handler for performing post simulation operations based on the batch configuration.

    Args:
        batchConfig (dict): The batch configuration dictionary.

    Returns:
        None
    """

    if not "postSimulationInfo" in batchConfig:
        return


    cluster_handler(batchConfig)
    endpoint_handler(batchConfig)
    directory_cleanup_handler(batchConfig)

    
######################################################################################################
def cluster_handler(batchConfig: Dict) -> None:
    postSimulationInfo = batchConfig["postSimulationInfo"]
    pathInfo = batchConfig["pathInfo"]
    if not "clusterInfo" in postSimulationInfo:
        return

    clusterInfo: Dict = postSimulationInfo["clusterInfo"]
    allClusterPdbs: list[Union[PathLike, str]] = drClusterizer.clustering_manager(pathInfo, clusterInfo)
    if "removeAtoms" in clusterInfo:
        removeAtomsSelections = clusterInfo["removeAtoms"]
        for removeAtomsSelection in removeAtomsSelections:
            remove_atoms_from_pdbs(allClusterPdbs, removeAtomsSelection)
    if "collate" in clusterInfo:
        if clusterInfo["collate"]:
            collate_pdbs(allClusterPdbs, pathInfo)
######################################################################################################
def endpoint_handler(batchConfig: Dict) -> None:
    postSimulationInfo = batchConfig["postSimulationInfo"]
    pathInfo = batchConfig["pathInfo"]

    if not "endPointInfo" in postSimulationInfo:
        return
    endpointInfo: Dict = postSimulationInfo["endPointInfo"]
    endpointPdbs: List[Union[PathLike, str]] = get_endpoint_pdbs(endpointInfo, pathInfo)
    if "removeAtoms" in endpointInfo:
        removeAtomsSelections = [sele["selection"] for sele in endpointInfo["removeAtoms"]]
        remove_atoms_from_pdbs(endpointPdbs, removeAtomsSelections)
    if "collate" in endpointInfo:
        if endpointInfo["collate"]:
            collate_pdbs(endpointPdbs, pathInfo)
######################################################################################################
def directory_cleanup_handler(batchConfig: dict) -> None:
    postSimulationInfo = batchConfig["postSimulationInfo"]
    print(postSimulationInfo)
    if not "directoryCleanUpInfo" in postSimulationInfo:
        return
    directoryCleanupInfo: Dict = postSimulationInfo["directoryCleanUpInfo"]
    if  "removeAllSimulationDirs" in directoryCleanupInfo:
        if directoryCleanupInfo["removeAllSimulationDirs"]:
            remove_siulation_directories(batchConfig)
            return
    if "stepsToRemove" in directoryCleanupInfo:
        stepsToRemove = directoryCleanupInfo["stepsToRemove"]
        remove_step_directories(batchConfig, stepsToRemove)

######################################################################################################
def remove_step_directories(batchConfig: dict, stepsToRemove: list) -> None:
    inputDir = batchConfig["pathInfo"]["inputDir"]
    outDir = batchConfig["pathInfo"]["outputDir"]

    dirsToRemove = [p.join(outDir,p.splitext(file)[0],stepName) for 
                    file in os.listdir(inputDir) 
                    if p.splitext(file)[1] == ".pdb" 
                    for stepName in stepsToRemove]

    for dir in dirsToRemove:
        rmtree(dir)
######################################################################################################
def remove_siulation_directories(batchConfig: dict) -> None:
    print("WONK")
    inputDir = batchConfig["pathInfo"]["inputDir"]
    outDir = batchConfig["pathInfo"]["outputDir"]

    dirsToRemove = [p.join(outDir,p.splitext(file)[0]) for file in os.listdir(inputDir) if p.splitext(file)[1] == ".pdb"]
    for dir in dirsToRemove:
        print(dir)
        rmtree(dir)

######################################################################################################
def collate_pdbs(pdbFiles, pathInfo) -> None:
    print(f"-->\tCollating {len(pdbFiles)} PDB files into per-step directories")
    for pdbFile in pdbFiles:
        stepName = p.basename(p.dirname(pdbFile))
        stepCollateDir = p.join(pathInfo["outputDir"], "00_collated_pdbs", stepName)
        os.makedirs(stepCollateDir,exist_ok=True)
        copy(pdbFile, stepCollateDir)


######################################################################################################
def get_endpoint_pdbs(endPointInfo: Dict, pathInfo: Dict ) -> List[Union[PathLike, str]]:

    print("-->\tGetting endpoint PDB files")
    outDir = pathInfo["outputDir"]

    notRunDirs = ["00_configs", "01_ligand_parameters", "00_collated_pdbs", "00_drMD_logs"]

    runDirs = [p.join(outDir, dir) for dir in os.listdir(outDir) if not dir in notRunDirs]
    stepDirs = [p.join(runDir,stepDir) for runDir in runDirs for stepDir in endPointInfo["stepNames"]]

    endpointPdbs = [p.join(stepDir,stepPdb) for stepDir in stepDirs 
                    for stepPdb in os.listdir(stepDir)
                      if p.splitext(stepPdb)[1] == ".pdb"]
    return endpointPdbs
######################################################################################################
def remove_atoms_from_pdbs(
    pdbFiles: List[Union[PathLike, str]],
    removeAtomsSelections: List[Dict]) -> None:  

    for pdbFile in pdbFiles:
        pdbDf = pdbUtils.pdb2df(pdbFile)
        for removeAtomsSelection in removeAtomsSelections:
            indexesToRemove = drSelector.get_atom_indexes(removeAtomsSelection, pdbFile)
            droppedDf = pdbDf.drop(indexesToRemove)
        pdbUtils.df2pdb(droppedDf, pdbFile)
######################################################################################################

