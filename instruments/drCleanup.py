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

    
    postSimulationInfo: Dict = batchConfig["postSimulationInfo"]
    pathInfo: Dict = batchConfig["pathInfo"]


    if "clusterInfo" in postSimulationInfo:
        clusterInfo: Dict = postSimulationInfo["clusterInfo"]
        allClusterPdbs: list[Union[PathLike, str]] = drClusterizer.clustering_manager(pathInfo, clusterInfo)
        if "removeAtoms" in clusterInfo:
            removeAtomsSelections = clusterInfo["removeAtoms"]
            for removeAtomsSelection in removeAtomsSelections:
                remove_atoms_from_pdbs(allClusterPdbs, removeAtomsSelection)
        if "collate" in clusterInfo:
            if clusterInfo["collate"]:
                collate_pdbs(allClusterPdbs, pathInfo)

    if "endPointInfo" in postSimulationInfo:
        endpointInfo: Dict = postSimulationInfo["endPointInfo"]
        endpointPdbs: List[Union[PathLike, str]] = get_endpoint_pdbs(endpointInfo, pathInfo)
        if "removeAtoms" in endpointInfo:
            removeAtomsSelections = endpointInfo["removeAtoms"]
            for removeAtomsSelection in removeAtomsSelections:
                remove_atoms_from_pdbs(endpointPdbs, removeAtomsSelection)
        if "collate" in endpointInfo:
            if endpointInfo["collate"]:
                collate_pdbs(endpointPdbs, pathInfo)
######################################################################################################
def collate_pdbs(pdbFiles, pathInfo) -> None:
    print(f"-->\tCollating {len(pdbFiles)} PDB files into per-step directories")
    for pdbFile in pdbFiles:
        stepName = p.basename(p.dirname(p.dirname(pdbFile)))
        stepCollateDir = p.join(pathInfo["outputDir"], "00_collated_pdbs", stepName)
        os.makedirs(stepCollateDir,exist_ok=True)
        copy(pdbFile, stepCollateDir)


######################################################################################################
def get_endpoint_pdbs(endPointInfo: Dict, pathInfo: Dict ) -> List[Union[PathLike, str]]:
    outDir = pathInfo["outputDir"]

    notRunDirs = ["00_configs", "01_ligand_parameters", "00_collated_pdbs"]

    runDirs = [p.join(outDir, dir) for dir in os.listdir(outDir) if not dir in notRunDirs]
    stepDirs = [p.join(runDir,stepDir) for runDir in runDirs for stepDir in endPointInfo["stepNames"]]

    endpointPdbs = [p.join(stepDir,stepPdb) for stepDir in stepDirs 
                    for stepPdb in os.listdir(stepDir)
                      if p.splitext(stepPdb)[1] == ".pdb"]
    return endpointPdbs
######################################################################################################
def remove_atoms_from_pdbs(
    pdbFiles: List[Union[PathLike, str]],
    removeAtomsSelection: List[Dict]) -> None:  
    """
    Uses selection syntax to remove atoms from a list of pdb files.
    Args:
        pdbFiles (List[Union[PathLike, str]]): List of pdb files.
        removeAtomsSelection (Dict): Dict containing either a keyword or custom selection syntax
    """
    
    for pdbFile in pdbFiles:
        indexesToRemove = drSelector.get_atom_indexes(removeAtomsSelection, pdbFile)
        pdbDf = pdbUtils.pdb2df(pdbFile)
        droppedDf = pdbDf.drop(indexesToRemove)
        pdbUtils.df2pdb(droppedDf, pdbFile)
######################################################################################################

