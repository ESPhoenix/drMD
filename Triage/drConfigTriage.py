## BASIC PYTHON LIBRARIES
import argpass
import yaml
import os
from os import path as p
from pathlib import Path
## PARALLELISATION LIBRARIES
import multiprocessing as mp

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSplash

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## CLEAN CODE
from typing import Tuple, Union, List, Dict
from os import PathLike
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath


#####################################################################################
def validate_config(config) -> FilePath:
    """
    Main script for drConfigInspector:
    1. Accepts --config flag arg with argpass, reads yaml file into a dict
    2. Checks paths in pathInfo to see if they are real
    3. Checks inputs of cpuInfo to see if they are correct
    4. Checks each dockingOrder in dockingOrders to see if they are correct

    Returns:
    - config (dict)
    """

    ## set up logging
    topDir: DirectoryPath = os.getcwd()
    configTriageLog: FilePath = p.join(topDir,"config_triage.log")
    drLogger.setup_logging(configTriageLog)
    drLogger.log_info(f"Checking config file...",True)

    ## get some defaults for config
    configDefaults = init_config_defaults(topDir)
    config = {**configDefaults, **config}
    ## check each major section in config file
    ## throw errors with a nice splash screen if something goes wrong
    noDisordersFound = True
    configDisorders = {}
    config, configDisorders["miscInfo"], noDisordersFound = check_miscInfo(config, configDefaults, noDisordersFound)

    config, configDisorders["pathInfo"], noDisordersFound = check_pathInfo(config, configDefaults, noDisordersFound)
    config, configDisorders["hardwareInfo"], noDisordersFound  = check_hardwareInfo(config, configDefaults, noDisordersFound)

    configDisorders["simulationInfo"], noDisordersFound  = check_simulationInfo(config, noDisordersFound)

    configDisorders["postSimulationInfo"], noDisordersFound  = check_postSimulationInfo(config,noDisordersFound)
  
    configDisorders["ligandInfo"], noDisordersFound  = check_ligandInfo(config, noDisordersFound)
                
    if noDisordersFound:
        drLogger.log_info("No disorders found in batch config, proceeding with simulations...", True)
        return config, configTriageLog
    else:
        drSplash.print_config_error(configDisorders)
       

#####################################################################################
def init_config_defaults(topDir: DirectoryPath) -> dict:
    configDefaults = {
        "pathInfo": {
            "inputDir": topDir,
            "outputDir": p.join(topDir, "outputs")
        },
        "hardwareInfo": {
            "parallelCPU": 1,
            "platform": "CPU",
            "subprocessCpus": 1
        },

        "miscInfo": {   
            "pH": 7,
            "firstAidMaxRetries": 10,
            "boxGeometry": "cubic",
            "writeMyMethodsSection": True,
            "skipPdbTriage": False,
            "trajectorySelections": [{"selection": {"keyword": 'all'}}]}
    }
    return configDefaults


#####################################################################################
def check_pathInfo(config: dict, configDefaults: dict, noDisorders: bool) -> Tuple[dict, dict, bool]:
    """
    Checks for pathInfo entry in config
    Checks paths in pathInfo to see if they are real
    Don't check outputDir, this will be made automatically
    """
    noDisorders = True
    pathInfoDisorders = {}
    ## log this check
    drLogger.log_info(f"Checking pathInfo...")

    ## check to see if pathInfo in config
    pathInfo = config.get("pathInfo", None)
    if pathInfo is None:
        config["pathInfo"] = configDefaults["pathInfo"]
        drLogger.log_info(f"No pathInfo specified, using defaults")
        pathInfoDisorders["inptDir"] = "Automatic Default Used!"
        pathInfoDisorders["outputDir"] = "Automatic Default Used!"

        return config, pathInfoDisorders, False

    ## validate inputDir
    inputDir = pathInfo.get("inputDir", None)
    if inputDir is None:
        config["pathInfo"]["inputDir"] = configDefaults["pathInfo"]["inputDir"]
        pathInfoDisorders["inptDir"] = "Automatic Default Used!"
        noDisorders = False
    else:
        inputDirPathProblem = validate_path(f"inputDir", inputDir)
        if inputDirPathProblem is not None:
            pathInfoDisorders["inptDir"] = inputDirPathProblem
            noDisorders = False
        else:
            pathInfoDisorders["inptDir"] = None

    ## validate outputDir
    outputDir = pathInfo.get("outputDir", None)
    if outputDir is None:
        config["pathInfo"]["outputDir"] = configDefaults["pathInfo"]["outputDir"]
        pathInfoDisorders["outputDir"] = "Automatic Default Used!"
        noDisorders = False
    else:
        if not isinstance(Path(outputDir), PathLike):
            pathInfoDisorders["outputDir"] = "outputDir must be a PathLike string!"
            noDisorders = False
        else:
            pathInfoDisorders["outputDir"] = None

    ## log that pathInfo is correct
    drLogger.log_info(f"pathInfo is correct...")

    return config, pathInfoDisorders, noDisorders




#####################################################################################
def check_hardwareInfo(config: dict, configDefaults: dict, noDisorders) -> Tuple[dict, dict, bool]:
    """
    Checks hardwareInfo in config
    Makes sure that CPU allocations are properly formatted
    Makes sure that the "Platform" specified is an allowed value 
    """
    ## log this check
    drLogger.log_info(f"Checking hardwareInfo...")
    noDisorders = True
    hardwareInfoDisorders = {}
    ## check if hardwareInfo in config
    hardwareInfo = config.get("hardwareInfo", None)
    if hardwareInfo is None:
        config["hardwareInfo"] = configDefaults["hardwareInfo"]
        for argName in ["parallelCPU", "platform", "subprocessCpus"]:
            hardwareInfoDisorders[argName] = "Automatic Default Used!"
        drLogger.log_info(f"No hardwareInfo specified, using defaults")
        return config, hardwareInfoDisorders, True

    ## validate parallelCPU
    parallelCPU = hardwareInfo.get("parallelCPU", None)
    if parallelCPU is None:
        ## use a default value
        config["hardwareInfo"]["parallelCPU"] = configDefaults["hardwareInfo"]["parallelCPU"]
        hardwareInfoDisorders["parallelCPU"] = "Automatic Default Used!"
    else:
        if not isinstance(parallelCPU, int):
            hardwareInfoDisorders["parallelCPU"] = "parallelCPU is not an int"
            noDisorders = False
        else:
            if parallelCPU < 1:
                hardwareInfoDisorders["parallelCPU"] = "parallelCPU is less than 1, this must be a positive integer"
                noDisorders = False
            else:
                hardwareInfoDisorders["parallelCPU"] = None


    ## validate subprocessCpus
    subprocessCpus = hardwareInfo.get("subprocessCpus", None)
    if subprocessCpus is None:
        # use a default value
        config["hardwareInfo"]["subprocessCpus"] = configDefaults["hardwareInfo"]["subprocessCpus"]
        hardwareInfoDisorders["subprocessCpus"] = "Automatic Default Used!"
    else:
        if not isinstance(subprocessCpus, int):
            hardwareInfoDisorders["subprocessCpus"] = "subprocessCpus is not an int"
            noDisorders = False
        else:
            if subprocessCpus < 1:
                hardwareInfoDisorders["subprocessCpus"] = "subprocessCpus is less than 1, this must be a positive integer"
                noDisorders = False
            else:
                hardwareInfoDisorders["subprocessCpus"] = None

    ## check to see if the number of CPU cores requested is less than the number of CPU cores available
    if isinstance(parallelCPU,int) and isinstance(subprocessCpus,int):
        systemCpus = mp.cpu_count()
        if parallelCPU * subprocessCpus > systemCpus:
            hardwareInfoDisorders["totalCpuUseage"] = "Number for CPU cores requested exceeds number of CPU cores available, change the values of parallelCPU and subprocessCpus"
            noDisorders = False
    ## validate platform
    platform = hardwareInfo.get("platform", None)
    if platform is None:
        ## use a default value
        config["hardwareInfo"]["platform"] = configDefaults["hardwareInfo"]["platform"]
        hardwareInfoDisorders["platform"] = "Automatic Default Used!"
    else:
        if platform not in ["CUDA", "OPENCL", "CPU"]:
            hardwareInfoDisorders["platform"] = "platform is not 'CUDA', 'OPENCL', or 'CPU'"
            noDisorders = False
        else:
            hardwareInfoDisorders["platform"] = None


    return config, hardwareInfoDisorders, noDisorders



#####################################################################################
def check_miscInfo(config, configDefaults, noDisorders) -> Tuple[dict,dict,bool]:
    ## log this check
    drLogger.log_info(f"Checking miscInfo...")
    noDisorders = True
    miscInfoDisorders = {}

    miscInfo = config.get("miscInfo", None)
    ## if no miscInfo in config, use defaults and report disorders
    if miscInfo is None:
        config["miscInfo"] = configDefaults["miscInfo"]
        for argName in ["pH",
                         "firstAidMaxRetries",
                           "boxGeometry",
                             "writeMyMethodsSection",
                               "skipPdbTriage",
                                 "trajectorySelections"]:
            miscInfoDisorders[argName] = "Automatic Default Used!"

        return config, miscInfoDisorders, True
    
    ## validate pH
    pH = miscInfo.get("pH", None)
    if pH is None:
        ## use a default value
        config["miscInfo"]["pH"] = configDefaults["miscInfo"]["pH"]
        miscInfoDisorders["pH"] = "No pH specified, using default of 7"
    else:
        if not isinstance(pH, (int, float)):
            miscInfoDisorders["pH"] = "pH must be an int or float between 0 and 14"
            noDisorders = False
        else:
            if pH < 0 or pH > 14:
                miscInfoDisorders["pH"] = "pH must be an int or float between 0 and 14"
                noDisorders = False
            else:
                miscInfoDisorders["pH"] = None

    ## validate firstAidMaxRetries
    firstAidMaxRetries = miscInfo.get("firstAidMaxRetries", None)
    if firstAidMaxRetries is None:
        ## use a default value
        config["miscInfo"]["firstAidMaxRetries"] = configDefaults["miscInfo"]["firstAidMaxRetries"]
        miscInfoDisorders["firstAidMaxRetries"] = "No firstAidMaxRetries specified, using default of 10"
        noDisorders = False
    else:
        if not isinstance(firstAidMaxRetries, int):
            miscInfoDisorders["firstAidMaxRetries"] = "firstAidMaxRetries must be an int greater than 0"
            noDisorders = False
        else:
            if firstAidMaxRetries < 1:
                miscInfoDisorders["firstAidMaxRetries"] = "firstAidMaxRetries must be an int greater than 0"
                noDisorders = False
            else:
                miscInfoDisorders["firstAidMaxRetries"] = None

    ## validate boxGeometry
    boxGeometry = miscInfo.get("boxGeometry", None)
    if boxGeometry is None:
        ## use a default value
        config["miscInfo"]["boxGeometry"] = configDefaults["miscInfo"]["boxGeometry"]
        miscInfoDisorders["boxGeometry"] = "No boxGeometry specified, using default of 'cubic'"
    else:
        if boxGeometry not in ["cubic", "ocatahedral"]:
            miscInfoDisorders["boxGeometry"] = "boxGeometry must be either 'cubic' or 'ocatahedral'"
            noDisorders = False
        else:
            miscInfoDisorders["boxGeometry"] = None

    ## validate skipPdbTriage
    skipPdbTriage = miscInfo.get("skipPdbTriage", None)
    if skipPdbTriage is None:
        ## use a default value
        config["miscInfo"]["skipPdbTriage"] = configDefaults["miscInfo"]["skipPdbTriage"]
        miscInfoDisorders["skipPdbTriage"] = "No skipPdbTriage specified, using default of False"
    else:
        if not isinstance(skipPdbTriage, bool):
            miscInfoDisorders["boxGeometry"] = "skipPdbTriage must be either True or False"
            noDisorders = False
        else:
            miscInfoDisorders["skipPdbTriage"] = None
    
    ## validate writeMyMethodsSection
    writeMyMethodsSection = miscInfo.get("writeMyMethodsSection", None)
    if writeMyMethodsSection is None:
        ## use a default value
        config["miscInfo"]["writeMyMethodsSection"] = configDefaults["miscInfo"]["writeMyMethodsSection"]
        miscInfoDisorders["writeMyMethodsSection"] = "No writeMyMethodsSection specified, using default of True"
    else:
        if not isinstance(writeMyMethodsSection, bool):
            miscInfoDisorders["writeMyMethodsSection"] = "writeMyMethodsSection must be either True or False"
            noDisorders = False
        else:
            miscInfoDisorders["writeMyMethodsSection"] = None

    ## validate trajectorySelections
    trajectorySelections = miscInfo.get("trajectorySelections", None)
    if trajectorySelections is None:
        ## use a default value  
        config["miscInfo"]["trajectorySelections"] = configDefaults["miscInfo"]["trajectorySelections"]
        miscInfoDisorders["trajectorySelections"] = "No trajectorySelections specified, all atoms will be selected by default"
    else:
        if not isinstance(trajectorySelections, list):
            miscInfoDisorders["trajectorySelections"] = "trajectorySelections must be a list of selection dictionaries (see README for syntax!)"
            noDisorders = False
        elif len(trajectorySelections) == 0:
            miscInfoDisorders["trajectorySelections"] = "trajectorySelections must be a list of selection dictionaries (see README for syntax!)"
            noDisorders = False
        else:
            for selection in trajectorySelections:
                selectionDisorder = check_selection(selection)
                if len(selectionDisorder) > 0:
                    miscInfoDisorders["trajectorySelections"] = selectionDisorder
                    noDisorders = False
                else:
                    miscInfoDisorders["trajectorySelections"] = None

    return config, miscInfoDisorders, noDisorders

#####################################################################################
def check_simulationInfo(config: dict, noDisorders) -> Tuple[dict, bool]:
    """
    Checks for simulationInfo in config
    Depending on the type of simulation, checks your parameters
    """
    ## log this check
    noDisorders = True
    simulationInfo = config.get("simulationInfo", None)


    ## check for problems with the simuationInfo dict
    if simulationInfo is None:
        return config,  "No simulationInfo found", False
    if not isinstance(simulationInfo, list):
        return config, "simulationInfo must be a list of dicts", False
    if len(simulationInfo) == 0:
        return config, "simulationInfo must have at least one entry", False
    
    simulationInfoDisorders = {}
    counter = 0
    for simulation in simulationInfo:
        counter += 1
        disorders = {}
        drLogger.log_info(f"Checking {simulation['stepName']}...")
        disorders, stepName, simulationType, noDisorders  = check_shared_simulation_options(simulation, disorders, noDisorders)
        ## if we don't have a simulationType or stepName, further checks wont work.
        if stepName is None:
            disorders["stepName"] = "stepName must be specified"
            simulationInfoDisorders[f"Unnamed_step_{str(counter)}"] = "stepName must be specified"
            noDisorders = False
            continue
        if simulationType is None:
            disorders["simulationType"] = "simulationType must be specified"
            simulationInfoDisorders[stepName] = disorders
            noDisorders = False
            continue

        if simulationType in ["NVT", "NPT"]:
            disorders, noDisorders = check_nvt_npt_options(simulation, stepName, disorders, noDisorders)
        elif simulationType == "META":
            disorders, noDisorders = check_metadynamics_options(simulation, stepName, disorders, noDisorders)

        restraintsInfo = simulation.get("restraintInfo", None)
        if restraintsInfo:
            disorders, noDisorders = check_restraintInfo(restraintsInfo, disorders, noDisorders)

        simulationInfoDisorders[stepName] = disorders

    return simulationInfoDisorders, noDisorders

#################################################################################################


def check_restraintInfo(restraintInfo: dict, disorders: dict, noDisorders: bool) -> Tuple[dict, bool]:
    if not isinstance(restraintInfo, list):
        disorders["restraintInfo"] = "restraintInfo must be a dictionary"
        return disorders, False
    if len(restraintInfo) == 0:
        disorders["restraintInfo"] = "restraintInfo must have at least one entry"
        return disorders, False
    
    ## check each entry in restraintInfo

    for restraintIndex, info in enumerate(restraintInfo):
        if not isinstance(info, dict):
            disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo must be a dictionary (see README for syntax!)"
            noDisorders = False
        else:
            ## check restraintType
            restraintType = info.get("restraintType", None)
            if restraintType is None:
                disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo must have a 'restraintType' key"
                noDisorders = False
            else:
                if not isinstance(restraintType, str):
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo['restraintType'] must be a string"
                    noDisorders = False
                if not restraintType in ["position", "distance", "angle", "torsion"]:
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo['restraintType'] must be one of 'position', 'distance', 'angle', 'torsion'"
                    noDisorders = False
            ## check selection for restraint to act upon
            restrantSelection = info.get("selection", None)
            if  restrantSelection is None:
                disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo must have a 'selection' key"
                noDisorders = False
            else:
                if not isinstance(restrantSelection, dict):
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo['selection'] must be a dictionary (see README for syntax!)"
                    noDisorders = False
                else:
                    selectionDisorder = check_selection({"selection": restrantSelection})
                    if len(selectionDisorder) > 0:
                        disorders["restraintInfo"][f"restraint_{restraintIndex}"] = selectionDisorder
                        noDisorders = False
            ## check parameters for restraint
            restraintParamers = info.get("parameters", None)
            if restraintParamers is None:
                disorders["restraintInfo"][f"restraint_{restraintIndex}"] = "each entry in restraintInfo must have a 'parameters' key"
                noDisorders = False
            else:
                restraintParamProblems = check_restraint_parameters(restraintType, restraintParamers)
                if len(restraintParamProblems) > 0:
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = restraintParamProblems
                    noDisorders = False
                else:
                    disorders["restraintInfo"] = {}
                    disorders["restraintInfo"][f"restraint_{restraintIndex}"] = None

    return disorders, noDisorders

    


#####################################################################################
def check_postSimulationInfo(config: dict, noDisorders) -> Tuple[dict,bool]:
    """
    Checks for postSimulationInfo in config
    """
    postSimulationInfoDisorders = {}
    ## check for postSimulationInfo in config
    postSimulationInfo = config.get("postSimulationInfo", None)
    if postSimulationInfo is None:
        return None
    
    endPointInfo = postSimulationInfo.get("endPointInfo", None)
    if endPointInfo is not None:
        endPointInfoDisorders, noDisorders = check_endPointInfo(endPointInfo, noDisorders)
        postSimulationInfoDisorders["endPointInfo"] = endPointInfoDisorders

    clusterInfo = postSimulationInfo.get("clusterInfo", None)
    if clusterInfo is not None:
        postSimulationInfoDisorders["clusterInfo"] = check_clusterInfo(clusterInfo, noDisorders)

    return postSimulationInfoDisorders, noDisorders



#####################################################################################
def check_endPointInfo(endPointInfo: dict, noDisorders) -> Tuple[dict, bool]:
    """
    Checks for endPointInfo in config
    """
    ## log this check
    drLogger.log_info(f"Checking endPointInfo...")

    endPointDisorders = {}
    ## check if endPointInfo is a dictionary
    if not isinstance(endPointInfo, dict):
        return "endPointInfo must be a dictionary"
    ## check stepNames
    stepNames = endPointInfo.get("stepNames", None)
    if stepNames is None:
        endPointDisorders["stepNames"] = "endPointInfo must have a 'stepNames' entry"
        noDisorders = False
    else:
        if not isinstance(stepNames, list):
            endPointDisorders["stepNames"] = "stepNames must be a list"
            noDisorders = False
        ## ensure that stepNames is a list of strings
        if not all(isinstance(stepName, str) for stepName in stepNames):
            endPointDisorders["stepNames"] = "stepNames must be a list of strings"
            noDisorders = False
        ## ensure that stepNames is not empty
        if  len(stepNames) == 0:
            endPointDisorders["stepNames"] = "stepNames must not be empty"
            noDisorders = False
        else:
            endPointDisorders["stepNames"] = None
    ## check removeAtoms
    removeAtoms = endPointInfo.get("removeAtoms", None)
    if removeAtoms is not None:
        if not isinstance(removeAtoms, list):
            endPointDisorders["removeAtoms"] = "removeAtoms must be a list"
            noDisorders = False
        if not all(isinstance(removeAtom, dict) for removeAtom in removeAtoms):
            endPointDisorders["removeAtoms"] = "removeAtoms must be a list of dictionaries"
            noDisorders = False
        for selection in removeAtoms:
            removeAtomSelectionDisorders = check_selection(selection)
            if len(removeAtomSelectionDisorders) > 0:
                endPointDisorders["removeAtoms"] = removeAtomSelectionDisorders
                noDisorders = False
            else:
                endPointDisorders["removeAtoms"] = None
        
    return endPointDisorders, noDisorders
#####################################################################################
def check_clusterInfo(clusterInfo: dict, noDisorders: bool) -> Tuple[dict,bool]:
    """
    Checks for clusterInfo in config
    """
    ## log this check
    drLogger.log_info(f"Checking clusterInfo...")
    clusterInfoDisorders = {}
    ## check if clusterInfo is a dictionary
    if not isinstance(clusterInfo, dict):
        return "clusterInfo must be a dictionary", False
    ## check stepNames
    stepNames = clusterInfo.get("stepNames", None)
    if stepNames is None:
        clusterInfoDisorders["stepNames"] = "clusterInfo must have a 'stepNames' entry"
        noDisorders = False
    else:
        if not isinstance(stepNames, list):
            clusterInfoDisorders["stepNames"] = "clusterInfo['stepNames'] must be a list"
            noDisorders = False
        ## ensure that stepNames is a list of strings
        elif not all(isinstance(stepName, str) for stepName in stepNames):
            clusterInfoDisorders["stepNames"] = "clusterInfo['stepNames'] must be a list of strings"
            noDisorders = False
        ## ensure that stepNames is not empty
        elif  len(stepNames) == 0:
            clusterInfoDisorders["stepNames"] = "clusterInfo['stepNames'] must not be empty"
            noDisorders = False
        else:
            clusterInfoDisorders["stepNames"] = None

    ## check nClusters
    nClusters = clusterInfo.get("nClusters", None)
    if nClusters is None:
        clusterInfoDisorders["nClusters"] = "nClusters must be specified as a positive int"
        noDisorders = False
    else:
        if not isinstance(nClusters, int):
            clusterInfoDisorders["nClusters"] = "nClusters must be an int"
            noDisorders = False
        elif nClusters < 1:
            clusterInfoDisorders["nClusters"] = "nClusters must be specified as a positive int"
            noDisorders = False
        else:
            clusterInfoDisorders["nClusters"] = None

    ## check clusterBy
    clusterBy = clusterInfo.get("clusterBy", None)
    if clusterBy is None:
        clusterInfoDisorders["clusterBy"] = "clusterInfo must have a 'clusterBy' entry"
        noDisorders = False
    else:
        if not isinstance(clusterBy, list):
            clusterInfoDisorders["clusterBy"] = "clusterBy must be a list"
            noDisorders = False
        elif not all(isinstance(clusterSelection, dict) for clusterSelection in clusterBy):
            clusterInfoDisorders["clusterBy"] = "clusterBy must be a list of dictionaries"
            noDisorders = False
        elif len(clusterBy) == 0:
            clusterInfoDisorders["clusterBy"] = "clusterBy must not be empty"
            noDisorders = False
        for clusterSelection in clusterBy:
            clusterSelectionDisorders = check_selection(clusterSelection)
            if len(clusterSelectionDisorders) > 0:
                clusterInfoDisorders["clusterBy"] = clusterSelectionDisorders
                noDisorders = False
            else:
                clusterInfoDisorders["clusterBy"] = None
    ## check removeAtoms
    removeAtoms = clusterInfo.get("removeAtoms", None)
    if removeAtoms is not None:
        if not isinstance(removeAtoms, list):
            clusterInfoDisorders["removeAtoms"] = "endPointInfo['removeAtoms'] must be a list"
            noDisorders = False
        elif not all(isinstance(removeAtom, dict) for removeAtom in removeAtoms):
            clusterInfoDisorders["removeAtoms"] = "endPointInfo['removeAtoms'] must be a list of dictionaries"
            noDisorders = False
        for selection in removeAtoms:
            removeAtomSelectionDisorders = check_selection(selection)
            if len(removeAtomSelectionDisorders) > 0:
                clusterInfoDisorders["removeAtoms"] = removeAtomSelectionDisorders
                noDisorders = False
            else:
                clusterInfoDisorders["removeAtoms"] = None

    return clusterInfoDisorders, noDisorders


#########################################################################
def check_ligandInfo(config: dict, noDisorders: bool) -> Tuple[dict, bool]:
    """
    Checks optional ligandInfo entry in config

    Args:
        config (dict): The main configuration dictionary

    Raises:
        TypeError: If ligandInfo is not a list of dictionaries, or if ligandName is not a string
        TypeError: If protons, charge, toppar, or mol2 is not a boolean
        ValueError: If ligandInfo is an empty list
        ValueError: If ligandInfo does not have at least one entry
    """
    ## log this check
    drLogger.log_info(f"Checking ligandInfo...")
    # Check if ligandInfo in config
    ligandInfo = config.get("ligandInfo", None)
    ## if there is no ligandInfo specified, not a problem, return None
    if ligandInfo is None:
        return None, noDisorders
    elif len(ligandInfo) == 0:
        return "ligandInfo must be a contain at least one ligand dictionary (see README for more info)", False
    
    ligandInfoDisorders = {}
    # Check each entry in ligandInfo
    for ligandIndex, ligand in enumerate(ligandInfo):
        # Check if ligand is a dictionary
        if not isinstance(ligand, dict):
            ligandInfoDisorders[f"ligand_{ligandIndex}"] = "ligand entry must be a dictionary"
            noDisorders = False 
            continue
        ## check ligandName
        ligandName = ligand.get("ligandName", None)
        if ligandName is None:
            ligandInfoDisorders[f"ligand_{ligandIndex}"]["ligandName"] = "each ligand entry must have a ligandName entry as a unique string"
            noDisorders = False
            ## set a temporary ligand name for reporting disorders
            ligandName = f"ligand_{ligandIndex}"
        else:
            ligandInfoDisorders[ligandName] = {}
        if not isinstance(ligandName, str):
            ligandInfoDisorders[ligandName]["ligandName"] = "each ligand entry must have a ligandName entry as a unique string"
            noDisorders = False
        else:
            ligandInfoDisorders[ligandName]["ligandName"] = None

        ## check boolean flags
        for argName in ["protons", "toppar", "mol2"]:
            argValue = ligand.get(argName, None)
            if argValue is None:
                ligandInfoDisorders[ligandName][argName] = f"each ligand entry must have a {argName} entry as a bool"
                noDisorders = False
            elif not isinstance(argValue, bool):
                ligandInfoDisorders[ligandName][argName] = f"{argName} must be a bool"
                noDisorders = False
            else:
                ligandInfoDisorders[ligandName][argName] = None
        ## check charge
        charge = ligand.get("charge", None)
        if charge is None:
            ligandInfoDisorders[ligandName]["charge"] = "each ligand entry must have a charge entry as an int"
            noDisorders = False
        elif not isinstance(charge, int):
            ligandInfoDisorders[ligandName]["charge"] = "charge must be an integer"
            noDisorders = False
        else:
            ligandInfoDisorders[ligandName]["charge"] = None
 
    return ligandInfoDisorders, noDisorders




#########################################################################
def check_cluster_trajectory_options(clusterTrajectory: dict, stepName: str) -> None:
    """
    Checks for options in clusterTrajectory.

    Args:
        clusterTrajectory (dict): A dictionary containing options for clustering a trajectory.
        stepName (str): The name of the step.

    Raises:
        TypeError: If clusterTrajectory is not a dictionary.
        ValueError: If clusterTrajectory is empty.
        TypeError: If nClusters is not an integer.
    """
    # Check if clusterTrajectory is a dictionary
    if not isinstance(clusterTrajectory, dict):
        raise TypeError(f"clusterTrajectory in {stepName} must be a dict")
    
    # Check if clusterTrajectory has at least one entry
    if len(clusterTrajectory) == 0:
        raise ValueError(f"clusterTrajectory in {stepName} must contain some entries")
    
    # Check for required args for clustering a trajectory
    nClusters, selection = check_info_for_args(clusterTrajectory, stepName, ["nClusters", "selection"], optional=False)
    
    # Check if nClusters is an integer
    if not isinstance(nClusters, int):
        raise TypeError(f"nClusters in {stepName} must be an int")
    
    # Check selection for clustering a trajectory
    check_selection(selection, stepName)



#########################################################################
def check_restraint_parameters(restraintType: str, parameters: dict) -> None:

    parameterProblems = []
    # Check for force constants
    forceConstant = parameters.get("k", None)
    if not forceConstant:
        parameterProblems.append("No force constat k provided for restraints")
    else:
        if not isinstance(forceConstant, (int, float)):
            parameterProblems.append("force constants must be a positive number")
        elif forceConstant <= 0:
            parameterProblems.append("force constants must be a positive number")

    ## deal with phi0 in torsions
    if restraintType.upper() == "TORSION":
        phi0 = parameters.get("phi0", None)
        if  phi0 is None:
            parameterProblems.append("phi0 parameter must be provided for torsion restraints")
        else:
            if not isinstance(phi0, (int, float)):
                parameterProblems.append("phi0 parameter must be a number for torsion restraints")
            if phi0 < -180 or phi0 > 180:
                parameterProblems.append("phi0 parameter must be between -180 and 180 for torsion restraints")
        
    elif restraintType.upper() == "DISTANCE":
        r0 = parameters.get("r0", None)
        if  r0 is None:
            parameterProblems.append("r0 parameter must be provided for distance restraints")
        else:
            if not isinstance(r0, (int, float)):
                parameterProblems.append("r0 parameter must be a number for distance restraints")
            if r0 < 0:
                parameterProblems.append("r0 parameter must be positive for distance restraints")

    elif restraintType.upper() == "ANGLE":
        theta0 = parameters.get("theta0", None)
        if  theta0 is None:
            parameterProblems.append("theta0 parameter must be provided for angle restraints")
        else:
            if not isinstance(theta0, (int, float)):
                parameterProblems.append("theta0 parameter must be a number for angle restraints")
            if theta0 < 0 or theta0 > 360:
                parameterProblems.append("theta0 parameter must be between 0 and 360 for angle restraints")


    return parameterProblems

#########################################################################
def check_selection(selection: dict) -> list:

    selectionDisorders = []
    ## check keyword
    subSelection = selection.get("selection", None)
    if subSelection is None:
        return ["No selection specified in selection"]

    keyword = subSelection.get("keyword", None)
    if keyword is None:
        return ["No keyword specified in selection"]

    if not isinstance(keyword, str):
        return [f"keyword must be a string, not {type(keyword)}"]
    
    if not keyword in ["all", "protein", "ligand", "water", "ions", "custom"]:
        return [f"selection keywords incorrect see README.md for more details"]
    
    ## check custom selection syntax
    if keyword == "custom":
        customSelection = subSelection.get("customSelection", None)
        if customSelection is None:
            selectionDisorders.append("No customSelction specified in selection")
            return selectionDisorders
        if not isinstance(customSelection, list):
            selectionDisorders.append("customSelection must be a list of selection dictionaries (see README for more details)")
            return selectionDisorders
        ## check each selection
        for customSelctionDict in customSelection:
            ## check chainId
            chainId = customSelctionDict.get("CHAIN_ID", None)
            if chainId is None:
                selectionDisorders.append("No CHAIN_ID specified in selection")
            ## check resName
            resName = customSelctionDict.get("RES_NAME", None)
            if resName is None:
                selectionDisorders.append("No RES_NAME specified in selection")
            elif resName == "all":
                pass
            elif not isinstance(resName, str):
                selectionDisorders.append("RES_NAME must be a three-letter string")
            elif len(resName) != 3:
                selectionDisorders.append("RES_NAME must be a three-letter string")
            ## check resId
            resId = customSelctionDict.get("RES_ID", None)
            if resId is None:
                selectionDisorders.append("No RES_ID specified in selection")
            elif resId == "all":
                pass
            elif not isinstance(resId, int):
                selectionDisorders.append("RES_ID must be an integer")
            ## check atomName
            atomName = customSelctionDict.get("ATOM_NAME", None)
            if atomName is None:
                selectionDisorders.append("No ATOM_NAME specified in selection")
            elif not isinstance(atomName, str):
                selectionDisorders.append("ATOM_NAME must be a string")
            elif len(atomName) > 4:
                selectionDisorders.append("ATOM_NAME must be a string less than 4 characters")

    return selectionDisorders


#########################################################################
def check_metadynamics_options(simulation: dict, stepName: str, disorders: dict, noDisorders: bool) -> Tuple[dict,bool]:
    drLogger.log_info(f"Checking metaDynamicsInfo for {stepName}...")
    ## check metaDynamicsInfo
    metaDynamicsInfo = simulation.get("metaDynamicsInfo", None)
    if metaDynamicsInfo is None:
        disorders["metaDynamicsInfo"] = "No metadynamicsInfo found in step with simluationType of META"
        return disorders, False
    if not isinstance(metaDynamicsInfo, dict):
        disorders["metaDynamicsInfo"] = "metaDynamicsInfo must be a dictionary"
        return disorders, False
    
    ## check height parameter
    height = metaDynamicsInfo.get("height", None)
    if height is not None:
        disorders["metaDynamicsInfo"]["height"] = "No height specified in metaDynamicsInfo"
        noDisorders = False
    else:
        if not isinstance(height, (int, float)):
            disorders["metaDynamicsInfo"]["height"] = "height must be a number"
            noDisorders = False
        elif height <= 0:
            disorders["metaDynamicsInfo"]["height"] = "height must be positive"
            noDisorders = False
        else:
            disorders["metaDynamicsInfo"]["height"] = None

    
    ## check biasFactor parameter
    biasFactor = metaDynamicsInfo.get("biasFactor", None)
    if biasFactor is None:
        disorders["metaDynamicsInfo"]["biasFactor"] = "No biasFactor specified in metaDynamicsInfo"
        noDisorders = False
    else:
        if not isinstance(biasFactor, (int, float)):
            disorders["metaDynamicsInfo"]["biasFactor"] = "biasFactor must be a number"
            noDisorders = False
        elif biasFactor <= 0:
            disorders["metaDynamicsInfo"]["biasFactor"] = "biasFactor must be positive"
            noDisorders = False
        else:
            disorders["metaDynamicsInfo"]["biasFactor"] = None
    
    ## check biases parameter
    biases = metaDynamicsInfo.get("biases", None)
    if biases is None:
        disorders["metaDynamicsInfo"]["biases"] = "No biases specified in metaDynamicsInfo"
        noDisorders = False
    else:
        if not isinstance(biases, list):
            disorders["metaDynamicsInfo"]["biases"] = "biases must be a list of biases (check README for more details)"
            noDisorders = False
        elif len(biases) == 0:
            disorders["metaDynamicsInfo"]["biases"] = "biases must contain at least one bias variable"
            noDisorders = False
        else:
            for biasCount, bias in enumerate(biases):
                if not isinstance(bias, dict):
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = "biases must be a list of biases (check README for more details)"
                    noDisorders = False
                biasVar = bias.get("biasVar", None)
                if biasVar is None:
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = "No biasVar specified in bias"
                    noDisorders = False
                else:
                    if not isinstance(biasVar, str):
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = "biasVar must be 'rmsd', 'torsion', 'distance', or 'angle'"
                        noDisorders = False
                    elif not biasVar.lower() in ["rmsd", "torsion", "distance", "angle"]:
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = "biasVar must be 'rmsd', 'torsion', 'distance', or 'angle'"
                        noDisorders = False
                biasSelection = bias.get("selection", None)
                if biasSelection is None:
                    disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = "No selection specified in bias"
                    noDisorders = False
                else:
                    if not isinstance(biasSelection, dict):
                        disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = "selection must be a dictionary"
                        noDisorders = False
                    else:
                        selectionDisorders = check_selection(biasSelection, stepName)
                        if len(selectionDisorders) > 0:
                            disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = selectionDisorders
                            noDisorders = False
                        else:
                            disorders["metaDynamicsInfo"]["biases"][f"bias_{biasCount}"] = None
                        
    return disorders, noDisorders
#########################################################################
def check_nvt_npt_options(simulation: dict, stepName: str, disorders: dict, noDisorders: bool) -> Tuple[dict,bool]:
    ## check for required args for a nvt or npt simulation

    ## check duration
    duration = simulation.get("duration", None)
    if duration is None:
        disorders["duration"] = "No duration specified in simulation"
        noDisorders = False
    else:
        timeCheckProblems = check_time_input(duration, "duration", stepName)
        if timeCheckProblems is not None:
            disorders["duration"] = timeCheckProblems
            noDisorders = False
        else:
            disorders["duration"] = None

    # check logInterval
    logInterval = simulation.get("logInterval", None)
    if logInterval is None:
        disorders["logInterval"] = "No logInterval specified in simulation"
        noDisorders = False
    else:
        timeCheckProblems = check_time_input(logInterval, "logInterval", stepName)
        if timeCheckProblems is not None:
            disorders["logInterval"] = timeCheckProblems
            noDisorders = False
        else:
            disorders["logInterval"] = None

    ## check heavyProtons
    heavyProtons = simulation.get("heavyProtons", None)
    if heavyProtons is None:
        disorders["heavyProtons"] = "No heavyProtons specified in simulation"
        noDisorders = False
    else:
        if not isinstance(heavyProtons, bool):
            disorders["heavyProtons"] = "heavyProtons must be a boolean"
            noDisorders = False
        else:
            disorders["heavyProtons"] = None
        
    return disorders, noDisorders

        

#########################################################################
def check_time_input(timeInputValue: str, timeInputName: str, stepName: str) -> None:
    problemText = f"{timeInputName} in {stepName} must be in the format \"10 ns\""

    if not isinstance(timeInputValue, str):
        return problemText
    timeInputData = timeInputValue.split()
    try:
        numberAsInt = int(timeInputData[0])
    except:
        return problemText
    if not timeInputData[1] in ["fs","ps","ns","ms"]:
        return problemText
    
    return None

#########################################################################
def check_shared_simulation_options(simulation: dict, disorders: dict, noDisorders) -> Tuple[dict, str, str, bool]:
    ## check simulation step name
    stepName = simulation.get("stepName", None)
    if stepName is None: 
        disorders["stepName"] = "No stepName specified in simulation"
        noDisorders = False
    elif not isinstance(stepName, str):
        disorders["stepName"] = "stepName must be a string"
        noDisorders = False
    elif " " in stepName:
        disorders["stepName"] = "No whitespace allowed in stepName"
        noDisorders = False
    else:
        disorders["stepName"] = None


    ## check simulationType
    simulationType = simulation.get("simulationType", None)
    if simulationType is None:
        disorders["simulationType"] = "No simulationType specified in simulation"
        noDisorders = False
    elif not simulationType.upper() in ["EM", "NVT", "NPT", "META"]:
        disorders["simulationType"] = "simulationType in simulation must be one of the following: 'EM', 'NVT', 'NPT', 'META'"
        noDisorders = False
    else:
        disorders["simulationType"] = None
    ## check for either temperature or temparatureRange in simulation
    temperature = simulation.get("temperature", None)
    tempRange = simulation.get("temperatureRange", None)
    ## check for both temperature and tempRange (this is not allowed!)
    if temperature and tempRange:
        disorders["temperature"] = "Cannot specify both temperature and temperatureRange in simulation"
        noDisorders = False
    ## check for neither temperature or tempRange (this is also not allowed!)
    elif not temperature and not tempRange:
        disorders["temperature"] = "Must specify either temperature or temperatureRange in simulation"
        noDisorders = False
    ## check temperature
    elif temperature:
        if not isinstance(temperature, int):
            disorders["temperature"] = "Temperature in simulation must be an int"
            noDisorders = False
        if temperature < 0:
            disorders["temperature"] = "Temperature in simulation must be a positive int"
            noDisorders = False
        else:
            disorders["temperature"] = None
    ## check tempRange
    elif tempRange:
        if not isinstance(tempRange, list):
            disorders["tempRange"] = "TemperatureRange in simulation must be a list of ints"
            noDisorders = False
        if len(tempRange) == 0:
            disorders["tempRange"] = "TemperatureRange in simulation must be a list of at least one int"
            noDisorders = False
        for temp in tempRange:
            if not isinstance(temp, int):
                disorders["tempRange"] = "TemperatureRange in simulation must be a list of ints"
                noDisorders = False
            if temp < 0:
                disorders["tempRange"] = "TemperatureRange in simulation must be a list of positive ints"
                noDisorders = False

    return disorders, stepName, simulationType, noDisorders

#########################################################################
def validate_path(argName: str, argPath: Union[PathLike, str]) -> str:
    """
    Check to see if a path variable is indeed the correct type
    Check to see if the path exists
    """
    if  not isinstance(Path(argPath), (PathLike, str)) :
        return f"The config argument {argName} = {argPath} is not a PathLike."
    # Check if the path exists
    if not p.exists(argPath):
        return f"The config argument {argName} = {argPath} does not exist."
    return None
#####################################################################################
def get_config_input_arg() -> FilePath:
    """
    Sets up argpass to read the config.yaml file from command line
    Reads a YAML file using the "--config" flag with argpass

    Returns:
    - configFile (FilePath)
    """
    # create an argpass parser, read config file,
    parser = argpass.ArgumentParser()
    parser.add_argument(f"--config")
    args = parser.parse_args()

    configFile: Union[PathLike, str] = args.config

    return configFile
#####################################################################################

def read_input_yaml(configFile: FilePath) -> dict:
    """
    Reads YAML file into a dict

    Returns:
    - config (dict)
    """
    # Read config.yaml into a dictionary
    try:
        with open(configFile, "r") as yamlFile:
            config: dict = yaml.safe_load(yamlFile)
            return config
    except FileNotFoundError:
        raise FileNotFoundError(f"config file {configFile} not found")
    except yaml.YAMLError as exc:
        raise yaml.YAMLError(f"Error parsing YAML file:", exc)



#####################################################################################
def read_config(configYaml: str) -> dict:
    """
    Reads a config.yaml file and returns its contents as a dictionary.
    This is performed by drOperator on automatically generated configs,

    Args:
        configYaml (str): The path to the config.yaml file.

    Returns:
        dict: The contents of the config.yaml file as a dictionary.
    """
    # Open the config.yaml file and read its contents into a dictionary
    with open(configYaml, "r") as yamlFile:
        config: dict = yaml.safe_load(yamlFile)

    return config

#####################################################################################
