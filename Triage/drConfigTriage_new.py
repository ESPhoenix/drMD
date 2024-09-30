## BASIC PYTHON LIBRARIES
import argpass
import yaml
import os
from os import path as p

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
    configDisorders = []
    for function in [check_pathInfo, check_hardwareInfo, check_miscInfo]:
        config, disorders = function(config, configDefaults)
        if len (disorders) > 0:
            configDisorders.append(disorders)



    simulationInfoDisorders = check_simulationInfo(config)
    if len(simulationInfoDisorders) > 0:
        configDisorders.append(simulationInfoDisorders)

    endPointDisorders, clusterDisorders = check_postSimulationInfo(config)


    if len(endPointDisorders) > 0:
        configDisorders.append(endPointDisorders)
    if len(clusterDisorders) > 0:
        configDisorders.append(clusterDisorders)

    for problem in configDisorders:
       print(problem)
    exit()

    drLogger.log_info(f"Config file is correct", True)
    drLogger.close_logging()
    return configTriageLog
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
            "skipPdbTriage": False
        }
    }
    return configDefaults


#####################################################################################
def check_pathInfo(config: dict, configDefaults: dict) -> None:
    """
    Checks for pathInfo entry in config
    Checks paths in pathInfo to see if they are real
    Don't check outputDir, this will be made automatically
    """

    pathInfoDisorders = []
    ## log this check
    drLogger.log_info(f"Checking pathInfo...")

    ## check to see if pathInfo in config
    pathInfo = config.get("pathInfo", None)
    if pathInfo is None:
        config["pathInfo"] = configDefaults["pathInfo"]
        drLogger.log_info(f"No pathInfo specified, using defaults")
        return config, None

    ## validate inputDir
    inputDir = pathInfo.get("inputDir", None)
    if inputDir is None:
        config["pathInfo"]["inputDir"] = configDefaults["pathInfo"]["inputDir"]
        drLogger.log_info(f"No inputDir specified, using current dir for input")
        # pathInfoDisorders.append("No inputDir specified, using current dir for input")
    else:
        inputDirPathProblem = validate_path(f"inputDir", inputDir)
        if inputDirPathProblem is not None:
            pathInfoDisorders.append(inputDirPathProblem)

    ## validate outputDir
    outputDir = pathInfo.get("outputDir", None)
    if outputDir is None:
        config["pathInfo"]["outputDir"] = configDefaults["pathInfo"]["outputDir"]
        drLogger.log_info(f"No outputDir specified, using current dir for outputDir")

    ## log that pathInfo is correct
    drLogger.log_info(f"pathInfo is correct...")

    if len(pathInfoDisorders) > 0:
        return config, {"pathInfo": pathInfoDisorders}
    else:
        return config, []



#####################################################################################
def check_hardwareInfo(config: dict, configDefaults: dict) -> None:
    """
    Checks hardwareInfo in config
    Makes sure that CPU allocations are properly formatted
    Makes sure that the "Platform" specified is an allowed value 
    """
    ## log this check
    drLogger.log_info(f"Checking hardwareInfo...")

    hardwareInfoDisorders = []
    ## check if hardwareInfo in config
    hardwareInfo = config.get("hardwareInfo", None)
    if hardwareInfo is None:
        config["hardwareInfo"] = configDefaults["hardwareInfo"]
        drLogger.log_info(f"No hardwareInfo specified, using defaults")
        return config, None

    ## validate parallelCPU
    parallelCPU = hardwareInfo.get("parallelCPU", None)
    if parallelCPU is None:
        ## use a default value
        config["hardwareInfo"]["parallelCPU"] = configDefaults["hardwareInfo"]["parallelCPU"]
        drLogger.log_info(f"No parallelCPU specified, using default of 1")
    else:
        if not isinstance(parallelCPU, int):
            hardwareInfoDisorders.append("parallelCPU is not an int")
        else:
            if parallelCPU < 1:
                hardwareInfoDisorders.append("parallelCPU is less than 1, this must be a positive integer")


    ## validate subprocessCpus
    subprocessCpus = hardwareInfo.get("subprocessCpus", None)
    if subprocessCpus is None:
        # use a default value
        config["hardwareInfo"]["subprocessCpus"] = configDefaults["hardwareInfo"]["subprocessCpus"]
        drLogger.log_info(f"No subprocessCpus specified, using default of 1")

    else:
        if not isinstance(subprocessCpus, int):
            hardwareInfoDisorders.append("subprocessCpus is not an int")
        else:
            if subprocessCpus < 1:
                hardwareInfoDisorders.append("subprocessCpus is less than 1, this must be a positive integer")

    ## check to see if the number of CPU cores requested is less than the number of CPU cores available
    if isinstance(parallelCPU,int) and isinstance(subprocessCpus,int):
        systemCpus = mp.cpu_count()
        if parallelCPU * subprocessCpus > systemCpus:
            hardwareInfoDisorders.append("Number for CPU cores requested exceeds number of CPU cores available, change the values of parallelCPU and subprocessCpus")

    ## validate platform
    platform = hardwareInfo.get("platform", None)
    if platform is None:
        ## use a default value
        config["hardwareInfo"]["platform"] = configDefaults["hardwareInfo"]["platform"]
        drLogger.log_info(f"No platform specified, using default of 'CPU'")
    else:
        if platform not in ["CUDA", "OPENCL", "CPU"]:
            hardwareInfoDisorders.append("platform is not 'CUDA', 'OPENCL', or 'CPU'")
    ## log that hardwareInfo is correct
    drLogger.log_info(f"hardwareInfo is correct...")

    if len(hardwareInfoDisorders) > 0:
        return config, {"hardwareInfo": hardwareInfoDisorders}
    else:
        return config, []


#####################################################################################
def check_miscInfo(config, configDefaults):
    ## log this check
    drLogger.log_info(f"Checking miscInfo...")

    miscInfoDisorders = []

    miscInfo = config.get("miscInfo", None)
    if miscInfo is None:
        drLogger.log_info(f"No miscInfo specified, using defaults")
        config["miscInfo"] = configDefaults["miscInfo"]
        return config, None
    
    ## validate pH
    pH = miscInfo.get("pH", None)
    if pH is None:
        ## use a default value
        drLogger.log_info(f"No pH specified, using default of 7")
        config["miscInfo"]["pH"] = configDefaults["miscInfo"]["pH"]
        
    else:
        if not isinstance(pH, (int, float)):
            miscInfoDisorders.append("pH must be an int or float between 0 and 14")
        else:
            if pH < 0 or pH > 14:
                miscInfoDisorders.append("pH must be an int or float between 0 and 14")

    ## validate firstAidMaxRetries
    firstAidMaxRetries = miscInfo.get("firstAidMaxRetries", None)
    if firstAidMaxRetries is None:
        ## use a default value
        drLogger.log_info(f"No firstAidMaxRetries specified, using default of 10")
        config["miscInfo"]["firstAidMaxRetries"] = configDefaults["miscInfo"]["firstAidMaxRetries"]
    else:
        if not isinstance(firstAidMaxRetries, int):
            miscInfoDisorders.append("firstAidMaxRetries must be an int greater than 0")
        else:
            if firstAidMaxRetries < 1:
                miscInfoDisorders.append("firstAidMaxRetries must be an int greater than 0")

    ## validate boxGeometry
    boxGeometry = miscInfo.get("boxGeometry", None)
    if boxGeometry is None:
        ## use a default value
        drLogger.log_info(f"No boxGeometry specified, using default of 'cubic'")
        config["miscInfo"]["boxGeometry"] = configDefaults["miscInfo"]["boxGeometry"]
    else:
        if boxGeometry not in ["cubic", "ocatahedral"]:
            miscInfoDisorders.append("boxGeometry must be either 'cubic' or 'ocatahedral'")

    ## validate skipPdbTriage
    skipPdbTriage = miscInfo.get("skipPdbTriage", None)
    if skipPdbTriage is None:
        ## use a default value
        drLogger.log_info(f"No skipPdbTriage specified, using default of False")
        config["miscInfo"]["skipPdbTriage"] = configDefaults["miscInfo"]["skipPdbTriage"]
    else:
        if not isinstance(skipPdbTriage, bool):
            miscInfoDisorders.append("skipPdbTriage must be either True or False")
    
    ## validate writeMyMethodsSection
    writeMyMethodsSection = miscInfo.get("writeMyMethodsSection", None)
    if writeMyMethodsSection is None:
        ## use a default value
        drLogger.log_info(f"No writeMyMethodsSection specified, using default of True")
        config["miscInfo"]["writeMyMethodsSection"] = configDefaults["miscInfo"]["writeMyMethodsSection"]
    else:
        if not isinstance(writeMyMethodsSection, bool):
            miscInfoDisorders.append("writeMyMethodsSection must be either True or False")

    ## validate trajectorySelections
    trajectorySelections = miscInfo.get("trajectorySelections", None)
    if trajectorySelections is None:
        ## use a default value  
        drLogger.log_info(f"No trajectorySelections specified, using defaults")
        config["miscInfo"]["trajectorySelections"] = configDefaults["miscInfo"]["trajectorySelections"]
    else:
        if not isinstance(trajectorySelections, list):
            miscInfoDisorders.append("trajectorySelections must be a list of selection dictionaries (see README for syntax!)")
        elif len(trajectorySelections) == 0:
            miscInfoDisorders.append("trajectorySelections must be a list of selection dictionaries (see README for syntax!)")
        else:
            for selection in trajectorySelections:
                selectionDisorder = check_selection(selection)
                if len(selectionDisorder) > 0:
                    miscInfoDisorders.append(selectionDisorder)

    if len(miscInfoDisorders) > 0:
        return config, {"miscInfo": miscInfoDisorders}
    else:
        return config, []
#####################################################################################
def check_simulationInfo(config: dict) -> None:
    """
    Checks for simulationInfo in config
    Depending on the type of simulation, checks your parameters
    """
    ## log this check
    drLogger.log_info(f"Checking simulationInfo...")
    simulationInfo = config.get("simulationInfo", None)

    if simulationInfo is None:
        return config, {"simluationInfo": "No simulationInfo found"}
    if not isinstance(simulationInfo, list):
        return config, {"simluationInfo": "simulationInfo must be a list of dicts"}
    if len(simulationInfo) == 0:
        return config, {"simluationInfo": "simulationInfo must have at least one entry"}
    
    allSimulationProblems = {}
    for simulation in simulationInfo:
        drLogger.log_info(f"Checking {simulation['stepName']}...")
        simulationProblems, stepName, simulationType  = check_shared_simulation_options(simulation)
        ## if we don't have a simulationType or stepName, further checks wont work.
        if simulationType is None or stepName is None:
            allSimulationProblems.append(simulationProblems)
            continue
        if simulationType in ["NVT", "NPT"]:
            optionProblems = check_nvt_npt_options(simulation, stepName)
            simulationProblems.extend(optionProblems)
        elif simulationType == "META":
            metaProblems = check_metadynamics_options(simulation, stepName)
            simulationProblems.extend(metaProblems)

        restraintsInfo = simulation.get("restraintInfo", None)
        if restraintsInfo:
            restraintsInfoProblems = check_restraintInfo(restraintsInfo)
            if len(restraintsInfoProblems) > 0:
                simulationProblems.extend(restraintsInfoProblems)
        allSimulationProblems[stepName] = simulationProblems



    if len(allSimulationProblems) > 0:
        return  {"simulationInfo": allSimulationProblems}
    else:
        ## log that simulationInfo is correct
        drLogger.log_info(f"simulationInfo is correct...")
        return []
#################################################################################################


def check_restraintInfo(restraintInfo: dict) -> list:
    if not isinstance(restraintInfo, list):
        return ["restraintInfo must be a dictionary"]
    if len(restraintInfo) == 0:
        return ["restraintInfo must have at least one entry"]
    ## check each entry in restraintInfo
    restraintsInfoDisorders = []
    for info in restraintInfo:
        if not isinstance(info, dict):
            restraintsInfoDisorders.append("each entry in restraintInfo must be a dictionary (see README for syntax!)")
        else:
            ## check restraintType
            restraintType = info.get("restraintType", None)
            if restraintType is None:
                restraintsInfoDisorders.append("each entry in restraintInfo must have a 'restraintType' key")
            else:
                if not isinstance(restraintType, str):
                    restraintsInfoDisorders.append("each entry in restraintInfo['restraintType'] must be a string")
                if not restraintType in ["position", "distance", "angle", "torsion"]:
                    restraintsInfoDisorders.append("each entry in restraintInfo['restraintType'] must be one of 'position', 'distance', 'angle', 'torsion'")
            ## check selection for restraint to act upon
            restrantSelection = info.get("selection", None)
            if  restrantSelection is None:
                restraintsInfoDisorders.append("each entry in restraintInfo must have a 'selection' key")
            else:
                if not isinstance(restrantSelection, dict):
                    restraintsInfoDisorders.append("each entry in restraintInfo['selection'] must be a dictionary (see README for syntax!)")
                else:
                    selectionDisorder = check_selection({"selection": restrantSelection})
                    if len(selectionDisorder) > 0:
                        restraintsInfoDisorders.append(selectionDisorder)
            ## check parameters for restraint
            restraintParamers = info.get("parameters", None)
            if restraintParamers is None:
                restraintsInfoDisorders.append("each entry in restraintInfo must have a 'parameters' key")
            else:
                restraintParamProblems = check_restraint_parameters(restraintType, restraintParamers)
                if len(restraintParamProblems) > 0:
                    restraintsInfoDisorders.append(restraintParamProblems)

    return restraintsInfoDisorders

    


#####################################################################################
def check_postSimulationInfo(config: dict) -> None:
    """
    Checks for postSimulationInfo in config
    """
    ## log this check
    drLogger.log_info(f"Checking postSimulationInfo...")
    ## check for postSimulationInfo in config
    postSimulationInfo = config.get("postSimulationInfo", None)
    if postSimulationInfo is None:
        return []
    
    endPointInfo = postSimulationInfo.get("endPointInfo", None)
    if endPointInfo is not None:
        endPointDisorders = check_endpointInfo(endPointInfo)

    clusterInfo = postSimulationInfo.get("clusterInfo", None)
    if clusterInfo is not None:
        clusterDisorders = check_clusterInfo(clusterInfo)

    return endPointDisorders, clusterDisorders



#####################################################################################
def check_endpointInfo(endPointInfo: dict) -> None:
    """
    Checks for endPointInfo in config
    """
    ## log this check
    drLogger.log_info(f"Checking endPointInfo...")

    endPointDisorders = []
    ## check if endPointInfo is a dictionary
    if not isinstance(endPointInfo, dict):
        endPointDisorders.append("endPointInfo must be a dictionary")
        return {"endPointInfo": endPointDisorders}
    ## ensure that stepNames is present in endPointInfo
    stepNames = endPointInfo.get("stepNames", None)
    if stepNames is None:
        endPointDisorders.append("endPointInfo must have a 'stepNames' entry")
    else:
        if not isinstance(stepNames, list):
            endPointDisorders.append("endPointInfo['stepNames'] must be a list")
        ## ensure that stepNames is a list of strings
        if not all(isinstance(stepName, str) for stepName in stepNames):
            endPointDisorders.append("endPointInfo['stepNames'] must be a list of strings")
        ## ensure that stepNames is not empty
        if  len(stepNames) == 0:
            endPointDisorders.append("endPointInfo['stepNames'] must not be empty")

    removeAtoms = endPointInfo.get("removeAtoms", None)
    if removeAtoms is not None:
        if not isinstance(removeAtoms, list):
            endPointDisorders.append("endPointInfo['removeAtoms'] must be a list")
        if not all(isinstance(removeAtom, dict) for removeAtom in removeAtoms):
            endPointDisorders.append("endPointInfo['removeAtoms'] must be a list of dictionaries")
        for selection in removeAtoms:
            removeAtomSelectionDisorders = check_selection(selection)
            if len(removeAtomSelectionDisorders) > 0:
                endPointDisorders.extend(removeAtomSelectionDisorders)

    if len(endPointDisorders) > 0:  
        return {"endPointInfo": endPointDisorders}
    else:
        ## log that endPointInfo is correct
        drLogger.log_info(f"endPointInfo is correct...")
        return []
#####################################################################################
def check_clusterInfo(clusterInfo: dict) -> None:
    """
    Checks for clusterInfo in config
    """
    ## log this check
    drLogger.log_info(f"Checking clusterInfo...")
    clusterInfoDisorders = []
    ## check if clusterInfo is a dictionary
    if not isinstance(clusterInfo, dict):
        clusterInfoDisorders.append("clusterInfo must be a dictionary")
        return {"clusterInfo": clusterInfoDisorders}
    ## ensure mandetory args are present in clusterInfo
    stepNames = clusterInfo.get("stepNames", None)
    if stepNames is None:
        clusterInfoDisorders.append("clusterInfo must have a 'stepNames' entry")
    else:
        if not isinstance(stepNames, list):
            clusterInfoDisorders.append("clusterInfo['stepNames'] must be a list")
        ## ensure that stepNames is a list of strings
        if not all(isinstance(stepName, str) for stepName in stepNames):
            clusterInfoDisorders.append("clusterInfo['stepNames'] must be a list of strings")
        ## ensure that stepNames is not empty
        if  len(stepNames) == 0:
            clusterInfoDisorders.append("clusterInfo['stepNames'] must not be empty")

    clusterBy = clusterInfo.get("clusterBy", None)
    if clusterBy is None:
        clusterInfoDisorders.append("clusterInfo must have a 'clusterBy' entry")
    else:
        if not isinstance(clusterBy, list):
            clusterInfoDisorders.append("clusterInfo['clusterBy'] must be a list")
        if not all(isinstance(clusterSelection, dict) for clusterSelection in clusterBy):
            clusterInfoDisorders.append("clusterInfo['clusterBy'] must be a list of dictionaries")
        if len(clusterBy) == 0:
            clusterInfoDisorders.append("clusterInfo['clusterBy'] must not be empty")
        for clusterSelection in clusterBy:
            clusterSelectionDisorders = check_selection(clusterSelection)
            if len(clusterSelectionDisorders) > 0:
                clusterInfoDisorders.extend(clusterSelectionDisorders)

    removeAtoms = clusterInfo.get("removeAtoms", None)
    if removeAtoms is not None:
        if not isinstance(removeAtoms, list):
            clusterInfoDisorders.append("endPointInfo['removeAtoms'] must be a list")
        if not all(isinstance(removeAtom, dict) for removeAtom in removeAtoms):
            clusterInfoDisorders.append("endPointInfo['removeAtoms'] must be a list of dictionaries")
        for selection in removeAtoms:
            removeAtomSelectionDisorders = check_selection(selection)
            if len(removeAtomSelectionDisorders) > 0:
                clusterInfoDisorders.extend(removeAtomSelectionDisorders)

    if len(clusterInfoDisorders) > 0:
        return {"clusterInfo": clusterInfoDisorders}
    else:
        ## log that clusterInfo is correct
        drLogger.log_info(f"clusterInfo is correct...")
        return []


#########################################################################
def check_ligandInfo(config: dict) -> None:
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
    ligandInfo,  = check_info_for_args(config, "config", ["ligandInfo"], optional=True)
    inputDir = config["pathInfo"]["inputDir"]

    if not ligandInfo:
        return

    # Check each entry in ligandInfo
    for ligand in ligandInfo:
        # Check if ligand is a dictionary
        if not isinstance(ligand, dict):
            raise TypeError(f"ligandInfo must be a list of dicts")
        # Check if ligand has at least one entry
        if len(ligand) == 0:
            raise ValueError(f"ligandInfo must have at least one entry")
        
        # Check each argument in the ligand dictionary
        ligandName, protons, charge, toppar, mol2  = check_info_for_args(
            ligand, "ligand", ["ligandName", "protons", "charge", "toppar", "mol2"], optional=True
        )
        
        # Check if ligandName is a string
        if not isinstance(ligandName, str):
            raise TypeError(f"ligandName must be a string")
        
        # Check if protons, charge, toppar, and mol2 are booleans
        for arg_value, arg_name in zip([protons, toppar, mol2], ["protons", "toppar", "mol2"]):
            if not isinstance(arg_value, bool):
                raise TypeError(f"{arg_name} must be a boolean")

        # Check if charge is an int
        if not isinstance(charge, int):
            raise TypeError(f"charge must be an int")
        
        # Check if ligand is in input pdb file
        inputPdbs = [file for file in os.listdir(inputDir) if file.endswith(f".pdb")]
        for inputPdb in inputPdbs:
            inputDf = pdbUtils.pdb2df(p.join(inputDir, inputPdb))
            if not ligandName in inputDf["RES_NAME"].values:
                raise ValueError(f"Ligand name \"{ligandName}\" not found in input file: {inputPdb}")
    ## log that ligandInfo is correct
    drLogger.log_info(f"ligandInfo is correct...")




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
def check_metadynamics_options(simulation: dict, stepName: str) -> None:
    drLogger.log_info(f"Checking metaDynamicsInfo for {stepName}...")
    ## check metaDynamicsInfo
    metaProblems = []
    metaDynamicsInfo = simulation.get("metaDynamicsInfo", None)
    if metaDynamicsInfo is None:
        return ["No metaDynamicsInfo specified in simulation"]
    if not isinstance(metaDynamicsInfo, dict):
        return ["metaDynamicsInfo must be a list"]
    
    ## check height parameter
    height = metaDynamicsInfo.get("height", None)
    if height is not None:
        metaProblems.append("No height specified in metaDynamicsInfo")
    else:
        if not isinstance(height, (int, float)):
            metaProblems.append("height must be a number")
        elif height <= 0:
            metaProblems.append("height must be positive")
    
    ## check biasFactor parameter
    biasFactor = metaDynamicsInfo.get("biasFactor", None)
    if biasFactor is None:
        metaProblems.append("No biasFactor specified in metaDynamicsInfo")
    else:
        if not isinstance(biasFactor, (int, float)):
            metaProblems.append("biasFactor must be a number")
        elif biasFactor <= 0:
            metaProblems.append("biasFactor must be positive")
    
    ## check biases parameter
    biases = metaDynamicsInfo.get("biases", None)
    if biases is None:
        metaProblems.append("No biases specified in metaDynamicsInfo")
    else:
        if not isinstance(biases, list):
            metaProblems.append("biases must be a list of biases (check README for more details)")
        elif len(biases) == 0:
            metaProblems.append("biases must contain at least one bias variable")
        else:
            for bias in biases:
                if not isinstance(bias, dict):
                    metaProblems.append("biases must be a list of biases (check README for more details)")
                biasVar = bias.get("biasVar", None)
                if biasVar is None:
                    metaProblems.append("No biasVar specified in bias")
                else:
                    if not isinstance(biasVar, str):
                        metaProblems.append("biasVar must be 'rmsd', 'torsion', 'distance', or 'angle'")
                    elif not biasVar.lower() in ["rmsd", "torsion", "distance", "angle"]:
                        metaProblems.append("biasVar must be 'rmsd', 'torsion', 'distance', or 'angle'")
                biasSelection = bias.get("selection", None)
                if biasSelection is None:
                    metaProblems.append("No selection specified in bias")
                else:
                    if not isinstance(biasSelection, dict):
                        metaProblems.append("selection must be a dictionary")
                    else:
                        selectionDisorders = check_selection(biasSelection, stepName)
                        if len(selectionDisorders) > 0:
                            metaProblems.extend(selectionDisorders)
                        
    if len(metaProblems) > 0:
        return metaProblems
    else:
        drLogger.log_info(f"MetaDynamicsInfo is correct...")
        return []

#########################################################################
def check_nvt_npt_options(simulation: dict, stepName: str) -> list:

    optionsProblems = []
    ## check for required args for a nvt or npt simulation
    duration = simulation.get("duration", None)
    if duration is None:
        optionsProblems.append("No duration specified in simulation")
    else:
        timeCheckProblems = check_time_input(duration, "duration", stepName)
        if timeCheckProblems is not None:
            optionsProblems.append(timeCheckProblems)
        
    logInterval = simulation.get("logInterval", None)
    if logInterval is None:
        optionsProblems.append("No logInterval specified in simulation")
    else:
        timeCheckProblems = check_time_input(logInterval, "logInterval", stepName)
        if timeCheckProblems is not None:
            optionsProblems.append(timeCheckProblems)


    heavyProtons = simulation.get("heavyProtons", None)
    if heavyProtons is None:
        optionsProblems.append("No heavyProtons specified in simulation")
    else:
        if not isinstance(heavyProtons, bool):
            optionsProblems.append("heavyProtons must be a boolean")
        
    return optionsProblems

        

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
def check_shared_simulation_options(simulation) -> Tuple[list, str, str]:
    sharedOptionProblems = []
    ## check simulation step name
    stepName = simulation.get("stepName", None)
    if stepName is None: 
        sharedOptionProblems.append("No stepName specified in simulation")
    elif not isinstance(stepName, str):
        sharedOptionProblems.append("stepName must be a string")
    elif " " in stepName:
        sharedOptionProblems.append("No whitespace allowed in stepName")


    ## check simulationType
    simulationType = simulation.get("simulationType", None)
    if simulationType is None:
        sharedOptionProblems.append("No simulationType specified in simulation")
    elif not simulationType in ["EM", "NVT", "NPT", "META"]:
        sharedOptionProblems.append(f"simulationType in simulation must be one of the following:\nEM, NVT, NPT, META")

    ## check for either temperature or temparatureRange in simulation
    temperature = simulation.get("temperature", None)
    tempRange = simulation.get("temperatureRange", None)

    if temperature and tempRange:
        sharedOptionProblems.append("Cannot specify both temperature and temperatureRange in simulation")
    elif not temperature and not tempRange:
        sharedOptionProblems.append("Must specify either temperature or temperatureRange in simulation")
    elif temperature:
        if not isinstance(temperature, int):
            sharedOptionProblems.append("Temperature in simulation must be an int")
        if temperature < 0:
            sharedOptionProblems.append("Temperature in simulation must be a positive int")
    elif tempRange:
        if not isinstance(tempRange, list):
            sharedOptionProblems.append("TemperatureRange in simulation must be a list of ints")
        if len(tempRange) == 0:
            sharedOptionProblems.append("TemperatureRange in simulation must be a list of at least one int")
        for temp in tempRange:
            if not isinstance(temp, int):
                sharedOptionProblems.append("TemperatureRange in simulation must be a list of ints")
            if temp < 0:
                sharedOptionProblems.append("TemperatureRange in simulation must be a list of positive ints")

    return sharedOptionProblems, stepName, simulationType

#########################################################################
def validate_path(argName: str, argPath: Union[PathLike, str]) -> str:
    """
    Check to see if a path variable is indeed the correct type
    Check to see if the path exists
    """
    if  not isinstance(argPath, (PathLike, str)) :
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
