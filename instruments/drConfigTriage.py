## basic libraries
import argpass
import yaml
import os
from os import path as p
import re
import multiprocessing as mp

## clean code
from typing import Tuple, Union
from os import PathLike
from instruments.drCustomClasses import FilePath, DirectoryPath
## drMD modules
from instruments import drSplash
from instruments import drLogger
## custom modules
from pdbUtils import pdbUtils
#####################################################################################
def read_and_validate_config() -> Tuple[dict, FilePath]:
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
    drLogger.log_info("-->\tChecking config file...",True)

    ## try to read config yaml file into a dict
    ## throw errors with a nice splash screen if something goes wrong
    try:
        config: dict = read_input_yaml()
    except FileNotFoundError as e:
        drSplash.print_config_error(e) 
    except yaml.YAMLError as e:
        drSplash.print_config_error(e) 
    except KeyError as e:
        drSplash.print_config_error(e) 
    except TypeError as e:
        drSplash.print_config_error(e) 
    except ValueError as e:
        drSplash.print_config_error(e) 
    
    ## check each major section in config file
    ## throw errors with a nice splash screen if something goes wrong
    for function in [check_pathInfo,
                      check_generalInfo,
                        check_ligandInfo,
                          check_simulationInfo,
                          check_postSimulationInfo]:
        try:
            function(config)
        except FileNotFoundError as e:
            drSplash.print_config_error(e) 
        except yaml.YAMLError as e:
            drSplash.print_config_error(e) 
        except KeyError as e:
            drSplash.print_config_error(e) 
        except TypeError as e:
            drSplash.print_config_error(e) 
        except ValueError as e:
            drSplash.print_config_error(e) 

    drLogger.log_info("-->\tConfig file is correct", True)
    drLogger.close_logging()
    return config, configTriageLog
#####################################################################################
def check_postSimulationInfo(config: dict) -> None:
    """
    Checks for postSimulationInfo in config
    """
    ## log this check
    drLogger.log_info("-->\tChecking postSimulationInfo...")
    ## check for postSimulationInfo in config
    postSimulationInfo, = check_info_for_args(config, "config", ["postSimulationInfo"], optional=True)
    if not postSimulationInfo:
        drLogger.log_info("-->\tOptional postSimulationInfo is not present")
        return
    optionalArgs = ["endPointInfo", "clusterInfo"]
    endPointInfo, clusterInfo = check_info_for_args(postSimulationInfo, "postSimulationInfo", optionalArgs, optional=True) 
    
    if endPointInfo:
        check_endpointInfo(endPointInfo)
    if clusterInfo:
        check_clusterInfo(clusterInfo)
#####################################################################################
def check_endpointInfo(endPointInfo: dict) -> None:
    """
    Checks for endPointInfo in config
    """
    ## log this check
    drLogger.log_info("-->\tChecking endPointInfo...")
    ## check if endPointInfo is a dictionary
    if not isinstance(endPointInfo, dict):
        raise TypeError("-->\tendPointInfo must be a dict")
    ## ensure that stepNames is present in endPointInfo
    mandatoryArgs = ["stepNames", "collate"]
    stepNames, collate = check_info_for_args(endPointInfo, "endPointInfo", mandatoryArgs, optional=False)
    if not isinstance(stepNames, list):
        raise TypeError("-->\tstepNames must be a list")
    ## ensure that stepNames is a list of strings
    if not all(isinstance(stepName, str) for stepName in stepNames):
        raise TypeError("-->\tstepNames must be a list of strings")
    ## ensure that stepNames is not empty
    if  len(stepNames) == 0:
        raise ValueError("-->\tstepNames must not be empty")
    ## check for optional removeAtoms arg in endPointInfo
    removeAtoms = check_info_for_args(endPointInfo, "endPointInfo", ["removeAtoms"], optional=True)


    if removeAtoms:
        ## check each selection in remveAtoms
        removeAtomsSelections = [sele for sele in removeAtoms[0]]
        if len(removeAtomsSelections) == 0:
            raise ValueError("-->\tremoveAtoms must not be empty")
        for removeAtomsSelection in removeAtomsSelections:
            check_selection(removeAtomsSelection["selection"], "removeAtoms")

    ## ensure collate is a bool
    if not isinstance(collate, bool):
        raise TypeError("-->\tcollate must be a bool")
    ## log that endPointInfo is correct
    drLogger.log_info("-->\tendPointInfo is correct...")
#####################################################################################
def check_clusterInfo(clusterInfo: dict) -> None:
    """
    Checks for clusterInfo in config
    """
    ## log this check
    drLogger.log_info("-->\tChecking clusterInfo...")
    ## check if clusterInfo is a dictionary
    if not isinstance(clusterInfo, dict):
        raise TypeError("-->\tclusterInfo must be a dict")
    ## ensure mandetory args are present in clusterInfo
    mandetoryArgs = ["stepNames", "clusterBy", "collate"]
    stepNames, clusterBy, collate = check_info_for_args(clusterInfo, "clusterInfo", mandetoryArgs, optional=False)
    if not isinstance(stepNames, list):
        raise TypeError("-->\tstepNames must be a list")
    ## ensure that stepNames is a list of strings
    if not all(isinstance(stepName, str) for stepName in stepNames):
        raise TypeError("-->\tstepNames must be a list of strings")
    ## ensure that stepNames is not empty
    if  len(stepNames) == 0:
        raise ValueError("-->\tstepNames must not be empty")
    ## ensure that clusterBy is a dict
    if not isinstance(clusterBy, dict):
        raise TypeError("-->\tclusterBy must be a string")
    clusterSelection = check_info_for_args(clusterBy, "clusterBy", ["selection"], optional=False)[0]
    check_selection(clusterSelection, "clusterBy")
    ## ensure collate is a bool
    if not isinstance(collate, bool):
        raise TypeError("-->\tcollate must be a bool")
    ## check for optional removeAtoms arg in clusterInfo
    removeAtoms = check_info_for_args(clusterInfo, "clusterInfo", ["removeAtoms"], optional=True)
    ## check each selection in removeAtoms
    if removeAtoms:
        for removeAtomsSelection in removeAtoms[0]:
            check_selection(removeAtomsSelection["selection"], "removeAtoms")

    ## log that clusterInfo is correct
    drLogger.log_info("-->\tclusterInfo is correct...")


    


#####################################################################################
def check_pathInfo(config: dict) -> None:
    """
    Checks for pathInfo entry in config
    Checks paths in pathInfo to see if they are real
    Don't check outputDir, this will be made automatically
    """
    ## log this check
    drLogger.log_info("-->\tChecking pathInfo...")
    ## check if pathInfo in config
    pathInfo, = check_info_for_args(config, "config", ["pathInfo"], optional=False)
    ## check for required args in pathInfo
    inputDir, outputDir = check_info_for_args(pathInfo, "pathInfo", ["inputDir", "outputDir"], optional=False)
    validate_path("inputDir", inputDir)
    ## log that pathInfo is correct
    drLogger.log_info("-->\tpathInfo is correct...")

#########################################################################
def check_generalInfo(config: dict) -> None:
    """
    Checks generalInfo in config
    Makes sure that CPU allocations are properly formatted
    Makes sure that the "Platform" specified is an allowed value 
    """
    ## log this check
    drLogger.log_info("-->\tChecking generalInfo...")
    ## check if generalInfo in config
    generalInfo, = check_info_for_args(config, "config", ["generalInfo"], optional= False)
    ## check for required args in generalInfo
    parallelCPU, platform, subprocessCpus = check_info_for_args(generalInfo, "generalInfo", ["parallelCPU", "platform", "subprocessCpus"], optional= False)
    ## check that generalInfo cpu arguments are int values and are positive
    for argValue, argName in zip([parallelCPU, subprocessCpus], ["parallelCPU", "subprocessCpus"]):
        if not isinstance(argValue, int):
            raise TypeError(f"-->\tThe config argument {argName} = {argValue} is not a an int type.")
        if argValue < 1:
            raise ValueError(f"-->\tThe config argument {argName} = {argValue} must be a int greater than 1")
    ## check that your computer has enough CPUs
    if parallelCPU * subprocessCpus > mp.cpu_count():
        raise ValueError("-->\ttotalCpuUseage argument exceeds your computers number of cores")
    ## ensure that compatable platform has been chosen 
    if not platform in ["CUDA", "OpenCL", "CPU"]:
        raise ValueError("-->\tPlatform must be CUDA, OpenCL, or CPU")
    ## log that generalInfo is correct
    drLogger.log_info("-->\tgeneralInfo is correct...")
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
    drLogger.log_info("-->\tChecking ligandInfo...")
    # Check if ligandInfo in config
    ligandInfo,  = check_info_for_args(config, "config", ["ligandInfo"], optional=True)
    inputDir = config["pathInfo"]["inputDir"]

    if not ligandInfo:
        return

    # Check each entry in ligandInfo
    for ligand in ligandInfo:
        # Check if ligand is a dictionary
        if not isinstance(ligand, dict):
            raise TypeError("-->\tligandInfo must be a list of dicts")
        # Check if ligand has at least one entry
        if len(ligand) == 0:
            raise ValueError("-->\tligandInfo must have at least one entry")
        
        # Check each argument in the ligand dictionary
        ligandName, protons, charge, toppar, mol2  = check_info_for_args(
            ligand, "ligand", ["ligandName", "protons", "charge", "toppar", "mol2"], optional=True
        )
        
        # Check if ligandName is a string
        if not isinstance(ligandName, str):
            raise TypeError("-->\tligandName must be a string")
        
        # Check if protons, charge, toppar, and mol2 are booleans
        for arg_value, arg_name in zip([protons, toppar, mol2], ["protons", "toppar", "mol2"]):
            if not isinstance(arg_value, bool):
                raise TypeError(f"-->\t{arg_name} must be a boolean")

        # Check if charge is an int
        if not isinstance(charge, int):
            raise TypeError("-->\tcharge must be an int")
        
        # Check if ligand is in input pdb file
        inputPdbs = [file for file in os.listdir(inputDir) if file.endswith(".pdb")]
        for inputPdb in inputPdbs:
            inputDf = pdbUtils.pdb2df(p.join(inputDir, inputPdb))
            if not ligandName in inputDf["RES_NAME"].values:
                raise ValueError(f"-->\tLigand name \"{ligandName}\" not found in input file: {inputPdb}")
    ## log that ligandInfo is correct
    drLogger.log_info("-->\tligandInfo is correct...")


#########################################################################
def check_simulationInfo(config: dict) -> None:
    """
    Checks for simulationInfo in config
    Depending on the type of simulation, checks your parameters
    """
    ## log this check
    drLogger.log_info("-->\tChecking simulationInfo...")
    ## check for simulationInfo in config
    simulationInfo,  = check_info_for_args(config, "config", ["simulationInfo"], optional=False)
    ## ensure that simulationInfo is a non-zero-length list
    if not isinstance(simulationInfo, list):
        raise TypeError("-->\tsimulationInfo must be a list containing dicts")
    if len(simulationInfo) == 0:
        raise ValueError("-->\tsimulationInfo must have at least one entry")
    ## look through each entry in simulationInfo
    for simulation in simulationInfo:
        ## log this check of this simulation step
        drLogger.log_info(f"-->\tChecking {simulation['stepName']}...")
        simulationType, stepName = check_shared_simulation_options(simulation)
        drLogger.log_info(f"\t\tsimulationType: {simulationType}")
        if simulationType.upper() in ["NVT", "NPT"]:
            check_nvt_npt_options(simulation, stepName)

        elif simulationType.upper() == "META":
            check_metadynamics_options(simulation, stepName)

        restraintInfo, = check_info_for_args(simulation, stepName, ["restraintInfo"], optional=True)
        if restraintInfo:
            drLogger.log_info("\t\twith restraints")
            check_restraints_options(restraintInfo, stepName)

        ## check for illigal args
        alowedArgsForSimulations = ["stepName", "simulationType", "duration", "maxIterations",
                                     "timestep", "temp", "logInterval", "restraintInfo",
                                       "metaDynamicsInfo", "clusterTrajectory"]
        check_info_for_illigal_args(simulation, stepName, alowedArgsForSimulations)

    ## log that simulationInfo is correct
    drLogger.log_info("-->\tsimulationInfo is correct...")

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
        raise TypeError("-->\tclusterTrajectory in {stepName} must be a dict")
    
    # Check if clusterTrajectory has at least one entry
    if len(clusterTrajectory) == 0:
        raise ValueError("-->\tclusterTrajectory in {stepName} must contain some entries")
    
    # Check for required args for clustering a trajectory
    nClusters, selection = check_info_for_args(clusterTrajectory, stepName, ["nClusters", "selection"], optional=False)
    
    # Check if nClusters is an integer
    if not isinstance(nClusters, int):
        raise TypeError(f"-->\tnClusters in {stepName} must be an int")
    
    # Check selection for clustering a trajectory
    check_selection(selection, stepName)

#########################################################################
def check_restraints_options(restraintInfo: dict, stepName: str) -> None:
    """
    Checks for options in restraintInfo.

    Args:
        restraintInfo (dict): A dictionary containing options for restraints.
        stepName (str): The name of the step.

    Raises:
        TypeError: If restraintInfo is not a list.
        ValueError: If restraintInfo is empty.
        TypeError: If restraintType is not one of the following: POSITION, TORSION, DISTANCE, ANGLE.
        ValueError: If restraintInfo does not contain any restraints.
    """
    # Check if restraintInfo is a list
    if not isinstance(restraintInfo, list):
        raise TypeError(f"-->\trestraintsInfo in {stepName} must be a list of restraints")
    
    # Check if restraintInfo has at least one entry
    if len(restraintInfo) == 0:
        raise ValueError(f"-->\trestraintsInfo is must contain some restraints")
    
    # Loop through each restraint in restraintInfo
    for restraint in restraintInfo:
        # Check for required args for restraints
        restraintType, selection, parameters = check_info_for_args(restraint, "restraint", ["restraintType", "selection", "parameters"], optional=False)
        
        # Check if restraintType is one of the following: POSITION, TORSION, DISTANCE, ANGLE
        if not restraintType.upper() in ["POSITION", "TORSION", "DISTANCE", "ANGLE"]:
            raise ValueError("-->\trestraintType must be one of the following:\n-->\tPOSITION, DISTANCE, ANGLE, TORSION")
        
        # Check selection for restraints
        check_selection(selection, stepName)
        
        # Check restraint parameters
        check_restraint_parameters(restraintType, parameters)
#########################################################################
def check_restraint_parameters(restraintType: str, parameters: dict) -> None:
    k, = check_info_for_args(parameters, "parameters", ["k"], optional=False)
    if not isinstance(k, (int, float)):
        raise TypeError(f"-->\k force constants parameters must be a number for all resraints")
    if k <= 0:
        raise ValueError(f"-->\tk force constants  parameters must be positive for all restraints")

    if restraintType.upper() == "TORSION":
        phi, = check_info_for_args(parameters, "parameters", ["phi"], optional=False)
        if not isinstance(phi, (int, float)):
            raise TypeError(f"-->\phi0 parameters must be a number for torsion restraints")
        if phi < -180 or phi > 180:
            raise ValueError(f"-->\tphi0 parameters must be between -180 and 180 for torsion restraints")
        
    elif restraintType.upper() == "DISTANCE":
        r0, = check_info_for_args(parameters, "parameters", ["r0"], optional=False)
        if not isinstance(r0, (int, float)):
            raise TypeError(f"-->\r0  parameters must be a number for distance restraints")
        if r0 < 0:
            raise ValueError(f"-->\r0 parameters must be positive for distance restraints")

    elif restraintType.upper() == "ANGLE":
        theta0, = check_info_for_args(parameters, "parameters", ["theta0"], optional=False)
        if not isinstance(theta0, (int, float)):
            raise TypeError(f"-->\theta0  parameters must be a number for angle restraints")
        if theta0 < 0 or theta0 > 180:
            raise ValueError(f"-->\theta0 parameters must be between 0 and 180 for angle restraints")

#########################################################################
def check_selection(selection: dict, stepName: str) -> None:
    keyword, = check_info_for_args(selection, "selection", ["keyword"], optional=False)
    if not isinstance(keyword, str):
        raise TypeError(f"-->\tselection keywords incorrect in {stepName}, see README.md for more details")
    if not keyword in ["all", "protein", "ligands", "water", "ions", "custom"]:
        raise ValueError(f"-->\tselection keywords incorrect in {stepName}, see README.md for more details")

    if keyword == "custom":
        customSelection, = check_info_for_args(selection, stepName, ["customSelection"], optional=False)
        if not isinstance(customSelection, list):
            raise TypeError(f"-->\tselectionSyntax in {stepName} must be a list of dcits")
        for selection in customSelection:
            chainId, residueName, residueNumber, atomName = check_info_for_args(selection, stepName, ["CHAIN_ID", "RES_NAME", "RES_ID", "ATOM_NAME"], optional=False)
            if not isinstance(chainId, (str,list)):
                raise TypeError(f"""-->\tError in custom selection syntax in {stepName}:
-->\tCHAIN_ID must be a string, a list of strings or the keyword 'all'
                                """)
            if not isinstance(residueName, (str,list)):
                raise TypeError(f"""-->\tError in custom selection syntax in {stepName}:
-->\tRES_NAME must be a string, a list of strings or the keyword 'all'
                                """)
            if not isinstance(residueNumber, (int,list,str)):
                raise TypeError(f"""-->\tError in custom selection syntax in {stepName}:
-->\tRES_ID must be an int or a list of ints or the keyword 'all'
                                """)
            if not isinstance(atomName, (str,list)):
                raise TypeError(f"""-->\tError in custom selection syntax in {stepName}:
-->\tATOM_NAME must be a string, a list of strings or the keyword 'all'
                                """)

#########################################################################
def check_metadynamics_options(simulation: dict, stepName: str) -> None:
    manditoryArgs = ["height", "biasFactor", "biases"]
    height, biasFactor, biases = check_info_for_args(simulation, stepName, manditoryArgs, optional=False)
    if not isinstance(height, (int, float)):
        raise TypeError(f"-->\tMetaDynamics height must be a number")
    if not isinstance(biasFactor, (int, float)):
        raise TypeError(f"-->\tMetaDynamics biasFactor must be a number")
    if not isinstance(biases, list):
        raise TypeError(f"-->\tMetaDynamics biases must be a bias variables")
    if len(biases) == 0:
        raise ValueError(f"-->\tMetaDynamics biases must contain at least one bias variable")
    
    for bias in biases:
        manditoryArgs = ["biasVar", "selection"]
        biasVar, selection = check_info_for_args(bias, stepName, manditoryArgs, optional=False)
        if not biasVar.lower() in ["rmsd", "torsion", "distance", "angle"]:
            raise ValueError(f"-->\tMetaDynamics biasVar must be one of the following:\n-->\tRMSD, TORSION, DISTANCE, ANGLE")
        check_selection(selection, stepName)


#########################################################################
def check_nvt_npt_options(simulation: dict, stepName: str) -> None:
    ## check for required args for a nvt or npt simulation
    duration, logInterval = check_info_for_args(simulation, stepName, ["duration", "logInterval"], optional= False)
    ## ensure that duration and log interval are correctly formatted
    for timeInputValue, timeInputName in zip([duration, logInterval],["duration", "logInterval"]):
        check_time_input(timeInputValue, timeInputName, stepName)

        

#########################################################################
def check_time_input(timeInputValue: str, timeInputName: str, stepName: str) -> None:
    if not isinstance(timeInputValue, str):
        raise TypeError(f"-->\t{timeInputName} in {stepName} must be in the format \"10 ns\"")
    timeInputData = timeInputValue.split()
    try:
        numberAsFloat = int(timeInputData[0])
    except:
        raise ValueError(f"-->\t{timeInputName} in {stepName} must be in the format \"10 ns\"")
    if not timeInputData[1] in ["fs","ps","ns","ms"]:
        raise ValueError(f"-->\t{timeInputName} in {stepName} must be in the format \"10 ns\"")

#########################################################################
def check_shared_simulation_options(simulation) -> str:
    ## check for shared required args in simulation
    stepName, simulationType, temp = check_info_for_args(simulation, "simulation", ["stepName", "simulationType", "temp"], optional=False)
    ## make sure stepName is correct
    if not isinstance(stepName, str):
        raise TypeError(f"-->\tstepName {str(stepName)} must be a a string")
    if " " in stepName:
        raise ValueError(f"-->\tNo whitespace allowed in stepName {stepName}")
    ## make sure simulationType is an allowed value
    if not simulationType.upper() in ["EM", "NVT", "NPT", "META"]:
        raise ValueError(f"-->\tsimulationType in {stepName} must be the following:\n-->\tEM, NVT, NPT, META")
    ## make sure temp is an int
    if not isinstance(temp, int):
        raise TypeError(f"-->\tTemperature in {stepName} must be an int")
    if temp < 0:
        raise ValueError(f"-->\tTemperature in {stepName} must be a positive int")

    return simulationType, stepName

#########################################################################
def check_info_for_args(info: dict, infoName: str,  argNames: list, optional: bool) -> list:
    """
    Simple check to see if list of keys is in a dict
    If optional is set to "False", raises a KeyError
    If all is as expected returns a list of values to be unpacked
    """
    ## init empty list to append to 
    unpackedDicts: list = []
    for  argName in argNames:
        isArg, argValue = check_dict_for_key(info, argName)
        ## if a required arg is not in the dict, raise a KeyError
        if not isArg and not optional:
            raise KeyError(f"-->\tArgument {argName} not found in {infoName}")
        unpackedDicts.append(argValue)
    return unpackedDicts

#########################################################################
def check_info_for_illigal_args(info: dict, infoName: str,  allowedKeys: list) -> None:
    """
    Simple check to see if list of keys is in a dict
    If optional is set to "False", raises a KeyError
    If all is as expected returns a list of values to be unpacked
    """
    keysInInfo = list(info.keys())
    illegalKeys = list(set(keysInInfo) - set(allowedKeys))
    if len(illegalKeys) > 0:
        raise ValueError(f"-->\t{infoName} contains illegal keys: {illegalKeys}")
#########################################################################
def check_dict_for_key(info: dict, key: any) -> Tuple[bool, any]:
    """
    Checks to see if a key is in a dict
    If it is, return "True" and the associated value 
    """
    if key in info:
        if info[key] == False:
            return True, False
        else:
            return True, info[key]
    return False, False
#########################################################################
def validate_path(argName: str, argPath: Union[PathLike, str]) -> None:
    """
    Check to see if a path variable is indeed the correct type
    Check to see if the path exists
    """
    if  not isinstance(argPath, (PathLike, str)) :
        raise TypeError(f"-->\tThe config argument {argName} = {argPath} is not a PathLike.")
    # Check if the path exists
    if not p.exists(argPath):
        raise FileNotFoundError(f"-->\tThe config argument {argName} = {argPath} does not exist.")
#####################################################################################
def read_input_yaml() -> dict:
    """
    Reads a YAML file using the "--config" flag with argpass
    Reads YAML file into a dict

    Returns:
    - config (dict)
    """
    # create an argpass parser, read config file,
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    configFile: Union[PathLike, str] = args.config
    # Read config.yaml into a dictionary
    try:
        with open(configFile, "r") as yamlFile:
            config: dict = yaml.safe_load(yamlFile)
            return config
    except FileNotFoundError:
        raise FileNotFoundError(f"-->\tconfig file {configFile} not found")
    except yaml.YAMLError as exc:
        raise yaml.YAMLError("-->\tError parsing YAML file:", exc)



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
