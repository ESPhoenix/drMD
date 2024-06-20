## basic libraries
import argpass
import yaml
from os import path as p
import re
## clean code
from typing import Tuple, Union
from os import PathLike
import multiprocessing as mp
#####################################################################################
def read_and_validate_config() -> dict:
    """
    Main script for drConfigInspector:
    1. Accepts --config flag arg with argpass, reads yaml file into a dict
    2. Checks paths in pathInfo to see if they are real
    3. Checks inputs of cpuInfo to see if they are correct
    4. Checks each dockingOrder in dockingOrders to see if they are correct

    Returns:
    - config (dict)
    """

    print("checking config file...")
    config: dict = read_input_yaml()
    check_pathInfo(config)
    check_generalInfo(config)
    check_ligandInfo(config)
    check_simulationInfo(config)
    print("... config file is correct")
    return config
#####################################################################################
def check_pathInfo(config: dict) -> None:
    """
    Checks for pathInfo entry in config
    Checks paths in pathInfo to see if they are real
    Don't check outputDir, this will be made automatically
    """
    ## check if pathInfo in config
    pathInfo, = check_info_for_args(config, "config", ["pathInfo"], optional=False)
    ## check for required args in pathInfo
    inputDir, outputDir = check_info_for_args(pathInfo, "pathInfo", ["inputDir", "outputDir"], optional=False)
    validate_path("inputDir", inputDir)

#########################################################################
def check_generalInfo(config: dict) -> None:
    """
    Checks generalInfo in config
    Makes sure that CPU allocations are properly formatted
    Makes sure that the "Platform" specified is an allowed value 
    """
    ## check if generalInfo in config
    generalInfo, = check_info_for_args(config, "config", ["generalInfo"], optional= False)
    ## check for required args in generalInfo
    parallelCPU, platform, subprocessCpus = check_info_for_args(generalInfo, "generalInfo", ["parallelCPU", "platform", "subprocessCpus"], optional= False)
    ## check that generalInfo cpu arguments are int values and are positive
    for argValue, argName in zip([parallelCPU, subprocessCpus], ["parallelCPU", "subprocessCpus"]):
        if not isinstance(argValue, int):
            raise TypeError(f"The config argument {argName} = {argValue} is not a an int type.")
        if argValue < 1:
            raise ValueError(f"The config argument {argName} = {argValue} must be a int greater than 1")
    ## check that your computer has enough CPUs
    if parallelCPU * subprocessCpus > mp.cpu_count():
        raise ValueError("totalCpuUseage argument exceeds your computers number of cores")
    ## ensure that compatable platform has been chosen 
    if not platform in ["CUDA", "OpenCL", "CPU"]:
        raise ValueError("Platform must be CUDA, OpenCL, or CPU")
#########################################################################
def check_ligandInfo(config: dict) -> None:
    ## TODO: come back to this!
    """
    Checks optional ligandInfo entry in config
    """
    ligandInfo,  = check_info_for_args(config, "config", ["ligandInfo"], optional=True)
    if not ligandInfo:
        return
    ##TODO: redo this ligandInfo - remove nLigands etc
#########################################################################
def check_simulationInfo(config: dict) -> None:
    """
    Checks for simulationInfo in config
    Depending on the type of simulation, checks your parameters
    """
    ## check for simulationInfo in config
    simulationInfo,  = check_info_for_args(config, "config", ["simulationInfo"], optional=False)
    ## ensure that simulationInfo is a non-zero-length list
    if not isinstance(simulationInfo, list):
        raise TypeError("simulationInfo must be a list containing dicts")
    if len(simulationInfo) == 0:
        raise ValueError("simulationInfo must have at least one entry")
    ## look through each entry in simulationInfo
    for simulation in simulationInfo:
        simulationType, stepName = check_shared_simulation_options(simulation)
        if simulationType.upper() in ["NVT", "NPT"]:
            check_nvt_npt_options(simulation, stepName)

        elif simulationType.upper() == "META":
            check_meta_options(simulation, stepName)

        restraintInfo, = check_info_for_args(simulation, stepName, ["restraintInfo"], optional=True)
        if restraintInfo:
            check_restraints_options(restraintInfo, stepName)
#########################################################################
def check_restraints_options(restraintInfo: dict, stepName: str) -> None:
    if not isinstance(restraintInfo, list):
        raise TypeError(f"restraintsInfo in {stepName} must be a list of restraints")
    if len(restraintInfo) == 0:
        raise ValueError(f"restraintsInfo is must contain some restraints")
    for restraint in restraintInfo:
        print(restraint)
        restraintType, selection, parameters = check_info_for_args(restraint, "restraint", ["restraintType", "selection", "parameters"], optional=False)
        if not restraintType.upper() in ["POSITION", "TORSION", "DISTANCE", "ANGLE"]:
            raise ValueError("restraintType must be one of the following:\n POSITION, DISTANCE, ANGLE, TORSION")
        check_selection(selection, stepName)
        check_restraint_parameters(restraintType, parameters)
#########################################################################
def check_restraint_parameters(restraintType: str, parameters: dict) -> None:
    k, = check_info_for_args(parameters, "parameters", ["k"], optional=False)
    if not isinstance(k, (int, float)):
        raise TypeError(f"restraint parameters must be a number")
    if k <= 0:
        raise ValueError(f"restraint parameters must be positive")

    if restraintType.upper() == "TORSION":
        phi, = check_info_for_args(parameters, "parameters", ["phi"], optional=False)
        if not isinstance(phi, (int, float)):
            raise TypeError(f"restraint parameters must be a number")
        if phi < -180 or phi > 180:
            raise ValueError(f"restraint parameters must be between -180 and 180")
        
    elif restraintType.upper() == "DISTANCE":
        r0, = check_info_for_args(parameters, "parameters", ["r0"], optional=False)
        if not isinstance(r0, (int, float)):
            raise TypeError(f"restraint parameters must be a number")
        if r0 < 0:
            raise ValueError(f"restraint parameters must be positive")

    elif restraintType.upper() == "ANGLE":
        theta0, = check_info_for_args(parameters, "parameters", ["theta0"], optional=False)
        if not isinstance(theta0, (int, float)):
            raise TypeError(f"restraint parameters must be a number")
        if theta0 < 0 or theta0 > 180:
            raise ValueError(f"restraint parameters must be between 0 and 180")

#########################################################################
def check_selection(selection: dict, stepName: str) -> None:
    keyword, = check_info_for_args(selection, "selection", ["keyword"], optional=False)
    if not isinstance(keyword, str):
        raise TypeError(f"selection keywords incorrect in {stepName}, see README.md for more details")
    if not keyword in ["all", "protein", "ligands", "water", "ions", "custom"]:
        raise ValueError(f"selection keywords incorrect in {stepName}, see README.md for more details")

    if keyword == "custom":
        customSelection, = check_info_for_args(selection, stepName, ["customSelection"], optional=False)
        if not isinstance(customSelection, list):
            raise TypeError(f"selectionSyntax in {stepName} must be a list of dcits")
        for selection in customSelection:
            chainId, residueName, residueNumber, atomName = check_info_for_args(selection, stepName, ["CHAIN_ID", "RES_NAME", "RES_ID", "ATOM_NAME"], optional=False)
            if not isinstance(chainId, (str,list)):
                raise TypeError("CHAIN_ID must be a string or list of strings")
            if not isinstance(residueName, (str,list)):
                raise TypeError("RES_NAME must be a string or list of strings")
            if not isinstance(residueNumber, (int,list)):
                raise TypeError("RES_ID must be an int or list of ints")
            if not isinstance(atomName, (str,list)):
                raise TypeError("ATOM_NAME must be a string or list of strings")
#########################################################################
def check_meta_options(simulation: dict, stepName: str) -> None:
    ### TODO: rework how metaDynamicsInfo works
    return


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
        raise TypeError(f"{timeInputName} in {stepName} must be in the format \"10 ns\"")
    timeInputData = timeInputValue.split()
    try:
        numberAsFloat = int(timeInputData[0])
    except:
        raise ValueError(f"{timeInputName} in {stepName} must be in the format \"10 ns\"")
    if not timeInputData[1] in ["fs","ps","ns","ms"]:
        raise ValueError(f"{timeInputName} in {stepName} must be in the format \"10 ns\"")

#########################################################################
def check_shared_simulation_options(simulation) -> str:
    ## check for shared required args in simulation
    stepName, simulationType, temp = check_info_for_args(simulation, "simulation", ["stepName", "simulationType", "temp"], optional=False)
    ## make sure stepName is correct
    if not isinstance(stepName, str):
        raise TypeError(f"stepName {str(stepName)} must be a a string")
    if " " in stepName:
        raise ValueError(f"No whitespace allowed in stepName {stepName}")
    ## make sure simulationType is an allowed value
    if not simulationType.upper() in ["EM", "NVT", "NPT", "META"]:
        raise ValueError(f"simulationType in {stepName} must be the following:\n EM, NVT, NPT, META")
    ## make sure temp is an int
    if not isinstance(temp, int):
        raise TypeError(f"Temperature in {stepName} must be an int")
    if temp < 0:
        raise ValueError(f"Temperature in {stepName} must be a positive int")

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
            raise KeyError(f"Argument {argName} not found in {infoName}")
        unpackedDicts.append(argValue)
    return unpackedDicts
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
        raise TypeError(f"The config argument {argName} = {argPath} is not a PathLike.")
    # Check if the path exists
    if not p.exists(argPath):
        raise FileNotFoundError(f"The config argument {argName} = {argPath} does not exist.")
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
        raise FileNotFoundError(f"config file {configFile} not found")
    except yaml.YAMLError as exc:
        raise yaml.YAMLError("Error parsing YAML file:", exc)



#####################################################################################
def read_config(configYaml: str) -> dict:
    """
    Reads a config.yaml file and returns its contents as a dictionary.

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
def validate_config(config: dict) -> None:
    """
    Validates the configuration dictionary.

    Args:
        config (dict): The configuration dictionary.

    Raises:
        SystemExit: If the number of proteins in the config does not match nProteins,
                    if the number of ligands in the config does not match nLigands,
                    or if the input PDB file does not exist.

    Returns:
        None
    """
    # Check if the number of proteins in the config matches nProteins
    if not config["proteinInfo"]["nProteins"] == len(config["proteinInfo"]["proteins"]):
        print("Number of proteins in config does not match nProteins")
        exit()

    # Check if the number of ligands in the config matches nLigands
    if "ligandInfo" in config and config["ligandInfo"] is not None:
        if not config["ligandInfo"]["nLigands"] == len(config["ligandInfo"]["ligands"]):
            print("Number of ligands in config does not match nLigands")
            exit()

    # Check if the input PDB file exists
    if not p.isfile(config["pathInfo"]["inputPdb"]):
        print("Input PDB does not exist")
        exit()

