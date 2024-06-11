# chech health of a config file 
import yaml
from os import path as p
import re
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


#####################################################################################
def validate_config(config: dict) -> None:
    """
    Validates the configuration dictionary.

    Args:
        config (dict): The configuration dictionary.

    Raises:
        AssertionError: If any of the required sections, keys, or values are missing or invalid.

    Returns:
        None
    """

    # Assert the presence of main sections
    required_sections = ['pathInfo', 'generalInfo','simulationInfo','cleanUpInfo']
    for section in required_sections:
        assert section in config, f"Missing required section: {section}"

    # Validate pathInfo
    assert isinstance(config['pathInfo']['inputDir'], str), "inputDir must be a string"
    assert isinstance(config['pathInfo']['outputDir'], str), "outputDir must be a string"
    assert p.exists(config['pathInfo']['inputDir']), f"Input directory does not exist: {config['pathInfo']['inputDir']}"
    assert p.exists(p.dirname(config['pathInfo']['outputDir'])), f"At least parent directory of outputDir should exist: {config['pathInfo']['outputDir']}"

    # Validate generalInfo
    assert isinstance(config['generalInfo']['parallelCPU'], int), "parallelCPU must be an integer"
    assert isinstance(config['generalInfo']['platform'], str), "platform must be a string"
    assert config['generalInfo']['platform'] in ["CPU", "CUDA", "OpenCL"], "Platform must be either 'CPU', 'OpenCL' or 'CUDA'"
    assert isinstance(config['generalInfo']['subprocessCpus'], int), "subprocessCpus must be an integer"

    # Validate ligandInfo
    if 'ligandInfo' in config:
        assert isinstance(config['ligandInfo']['nLigands'], int), "nLigands must be an integer"
        assert isinstance(config['ligandInfo']['ligands'], list), "ligands must be a list"
        for ligand in config['ligandInfo']['ligands']:
            assert isinstance(ligand['ligandName'], str), "ligandName must be a string"
            assert isinstance(ligand['protons'], bool), "protons must be a boolean"
            assert isinstance(ligand['charge'], int), "charge must be an integer"

    # Validate simulationInfo
    filename_invalid_chars_pattern = re.compile(r'[^\w-]') # Allows letters, numbers, underscores, and hyphens
    permittedSimulationTypes = ["EM","NVT","NPT","META"]
    for step in config['simulationInfo']:
        
        assert isinstance(step['stepName'], str), "stepName must be a string"
        assert not filename_invalid_chars_pattern.search(step['stepName']), f"stepName contains invalid characters for a filename: {step['stepName']}. Allows letters, numbers, underscores, and hyphens"

        assert isinstance(step['type'], str), "type must be a string"
        assert step['type'].upper() in permittedSimulationTypes, f"simulation type can only be one of {permittedSimulationTypes}. You entered: {step['type']}"
        assert isinstance(step['temp'], int), "temp must be an integer"
        if 'maxIterations' in step:
            assert isinstance(step['maxIterations'], int), "maxIterations must be an integer"
            assert step['maxIterations']>=0, "maxIterations must be positive integer or 0 (for complete EM step)"
        
        # assert isinstance(step["timestep"], str), "timestep must be a string, eg \"2 fs\""
        # assert isinstance(step["duration"], str), "duration must be a string, eg \"10 ns\""


    # Validate cleanUpInfo
    assert isinstance(config['cleanUpInfo']['getEndpointPdbs'], bool), "getEndpointPdbs must be a boolean"
    assert isinstance(config['cleanUpInfo']['removeWaters'], bool), "removeWaters must be a boolean"
    assert isinstance(config['cleanUpInfo']['removeIons'], bool), "removeIons must be a boolean"
    assert isinstance(config['cleanUpInfo']['keepFileNames'], bool), "keepFileNames must be a boolean"

    print("Configuration is valid.")

# Example usage
# Load your config dictionary from the YAML file first
# Then pass it to the validation function
# validate_config(config)
