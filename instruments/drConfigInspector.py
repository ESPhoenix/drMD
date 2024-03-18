# chech health of a config file 
import yaml
from os import path as p
import re
#####################################################################################
def read_config(configYaml):

    ## Read config.yaml into a dictionary
    with open(configYaml,"r") as yamlFile:
        config = yaml.safe_load(yamlFile) 

    return config

#####################################################################################
def validate_config(config):
    if not config["proteinInfo"]["nProteins"] == len(config["proteinInfo"]["proteins"]):
        print("Number of proteins in config does not match nProteins")
        exit()
    if "ligandInfo" in config and not "ligandInfo" == None:
        if not config["ligandInfo"]["nLigands"] == len(config["ligandInfo"]["ligands"]):
            print("Number of ligands in config does not match nLigands")
            exit()
    if not p.isfile(config["pathInfo"]["inputPdb"]):
        print("Input PDB does not exist")
        exit()   


#####################################################################################
def validate_config(config):
    # Assert the presence of main sections
    required_sections = ['pathInfo', 'generalInfo','simulationInfo','cleanUpInfo']
    for section in required_sections:
        assert section in config, f"Missing required section: {section}"

    # Validate pathInfo
    assert 'inputDir' in config['pathInfo'], "Missing inputDir in pathInfo"
    assert 'outputDir' in config['pathInfo'], "Missing outputDir in pathInfo"
    assert p.exists(config['pathInfo']['inputDir']), f"Input directory does not exist: {config['pathInfo']['inputDir']}"
    assert p.exists(p.dirname(config['pathInfo']['outputDir'])), f"At least parent directory of outputDir should exist: {config['pathInfo']['outputDir']}"

    # Validate generalInfo
    assert 'parallelCPU' in config['generalInfo'], "Missing parallelCPU in generalInfo"
    assert isinstance(config['generalInfo']['parallelCPU'], int), "parallelCPU must be an integer"
    assert 'platform' in config['generalInfo'], "Missing platform in generalInfo"
    assert config['generalInfo']['platform'] in ["CPU", "CUDA", "OpenCl"], "Platform must be either 'CPU', 'OpenCl' or 'CUDA'"
    assert 'subprocessCpus' in config['generalInfo'], "Missing subprocessCpus in generalInfo"
    assert isinstance(config['generalInfo']['subprocessCpus'], int), "subprocessCpus must be an integer"

    # Validate ligandInfo
    if 'ligandInfo' in config:
        assert 'nLigands' in config['ligandInfo'], "Missing nLigands in ligandInfo"
        assert isinstance(config['ligandInfo']['nLigands'], int), "nLigands must be an integer"
        assert 'ligands' in config['ligandInfo'], "Missing ligands in ligandInfo"
        for ligand in config['ligandInfo']['ligands']:
            assert 'ligandName' in ligand, "Missing ligandName in ligand"
            assert 'protons' in ligand, "Missing protons in ligand"
            assert isinstance(ligand['protons'], bool), "protons must be a boolean"
            assert 'charge' in ligand, "Missing charge in ligand"
            assert isinstance(ligand['charge'], int), "charge must be an integer"

    # Validate simulationInfo
    filename_invalid_chars_pattern = re.compile(r'[^\w-]') # Allows letters, numbers, underscores, and hyphens
    required_keys = ['stepName', 'type', 'temp']
    permittedSimulationTypes = ["EM","NVT","NPT"]
    for step in config['simulationInfo']:
        
        for key in required_keys:
            assert key in step, f"Missing {key} in simulation step"
        
        assert not filename_invalid_chars_pattern.search(step['stepName']), f"stepName contains invalid characters for a filename: {step['stepName']}. Allows letters, numbers, underscores, and hyphens"

        assert step['type'].upper() in permittedSimulationTypes, f"simulation type can only be one of {permittedSimulationTypes}. You entered: {step['type']}"
        assert isinstance(step['temp'], int), "temp must be an integer"
        if 'maxIterations' in step:
            assert isinstance(step['maxIterations'], int), "maxIterations must be an integer"
            assert step['maxIterations']>=0, "maxIterations must be positive integer or 0 (for complete EM step)"
    # Validate cleanUpInfo
    required_keys = ['getEndpointPdbs', 'removeWaters', 'removeIons', 'keepFileNames']
    for key in required_keys:
        assert key in config['cleanUpInfo'], f"Missing {key} in cleanUpInfo"
        assert isinstance(config['cleanUpInfo'][key], bool), f"{key} must be a boolean"

    print("Configuration is valid.")

# Example usage
# Load your config dictionary from the YAML file first
# Then pass it to the validation function
# validate_config(config)
