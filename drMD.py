## BASIC LIBS
import os
from os import path as p
## INPUT LIBS
import yaml
import argpass
## CUSTOM DR MD MODULES
from instruments.pdbUtils import pdb2df
import instruments.drPrep as drPrep
import instruments.drCleanup as drCleanup
import instruments.drOperator as drOperator
import instruments.drConfigInspector as drConfigInspector

## Multiprocessing
import concurrent.futures as cf
from tqdm import tqdm
from subprocess import run
#####################################################################################################
def read_inputs():
    ## create an argpass parser, read config file, snip off ".py" if on the end of file
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    batchConfig=args.config
    ## Read config.yaml into a dictionary
    with open(batchConfig,"r") as yamlFile:
        batchConfig = yaml.safe_load(yamlFile) 
    return batchConfig
######################################################################################################
def main():
    ## SORT OUT DIRECTORIES
    topDir = os.getcwd()
    batchConfig = read_inputs()
    drConfigInspector.validate_config(batchConfig)
    outDir = batchConfig["pathInfo"]["outputDir"]
    yamlDir = p.join(outDir,"00_configs")
    pdbDir = batchConfig["pathInfo"]["inputDir"]
    simInfo = batchConfig["simulationInfo"]
    parallelCPU = batchConfig["generalInfo"]["parallelCPU"]
    subprocessCpus = batchConfig["generalInfo"]["subprocessCpus"]
    
    manage_cpu_usage_for_subprocesses("ON",subprocessCpus)

    os.makedirs(yamlDir,exist_ok=True)
    if parallelCPU == 1:
        run_serial(batchConfig, pdbDir, outDir, yamlDir, simInfo)
    elif parallelCPU > 1:
        run_parallel(parallelCPU, batchConfig, pdbDir, outDir, yamlDir, simInfo)


    manage_cpu_usage_for_subprocesses("OFF")

######################################################################################################
def manage_cpu_usage_for_subprocesses(mode , subprocessCpus = None):
    if mode == "ON":
        # Set environment variables directly in Python
        os.environ['OMP_NUM_THREADS'] = str(subprocessCpus)
        os.environ['OPENMM_CPU_THREADS'] = str(subprocessCpus)

        omp_num_threads = os.environ.get('OMP_NUM_THREADS')
        print("OMP_NUM_THREADS:", omp_num_threads)
    elif mode == "OFF":
        # remove cpu useage limits
        run(f"unset OMP_NUM_THREADS", shell=True)
        run(f"unset OPENMM_CPU_THREADS", shell=True)

#####################################################################################################
def process_pdb_file(pdbFile, pdbDir, outDir, yamlDir, simInfo, batchConfig):
    # Skip if not a PDB file
    fileData = p.splitext(pdbFile)
    if not fileData[1] == ".pdb":
        return
    # Extract basic info, make dirs
    protName = fileData[0]
    pdbPath = p.join(pdbDir, pdbFile)
    runDir = p.join(outDir, protName)
    os.makedirs(runDir, exist_ok=True)

    # Convert to DataFrame, extract rest of info
    pdbDf = pdb2df(pdbPath)
    proteinInfo, ligandInfo, generalInfo = extract_info(pdbDf, pdbDir, protName, yamlDir, batchConfig)

    # Get path info
    pathInfo = {
        "inputDir": pdbDir,
        "inputPdb": pdbPath,
        "outputDir": runDir,
        "outputName": protName
    }

    # Combine infos into one dict (ignore ligInfo if empty)
    if ligandInfo is None:
        config = {
            "pathInfo": pathInfo,
            "generalInfo": generalInfo,
            "proteinInfo": proteinInfo,
            "simulationInfo": simInfo
        }
    else:
        config = {
            "pathInfo": pathInfo,
            "generalInfo": generalInfo,
            "proteinInfo": proteinInfo,
            "ligandInfo": ligandInfo,
            "simulationInfo": simInfo
        }
    # Write config to YAML
    configYaml = p.join(yamlDir, f"{protName}_config.yaml")
    with open(configYaml, "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    drOperator.drMD_protocol(configYaml)

###################################################################################################### 
def run_serial(batchConfig, pdbDir, outDir, yamlDir, simInfo):
    for pdbFile in os.listdir(pdbDir):
        fileData = p.splitext(pdbFile)
        if not fileData[1] == ".pdb":
            continue  
        process_pdb_file(pdbFile, pdbDir, outDir, yamlDir, simInfo, batchConfig)
    ## CLEAN UP
    drCleanup.clean_up_handler(batchConfig)
######################################################################################################
def run_parallel(parallelCPU, batchConfig, pdbDir, outDir, yamlDir, simInfo):

    pdbFiles = [pdbFile for pdbFile in os.listdir(pdbDir) if p.splitext(pdbFile)[1] == ".pdb"]

    # Create a ThreadPoolExecutor with the desired number of worker threads
    with cf.ThreadPoolExecutor(max_workers=parallelCPU) as executor:
        # Prepare the arguments for process_pdb_file
        inputArgs = [(pdbFile, pdbDir, outDir, yamlDir, simInfo, batchConfig) for pdbFile in pdbFiles]

        # Use starmap to apply process_pdb_file to each PDB file
        list(tqdm(executor.map(lambda args: process_pdb_file(*args), inputArgs), total=len(inputArgs)))

    # CLEAN UP
    drCleanup.clean_up_handler(batchConfig)
######################################################################################################
def extract_info(pdbDf,pdbDir,protName,yamlDir,batchConfig): ## gets info from pdb file, writes a config file
    generalInfo = batchConfig["generalInfo"]
    
    ## GET PROTEIN INFORMATION
    aminoAcids =   ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                    'SER', 'THR', 'TRP', 'TYR', 'VAL']
    protDf = pdbDf[pdbDf["RES_NAME"].isin(aminoAcids)]

    protH = False
    if (protDf["ELEMENT"] == "H").any():
        protH = True

    proteinInfo = {"nProteins":1,
                   "proteins":[{"proteinName":f"{protName}",
                                "protons":protH}]}   
    ## GET LIGAND INFORMATION 
    ligsDf = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcids)]
    ligNames = ligsDf["RES_NAME"].unique().tolist()

    ## SKIP IF NOT LIGAND
    if len(ligNames) == 0:
        ligandInfo = None
        return proteinInfo, ligandInfo, generalInfo 
    
    ## USE ligandInfo IF SUPPLIED IN BATCH CONFIG
    if "ligandInfo" in batchConfig:
        ligandInfo = batchConfig["ligandInfo"]
        return proteinInfo, ligandInfo, generalInfo
    ## CREATE ligandInfo AUTOMATICALLY (WORKS FOR SIMPLE LIGANDS)
    else:
        ligList = []
        for ligName in ligNames:
            ligDf = ligsDf[ligsDf["RES_NAME"] == ligName]
            # deal with protonation
            ligH = False
            if (ligDf["ELEMENT"] == "H").any():
                ligH = True

            # deal with mol2
            ligMol2 = p.join(pdbDir,f"{ligName}.mol2")
            mol2 = False
            if p.isfile(ligMol2):
                mol2 = True
            # deal with frcmod
            ligFrcmod = p.join(pdbDir,f"{ligName}.frcmod")
            frcmod = False
            if p.isfile(ligFrcmod):
                frcmod = True  
            # deal with charge
            charge = drPrep.find_ligand_charge(ligDf,ligName,yamlDir,pH=7.4)
            # write to temporary dict, then to ligandInfo for config
            tmpDict = {"ligandName":ligName,
                    "protons":   ligH,
                    "mol2":      mol2,
                    "toppar":    frcmod,
                    "charge":    charge}
            ligList.append(tmpDict)
        nLigands = len(ligList)
        ligandInfo = {"nLigands":nLigands,
                    "ligands":ligList}

        return proteinInfo, ligandInfo, generalInfo
######################################################################################################
main()