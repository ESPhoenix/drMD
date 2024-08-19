import os
from os import path as p
import sys
import yaml
from glob import glob
import numpy as np
## drMD LIBS
# Get the parent directory of the current script
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from typing import Optional, Dict, List, Set
from instruments.drCustomClasses import FilePath, DirectoryPath

from pdbUtils import pdbUtils




##########################################################################################
def main(batchConfigYaml: Dict, configDir: DirectoryPath, outDir: DirectoryPath) -> None:
    ## find all config files
    configDicts = get_config_dicts(configDir)
    ## read batch config
    with open(batchConfigYaml, 'r') as f:
        batchConfig = yaml.safe_load(f)

    
    ## make outDir if it doesn't exist
    os.makedirs(outDir, exist_ok=True)
    ## create methods file
    methodsFile = p.join(outDir, "drMD_AutoMethods.md")
    ## write a header to the methods file
    with open(methodsFile, "w") as f:
        f.write("This methods section was automatically generated my drMD\n\n")
    ## write ligand-prep related methods
    write_ligand_parameterisation_methods(configDicts, methodsFile)
    ## write protein-prep related methods
    write_protein_preparation_methods(configDicts, methodsFile)
    ## write solvation related methods
    write_solvation_charge_balence_methods(batchConfig, configDicts, methodsFile)
##########################################################################################

def get_config_dicts(configDir: DirectoryPath) -> List[Dict]:

    configFiles = [p.join(configDir, file) for file in os.listdir(configDir) if p.splitext(file)[1] == ".yaml"]

    configDicts = []
    for configFile in configFiles:
        with open(configFile, 'r') as f:
            configDicts.append(yaml.safe_load(f))

    return configDicts


    

##########################################################################################
def write_ligand_parameterisation_methods(configFiles, methodsFile) -> None:
    ## find all ligands in config files, see what drMD has done with them
    obabelProtonatedLigands = set()
    antechamberChargesLigands = set()
    parmchkParamsLigands = set()

    allLigandNames = set()
    for configFile in configFiles:
        ligandInfo = configFile.get("ligandInfo", False)
        if not ligandInfo:
            continue
        for ligand in ligandInfo:
            allLigandNames.add(ligand["ligandName"])
            if ligand["protons"]:
                obabelProtonatedLigands.add(ligand["ligandName"])
            if ligand["mol2"]:
                antechamberChargesLigands.add(ligand["ligandName"])
            if ligand["toppar"]:
                parmchkParamsLigands.add(ligand["ligandName"])


    ## if no ligands in ligandInfo, then return an empty string 
    if len(allLigandNames) == 0:
        return 
    ## open methods file and append to it
    with open (methodsFile, "a") as methods:
        ## when we have ligands, but none have been processed by drMD, write a warning to fill this in manually
        if len(obabelProtonatedLigands) > 0 and len(antechamberChargesLigands) > 0 and len(parmchkParamsLigands) > 0:
            methods.write(f"\n\n**WARNING** drMD did not run any automated procedures to parameterise your ligands. \
                          You will need to fill in this section manually.\n\n")
            return

        methods.write(f"Ligand parameter generation was performed using the following procedure: ")

        if obabelProtonatedLigands == antechamberChargesLigands == parmchkParamsLigands:
            methods.write(f"All ligands were protonated using OpenBabel [Ref. {cite('obabel')}], newly added hydrogen atoms were then renamed to ensure compatability with the AMBER forcefeild. ")
            methods.write(f"Partial charges of all ligands were calculated, and atom types were assigned using antechamber [Ref. {cite('antechamber')}], and the parameters for the ligands were generated using parmchk [Ref. {cite('parmchk')}].")

        else:
            if len(obabelProtonatedLigands) > 0:
                ligands = format_list(list(obabelProtonatedLigands))
                methods.write(f"The ligands {ligands} were protonated using OpenBabel [Ref. {cite('obabel')}], novel hydrogen atoms were then renamed to ensure compatability with the AMBER forcefeild")
            if len(antechamberChargesLigands) > 0:
                ligands = format_list(list(antechamberChargesLigands))
                methods.write(f"For the ligands {ligands} partial charges were calculated and atom types were assigned using antechamber [Ref. {cite('antechamber')}]")
            if len(parmchkParamsLigands) > 0:
                ligands = format_list(list(parmchkParamsLigands))
                methods.write(f"The parameters for the ligands {ligands} were generated using parmchk [Ref. {cite('parmchk')}]")

        methods.write("\n")
##########################################################################################


def write_protein_preparation_methods(configDicts, methodsFile) -> None:
    ## work out which proteins were protonated, collect in a dict
    proteinsProtonated = {}
    for config in configDicts:
        proteinInfo = config.get("proteinInfo", False)
        isProteinProtonated = proteinInfo.get("protons", False)
        inputPdbName = p.splitext(p.basename(config["pathInfo"]["inputPdb"]))[0]
        proteinsProtonated[inputPdbName] = isProteinProtonated
    ## get the pH of simulations
    pH = configDicts[0]["miscInfo"]["pH"]

    ## open methods file and append to it
    with open (methodsFile, "a") as methods:
        ## if all of the proteins were protonated before drMD, write a warning
        if all(proteinsProtonated.values()):
            methods.write(f"\n\n**WARNING** drMD did not run any automated procedures to prepare your proteins.\
                          You will need to fill in this section manually.\n\n")
        ## if none of the proteins were protonated before drMD, write automated protonation procedure
        elif not any(proteinsProtonated.values()):
            methods.write(f"\nAll proteins were protonated using software pdb2pqr [Ref. {cite('pdb2pqr')}] \
                            which uses ProPKA to calculate per-residue proton affinities [Ref. {cite('propka')}]. ")
            methods.write(f"Proteins were protonated using the following pH: {pH}. ")
            methods.write("This process also automatically creates disulfide bonds as appropriate.")
        ## if some, but not all of the proteins were protonated before drMD, write automated protonation procedure
        else:
            protonatedProteinNames = [protName for protName in proteinsProtonated if proteinsProtonated[protName]]
            protonatedText = format_list(protonatedProteinNames)
            nonProtonatedProteinNames = [protName for protName in proteinsProtonated if not proteinsProtonated[protName]]
            nonProtonatedText = format_list(nonProtonatedProteinNames)


            methods.write(f"\nThe proteins {nonProtonatedText} were protonated using software pdb2pqr [Ref. {cite('pdb2pqr')}]")
            methods.write(f"which uses ProPKA to calculate per-residue proton affinities [Ref. {cite('propka')}].")
            methods.write(f"Proteins were protonated using the following pH: {pH}")
            methods.write("This process also automatically creates disulfide bonds as appropriate.")
            methods.write(f"\n\n**WARNING** drMD did not protonate proteins {protonatedText}. ")
            methods.write("You will need to fill in this secion manually.\n\n")


        methods.write("\n")

def count_waters(pdbFile: FilePath) -> int:
    pdbDf = pdbUtils.pdb2df(pdbFile)

    waterDf = pdbDf[pdbDf["RES_NAME"] == "WAT"]

    return len(waterDf) / 3


def write_solvation_charge_balence_methods(batchConfig, configDicts, methodsFile):

    boxGeometry = configDicts[0]["miscInfo"]["boxGeometry"]
    approxWaterCount = get_approximate_water_count(configDicts, batchConfig)

    with open(methodsFile, "a") as methods:
        methods.write(f"All proteins were placed in a {boxGeometry} solvation box with a 10 Å buffer between")
        methods.write(f" the protein and the nearest edge of the box. ")


        methods.write(f"Approximately {approxWaterCount} TIP3P water molecules were added to this box. ")


    ## TODO: amount of CL- / Na+ added


def get_approximate_water_count(configDicts, batchConfig):
    waterCounts = []
    outDir = batchConfig["pathInfo"]["outputDir"]
    for config in configDicts:
        proteinInfo = config.get("proteinInfo", False)
        inputPdbName = p.splitext(p.basename(config["pathInfo"]["inputPdb"]))[0]
        prepDir = p.join(outDir, inputPdbName, "00_prep")

        if "ligandInfo" in config:
            solvationDir = p.join(prepDir,"WHOLE")
        else:
            solvationDir = p.join(prepDir,"PROT")
        solvatedPdb = [p.join(solvationDir, file) for file in os.listdir(solvationDir) if file.endswith("solvated.pdb")][0]
        print(solvatedPdb)

        waterCount = count_waters(solvatedPdb)
        waterCounts.append(waterCount)
    
    waterCounts = np.array(waterCounts)
    averageWaterCount = np.mean(waterCounts)
    approximateAverageWaterCount = int(round(averageWaterCount, 2 - len(str(int(averageWaterCount)))))


    return approximateAverageWaterCount
    

##########################################################################################
def format_list(ligandList: List[str]) -> str:
    if len(ligandList) == 1:
        return ligandList[0]
    else:
        return ", ".join(ligandList[:-1]) + " and " + ligandList[-1]
##########################################################################################





    
    # During this step, any disulphide bonds present in the proteins were created based on the interatomic distances of cysteine gamma sulfur atoms.
    # """


# def write_solvation_charge_balence_methods(config):
#     f"""
#     1. Solvation:
#         - Protein was placed in a solvent box with dimensions {"DIMS"}
#         - Approximately {approxWaterCount} TIP3P waters were added to this box


#     2. Charge Balancing:
#         - To balence the charge of our system, {nCls} Chloride ions were added // {nNas} sodium ions were added.
#     """


def write_per_step_simulation_methods(simulationInfo):
    pass

def  write_simulation_methods(simulationInfo):
    pass


def cite(key) -> str:

    doiDict = {
        "pdb2pqr": ["10.1093/nar/gkm276", "10.1093/nar/gkh381"],
        "propka": ["10.1021/ct200133y", "10.1021/ct100578z"],
        "obabel": ["10.1186/1758-2946-3-33"],
        "antechamber": ["10.1002/jcc.20035", "10.1016/j.jmgm.2005.12.005"],
        "parmchk": ["10.1002/jcc.20035", "10.1016/j.jmgm.2005.12.005"]
    }


    citation = doiDict.get(key, "CITATION NOT FOUND")


    return format_list(citation)

if __name__ == "__main__":
    configDir = "/home/esp/scriptDevelopment/drMD/03_outputs/00_configs"
    outDir = "/home/esp/scriptDevelopment/drMD/03_outputs/00_methods"
    batchConfigYaml = "/home/esp/scriptDevelopment/drMD/prescriptions/em_config.yaml"
    main(batchConfigYaml, configDir, outDir)

