import os
from os import path as p
import sys
import yaml
from glob import glob
import numpy as np
import inflect
## drMD LIBS
# Get the parent directory of the current script
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from typing import Optional, Dict, List, Set, Union
from instruments.drCustomClasses import FilePath, DirectoryPath

from pdbUtils import pdbUtils

##########################################################################################
def main(batchConfigYaml: Dict, configDir: DirectoryPath, outDir: DirectoryPath) -> None:
    global inflecter
    inflecter = inflect.engine()

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
    ## write simulation related methods
    simulationInfo = configDicts[0]["simulationInfo"]
    write_simulation_methods(methodsFile, simulationInfo)
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
            methods.write(f"Proteins were protonated using the following pH: {pH}. ")
            methods.write("This process also automatically creates disulfide bonds as appropriate.")
            methods.write(f"\n\n**WARNING** drMD did not protonate proteins {protonatedText}. ")
            methods.write("You will need to fill in this section manually.\n\n")


        methods.write("\n")

def count_waters(pdbFile: FilePath) -> int:
    pdbDf = pdbUtils.pdb2df(pdbFile)

    waterDf = pdbDf[pdbDf["RES_NAME"] == "WAT"]

    return len(waterDf) / 3


def write_solvation_charge_balence_methods(batchConfig, configDicts, methodsFile):
    boxGeometry = configDicts[0]["miscInfo"]["boxGeometry"]
    approxWaterCount, counterIonCounts = get_solvation_atom_counts(configDicts, batchConfig)

    with open(methodsFile, "a") as methods:
        ## info on solvation box
        methods.write(f"All proteins were placed in a {boxGeometry} solvation box with a 10 Å buffer between")
        methods.write(f" the protein and the nearest edge of the box. ")
        ## average water count
        methods.write(f"Approximately {approxWaterCount} TIP3P water molecules were added to this box. ")
        ## info on counter ions
        methods.write(f"Sodium and Chloride ions were added to the box to balence the charge of the system.\n")
        methods.write(f"A table showing the counts of counter ions are provided below:\n\n")
        methods.write("|\t Protein Name\t| \tSodium Ions\t| \tChloride Ions\t|\n")
        methods.write("|\t------------\t| \t------------\t| \t------------\t|\n")
        for protName in counterIonCounts:
            methods.write(f"|\t{protName}|{counterIonCounts[protName]['Sodium']}| {counterIonCounts[protName]['Chloride']}|\n")

        methods.write("\n")
################################################################################
def get_solvation_atom_counts(configDicts, batchConfig):
    waterCounts = []
    counterIonCounts = {}
    outDir = batchConfig["pathInfo"]["outputDir"]
    for config in configDicts:
        protName = config["proteinInfo"]["proteinName"]
        ## find a the solvated pdb file
        inputPdbName = p.splitext(p.basename(config["pathInfo"]["inputPdb"]))[0]
        prepDir = p.join(outDir, inputPdbName, "00_prep")

        if "ligandInfo" in config:
            solvationDir = p.join(prepDir,"WHOLE")
        else:
            solvationDir = p.join(prepDir,"PROT")
        solvatedPdb = [p.join(solvationDir, file) for file in os.listdir(solvationDir) if file.endswith("solvated.pdb")][0]

        waterCount = count_waters(solvatedPdb)
        waterCounts.append(waterCount)

        nNa, nCl = count_ions(solvatedPdb)

        counterIonCounts[protName] = {"Sodium": nNa, "Chloride": nCl}

    
    waterCounts = np.array(waterCounts)
    averageWaterCount = np.mean(waterCounts)
    approximateAverageWaterCount = int(round(averageWaterCount, 2 - len(str(int(averageWaterCount)))))


    return approximateAverageWaterCount, counterIonCounts
##########################################################################################

def count_ions(pdbFile):
    pdbDf = pdbUtils.pdb2df(pdbFile)

    naDf = pdbDf[pdbDf["RES_NAME"] == "Na+"]
    clDf = pdbDf[pdbDf["RES_NAME"] == "Cl-"]

    return len(naDf) , len(clDf)

##########################################################################################
def format_list(inputList: List[str]) -> str:
    if len(inputList) == 1:
        return inputList[0]
    else:
        return ", ".join(inputList[:-1]) + ", and " + inputList[-1]
##########################################################################################

def get_progression_word(stepIndex: int, maxSteps: int) -> str:
    if stepIndex == 0:
        return "Initially,"
    elif stepIndex == maxSteps - 1:
        return "Finally,"
    else:
        return "Next,"
##########################################################################################
def get_simulation_type_text(sim: Dict) -> str:
    simulationType = sim["simulationType"]
    if simulationType.upper() == "NPT":
        return "a simulation was performed using the isothermal-isobaric (NpT) ensemble"
    elif simulationType.upper() == "NVT":
        return "a simulation was performed using using the canonical (NVT) ensemble"
    elif simulationType.upper() == "EM":
        return "an energy mimimisation step was performed using the steepest descent method"
    elif simulationType == "META":
        return "A metadynamics simulation was performed"
    
##########################################################################################
def get_restraints_methods_text(sim: Dict) -> str:
    
    if not "restraintInfo" in sim:
        return ""
    
    restraintInfo = sim["restraintInfo"]
    text = ""
    for restraint in restraintInfo:
        text += f"{inflecter.a(restraint['restraintType']).capitalize()} restraint"
        text += f" with a force constant of {restraint['parameters']['k']} {get_force_constant_units(restraint['restraintType'])} "
        text += f" {get_restraint_target(restraint)} "
        text += f"was applied to {selection_to_text(restraint['selection'])}. "
    return text

##########################################################################################
def get_force_constant_units(restraintType: str) -> str:
    if restraintType == "position":
        return "kJ mol<sup>-1</sup> nm<sup>-2</sup>"
    elif restraintType == "distance":
        return "kJ mol<sup>-1</sup> nm<sup>-2</sup>"
    elif restraintType == "angle":
        return "kJ mol<sup>-1</sup> rad<sup>-2</sup>"
    elif restraintType == "torsion":
        return "kJ mol<sup>-1</sup> rad<sup>-2</sup>"
##########################################################################################
def selection_to_text(selection: Dict) -> str:
    keyword = selection["keyword"]
    if not keyword == "custom":
        return f"all {keyword} atoms in the system"  

    text = "the following atoms: "

    selectionTexts = []
    customSelections = selection["customSelection"]
    for customSelection in customSelections:
        selectionText = ""
        ## deal with atoms
        atomName = customSelection["ATOM_NAME"]
        if  atomName == "all":
            selectionText += "all atoms in"
        else:
            selectionText += f"atom{identifier_list_to_str(atomName)}"

        ## deal with resId and resName together
        residueId = customSelection["RES_ID"]
        residueName = customSelection["RES_NAME"]
        if not residueId == "all" and not residueName == "all":
            if isinstance(residueId, str) and isinstance(residueName, str):
                selectionText += f" in residue {residueName}{residueId}"
            else:
                selectionText += f" in residue{identifier_list_to_str(residueId)}{identifier_list_to_str(residueName)}"
        elif not residueId == "all":
            selectionText += f" in residue{identifier_list_to_str(residueId)}"
        elif not residueName == "all":
            selectionText += f" in residue{identifier_list_to_str(residueName)}"
        ## deal with chain
        chainId = customSelection["CHAIN_ID"]
        if not chainId == "all":
            selectionText += f" in chain{identifier_list_to_str(chainId)}"

        selectionTexts.append(selectionText)


    selectionTexts = format_list(selectionTexts)

    text += selectionTexts

    return text

##########################################################################################
def  get_restraint_target(restraint):
    restraintType = restraint["restraintType"]

    if restraintType == "position":
        return ""
    elif restraintType == "distance":
        return f" and an equilibrium distance of {restraint['parameters']['r0']} Å"
    elif restraintType == "angle":
        return f" and an equilibrium angle of {restraint['parameters']['theta0']} degrees"    
    elif restraintType == "torsion":
        return f" and an equilibrium dihedral angle of {restraint['parameters']['phi0']} degrees"

##########################################################################################
def identifier_list_to_str(identifier: Union[str, list]) -> str:
    if isinstance(identifier, str):
        return " " + identifier
    else:
        return f"s {format_list(identifier)}"

##########################################################################################
def write_per_step_simulation_methods(methodsFile, sim, stepIndex, maxSteps):
    with open(methodsFile, "a") as methods:
        ## progression word
        methods.write(f"{get_progression_word(stepIndex, maxSteps)} ")
        ## simulation type
        methods.write(f"{get_simulation_type_text(sim)}.\n")
        ## deal with EM and maxIterations
        if sim["simulationType"] == "EM":
            if sim["maxIterations"] == "-1":
                methods.write(f"This energy minimisation was performed until it reached convergence.\n")
            else:
                methods.write(f"This energy minimisation was performed for {sim['maxIterations']} steps, or until it reached convergence.\n")
        ## deal with NPT, NVT, META
        else:
            methods.write(f"This simulation was performed for {sim['duration']}") 
            methods.write(f" with a timestep of {sim['timestep']}. ")
            if "temperature" in sim:
                methods.write(f" at {sim['temperature']} K") 
            elif "temperatureRange" in sim:
                methods.write(f"The temperature of this simulation was incremented through the range")
                methods.write(f" {format_list([str(temp) + ' K' for temp in sim['temperatureRange']])} in equal sized steps. ")

        ## deal with restraints
        methods.write(f"{get_restraints_methods_text(sim)}\n")
    
        methods.write("\n\n")
def  write_simulation_methods(methodsFile, simulationInfo):
    for stepIndex, sim in enumerate(simulationInfo):
        write_per_step_simulation_methods(methodsFile, sim, stepIndex, len(simulationInfo))






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
    batchConfigYaml = "/home/esp/scriptDevelopment/drMD/prescriptions/standard_MD_config.yaml"
    main(batchConfigYaml, configDir, outDir)

