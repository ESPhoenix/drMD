## BASIC LIBS
import os
from os import path as p
import subprocess
from subprocess import run
import string
from shutil import copy
import time
import pandas as pd
## drMD UTILS
from pdbUtils import pdbUtils
from instruments import drFixer 
## CLEAN CODE
from typing import Tuple, List, Dict, Union, Any, Optional
#####################################################################################

def prep_protocol(config: dict) -> Tuple[str, str, str]:
    """
    Prepares the protocol for the simulation.

    Args:
        config (dict): The configuration dictionary containing the necessary information.

    Returns:
        Tuple[str, str, str]: A tuple containing the merged PDB file, input coordinates file, and Amber parameter files.

    This function prepares the protocol for the simulation by performing the following steps:
    1. Retrieves the output directory from the configuration dictionary.
    2. Creates the preparation directory if it doesn't exist.
    3. Checks if the preparation is already complete by checking for the presence of the necessary files.
    4. If the preparation is already complete, returns the relevant files.
    5. If the preparation is not complete, proceeds with the main preparation protocol.
    6. Splits the input PDB file into protein and ligand files.
    7. Prepares the ligand parameters and outputs the ligand PDB files.
    8. Prepares the protein structure.
    9. Combines the protein and ligand PDB files.
    10. Merges the protein PDB files.
    11. Makes Amber parameter files using TLEAP.
    12. Returns the merged PDB file, input coordinates file, and Amber parameter files.

    If the configuration dictionary does not contain ligand information, the function proceeds with the following steps:
    1. Prepares the protein structure.
    2. Merges the protein PDB files.
    3. Makes Amber parameter files using TLEAP.
    4. Returns the merged PDB file, input coordinates file, and Amber parameter files.

    Note:
        The function assumes that the necessary directories and files are present in the configuration dictionary.
    """
    outDir: str = config["pathInfo"]["outputDir"]

    protName: str = config["proteinInfo"]["proteinName"]

    prepDir: str = p.join(outDir,"00_prep")
    prepLog: str = p.join(prepDir,"prep.log")

    os.makedirs(prepDir,exist_ok=True)

    ###### skip prep if complete ######
    amberParams: str = False
    inputCoords: str = False
    wholeDir: str = p.join(prepDir,"WHOLE")
    if p.isdir(p.join(wholeDir)):
        for file in os.listdir(wholeDir):
            if p.splitext(file)[1] == ".prmtop":
                amberParams: str = p.join(wholeDir,file)
            elif p.splitext(file)[1] == ".inpcrd":
                inputCoords: str = p.join(wholeDir,file)
            if p.splitext(file)[1] == ".pdb" and not file == f"{protName}.pdb":
                pdbFile: str = p.join(wholeDir, file)
            
    if amberParams and inputCoords:
        return pdbFile, inputCoords, amberParams
    
    ### MAIN PREP PROTOCOL ###
    if "ligandInfo" in config:
        ## SPLIT INPUT PDB INTO PROT AND ONE FILE PER LIGAND
        inputPdb: str = config["pathInfo"]["inputPdb"]
        split_input_pdb(inputPdb =inputPdb,
                        config = config,
                        outDir=prepDir)
        ## PREPARE LIGAND PARAMETERS, OUTPUT LIGAND PDBS
        ligandPdbs, ligandFileDict = prepare_ligand_parameters(config = config, outDir = prepDir, prepLog = prepLog)
        ## PREPARE PROTEIN STRUCTURE
        protPdb = prepare_protein_structure(config=config, outDir = prepDir, prepLog=prepLog)
        ## RE-COMBINE PROTEIN AND LIGAND PDB FILES
        wholePrepDir: str = p.join(prepDir,"WHOLE")
        os.makedirs(wholePrepDir,exist_ok=True)
        allPdbs: List[str] = [protPdb] + ligandPdbs
        outName: str = config["pathInfo"]["outputName"]
        mergedPdb: str = p.join(wholePrepDir,f"{protName}.pdb")
        pdbUtils.mergePdbs(pdbList=allPdbs, outFile = mergedPdb)

        mergedPdb = drFixer.reset_atom_numbers(pdbFile = mergedPdb)
        ## MAKE AMBER PARAMETER FILES WITH TLEAP
        inputCoords, amberParams, solvatedPdb  = make_amber_params(outDir = wholePrepDir,
                            ligandFileDict=ligandFileDict,
                            pdbFile= mergedPdb,
                            outName= outName,
                            prepLog= prepLog)
        # solvatedPdb = drFixer.reset_chains_residues(goodPdb=mergedPdb, badPdb=solvatedPdb)

        return solvatedPdb, inputCoords, amberParams
    
    else:
        ## PREPARE PROTEIN STRUCTURE
        protPdb = prepare_protein_structure(config=config, outDir = prepDir, prepLog = prepLog)  
        ## MERGE PROTEIN PDBS
        outName: str = config["pathInfo"]["outputName"]
        ## MAKE AMBER PARAMETER FILES WITH TLEAP
        inputCoords, amberParams, solvatedPdb = make_amber_params(outDir = p.join(prepDir,"PROT"),
                                                        pdbFile= protPdb,
                                                        outName= outName,
                                                        prepLog = prepLog)
        return solvatedPdb, inputCoords, amberParams

#####################################################################################
def find_ligand_charge(
    ligDf: pd.DataFrame, ligName: str, outDir: str, pH: float
) -> int:
    """
    This function uses propka to identify charges on a ligand.
    
    Args:
        ligDf (pd.DataFrame): Ligand DataFrame.
        ligName (str): Ligand name.
        outDir (str): Output directory.
        pH (float): pH value.
        
    Returns:
        int: Total charge of the ligand.
    """
    
    # Change working directory to the output directory
    os.chdir(outDir)
    
    # Fix atom names in the ligand dataframe
    ligDf = pdbUtils.fix_atom_names(ligDf)
    ligDf = ligDf[ligDf["ELEMENT"] != "H"]
    # Create a temporary pdb file from the ligand dataframe
    tmpPdb = p.join(outDir, f"{ligName}.pdb")
    pdbUtils.df2pdb(ligDf, tmpPdb, chain=False)
    
    # Run propka to predict charges on the ligand
    proPkaCommand = f"propka3 {tmpPdb}"
    run_with_log(proPkaCommand, False, None)
    
    # Read the propka output to extract charge at the specified pH
    proPkaFile = f"{ligName}.pka"
    with open(proPkaFile, "r") as f:
        pkaPredictions: List[List[str]] = []
        extract = False
        for line in f:
            # Extract predictions from the propka output
            if line.startswith("SUMMARY OF THIS PREDICTION"):
                extract = True
            if extract and line.startswith("-"):
                break
            if extract and not line == "\n":
                pkaPredictions.append(line.split())
        pkaPredictions = pkaPredictions[2:]
        
        # Calculate total charge
        totalCharge = 0
        for pred in pkaPredictions:
            atomType: str = pred[-1][0]
            atomProb: float = float(pred[-2])
            if atomType == "O":
                if atomProb < pH:
                    totalCharge -= 1
            elif atomType == "N":
                if atomProb > pH:
                    totalCharge += 1
    
    # Clean up temporary files
    [os.remove(p.join(outDir, file)) for file in os.listdir(outDir) if not p.splitext(file)[1] == ".yaml"]
    
    return totalCharge

#####################################################################################
def split_input_pdb(inputPdb: str, config: dict, outDir: str) -> None:
    """
    Split an input PDB file into separate PDB files for each ligand and a protein-only PDB file.

    Args:
        inputPdb (str): The path to the input PDB file.
        config (dict): The configuration dictionary containing information about the ligands.
        outDir (str): The output directory where the split PDB files will be saved.

    Returns:
        None
    """
    # Read whole pdb into a df
    pdbDf = pdbUtils.pdb2df(inputPdb)
    # Write each ligand to a separate pdb file
    ligandsInfo: List[dict] = config["ligandInfo"]
    for ligand in ligandsInfo:
        ligandName: str = ligand["ligandName"]
        ligPrepDir: str = p.join(outDir, ligandName)
        os.makedirs(ligPrepDir, exist_ok=True)
        ligDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"] == ligandName]
        pdbUtils.df2pdb(ligDf, p.join(ligPrepDir, f"{ligandName}.pdb"), chain=False)
        pdbDf.drop(pdbDf[pdbDf["RES_NAME"] == ligandName].index, inplace=True)
    # Write protein only to pdb file
    protPrepDir: str = p.join(outDir, "PROT")
    os.makedirs(protPrepDir, exist_ok=True)
    pdbUtils.df2pdb(pdbDf, p.join(protPrepDir, "PROT.pdb"))
#############################  PROTONATION & PDB CREATION ###############################


def ligand_protonation(
    ligand: Dict[str, Union[bool, int, str]],
    ligPrepDir: str,
    ligandName: str,
    ligandPdbs: List[str],
    prepLog: str,
) -> Tuple[str, List[str]]:
    """
    Protonates a ligand if specified in the configuration, otherwise performs protonation using Open Babel and pdb4amber.

    Args:
        ligand (Dict[str, Union[bool, int, str]]): A dictionary containing information about the ligand.
                                                   The dictionary should contain the following keys:
                                                   - "protons": Whether to protonate the ligand (bool).
                                                   - "charge": The charge of the ligand (int).
                                                   - "ligandName": The name of the ligand (str).
        ligPrepDir (str): The directory where the ligand files are prepared (str).
        ligandName (str): The name of the ligand (str).
        ligandPdbs (List[str]): A list of ligand pdb files (List[str]).
        prepLog (str): The log file for preparation (str).

    Returns:
        Tuple[str, List[str]]: A tuple containing the protonated ligand pdb file (str) and the updated list of ligand pdb files (List[str]).
    """

    # If the ligand is already protonated, return the ligand pdb file
    if ligand["protons"]:
        ligPdb: str = p.join(ligPrepDir, f"{ligandName}.pdb")
        ligDf: pd.DataFrame = pdbUtils.pdb2df(ligPdb)
        ligDf: pd.DataFrame = pdbUtils.fix_atom_names(ligDf)
        pdbUtils.df2pdb(ligDf, ligPdb)
        # rename_hydrogens(ligPdb, ligPdb)
        ligandPdbs.append(ligPdb)
        return ligPdb, ligandPdbs

    # If the ligand is not protonated, perform protonation using Open Babel and pdb4amber
    else:
        # # find pdb ligand pdb file
        ligPdb: str = p.join(ligPrepDir, f"{ligandName}.pdb")
        ligPdb_H: str = p.join(ligPrepDir, f"{ligandName}_H.pdb")

        # Protonate the ligand using Open Babel
        obabelCommand: str = f"obabel {ligPdb} -O {ligPdb_H} -h"
        run_with_log(obabelCommand, prepLog, ligPdb_H)

        ligPdb_newH: str = p.join(ligPrepDir, f"{ligandName}_newH.pdb")

        # Rename the hydrogens in the ligand pdb file
        rename_hydrogens(ligPdb_H, ligPdb_newH)

        # Run pdb4amber to get compatible types and fix atom numbering
        ligPdb_amber: str = p.join(ligPrepDir, f"{ligandName}_amber.pdb")
        pdb4amberCommand: str = f"pdb4amber -i {ligPdb_newH} -o {ligPdb_amber}"
        run_with_log(pdb4amberCommand, prepLog, ligPdb_amber)

        ligPdb = ligPdb_amber

        ligandPdbs.append(ligPdb)
        return ligPdb, ligandPdbs

###############################  MOL2 CREATION #####################################
def ligand_mol2(
    ligand: Dict[str, Any],
    input_dir: str,
    ligand_name: str,
    lig_param_dir: str,
    lig_prep_dir: str,
    lig_pdb: str,
    lig_file_dict: Dict[str, str],
    prep_log: str,
) -> Tuple[str, Dict[str, str]]:
    """
    Create a mol2 file for the ligand.

    This function looks for a mol2 file from the configuration, then in the ligParamDir. If not found, it creates a new
    mol2 file using antechamber.

    Args:
        ligand (Dict[str, Any]): Ligand information dictionary.
        input_dir (str): Path to the input directory.
        ligand_name (str): Name of the ligand.
        lig_param_dir (str): Path to the ligand parameter directory.
        lig_prep_dir (str): Path to the ligand preparation directory.
        lig_pdb (str): Path to the ligand pdb file.
        lig_file_dict (Dict[str, str]): Dictionary containing ligand file information.
        prep_log (str): Path to the preparation log file.

    Returns:
        Tuple[str, Dict[str, str]]: A tuple containing the path to the mol2 file and the updated lig_file_dict.
    """

    # Look for mol2 from config, then in ligParamDir, if not found, create new mol2
    if ligand.get("mol2"):  # Look in config
        lig_mol2: str = p.join(input_dir, f"{ligand_name}.mol2")

    elif p.isfile(p.join(lig_param_dir, f"{ligand_name}.mol2")):  # Look in ligParamDir
        lig_mol2: str = p.join(lig_param_dir, f"{ligand_name}.mol2")
    else:  # Convert to mol2 with antechamber
        charge: int = int(ligand["charge"])
        lig_mol2: str = p.join(lig_prep_dir, f"{ligand_name}.mol2")
        antechamber_command = (
            f"antechamber -i {lig_pdb} -fi pdb -o {lig_mol2} -fo mol2 -c bcc -s 2 -nc {charge}"
        )
        run_with_log(antechamber_command, prep_log, lig_mol2)
        # Copy to ligParamDir for future use
        copy(lig_mol2, p.join(lig_param_dir, f"{ligand_name}.mol2"))
    # Add mol2 path to lig_file_dict
    lig_file_dict.update({"mol2": lig_mol2})
    return lig_mol2, lig_file_dict
######################### TOPPAR CREATION ##########################################
def ligand_toppar(ligand: dict, inputDir: str, ligandName: str, ligParamDir: str, ligPrepDir: str, ligMol2: str, ligFileDict: dict, prepLog: str) -> dict:
    """
    This function looks for an frcmod file from the configuration, then in the ligParamDir. If not found, it creates a new frcmod file using parmchk2.

    Args:
        ligand (dict): Ligand information dictionary.
        inputDir (str): Path to the input directory.
        ligandName (str): Name of the ligand.
        ligParamDir (str): Path to the ligand parameter directory.
        ligPrepDir (str): Path to the ligand preparation directory.
        ligMol2 (str): Path to the ligand mol2 file.
        ligFileDict (dict): Dictionary containing ligand file information.
        prepLog (str): Path to the preparation log file.

    Returns:
        dict: The updated ligFileDict with the path to the frcmod file.
    """

    # Look for frcmod from config, then in ligParamDir, if not found, create new frcmod
    if ligand.get("toppar"):  # Look in config
        ligFrcmod: str = p.join(inputDir, f"{ligandName}.frcmod")

    elif p.isfile(p.join(ligParamDir, f"{ligandName}.frcmod")):  # Look in ligParamDir
        ligFrcmod: str = p.join(ligParamDir, f"{ligandName}.frcmod")

    else:  # Create new frcmod using parmchk2
        ligFrcmod: str = p.join(ligPrepDir, f"{ligandName}.frcmod")
        parmchk2Command: str = f"parmchk2 -i {ligMol2} -f mol2 -o {ligFrcmod}"
        run_with_log(parmchk2Command, prepLog, ligFrcmod)
        copy(ligFrcmod, p.join(ligParamDir, f"{ligandName}.frcmod"))

    # Add frcmod path to ligFileDict
    ligFileDict.update({"frcmod": ligFrcmod})

    return ligFileDict
#####################################################################################
def prepare_ligand_parameters(config: dict, outDir: str, prepLog: str = None) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Prepares the ligand parameters for a given configuration file.

    Args:
        config (dict): The configuration file containing the ligand information.
        outDir (str): The output directory where the ligand parameters will be saved.
        prepLog (str, optional): The path to the preparation log file. Defaults to None.

    Returns:
        Tuple[List[str], Dict[str, Dict[str, str]]]: A tuple containing two elements:
            - ligandPdbs (List[str]): A list of ligand PDB files.
            - ligandFileDict (Dict[str, Dict[str, str]]): A dictionary containing the ligand file information.

    This function reads the inputs from the configuration file and prepares the ligand parameters. It creates a directory to save the parameter files and initializes a list to store the PDB files and a dictionary to store all the information. For each ligand in the configuration file, it finds the files and directories, performs ligand protonation, generates the ligand mol2 file, and creates the ligand toppar file. The function returns the list of ligand PDB files and the dictionary containing the ligand file information.
    """
    # read inputs from config file
    ligandsInfo: dict = config["ligandInfo"]
    inputDir: str = config["pathInfo"]["inputDir"]
    mainDir: str = p.dirname(config["pathInfo"]["outputDir"])
    # create a dir to save parameter files in (saves re-running on subsequent runs)
    ligParamDir: str = p.join(mainDir,"01_ligand_parameters")
    os.makedirs(ligParamDir,exist_ok=True)
    # initialise list to store pdb files and dict to store all infod
    ligandPdbs: list = []
    ligandFileDict: dict = {}
    # for each ligand in config
    for ligand in ligandsInfo:
        ligFileDict: dict = {}
        # find files and directories
        ligandName: str = ligand["ligandName"]
        ligPrepDir: str = p.join(outDir,ligandName)
        os.chdir(ligPrepDir)

        ligPdb, ligandPdbs       = ligand_protonation(ligand,ligPrepDir,ligandName,ligandPdbs,prepLog)  

        ligMol2, ligFileDict    = ligand_mol2(ligand,inputDir,ligandName,ligParamDir,
                                              ligPrepDir,ligPdb,ligFileDict,prepLog)
        
        ligFileDict             =       ligand_toppar(ligand,inputDir,ligandName,ligParamDir,
                                                      ligPrepDir,ligMol2,ligFileDict,prepLog)

        ligandFileDict.update({ligandName:ligFileDict})
    return ligandPdbs, ligandFileDict
#####################################################################################
def rename_hydrogens(pdbFile: str, outFile: str) -> None:
    """
    Renames all hydrogens in a PDB file to new names.

    Args:
        pdbFile (str): Path to the input PDB file.
        outFile (str): Path to the output PDB file.

    Returns:
        None
    """
    # Read PDB file and get the DataFrame of hydrogens
    pdbDf: pd.DataFrame = pdbUtils.pdb2df(pdbFile)
    hDf: pd.DataFrame = pdbDf[pdbDf["ELEMENT"]=="H"]

    # Generate new names for hydrogens
    letters: List[str] = list(string.ascii_uppercase)
    numbers: List[str] = [str(i) for i in range(1,10)]
    newNameHs: List[str] = ["H"+letter+number for letter in letters for number in numbers]

    # Rename hydrogens in the DataFrame
    count: int = 0
    for index, _ in hDf.iterrows():
        pdbDf.loc[index,"ATOM_NAME"] = newNameHs[count]
        count += 1

    # Write the modified DataFrame back to the PDB file
    pdbUtils.df2pdb(pdbDf, outFile, chain=False)
#####################################################################################
def prepare_protein_structure(
    config: dict, outDir: str, prepLog: str
) -> List[str]:
    """
    Prepare the protein structure for simulations.

    This function prepares the protein structure by copying the input PDB file
    to the output directory if it doesn't exist.

    Args:
        config (dict): The configuration dictionary.
        outDir (str): The output directory.
        prepLog (str): The path to the preparation log file.

    Returns:
       Union[os.PathLike, str]: The PATH to a PDB file.
    """

    # Find files and directories
    protPrepDir: str = p.join(outDir, "PROT")  # Directory to prepare the protein
    os.makedirs(protPrepDir, exist_ok=True)
    os.chdir(protPrepDir)
    protPdb = p.join(protPrepDir, "PROT.pdb")  # Path to the protein PDB file
    # Check for PROT.pdb in protPrepDir (won't be there if noLigand)
    if not p.isfile(protPdb):
        # Copy the input PDB file to the output directory
        copy(config["pathInfo"]["inputPdb"], protPdb)

    return protPdb


#####################################################################################
def make_amber_params(
    outDir: str,
    pdbFile: str,
    outName: str,
    prepLog: str,
    ligandFileDict: Optional[Dict[str, Dict[str, str]]] = None
) -> Tuple[str, str]:
    """
    Prepare the protein structure for simulations using Amber.

    Args:
        outDir (str): The output directory.
        pdbFile (str): The path to the input PDB file.
        outName (str): The name of the prepared protein.
        prepLog (str): The path to the preparation log file.
        ligandFileDict (Optional[Dict[str, Dict[str, str]]]): A dictionary containing the paths to
            the mol2 and frcmod files for each ligand. Defaults to None.

    Returns:
        Tuple[str, str]: The paths to the input coordinate file and the Amber
            parameter files.
    """
    # Change the working directory to the output directory
    os.chdir(outDir)

    # Write the TLEAP input file
    tleapInput: str = p.join(outDir, "TLEAP.in")
    with open(tleapInput, "w") as f:
        # Load Amber force fields and TIP3P water model
        f.write("source oldff/leaprc.ff14SB\n")
        f.write("source leaprc.gaff2\n")
        f.write("source leaprc.water.tip3p\n\n")

        # Load Amber force fields for ions
        f.write("loadamberparams frcmod.ions1lm_126_tip3p\n")
        f.write("loadamberparams frcmod.ions234lm_126_tip3p\n")

        # Load molecular mechanics parameters for ligands
        if ligandFileDict:
            for ligandName, ligandInfo in ligandFileDict.items():
                ligMol2 = ligandInfo["mol2"]
                ligFrcmod = ligandInfo["frcmod"]
                f.write(f"{ligandName} = loadmol2 {ligMol2}\n")
                f.write(f"loadamberparams {ligFrcmod}\n\n")

        # Load the protein structure
        f.write(f"mol = loadpdb {pdbFile}\n")

        # Solvate the protein and add ions
        f.write("solvatebox mol TIP3PBOX 10.0\n")
        f.write("addions mol Na+ 0\n")
        f.write("addions mol Cl- 0\n")

        # Save the solvated protein structure and parameter files
        solvatedPdb: str = f"{outName}_solvated.pdb"
        f.write(f"savepdb mol {solvatedPdb}\n")
        prmTop: str = p.join(outDir, f"{outName}.prmtop")
        inputCoords: str = p.join(outDir, f"{outName}.inpcrd")
        f.write(f"saveamberparm mol {prmTop} {inputCoords}\n")
        f.write("quit\n")

    # Execute TLEAP and log the output
    tleapOutput: str = p.join(outDir, "TLEAP.out")
    amberParams: str = p.join(outDir, f"{outName}.prmtop")
    tleapCommand: str = f"tleap -f {tleapInput} > {tleapOutput}"
    run_with_log(tleapCommand, prepLog, amberParams)

    # Reset chain and residue IDs in Amber PDB file
    solvatedPdb: str = p.join(outDir, solvatedPdb)
    inputCoords: str = p.join(outDir, f"{outName}.inpcrd")
    ## reset chain and residue IDs in amber PDB
    solvatedPdb: str = p.join(outDir, solvatedPdb)
    # drFixer.reset_chains_residues(pdbFile, solvatedPdb)
    return inputCoords, amberParams, solvatedPdb
#####################################################################################
def run_with_log(
    command: str,
    prepLog: Optional[str] = None,
    expectedOutput: Optional[str] = None
) -> None:
    """
    Execute a command and log its output.

    Args:
        command (str): The command to execute.
        prepLog (Optional[str], optional): The path to the log file. Defaults to None.
        expectedOutput (Optional[str], optional): The path to the expected output file. Defaults to None.

    Returns:
        None
    """
    # Retry the command until it succeeds or reaches the maximum number of retries
    # Execute the command and capture its output
    process: subprocess.CompletedProcess[str] = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True, 
        env = os.environ
    )


    # Check if the expected output file exists
    if expectedOutput is None or p.isfile(expectedOutput):
        return
    else:
        raise FileNotFoundError(f"Expected output:\n\t {expectedOutput}\n\t\t not found.")


#####################################################################################