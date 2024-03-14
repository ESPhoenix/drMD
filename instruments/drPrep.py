## BASIC LIBS
import os
from os import path as p
import subprocess
from subprocess import run
import string
from shutil import copy
## drMD UTILS
from instruments.pdbUtils import pdb2df, df2pdb, fix_atom_names

#####################################################################################
def find_ligand_charge(ligDf,ligName,outDir,pH,prepLog):
    ## uses propka to identify charges on a ligand
    #make a temportaty pdb file from the ligand dataframe
    os.chdir(outDir)
    ligDf = fix_atom_names(ligDf)
    tmpPdb = p.join(outDir,f"{ligName}.pdb")
    df2pdb(ligDf,tmpPdb,chain=False)
    # run propka 
    proPkaCommand = f"propka3 {tmpPdb}"
    run_with_log(proPkaCommand,prepLog)

    proPkaFile = f"{ligName}.pka"
    # read propka output to extract charge at specified pH
    with open(proPkaFile,"r") as f:
        pkaPredictions = []
        extract = False
        for line in f:
            if line.startswith("SUMMARY OF THIS PREDICTION"):
                extract = True
            if extract and line.startswith("-"):
                break
            if extract and not line =="\n":
                pkaPredictions.append(line.split())
        pkaPredictions = pkaPredictions[2:]
        totalCharge = 0
        for pred in pkaPredictions:
            if pred[-1][0] == "O":
                if float(pred[-2]) < pH:
                    totalCharge += -1
            elif pred[-1][0] == "N":
                if float(pred[-2]) > pH:
                    totalCharge += 1
    # clean up
    os.remove(tmpPdb)
    os.remove(proPkaFile)
    return totalCharge

#####################################################################################
def split_input_pdb(inputPdb,config,outDir):
    # read whole pdb into a df
    pdbDf = pdb2df(inputPdb)
    # write each ligand to a separate pdb file
    ligandsDict = config["ligandInfo"]["ligands"]
    for ligand in ligandsDict:
        ligandName = ligand["ligandName"]
        ligPrepDir = p.join(outDir,ligandName)
        os.makedirs(ligPrepDir,exist_ok=True)
        ligDf = pdbDf[pdbDf["RES_NAME"]==ligandName]
        df2pdb(ligDf,p.join(ligPrepDir,f"{ligandName}.pdb"),chain=False)
        pdbDf.drop(pdbDf[pdbDf["RES_NAME"]==ligandName].index,inplace=True)
    # write protein only to pdb file
    protPrepDir = p.join(outDir,"PROT")
    os.makedirs(protPrepDir,exist_ok=True)
    df2pdb(pdbDf,p.join(protPrepDir,"PROT.pdb"))
#############################  PROTONATION & PDB CREATION ###############################
def ligand_protonation(ligand,ligPrepDir,ligandName,ligandPdbs, prepLog):
    if ligand["protons"]:
        ligPdb = p.join(ligPrepDir,f"{ligandName}.pdb")
        ligDf = pdb2df(ligPdb)
        ligDf = fix_atom_names(ligDf)
        df2pdb(ligDf, ligPdb)
        rename_hydrogens(ligPdb, ligPdb)
        ligandPdbs.append(ligPdb)
        return ligPdb, ligandPdbs
    else:
        # # find pdb ligand pdb file
        ligPdb = p.join(ligPrepDir,f"{ligandName}.pdb")
        # # fix atom names 
        # ligDf = pdb2df(ligPdb)
        # ligDf = fix_atom_names(ligDf)
        # ligPdb_fixNames = p.join(ligPrepDir,f"{ligandName}_fixNames.pdb")
        # df2pdb(ligDf, ligPdb_fixNames,chain = False)
        # # remove pre-existing protons with reduce
        # ligPdb_noH = p.join(ligPrepDir,f"{ligandName}_noH.pdb")
        # trimCommand = f"reduce -Trim {ligPdb_fixNames} > {ligPdb_noH}"
        # run_with_log(trimCommand,prepLog)
        # # add protons with reduce
        # ligPdb_h = p.join(ligPrepDir,f"{ligandName}_H.pdb")
        # reduceCommand = f"reduce {ligPdb_noH} > {ligPdb_h}"
        # run_with_log(reduceCommand,prepLog)
        # rename all new hydrogens H1, H2, H3 ... (fixes 4-character names)
        # ligPdb_newH = p.join(ligPrepDir,f"{ligandName}_newH.pdb")
        ligPdb_H = p.join(ligPrepDir,f"{ligandName}_H.pdb")
        obabelCommand = f"obabel {ligPdb} -O {ligPdb_H} -h"
        run_with_log(obabelCommand, prepLog)
        ligPdb_newH = p.join(ligPrepDir,f"{ligandName}_newH.pdb")
        rename_hydrogens(ligPdb_H,ligPdb_newH)
        # run pdb4amber to get compatable types and fix atom numbering
        ligPdb_amber = p.join(ligPrepDir,f"{ligandName}_amber.pdb")
        pdb4amberCommand = f"pdb4amber -i {ligPdb_newH} -o {ligPdb_amber}"
        run_with_log(pdb4amberCommand,prepLog)
        ligPdb = ligPdb_amber
        ligandPdbs.append(ligPdb)
        return ligPdb, ligandPdbs
###############################  MOL2 CREATION #####################################
def  ligand_mol2(ligand,inputDir,ligandName,ligParamDir,ligPrepDir,ligPdb,ligFileDict, prepLog):
    ####  MOL2 CREATION ####
    # look for mol2 from config, then in ligParamDir, if not found, create new mol2
    if ligand["mol2"]:  # look in config
        ligMol2 = p.join(inputDir,f"{ligandName}.mol2")

    elif p.isfile(p.join(ligParamDir,f"{ligandName}.mol2")): # look in ligParamDir
        ligMol2 = p.join(ligParamDir,f"{ligandName}.mol2")
    else:  # convert to mol2 with antechamber
        charge = ligand["charge"]
        ligMol2 = p.join(ligPrepDir,f"{ligandName}.mol2")
        antechamberCommand = f"antechamber -i {ligPdb} -fi pdb -o {ligMol2} -fo mol2 -c bcc -s 2 -nc {charge}"
        run_with_log(antechamberCommand,prepLog)
        # copy to ligParamDir for future use
        copy(ligMol2,p.join(ligParamDir,f"{ligandName}.mol2"))
    # add mol2  path to ligFileDict 
    ligFileDict.update({"mol2":ligMol2})
    return ligMol2, ligFileDict
######################### TOPPAR CREATION ##########################################

def ligand_toppar(ligand,inputDir,ligandName,ligParamDir,ligPrepDir,ligMol2,ligFileDict, prepLog):
    # look for frcmod from config, then in ligParamDir, if not found, create new frcmod
    if ligand["toppar"]: # look in config
        ligFrcmod = p.join(inputDir,f"{ligandName}.frcmod") 
    elif p.isfile(p.join(ligParamDir,f"{ligandName}.frcmod")): # look in ligParamDir
        ligFrcmod = p.join(ligParamDir,f"{ligandName}.frcmod")    
    else :    # use mol2 to generate amber parameters with parmchk2
        ligFrcmod = p.join(ligPrepDir,f"{ligandName}.frcmod")
        parmchk2Command = f"parmchk2 -i {ligMol2} -f mol2 -o {ligFrcmod}"
        run_with_log(parmchk2Command,prepLog)
        copy(ligFrcmod,p.join(ligParamDir,f"{ligandName}.frcmod"))
    # add frcmod path to ligFileDict
    ligFileDict.update({"frcmod":ligFrcmod})

    return ligFileDict


#####################################################################################
def prepare_ligand_parameters(config, outDir, prepLog):
    # read inputs from config file
    ligandsDict = config["ligandInfo"]["ligands"]
    inputDir = config["pathInfo"]["inputDir"]
    mainDir = p.dirname(config["pathInfo"]["outputDir"])
    # create a dir to save parameter files in (saves re-running on subsequent runs)
    ligParamDir = p.join(mainDir,"01_ligand_parameters")
    os.makedirs(ligParamDir,exist_ok=True)
    # initialise list to store pdb files and dict to store all info
    ligandPdbs = []
    ligandFileDict = {}
    # for each ligand in config
    for ligand in ligandsDict:
        ligFileDict = {}
        # find files and directories
        ligandName = ligand["ligandName"]
        ligPrepDir = p.join(outDir,ligandName)
        os.chdir(ligPrepDir)

        ligPdb, ligandPdbs       = ligand_protonation(ligand,ligPrepDir,ligandName,ligandPdbs,prepLog)  

        ligMol2, ligFileDict    = ligand_mol2(ligand,inputDir,ligandName,ligParamDir,
                                              ligPrepDir,ligPdb,ligFileDict,prepLog)
        
        ligFileDict             =       ligand_toppar(ligand,inputDir,ligandName,ligParamDir,
                                                      ligPrepDir,ligMol2,ligFileDict,prepLog)

        ligandFileDict.update({ligandName:ligFileDict})
    return ligandPdbs, ligandFileDict
#####################################################################################
def rename_hydrogens(pdbFile,outFile):
    pdbDf = pdb2df(pdbFile)
    hDf = pdbDf[pdbDf["ELEMENT"]=="H"]
    letters = list(string.ascii_uppercase)
    numbers = [str(i) for i in range(1,10)]
    newNameHs = []
    for letter in letters:
        for number in numbers:
            newNameHs.append("H"+letter+number)
    count = 0
    for index, row in hDf.iterrows():
        pdbDf.loc[index,"ATOM_NAME"] = newNameHs[count]
        count += 1
    df2pdb(pdbDf,outFile,chain=False)
#####################################################################################
def prepare_protein_structure(config, outDir, prepLog):
    proteinDict = config["proteinInfo"]["proteins"]
    proteinPdbs = []
    # for each protein in config
    for protein in proteinDict:
        # find files and directories
        protPrepDir = p.join(outDir,"PROT")
        os.makedirs(protPrepDir,exist_ok=True)
        os.chdir(protPrepDir)
        protPdb = p.join(protPrepDir,"PROT.pdb")
        # check for PROT.pdb in protPrepDir (won't be there if noLigand)
        if not p.isfile(protPdb):
            copy(config["pathInfo"]["inputPdb"],protPdb)

        if not protein["protons"]:
            # add protons with reduce
            protPdb_h = p.join(protPrepDir,"PROT_h.pdb")
            reduceCommand = f"reduce {protPdb} > {protPdb_h}"
            run_with_log(reduceCommand,prepLog)
          
             #run pdb4amber to get compatable types and fix atom numbering
            protPdb_amber = p.join(protPrepDir,"PROT_amber.pdb")
            pdb4amberCommand = f"pdb4amber -i {protPdb_h} -o {protPdb_amber}"
            run_with_log(pdb4amberCommand,prepLog)
 
            proteinPdbs.append(protPdb)
    return proteinPdbs


#####################################################################################
def make_amber_params(outDir, pdbFile, outName,prepLog,ligandFileDict=False):
    os.chdir(outDir)
    # write tleap input file
    tleapInput = p.join(outDir, "TLEAP.in")
    with open(tleapInput,"w") as f:
        # amber ff and tip3p ff
        f.write("source oldff/leaprc.ff14SB\n")
        f.write("source leaprc.gaff2\n")
        f.write("source leaprc.water.tip3p\n\n")
        if ligandFileDict:
            ## ligand mol2 and ff
            for ligandName in ligandFileDict:
                ligMol2 = ligandFileDict[ligandName]["mol2"]
                ligFrcmod = ligandFileDict[ligandName]["frcmod"]
                f.write(f"{ligandName} = loadmol2 {ligMol2}\n")
                f.write(f"loadamberparams {ligFrcmod}\n\n")
        ## ions ff
        f.write("loadamberparams frcmod.ions1lm_126_tip3p\n")
        ## whole protein structure
        f.write(f"mol = loadpdb {pdbFile}\n")
        # solvation and ions
        f.write("solvatebox mol TIP3PBOX 10.0\n")
        f.write("addions mol Na+ 0\n")
        f.write("addions mol Cl- 0\n")
        # save solvated pdb file
        f.write(f"savepdb mol {outName}.pdb\n")
        # save parameter files
        f.write(f"saveamberparm mol {outName}.prmtop {outName}.inpcrd\n")
        f.write("quit\n")
    tleapOutput = p.join(outDir,"TLEAP.out")
    tleapCommand = f"tleap -f {tleapInput} > {tleapOutput}"
    run_with_log(tleapCommand,prepLog)

    inputCoords = p.join(outDir, f"{outName}.inpcrd")
    amberParams = p.join(outDir, f"{outName}.prmtop")

    return inputCoords, amberParams
#####################################################################################
def run_with_log(command, prepLog):
    process = run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(prepLog,"a") as log:
        log.write(process.stdout)
#####################################################################################