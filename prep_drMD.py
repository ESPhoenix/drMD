## BASIC LIBS
import os
from os import path as p
import pandas as pd
from subprocess import run
import string
from shutil import copy
## drMD UTILS
from pdbUtils import *

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
        df2Pdb(ligDf,p.join(ligPrepDir,f"{ligandName}.pdb"))
        pdbDf.drop(pdbDf[pdbDf["RES_NAME"]==ligandName].index,inplace=True)
    # write protein only to pdb file
    protPrepDir = p.join(outDir,"PROT")
    os.makedirs(protPrepDir,exist_ok=True)
    df2Pdb(pdbDf,p.join(protPrepDir,"PROT.pdb"))
#####################################################################################
def prepare_ligand_parameters(config, outDir):
    ligandsDict = config["ligandInfo"]["ligands"]
    inputDir = config["pathInfo"]["inputDir"]
    ligandPdbs = []
    ligandFileDict = {}
    # for each ligand in config
    for ligand in ligandsDict:
        ligFileDict = {}
        # find files and directories
        ligandName = ligand["ligandName"]
        ligPrepDir = p.join(outDir,ligandName)
        os.chdir(ligPrepDir)
        ligPdb = p.join(ligPrepDir,f"{ligandName}.pdb")
        # add protons if not already present
        if not ligand["protons"]:
            # add protons with reduce
            ligPdb_h = p.join(ligPrepDir,f"{ligandName}_h.pdb")
            reduceCommand = f"reduce {ligPdb} > {ligPdb_h}"
            run(reduceCommand,shell=True)
            # rename all new hydrogens H1, H2, H3 ... (fixes 4-character names)
            ligPdb_newH = p.join(ligPrepDir,f"{ligandName}_newH.pdb")
            rename_hydrogens(ligPdb_h,ligPdb_newH)
            # run pdb4amber to get compatable types and fix atom numbering
            ligPdb_amber = p.join(ligPrepDir,f"{ligandName}_amber.pdb")
            pdb4amberCommand = f"pdb4amber -i {ligPdb_newH} -o {ligPdb_amber}"
            run(pdb4amberCommand,shell=True)
            ligPdb = ligPdb_amber
        ligandPdbs.append(ligPdb)

        if ligand["mol2"]:
            ligMol2 = p.join(inputDir,f"{ligandName}.mol2")
        else:
            # convert to mol2 with antechamber
            charge = ligand["charge"]
            ligMol2 = p.join(ligPrepDir,f"{ligandName}.mol2")
            antechamberCommand = f"antechamber -i {ligPdb} -fi pdb -o {ligMol2} -fo mol2 -c bcc -s 2 -nc {charge}"
            run(antechamberCommand, shell=True)
        # add mol2  path to config 
        ligFileDict.update({"mol2":ligMol2})

        # if frcmod file is provided, don't run param generation (expensive!)
        if ligand["toppar"]:
            ligFrcmod = p.join(inputDir,f"{ligandName}.frcmod")
        else :
            # use mol2 to generate amber parameters with parmchk2
            ligFrcmod = p.join(ligPrepDir,f"{ligandName}.frcmod")
            parmchk2Command = f"parmchk2 -i {ligMol2} -f mol2 -o {ligFrcmod}"
            run(parmchk2Command,shell=True)
        # add frcmod path to ligFileDict
        ligFileDict.update({"frcmod":ligFrcmod})
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
    df2Pdb(pdbDf,outFile)
#####################################################################################
def prepare_protein_structure(config, outDir):
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
            run(reduceCommand, shell=True)
            #run pdb4amber to get compatable types and fix atom numbering
            protPdb_amber = p.join(protPrepDir,"PROT_amber.pdb")
            pdb4amberCommand = f"pdb4amber -i {protPdb_h} -o {protPdb_amber}"
            run(pdb4amberCommand, shell = True)
            protPdb = protPdb
            proteinPdbs.append(protPdb)
    return proteinPdbs


#####################################################################################
def make_amber_params(outDir, pdbFile, outName,ligandFileDict=False):
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
    run(tleapCommand,shell=True)
    inputCoords = p.join(outDir, f"{outName}.inpcrd")
    amberParams = p.join(outDir, f"{outName}.prmtop")

    return inputCoords, amberParams