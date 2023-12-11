## BASIC LIBS
import os
from os import path as p
import pandas as pd
from subprocess import run
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
    ligandPdbs = []
    ligandFileDict = {}
    # for each ligand in config
    for ligand in ligandsDict:
        fileDict = {}
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
            # run pdb4amber to get compatable types and fix atom numbering
            ligPdb_amber = p.join(ligPrepDir,f"{ligandName}_amber.pdb")
            pdb4amberCommand = f"pdb4amber -i {ligPdb_h} -o {ligPdb_amber}"
            run(pdb4amberCommand,shell=True)
            ligPdb = ligPdb_amber
        ligandPdbs.append(ligPdb)
        # convert to mol2 with antechamber
        charge = ligand["charge"]
        ligMol2 = f"{ligandName}.mol2"
        antechamberCommand = f"antechamber -i {ligPdb} -fi pdb -o {ligMol2} -fo mol2 -c bcc -s 2 -nc {charge}"
        run(antechamberCommand, shell=True)
        # add mol2  path to config 
        fileDict.update({"mol2":ligMol2})
        # use mol2 to generate amber parameters with parmchk2
        ligFrcmod = f"{ligandName}.frcmod"
        parmchk2Command = f"parmchk2 -i {ligMol2} -f mol2 -o {ligFrcmod}"
        run(parmchk2Command,shell=True)
        # add frcmod path to config
        fileDict.update({"frcmod":ligFrcmod})
        ligandFileDict.update({ligandName:fileDict})
    return ligandPdbs, ligandFileDict
#####################################################################################
def prepare_protein_structure(config, outDir):
    proteinDict = config["proteinInfo"]["proteins"]
    proteinPdbs = []
    # for each protein in config
    for protein in proteinDict:
        # find files and directories
        protName = protein["proteinName"]
        protPrepDir = p.join(outDir,"PROT")
        os.chdir(protPrepDir)
        protPdb = p.join(protPrepDir,"PROT.pdb")
        if not protein["protons"]:
            # add protons with reduce
            protPdb_h = p.join(protPrepDir,"PROT_h.pdb")
            reduceCommand = f"reduce {protPdb} > {protPdb_h}"
            run(reduceCommand, shell=True)
            # run pdb4amber to get compatable types and fix atom numbering
            protPdb_amber = p.join(protPrepDir,"PROT_amber.pdb")
            pdb4amberCommand = f"pdb4amber -i {protPdb_h} -o {protPdb_amber}"
            run(pdb4amberCommand, shell = True)
            protPdb = protPdb_amber
            proteinPdbs.append(protPdb)
    return proteinPdbs
#####################################################################################
def make_amber_params(outDir,ligandFileDict, pdbFile, outName):
    os.chdir(outDir)
    # write tleap input file
    tleapInput = p.join(outDir, "TLEAP.in")
    with open(tleapInput,"w") as f:
        # amber ff and tip3p ff
        f.write("source oldff/leaprc.ff14SB\n")
        f.write("source leaprc.gaff2\n")
        f.write("source leaprc.water.tip3p\n\n")
        ## ligand mol2 and ff
        for ligandName in ligandFileDict:
            ligMol2 = ligandFileDict[ligandName]["mol2"]
            ligFrcmod = ligandFileDict[ligandName]["frcmod"]
            f.write(f"{ligandName} = loadmol2 {ligMol2}\n")
            f.write(f"loadamberparams {ligFrcmod}\n\n")
        ## ions ff
        f.write("loadamberparams frcmod.ions1lm_126_tip3p\n")
        ## whole protein structure
        f.write(f"loadpdb {pdbFile}\n")
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
    tleapCommand = f"tleap -s {tleapInput} > {tleapOutput}"
    run(tleapCommand,shell=True)