import os
from os import path as p 
from shutil import copy, rmtree
from pdbUtils import *

######################################################################################################
def  get_endpoint_pdbs(simulationInfo, outDir):                 
    excudeDirNames = ["00_configs","01_ligand_parameters"]
    for sim in simulationInfo:
        stepName = sim["stepName"]
        tag = sim["type"]
        collateDir = p.join(outDir,stepName)
        os.makedirs(collateDir,exist_ok=True)
        excudeDirNames.append(stepName)
        for protName in os.listdir(outDir):
            if protName in excudeDirNames:
                continue
            stepDir = p.join(outDir,protName,stepName)
            for file in os.listdir(stepDir):
                if p.splitext(file)[1] == ".pdb":
                    copy(p.join(stepDir,file),p.join(collateDir,f"{protName}_{tag}.pdb"))
######################################################################################################
def remove_atoms_from_pdb(simulationInfo, cleanUpInfo, outDir):
    for sim in simulationInfo:
        stepName = sim["stepName"]
        collateDir = p.join(outDir,stepName)
        for file in os.listdir(collateDir):
            if not p.splitext(file)[1] == ".pdb":
                continue
            pdbFile = p.join(collateDir,file)
            pdbDf = pdb2df(pdbFile)
            if "removeWaters" in cleanUpInfo:
                if cleanUpInfo["removeWaters"]:
                    pdbDf = pdbDf[~pdbDf["RES_NAME"].isin(["HOH"])].copy()
            if "removeIons" in cleanUpInfo:
                if cleanUpInfo["removeIons"]: ## PLACEHOLDER ION NAMES
                    pdbDf = pdbDf[~pdbDf["RES_NAME"].isin(["Na+","Cl-","Mg2+","F-"])].copy()
            df2Pdb(pdbDf, pdbFile)
######################################################################################################

