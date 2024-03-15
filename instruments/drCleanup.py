import os
from os import path as p 
from shutil import copy, rmtree
from instruments.pdbUtils import pdb2df, df2pdb
######################################################################################################
def clean_up_handler(batchConfig):
    # READ FROM batchConfig
    simulationInfo = batchConfig["simulationInfo"]
    pathInfo = batchConfig["pathInfo"]
    # RUN THROUGH OPTIONS IN cleanUpInfo
    if  not "cleanUpInfo" in batchConfig:
        return 
    cleanUpInfo = batchConfig["cleanUpInfo"]
    if "getEndpointPdbs" in cleanUpInfo:
        if cleanUpInfo["getEndpointPdbs"]:
            get_endpoint_pdbs(simulationInfo, pathInfo)
            if any(key in cleanUpInfo for key in ["removeWaters","removeIons"]):
                remove_atoms_from_pdb(simulationInfo, cleanUpInfo, pathInfo)
######################################################################################################
def get_endpoint_pdbs(simulationInfo, pathInfo):
    inputDir = pathInfo["inputDir"]
    outDir = pathInfo["outputDir"]

    inputNames = [p.splitext(pdbFile)[0] for pdbFile in os.listdir(inputDir) if p.splitext(pdbFile)[1] == ".pdb"]

    # make a new dir to collate pdbs to 
    collatedPdbDir = p.join(outDir, "collatedPdbs")
    os.makedirs(collatedPdbDir,exist_ok=True)
    excudeDirNames = ["00_configs","01_ligand_parameters"]
    for sim in simulationInfo:
        stepName = sim["stepName"]
        tag = sim["type"]
        collateSubDir = p.join(collatedPdbDir,stepName)
        os.makedirs(collateSubDir,exist_ok=True)
        excudeDirNames.append(stepName)
        for inputName in inputNames:
            if inputName in excudeDirNames:
                continue
            stepDir = p.join(outDir,inputName,stepName)
            for file in os.listdir(stepDir):
                if p.splitext(file)[1] == ".pdb":
                    copy(p.join(stepDir,file),p.join(collateSubDir,f"{inputName}_{tag}.pdb"))
######################################################################################################
def remove_atoms_from_pdb(simulationInfo, cleanUpInfo, pathInfo):
    outDir = pathInfo["outputDir"]
    collatedPdbDir = p.join(outDir, "collatedPdbs")

    for sim in simulationInfo:
        stepName = sim["stepName"]
        collateSubDir = p.join(collatedPdbDir,stepName)
        for file in os.listdir(collateSubDir):
            if not p.splitext(file)[1] == ".pdb":
                continue
            pdbFile = p.join(collateSubDir,file)
            pdbDf = pdb2df(pdbFile)
            if "removeWaters" in cleanUpInfo:
                if cleanUpInfo["removeWaters"]:
                    pdbDf = pdbDf[~pdbDf["RES_NAME"].isin(["HOH"])].copy()
            if "removeIons" in cleanUpInfo:
                if cleanUpInfo["removeIons"]: ## PLACEHOLDER ION NAMES
                    pdbDf = pdbDf[~pdbDf["RES_NAME"].isin(["Na+","Cl-","Mg2+","F-"])].copy()
            df2pdb(pdbDf, pdbFile)
######################################################################################################

