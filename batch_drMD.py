## BASIC LIBS
import os
from os import path as p
from shutil import copy, rmtree
## INPUT LIBS
import yaml
import argpass
## CUSTOM LIBS
from pdbUtils import *
from prep_drMD import *
#####################################################################################################
def read_inputs():
    ## create an argpass parser, read config file, snip off ".py" if on the end of file
    parser = argpass.ArgumentParser()
    parser.add_argument("--bc")
    args = parser.parse_args()

    batchConfig=args.bc
    ## Read config.yaml into a dictionary
    with open(batchConfig,"r") as yamlFile:
        batchConfig = yaml.safe_load(yamlFile) 
    return batchConfig
######################################################################################################
def main():
    ## SORT OUT DIRECTORIES
    topDir = os.getcwd()
    batchConfig = read_inputs()
    outDir = batchConfig["pathInfo"]["outputDir"]
    yamlDir = p.join(outDir,"00_configs")
    pdbDir = batchConfig["pathInfo"]["inputDir"]
    simInfo = batchConfig["simulationInfo"]
    os.makedirs(yamlDir,exist_ok=True)
    ## LOOP THROUGH PDB FILES IN INPUT DIR
    for pdbFile in os.listdir(pdbDir):
        ## SKIP IF NOT A PDB FILE
        fileData = p.splitext(pdbFile)
        if not fileData[1] == ".pdb":
            continue
        ## EXTRACT BASIC INFO, MAKE DIRS
        protName = fileData[0]
        pdbPath = p.join(pdbDir,pdbFile)
        runDir = p.join(outDir,protName)
        ## COVERT TO DF, EXTRACT REST OF INFO
        pdbDf = pdb2df(pdbPath)
        proteinInfo, ligandInfo = extract_info(pdbDf,pdbDir,protName,yamlDir)
        ## GET PATH INFO
        pathInfo = {"inputDir":pdbDir,
                    "inputPdb":pdbPath,
                    "outputDir":runDir,
                    "outputName":protName}
        ## COMBINE INFOS INTO ONE DICT (IGNORE LIGINFO IF EMPTY)
        if ligandInfo == None:
            config = {"pathInfo":pathInfo,
                  "proteinInfo":proteinInfo,
                  "simulationInfo":simInfo}
        else:
            config = {"pathInfo":pathInfo,
                    "proteinInfo":proteinInfo,
                    "ligandInfo":ligandInfo,
                    "simulationInfo":simInfo}
        ## WRITE CONFIG TO YAML
        yamlFile = p.join(yamlDir,f"{protName}_config.yaml")
        with open(yamlFile,"w") as f:
            yaml.dump(config,f,default_flow_style=False)

        ## RUN drMD MD PROTOCOL
        os.chdir(topDir)
        drMdCommand = f"python drMD.py --config {yamlFile}"
        run(drMdCommand,shell=True)
        
   ## CLEAN UP
    clean_up(batchConfig,outDir)
 
######################################################################################################
def clean_up(batchConfig,outDir):
    simulationInfo = batchConfig["simulationInfo"]
    cleanUpInfo = batchConfig["cleanUpInfo"]
    excudeDirNames = ["00_configs"]
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
                    
    if cleanUpInfo["pdbOnly"]:
        for dir in os.listdir(outDir):
            if dir in excudeDirNames:
                continue
            dirPath = p.join(outDir,dir)
            rmtree(dirPath)







######################################################################################################
def extract_info(pdbDf,pdbDir,protName,yamlDir): ## gets info from pdb file, writes a config file
    ## GET PROTEIN INFORMATION
    aminoAcids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
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
        return proteinInfo, ligandInfo
    
    ligList = []
    for ligName in ligNames:
        ligDf = ligsDf[ligsDf["RES_NAME"] == ligName]
        ligH = False
        if (ligDf["ELEMENT"] == "H").any():
            ligH = True
        ligMol2 = p.join(pdbDir,f"{ligName}.mol2")
        mol2 = False
        if p.isfile(ligMol2):
            mol2 = True
        ligFrcmod = p.join(pdbDir,f"{ligName}.frcmod")
        frcmod = False
        if p.isfile(ligFrcmod):
            frcmod = True       

        charge = find_ligand_charge(ligDf,ligName,yamlDir,pH=7.4)
        tmpDict = {"ligandName":ligName,
                   "protons":   ligH,
                   "mol2":      mol2,
                   "toppar":    frcmod,
                   "charge":    charge}
        ligList.append(tmpDict)
    nLigands = len(ligList)
    ligandInfo = {"nLigands":nLigands,
                  "ligands":ligList}

    return proteinInfo, ligandInfo
######################################################################################################
main()