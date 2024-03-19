import os
from os import path as p
import pandas as pd
import mdtraj as md
from scipy.signal import find_peaks
import yaml
#############################################################################################
def compute_delta_RMSF(analDir, referenceSystem = False):
    resultDfsToConcat =[]
    sysNames = [name for name in os.listdir(analDir) if p.isdir(p.join(analDir,name))]
    for sysName in sysNames:
        sysAnalDir = p.join(analDir, sysName)
        ## load contact dfs from sysAnalDir | concat into one big df
        dfsToConcat = []
        for file in os.listdir(sysAnalDir):
            if not p.splitext(file)[1] == ".csv":
                continue
            if  file.startswith("RMSF"):
                runDf = pd.read_csv(p.join(sysAnalDir,file),index_col="RES_ID")
                dfsToConcat.append(runDf)
        if len(dfsToConcat) == 0:
            return
        df = pd.concat(dfsToConcat, axis = 1, ignore_index=False) 
        meanDf = df.mean(axis=1)
        meanDf.name = sysName
        resultDfsToConcat.append(meanDf)
    meanDf = pd.concat(resultDfsToConcat, axis=1)
    colNames = []
    dfsToConcat = []
    if  referenceSystem:
        for i in range(0,len(sysNames)):
            sysName_A = sysNames[i]
            sysName_B = referenceSystem
            if sysName_A == sysName_B:
                continue
            colName = f"{sysName_A} - {sysName_B}"
            colNames.append(colName)
            delta_AB = meanDf[sysName_A] - meanDf[sysName_B]
            dfsToConcat.append(delta_AB)    
    else:
        for i in range(0,len(sysNames)):
            for j in range(0,len(sysNames)):
                if i == j:
                    continue
                sysName_A = sysNames[i]
                sysName_B = sysNames[j]
                colName = f"{sysName_A} - {sysName_B}"
                colNames.append(colName)
                delta_AB = meanDf[sysName_A] - meanDf[sysName_B]
                dfsToConcat.append(delta_AB)    

    deltaDf = pd.concat(dfsToConcat, axis=1)
    deltaDf.columns = colNames  
    deltaCsv = p.join(analDir, "delta_RMSF.csv")
    deltaDf.to_csv(deltaCsv)
    #### peak detection
    rmsfPeaks = {}
    for col in deltaDf:
        absCol = deltaDf[col].abs()
        peaks, _ = find_peaks(absCol, distance = 3, height = 0.4)
        rmsfPeaks[col] = peaks.tolist()

    rmsfPeaksYaml = p.join(analDir, "delta_RMSF_peaks.yaml")
    with open(rmsfPeaksYaml,"w") as file:
        yaml.dump(rmsfPeaks,file)

#############################################################################################
def check_RMSD(traj, outDir, inputDirName):
    rmsd = md.rmsd(traj, traj, 0)
    rmsdDf = pd.DataFrame(rmsd, columns = [f"RMSD_{inputDirName}"])
    rmsdCsv = p.join(outDir, f"RMSD_{inputDirName}.csv")
    rmsdDf.to_csv(rmsdCsv)
#############################################################################################
def check_RMSF(traj, pdbDf, outDir,inputDirName):
    rmsf = md.rmsf(traj, traj, 0)
    rmsfDf = pd.DataFrame(rmsf, columns = ["RMSF"])
    pdbDf["RMSF"] = rmsfDf
    perResidueRMSF = pdbDf.groupby("RES_ID")["RMSF"].mean()
    perResidueRMSF = perResidueRMSF.to_frame()
    perResidueRMSF.columns = [f"RMSF_{inputDirName}"]
    rmsfCsv = p.join(outDir, f"RMSF_{inputDirName}.csv")
    perResidueRMSF.to_csv(rmsfCsv)
#############################################################################################
def init_atoms_of_interest():
    interestingAtoms = {
        "ALA": ["CB"],
        "LEU": ["CB", "CG", "CD1", "CD2"],
        "GLY": [],
        "ILE": ["CB", "CG1", "CG2", "CD1"],
        "VAL": ["CB", "CG1", "CG2"],
        "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
        "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
        "ASP": ["CB", "CG", "OD1", "OD2"],
        "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
        "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
        "HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
        "LYS": ["CB", "CG", "CD", "CE", "NZ"],
        "SER": ["CB", "OG"],
        "THR": ["CB", "OG1", "CG2"],
        "ASN": ["CB", "CG", "OD1", "ND2"],
        "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
        "CYS": ["CB", "SG"],
        "MET": ["CB", "CG", "SD", "CE"],
        "PRO": ["CB", "CG", "CD"]
    }
    return interestingAtoms
#############################################################################################
def compute_contact_distances(traj, residuePairs, sysAnalDir, inputDirName):
    print("---->\tCalculating contact distances!")
    for resTag in residuePairs:
        print(f"-------->{resTag}")
        print(resTag)
        print(residuePairs[resTag].keys())
        # get contacts and labels from residuePairs dict
        pairwiseContacts = residuePairs[resTag]["contacts"]
        plotLabels = residuePairs[resTag]["labels"]
        ## calculate contact distances
        contactDistances, residueIds = md.compute_contacts(traj, pairwiseContacts)
        contactIds = residueIds[:, 1].tolist()
        contactIds = [x+1 for x in contactIds] # un-Zero index
        contactDf = pd.DataFrame(contactDistances, columns = plotLabels)
        outCsv = p.join(sysAnalDir, f"contacts_{resTag}_{inputDirName}.csv")
        contactDf.to_csv(outCsv)
#############################################################################################
def find_pairwise_contacts(traj, pdbDf, keyResidues):
    print("---->\t Finding pairwise contacts!")
    interestingAtoms = init_atoms_of_interest()
    residuePairs = {}
    for resTag in keyResidues:
        print(f"--------> {resTag}")

        ## find nearby residues 
        targetResId = int(keyResidues[resTag]["RES_ID"])
        targetResName = keyResidues[resTag]["RES_NAME"]
        resName = keyResidues[resTag]["RES_NAME"]
        resDf = pdbDf[(pdbDf["RES_ID"] == targetResId) &
                       (pdbDf["ATOM_NAME"].isin(interestingAtoms[resName])) &
                       (pdbDf["RES_NAME"] == targetResName)]
        if len(resDf) == 0:
            print(f"--------> No {resTag} in pdb")
            continue
        queryIndecies = resDf["ATOM_ID"].values -  1 ## zero-index
        neighbors = md.compute_neighbors(traj, 0.5, queryIndecies)
        uniqueNeighborAtoms = list(set([item + 1 for sublist in neighbors for item in sublist])) ## one-index!
        # exclue query residue and 2 residues adjacent
        exclueResidues = range(targetResId -2, targetResId +3)#[resId - 2, resId, resId + 2]
        neighborDf = pdbDf[(pdbDf["ATOM_ID"].isin(uniqueNeighborAtoms)) & (~pdbDf["RES_ID"].isin(exclueResidues))]
        nearbyResidues = neighborDf["RES_ID"].unique()
        # drop query residue from neighbors
        neighborDf = neighborDf.copy()
        neighborDf.loc[:,"plotLabels"] = neighborDf["RES_NAME"] + "_" + neighborDf["RES_ID"].astype(str)
        plotLabels = neighborDf["plotLabels"].unique().tolist()
        ## compute contacts, returns a list of [(X, Y), (X, Z)] contacts
        pairwiseContacts = [(targetResId -1, targetIndex - 1) for targetIndex in nearbyResidues] ## zero-index!
        residuePairs.update({resTag:{"contacts": pairwiseContacts,
                                             "labels": plotLabels}})
    return residuePairs