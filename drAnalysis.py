import os
from os import path as p
import pandas as pd
import mdtraj as md

import numpy as np
from scipy.spatial import distance
import argpass
import yaml
import math
### drMD modules ###
from module_pdbUtils import pdb2df
import module_drPlot as drPlot
#############################################################################################
def check_RMSD(traj, outDir):
    rmsd = md.rmsd(traj, traj, 0)
    rmsdDf = pd.DataFrame(rmsd, columns = ["RMSD"])
    drPlot.plot_rmsd(rmsdDf, outDir)
    return rmsdDf
#############################################################################################
def check_RMSF(traj, pdbDf, outDir):
    rmsf = md.rmsf(traj, traj, 0)
    rmsfDf = pd.DataFrame(rmsf, columns = ["RMSF"])
    pdbDf["RMSF"] = rmsfDf
    perResidueRMSF = pdbDf.groupby("RES_ID")["RMSF"].mean()
    perResidueRMSF = perResidueRMSF.to_frame()
    perResidueRMSF.columns = ["RMSF"]
    print(perResidueRMSF)
    drPlot.plot_rmsf(perResidueRMSF,  outDir)
    return rmsfDf

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
def find_contacts(traj, pdbDf, keyResidues, outDir):
    interestingAtoms = init_atoms_of_interest()
   # contactsPerTarget = find_nearby_residues(keyResidues, pdbFile)
    for resTag in keyResidues:
        print(f"-->{resTag}")
        ## find nearby residues 
        resId = keyResidues[resTag]["RES_ID"] 
        resName = keyResidues[resTag]["RES_NAME"]
        resDf = pdbDf[(pdbDf["RES_ID"] == resId) & (pdbDf["ATOM_NAME"].isin(interestingAtoms[resName]))]
        queryIndecies = resDf["ATOM_ID"].values -  1 ## zero-index
        neighbors = md.compute_neighbors(traj, 0.5, queryIndecies)
        uniqueNeighborAtoms = list(set([item + 1 for sublist in neighbors for item in sublist])) ## one-index!
        # exclue query residue and 2 residues adjacent
        exclueResidues = range(resId -2, resId +3)#[resId - 1, resId, resId + 1]
        neighborDf = pdbDf[(pdbDf["ATOM_ID"].isin(uniqueNeighborAtoms)) & (~pdbDf["RES_ID"].isin(exclueResidues))]
        nearbyResidues = neighborDf["RES_ID"].unique()
        # drop query residue from neighbors
        neighborDf = neighborDf.copy()
        neighborDf.loc[:,"plotLabels"] = neighborDf["RES_NAME"] + "_" + neighborDf["RES_ID"].astype(str)
        plotLabels = neighborDf["plotLabels"].unique().tolist()
        ## compute contacts
        pairwiseContacts = [(resId -1, targetIndex - 1) for targetIndex in nearbyResidues] ## zero-index!
        contactDistances, residuePairs = md.compute_contacts(traj, pairwiseContacts)
        contactIds = residuePairs[:, 1].tolist()
        contactIds = [x+1 for x in contactIds] # un-Zero index
        contactDf = pd.DataFrame(contactDistances, columns = plotLabels)
        drPlot.plot_distance_hist(contactDf, outDir, resTag, plotLabels)        
        continue
        ### OFF FOR NOW - DON'T UNDERSTAND RDFs
        ## calculate radial distribution function for pairs
        # radii, rdf = md.compute_rdf(traj, pairwiseContacts, bin_width = 0.05)
        # rdfDf = pd.DataFrame({'Radii': radii, 'RDF': rdf})
        # # peak detection with a cutoff of 500
        # peak_distances = [rdfDf['Radii'][i] for i in range(len(rdfDf['RDF'])) if rdfDf['RDF'][i] > 500]
        # for i, pair in enumerate(pairwiseContacts):
        #     pair_distances = contactDistances[i, :]
        #     for peak_distance in peak_distances:
        #         # Check if this pair often has a distance close to the peak distance
        #         close_to_peak = np.isclose(pair_distances, peak_distance, atol=0.01)
        #         if close_to_peak.mean() > 0.05:
        #             print(f"Pair {pair} is associated with peak at distance {round(peak_distance,4)}")
        # plot_rdf(rdfDf, outDir, tag = resTag)
#####################################################################################################
def read_inputs():
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    config=args.config
    ## Read config.yaml into a dictionary
    with open(config,"r") as yamlFile:
        config = yaml.safe_load(yamlFile) 
    return config

#############################################################################################
def main():
    # read and unpack config file
    config = read_inputs()
    pathInfo = config["pathInfo"]
    analDir = pathInfo["analDir"]
    analMenu = config["analysisMenu"]
    os.makedirs(analDir,exist_ok=True)
    keyResidues = config["keyResidues"]

    # reconstruct file structure from config
    mdDir = pathInfo["mdDir"]
    stepName = pathInfo["stepName"]
    repeats = pathInfo["repeats"]
    for sysTag in repeats:
        sysAnalDir = p.join(analDir, sysTag)
        os.makedirs(sysAnalDir,exist_ok=True)

        for inputDirName in repeats[sysTag]:
            simDir = p.join(mdDir,inputDirName,stepName)
            analysis_protocol(simDir, analMenu, keyResidues, sysAnalDir, inputDirName)
        plotting_protocol(analMenu, sysAnalDir, keyResidues)

        
def plotting_protocol(analMenu, sysAnalDir, keyResidues):
    ############ read analysis menu and do  plotting that has been ordered ############
    keyResiAnal = analMenu["keyResidueAnalysis"]
    wholeAnal = analMenu["wholeProteinAnalysis"]
    if keyResiAnal["contactDistances"]:
        for resTag in keyResidues: 
            drPlot.plot_distance_hist(sysAnalDir, resTag)



########################################################################
def analysis_protocol(simDir, analMenu, keyResidues, sysAnalDir, inputDirName):
    print(f"--> {inputDirName}")
    ## find files in simDir
    pdbFile, dcdFile = False, False
    for file in os.listdir(simDir):
        fileData = p.splitext(file)
        if fileData[1] == ".pdb":
            pdbFile = p.join(simDir,file)
        elif fileData[1] == ".dcd":
            dcdFile = p.join(simDir, file)

    ## skip if files not found
    if not pdbFile:
        print(f"-->\tNo PDB file found in {simDir}! EXITING")
        return
    elif not dcdFile:
        print(f"-->\tNo DCD file found in {simDir}! EXITING")
        return    

    ## load trajectory | remove solvent molecules
    traj = md.load_dcd(dcdFile, top = pdbFile)
    traj.remove_solvent([],True)
    ## load pdb file as a dataframe | remove solvent and ions
    pdbDf = pdb2df(pdbFile)
    pdbDf = pdbDf[~pdbDf["RES_NAME"].isin(["Cl-","Na+","HOH"])].copy()

    ############ read analysis menu and do ordered analysis ############
    keyResiAnal = analMenu["keyResidueAnalysis"]
    wholeAnal = analMenu["wholeProteinAnalysis"]
    if any(job for job in keyResiAnal.values() if job):
        residuePairs = find_pairwise_contacts(traj, pdbDf, keyResidues)
        if keyResiAnal["contactDistances"]:
            compute_contact_distances(traj, residuePairs, sysAnalDir, inputDirName)

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

#############################################################################################
if  __name__ == "__main__":
    main()