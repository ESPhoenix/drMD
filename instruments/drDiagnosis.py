import os
from os import path as p
import pandas as pd
import mdtraj as md
from scipy.signal import find_peaks
import yaml
import numpy as np

#############################################################################################
def get_residue_id_mapping(traj, pdbDf):
    trajRes = [residue.index for residue in traj.topology.residues]


    uniqueChains = pdbDf["CHAIN_ID"].unique().tolist()
    pdbRes = []
    chainIds = []
    for uniqueChain in uniqueChains:
        chainDf = pdbDf[pdbDf["CHAIN_ID"] == uniqueChain]
        chainRes = chainDf["RES_ID"].unique().tolist()
        pdbRes += chainRes
        chainIds += [uniqueChain for _ in chainRes]

    resMapping = {trajResi : {"RES_ID":pdbResi,
                              "CHAIN_ID":chainId}  
                              for trajResi, pdbResi, chainId  in zip(trajRes, pdbRes, chainIds) }


    return resMapping
#############################################################################################
def pdbResToTrajRes(resMap,  resId, chainId):
    return next((key for key, value in resMap.items() if value['RES_ID'] == resId and value['CHAIN_ID'] == chainId), None)


#############################################################################################
def compute_radial_distribution(traj, residuePairs, contactDf, sysAnalDir, inputDirName):
    print("---->\tCalculating radial basis function!")
    for resTag in residuePairs:
        print(f"--------> RDF for {resTag}")
        # unpack contacts from residuePairs
        pairwiseContacts = residuePairs[resTag]["contacts"]
        colNames = residuePairs[resTag]["labels"]
        # calculate radial distribution function for pairs
        radii, rdf = md.compute_rdf(traj, pairwiseContacts, bin_width = 0.01)
        rdfDf = pd.DataFrame({'Radii': radii, 'RDF': rdf})
        # peak detection with a cutoff of 500
        rdfPeaks, _ = find_peaks(rdfDf["RDF"], distance = 3)
        rdfPeaks = rdfPeaks.tolist()
        rdfInteractions = {}
        for rdfPeak in rdfPeaks:
            interactingResidues = []
            for colName in contactDf:
                peakDistance = rdfDf.loc[rdfPeak,"Radii"]
                # Check if this pair often has a distance close to the peak distance
                close_to_peak = np.isclose(contactDf[colName],peakDistance , atol=0.1)
                if close_to_peak.mean() > 0.30:
                    interactingResidues.append(colName)
            rdfInteractions[str(peakDistance)] = interactingResidues
        ## save rdf as csv 
        rdfCsv = p.join(sysAnalDir,f"RDF_{resTag}_{inputDirName}.csv")
        rdfDf.to_csv(rdfCsv)
        ## save interactions as yaml
        rdfYaml = p.join(sysAnalDir, f"RDF_interactions_{resTag}_{inputDirName}.yaml")
        with open(rdfYaml,"w") as file:
            yaml.dump(rdfInteractions,file, default_flow_style=False)
        return rdfDf
#############################################################################################
def compute_delta_RMSF(analDir, referenceSystem = False):
    print("---->\tCalculating delta RMSF!")
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
        yaml.dump(rmsfPeaks,file,default_flow_style=False)

#############################################################################################
def check_RMSD(traj, outDir, inputDirName):
    print("---->\tCalculating RMSD!")

    rmsd = md.rmsd(traj, traj, 0)
    rmsdDf = pd.DataFrame(rmsd, columns = [f"RMSD_{inputDirName}"])
    rmsdCsv = p.join(outDir, f"RMSD_{inputDirName}.csv")
    rmsdDf.to_csv(rmsdCsv)
#############################################################################################
def check_RMSF(traj, pdbDf, outDir,inputDirName):
    print("---->\tCalculating per residue RMSF!")

    rmsf = md.rmsf(traj, traj, 0)
    rmsfDf = pd.DataFrame(rmsf, columns = ["RMSF"])
    pdbDf["RMSF"] = rmsfDf
    perResidueRMSF = pdbDf.groupby("RES_ID")["RMSF"].mean()
    perResidueRMSF = perResidueRMSF.to_frame()
    perResidueRMSF.columns = [f"RMSF_{inputDirName}"]
    rmsfCsv = p.join(outDir, f"RMSF_{inputDirName}.csv")
    perResidueRMSF.to_csv(rmsfCsv)
#############################################################################################
def compute_atomic_distances(traj, atomPairs, sysAnalDir, inputDirName, tag):
    print("---->\tCalculating atomic distances!")
    distanceDfs = {}
    for atomTag in atomPairs:
        print(f"-------->{tag} : {atomTag}")
        # get contacts and labels from residuePairs dict
        pairwiseContacts = atomPairs[atomTag]["contacts"]
        if len(pairwiseContacts) == 0:
            continue
        plotLabels = atomPairs[atomTag]["labels"]
        ## calculate contact distances
        contactDistances  = md.compute_distances(traj, pairwiseContacts)
        distanceDf = pd.DataFrame(contactDistances, columns = plotLabels)
        outCsv = p.join(sysAnalDir, f"{tag}_{atomTag}_{inputDirName}.csv")
        distanceDf.to_csv(outCsv)
        distanceDfs.update({atomTag: distanceDfs})
    return distanceDfs
#############################################################################################
def compute_contact_distances(traj, residuePairs, sysAnalDir, inputDirName):
    print("---->\tCalculating contact distances!")
    contactDfs = {}
    for resTag in residuePairs:
        print(f"-------->{resTag}")
        # get contacts and labels from residuePairs dict
        pairwiseContacts = residuePairs[resTag]["contacts"]
        if len(pairwiseContacts) == 0:
            continue
        plotLabels = residuePairs[resTag]["labels"]

        ## calculate contact distances
        contactDistances, residueIds = md.compute_contacts(traj, pairwiseContacts)
        # contactIds = residueIds[:, 1].tolist()
        # contactIds = [x+1 for x in contactIds] # un-Zero index
        contactDf = pd.DataFrame(contactDistances, columns = plotLabels)
        outCsv = p.join(sysAnalDir, f"contacts_{resTag}_{inputDirName}.csv")
        contactDf.to_csv(outCsv)
        contactDfs.update({resTag:contactDf})
    return contactDfs
#############################################################################################
def find_pairwise_residue_contacts(traj, pdbDf, keyResidues):
    print("---->\t Finding pairwise contacts!")

    residueMapping = get_residue_id_mapping(traj, pdbDf)

    interestingAtoms, _, _ = init_atoms_of_interest()
    residuePairs = {}
    for resTag in keyResidues:
        print(f"--------> {resTag}")
        ## find nearby residues 
        targetResId = int(keyResidues[resTag]["RES_ID"])
        targetResName = keyResidues[resTag]["RES_NAME"]
        targetResChain = keyResidues[resTag]["CHAIN_ID"]

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

        neighborResAndChains= neighborDf[["RES_ID", "CHAIN_ID"]].drop_duplicates()
        neighborResidues = neighborResAndChains["RES_ID"].tolist()
        neighborChains = neighborResAndChains["CHAIN_ID"].tolist()

        ## compute contacts, returns a list of [(X, Y), (X, Z)] contacts
        ## translate pdb residue IDs to mdTraj residue Ids

        targetResId_traj = pdbResToTrajRes(residueMapping, targetResId, targetResChain)
        nearbyResIds_traj = [pdbResToTrajRes(residueMapping, neighborResId, neighborChain)
                              for neighborResId, neighborChain in zip(neighborResidues, neighborChains)]


        pairwiseContacts = [(targetResId_traj, nearbyIndex) for nearbyIndex in nearbyResIds_traj] 
 
        
        # make plot labels
        neighborDf = neighborDf.copy()
        neighborDf.loc[:,"plotLabels"] = resTag + "--" + neighborDf["RES_NAME"] + ":" + neighborDf["RES_ID"].astype(str)
        plotLabels = neighborDf["plotLabels"].unique().tolist()

        residuePairs.update({resTag:{"contacts": pairwiseContacts,
                                             "labels": plotLabels}})
        

    return residuePairs
#############################################################################################
def find_hydrogen_bonds(traj, pdbDf, keyResidues):
    print("---->\t Finding hydrogen bonds!")
    _, hydrogenBondDonors, hydrogenBondAcceptors = init_atoms_of_interest()
    for resTag in keyResidues:
        print(f"--------> {resTag}")
        ## get dataframes with residue of interest's hydrogen bond donors and acceptors
        targetResId = int(keyResidues[resTag]["RES_ID"])
        targetResName = keyResidues[resTag]["RES_NAME"]
        resName = keyResidues[resTag]["RES_NAME"]
        hydrogenDonorDf = pdbDf[(pdbDf["RES_ID"] == targetResId) &
                                (pdbDf["ATOM_NAME"].isin(hydrogenBondDonors[resName])) &
                                (pdbDf["RES_NAME"] == targetResName)]
        
        hydrogenAcceptorDf = pdbDf[(pdbDf["RES_ID"] == targetResId) &
                                    (pdbDf["ATOM_NAME"].isin(hydrogenBondAcceptors[resName])) &
                                    (pdbDf["RES_NAME"] == targetResName)]
        
        ##########################################################
        ## find neighbors for HBA's and HBD's in target residue
        donorAcceptorPairs = {}
        if len(hydrogenDonorDf) > 0:
            for idx, donorRow in hydrogenDonorDf.iterrows():
                ## compute neighboring atoms - get a unique list of them - produce a df 
                queryIndex = np.array([donorRow["ATOM_ID"] -  1]) ## zero-index
                neighbors = md.compute_neighbors(traj, 0.5, queryIndex)
                uniqueNeighborAtoms = list(set([item + 1 for sublist in neighbors for item in sublist])) ## one-index!
                neighborDf = pdbDf[(pdbDf["ATOM_ID"].isin(uniqueNeighborAtoms))]

                ## drop non-hydrogen bond acceptor atoms
                idxToDrop = []
                for idx, acceptorRow in neighborDf.iterrows():
                    if not acceptorRow["ATOM_NAME"] in hydrogenBondAcceptors[acceptorRow["RES_NAME"]]:
                        idxToDrop.append(idx)
                acceptorDf = neighborDf.drop(index=idxToDrop)

                ## make plot labels
                targetAtomTag = "-".join([resTag, donorRow["ATOM_NAME"]])
                acceptorDf["plotLabels"] =  acceptorDf["CHAIN_ID"] +":"+ acceptorDf["RES_NAME"] +":"+ acceptorDf["RES_ID"].astype(str) +":"+ acceptorDf["ATOM_NAME"]
                acceptorDf["plotLabels"] = targetAtomTag + " -- H -- " + acceptorDf["plotLabels"]
                plotLabels = acceptorDf["plotLabels"].unique().tolist()

                ## create outputs to feed to compute_distances, returns a list of [(X, Y), (X, Z)] contacts
                pairwiseContacts = [(donorRow["ATOM_ID"] -  1, acceptorIndex - 1) for acceptorIndex in acceptorDf["ATOM_ID"].to_list()] ## zero-index!
                donorAcceptorPairs.update({targetAtomTag:{"contacts": pairwiseContacts,
                                                    "labels": plotLabels}})
                
        ##########################################################
        ## find neighbors for HBA's and HBD's in target residue
        acceptorDonorPairs = {}
        if len(hydrogenAcceptorDf) > 0:
            for idx, acceptorRow in hydrogenAcceptorDf.iterrows():
                queryIndex = np.array([acceptorRow["ATOM_ID"] -  1]) ## zero-index
                neighbors = md.compute_neighbors(traj, 0.5, queryIndex)
                uniqueNeighborAtoms = list(set([item + 1 for sublist in neighbors for item in sublist])) ## one-index!
                neighborDf = pdbDf[(pdbDf["ATOM_ID"].isin(uniqueNeighborAtoms))]
                ## drop non-hydrogen bond acceptor atoms
                idxToDrop = []
                for idx, donorRow in neighborDf.iterrows():
                    if not donorRow["ATOM_NAME"] in hydrogenBondDonors[donorRow["RES_NAME"]]:
                        idxToDrop.append(idx)
                donorDf = neighborDf.drop(index=idxToDrop)
                ## make plot labels
                targetAtomTag = "_".join([resTag, acceptorRow["ATOM_NAME"]])
                donorDf["plotLabels"] =  donorDf["CHAIN_ID"] +":"+ donorDf["RES_NAME"] +":"+ donorDf["RES_ID"].astype(str) +":"+ donorDf["ATOM_NAME"]
                donorDf["plotLabels"] = targetAtomTag + " -- H -- " + donorDf["plotLabels"]

                plotLabels = donorDf["plotLabels"].unique().tolist()
                ## compute contacts, returns a list of [(X, Y), (X, Z)] contacts
                pairwiseContacts = [(acceptorRow["ATOM_ID"] -  1, donorIndex - 1) for donorIndex in donorDf["ATOM_ID"].to_list()] ## zero-index!
                acceptorDonorPairs.update({targetAtomTag:{"contacts": pairwiseContacts,
                                                    "labels": plotLabels}})

    return donorAcceptorPairs, acceptorDonorPairs


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

    hydrogenBondDonors = {
            "ALA": [],
        "LEU": ["N"],
        "GLY": ["N"],
        "ILE": ["N"],
        "VAL": ["N"],
        "PHE": ["N"],
        "TRP": ["N","NE1"],
        "TYR": ["N","OH"],
        "ASP": ["N"],
        "GLU": ["N"],
        "ARG": ["N","NH1", "NH2"],
        "HIS": ["N","ND1", "NE2"],
        "LYS": ["N","NZ"],
        "SER": ["N", "OG"],
        "THR": ["N","OG1"],
        "ASN": ["N","ND2"],
        "GLN": ["N","NE2"],
        "CYS": ["N","SG"],
        "MET": ["N"],
        "PRO": ["N"]
    }

    hydrogenBondAcceptors = {
        "ALA": ["O"],
        "LEU": ["O"],
        "GLY": ["O"],
        "ILE": ["O"],
        "VAL": ["O"],
        "PHE": ["O"],
        "TRP": ["O"],
        "TYR": ["O"],
        "ASP": ["O", "OD1", "OD2"],
        "GLU": ["O", "OE1", "OE2"],
        "ARG": ["O"],
        "HIS": ["O", "ND1", "NE2"],
        "LYS": ["O"],
        "SER": ["O", "OG"],
        "THR": ["O", "OG1"],
        "ASN": ["O", "OD1"],
        "GLN": ["O", "OE1"],
        "CYS": ["O", "SG"],
        "MET": ["O", "SD"],
        "PRO": ["O"]
    }

    return interestingAtoms, hydrogenBondDonors, hydrogenBondAcceptors