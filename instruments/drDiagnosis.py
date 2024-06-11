import os
from os import path as p
import pandas as pd
import mdtraj as md
from scipy.signal import find_peaks
import yaml
import numpy as np
from typing import List, Dict, Union, Any, Tuple
#############################################################################################
def get_residue_id_mapping(
    traj: md.Trajectory,
    pdb_df: pd.DataFrame
) -> dict:
    """
    Get a mapping between residue indices in the trajectory and residue IDs and chain IDs in the PDB file.

    Args:
        traj (mdtraj.Trajectory): The trajectory object containing the residue indices.
        pdb_df (pandas.DataFrame): The DataFrame containing the residue IDs and chain IDs in the PDB file.

    Returns:
        dict: A dictionary mapping residue indices in the trajectory to a dictionary containing the residue ID and chain ID.
    """
    # Get the residue indices in the trajectory
    traj_res: List[int] = [residue.index for residue in traj.topology.residues]

    # Get the unique chain IDs in the PDB file
    unique_chains: List[str] = pdb_df["CHAIN_ID"].unique().tolist()

    # Initialize lists to store the residue IDs and chain IDs
    pdb_res: List[int] = []
    chain_ids: List[str] = []

    # Iterate over each unique chain ID
    for unique_chain in unique_chains:
        # Get the DataFrame containing the residue IDs and chain IDs for the current chain
        chain_df: pd.DataFrame = pdb_df[pdb_df["CHAIN_ID"] == unique_chain]
        # Get the unique residue IDs for the current chain
        chain_res: List[int] = chain_df["RES_ID"].unique().tolist()
        # Add the residue IDs to the list
        pdb_res += chain_res
        # Add the chain IDs to the list
        chain_ids += [unique_chain for _ in chain_res]

    # Create a dictionary mapping residue indices in the trajectory to a dictionary containing the residue ID and chain ID
    res_mapping: Dict[int, Dict[str, Union[int, str]]] = {
        traj_resi: {"RES_ID": pdb_resi, "CHAIN_ID": chain_id}
        for traj_resi, pdb_resi, chain_id in zip(traj_res, pdb_res, chain_ids)
    }

    return res_mapping
#############################################################################################
def pdbResToTrajRes(resMap: Dict[int, Dict[str, Union[int, str]]], resId: int, chainId: str) -> Optional[int]:
    """
    Translates a residue ID and chain ID from a PDB file to its corresponding residue ID in a trajectory.

    Args:
        resMap (Dict[int, Dict[str, Union[int, str]]]): A dictionary mapping residue indices in the trajectory to a dictionary containing the residue ID and chain ID.
        resId (int): The residue ID from the PDB file.
        chainId (str): The chain ID from the PDB file.

    Returns:
        Optional[int]: The residue ID in the trajectory corresponding to the input residue ID and chain ID, or None if no match is found.
    """
    # Find the key in resMap whose value matches the input residue ID and chain ID
    return next((key for key, value in resMap.items() if value['RES_ID'] == resId and value['CHAIN_ID'] == chainId), None)



#############################################################################################
def compute_radial_distribution(
    traj: md.Trajectory,
    residuePairs: Dict[str, Dict[str, Any]],
    contactDf: pd.DataFrame,
    sysAnalDir: str,
    inputDirName: str,
) -> pd.DataFrame:
    """
    Compute the radial distribution function (RDF) for each pair of residues in the given trajectory
    and save the results as CSV and YAML files.

    Args:
        traj (mdtraj.Trajectory): The trajectory containing the residues.
        residuePairs (Dict[str, Dict[str, Any]]): A dictionary containing the pairs of residues for which to compute the RDF.
        contactDf (pd.DataFrame): A DataFrame containing the contacts between residues.
        sysAnalDir (str): The directory where the analysis files will be saved.
        inputDirName (str): The name of the input directory.

    Returns:
        pd.DataFrame: A DataFrame containing the RDF values for each pair of residues.
    """
    print("---->\tCalculating radial basis function!")

    # Iterate over each pair of residues
    for resTag in residuePairs:
        print(f"--------> RDF for {resTag}")

        # Unpack contacts from residuePairs
        pairwiseContacts: List[Tuple[int, int]] = residuePairs[resTag]["contacts"]
        colNames: List[str] = residuePairs[resTag]["labels"]

        # Calculate the RDF for the pair of residues
        radii, rdf = md.compute_rdf(traj, pairwiseContacts, bin_width=0.01)
        rdfDf: pd.DataFrame = pd.DataFrame({'Radii': radii, 'RDF': rdf})

        # Find peaks in the RDF
        rdfPeaks: np.ndarray = find_peaks(rdfDf["RDF"], distance=3)[0]

        # Find the residues that interact with the peak distance
        rdfInteractions: Dict[str, List[str]] = {}
        for rdfPeak in rdfPeaks:
            peakDistance: float = rdfDf.loc[rdfPeak, "Radii"]
            interactingResidues: List[str] = []
            for colName in contactDf:
                close_to_peak: np.ndarray = np.isclose(contactDf[colName], peakDistance, atol=0.1)
                if close_to_peak.mean() > 0.30:
                    interactingResidues.append(colName)
            rdfInteractions[str(peakDistance)] = interactingResidues

        # Save the RDF as a CSV file
        rdfCsv: str = p.join(sysAnalDir, f"RDF_{resTag}_{inputDirName}.csv")
        rdfDf.to_csv(rdfCsv)

        # Save the interactions as a YAML file
        rdfYaml: str = p.join(sysAnalDir, f"RDF_interactions_{resTag}_{inputDirName}.yaml")
        with open(rdfYaml, "w") as file:
            yaml.dump(rdfInteractions, file, default_flow_style=False)

    return rdfDf
#############################################################################################


def compute_delta_RMSF(analDir: str, referenceSystem: bool = False) -> None:
    """
    Calculate the delta RMSF for each pair of systems in the analysis directory.

    Args:
        analDir (str): The directory containing the analysis files.
        referenceSystem (bool, optional): If True, calculate the delta RMSF between
            each system and a reference system. Default is False.

    Returns:
        None
    """
    print("---->\tCalculating delta RMSF!")

    # Initialize lists to store dataframes and column names
    resultDfsToConcat: List[pd.DataFrame] = []
    colNames: List[str] = []
    dfsToConcat: List[pd.DataFrame] = []

    # Get the names of the systems in the analysis directory
    sysNames: List[str] = [
        name
        for name in os.listdir(analDir)
        if p.isdir(p.join(analDir, name))
    ]

    # Calculate the delta RMSF for each pair of systems
    for sysName in sysNames:
        sysAnalDir: str = p.join(analDir, sysName)

        # Load contact dfs from sysAnalDir | concat into one big df
        dfsToConcat = []
        for file in os.listdir(sysAnalDir):
            if not p.splitext(file)[1] == ".csv":
                continue
            if file.startswith("RMSF"):
                runDf: pd.DataFrame = pd.read_csv(
                    p.join(sysAnalDir, file), index_col="RES_ID"
                )
                dfsToConcat.append(runDf)
        if len(dfsToConcat) == 0:
            return
        df: pd.DataFrame = pd.concat(dfsToConcat, axis=1, ignore_index=False)
        meanDf: pd.DataFrame = df.mean(axis=1)
        meanDf.name = sysName
        resultDfsToConcat.append(meanDf)

    # Concatenate the mean RMSF values for each system
    meanDf: pd.DataFrame = pd.concat(resultDfsToConcat, axis=1)

    # Calculate the delta RMSF for each pair of systems
    if referenceSystem:
        for i in range(len(sysNames)):
            sysName_A: str = sysNames[i]
            sysName_B: str = referenceSystem
            if sysName_A == sysName_B:
                continue
            colName: str = f"{sysName_A} - {sysName_B}"
            colNames.append(colName)
            delta_AB: pd.DataFrame = meanDf[sysName_A] - meanDf[sysName_B]
            dfsToConcat.append(delta_AB)
    else:
        for i in range(len(sysNames)):
            for j in range(len(sysNames)):
                if i == j:
                    continue
                sysName_A: str = sysNames[i]
                sysName_B: str = sysNames[j]
                colName: str = f"{sysName_A} - {sysName_B}"
                colNames.append(colName)
                delta_AB: pd.DataFrame = meanDf[sysName_A] - meanDf[sysName_B]
                dfsToConcat.append(delta_AB)

    # Concatenate the delta RMSF values for each pair of systems
    deltaDf: pd.DataFrame = pd.concat(dfsToConcat, axis=1)
    deltaDf.columns = colNames

    # Save the delta RMSF as a CSV file
    deltaCsv: str = p.join(analDir, "delta_RMSF.csv")
    deltaDf.to_csv(deltaCsv)

    # Perform peak detection on the delta RMSF
    rmsfPeaks: Dict[str, List[int]] = {}
    for col in deltaDf:
        absCol: pd.Series = deltaDf[col].abs()
        peaks, _ = find_peaks(absCol, distance=3, height=0.4)
        rmsfPeaks[col] = peaks.tolist()

    # Save the peak information as a YAML file
    rmsfPeaksYaml: str = p.join(analDir, "delta_RMSF_peaks.yaml")
    with open(rmsfPeaksYaml, "w") as file:
        yaml.dump(rmsfPeaks, file, default_flow_style=False)

#############################################################################################
def check_RMSD(traj: md.Trajectory, outDir: str, inputDirName: str) -> None:
    """
    Calculate the root mean square deviation (RMSD) for a given trajectory and save it as a CSV file.

    Args:
        traj (mdtraj.Trajectory): The trajectory for which to calculate the RMSD.
        outDir (str): The directory where the CSV file will be saved.
        inputDirName (str): The name of the input directory.

    Returns:
        None
    """
    # Print a message indicating that RMSD calculation has started
    print("---->\tCalculating RMSD!")

    # Calculate the RMSD for the trajectory
    rmsd: np.ndarray = md.rmsd(traj, traj, 0)

    # Create a DataFrame with the RMSD values
    rmsdDf: pd.DataFrame = pd.DataFrame(rmsd, columns=[f"RMSD_{inputDirName}"])

    # Create the file path for the CSV file
    rmsdCsv: str = p.join(outDir, f"RMSD_{inputDirName}.csv")

    # Save the RMSD DataFrame as a CSV file
    rmsdDf.to_csv(rmsdCsv)
#############################################################################################
def check_RMSF(traj: md.Trajectory, pdbDf: pd.DataFrame, outDir: str ,inputDirName: str):
    """
    Calculate the Root Mean Square Fluctuations (RMSF) for each residue in the
    given trajectory and save it as a CSV file.

    Args:
        traj (mdtraj.Trajectory): The trajectory for which to calculate the RMSF.
        pdbDf (pandas.DataFrame): The DataFrame containing information about the
            residues in the PDB file.
        outDir (str): The directory where the CSV file will be saved.
        inputDirName (str): The name of the input directory.
    """
    # Print a message indicating that RMSF calculation has started
    print("---->\tCalculating per residue RMSF!")

    # Calculate the RMSF for each residue in the trajectory
    rmsf: md.RMSF = md.rmsf(traj, traj, 0)

    # Create a DataFrame with the RMSF values
    rmsfDf: pd.DataFrame = pd.DataFrame(rmsf, columns=["RMSF"])

    # Add the RMSF values to the DataFrame containing information about the residues
    pdbDf["RMSF"] = rmsfDf

    # Calculate the mean RMSF for each residue
    perResidueRMSF: pd.Series = pdbDf.groupby("RES_ID")["RMSF"].mean()

    # Convert the mean RMSF values to a DataFrame
    perResidueRMSF: pd.DataFrame = perResidueRMSF.to_frame()

    # Rename the columns of the DataFrame
    perResidueRMSF.columns  = [f"RMSF_{inputDirName}"]

    # Create the file path for the CSV file
    rmsfCsv: str = p.join(outDir, f"RMSF_{inputDirName}.csv")

    # Save the DataFrame as a CSV file
    perResidueRMSF.to_csv(rmsfCsv)
#############################################################################################
def compute_atomic_distances(
    traj: md.Trajectory,
    atomPairs: Dict[str, Dict[str, Any]],
    sysAnalDir: str,
    inputDirName: str,
    tag: str,
) -> Dict[str, pd.DataFrame]:
    """
    Calculate atomic distances for each pair of atoms in the given trajectory and
    save the results as CSV files.

    Args:
        traj (mdtraj.Trajectory): The trajectory containing the atoms.
        atomPairs (Dict[str, Dict[str, Any]]): A dictionary containing the pairs of atoms for which to compute the distances.
        sysAnalDir (str): The directory where the analysis files will be saved.
        inputDirName (str): The name of the input directory.
        tag (str): A tag to use in the output file names.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary containing the DataFrames with the atomic distances for each pair of atoms.
    """
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
        # Compute the distances between the specified atoms in the trajectory
        contactDistances  = md.compute_distances(traj, pairwiseContacts)

        # Create a DataFrame with the calculated distances
        distanceDf = pd.DataFrame(contactDistances, columns=plotLabels)

        # Create the file path for the CSV file
        outCsv = p.join(sysAnalDir, f"{tag}_{atomTag}_{inputDirName}.csv")

        # Save the DataFrame as a CSV file
        distanceDf.to_csv(outCsv)

        # Add the DataFrame to the dictionary of DataFrames
        distanceDfs.update({atomTag: distanceDf})

    return distanceDfs
#############################################################################################
def compute_contact_distances(
    traj: md.Trajectory,
    residuePairs: Dict[str, Dict[str, Any]],
    sysAnalDir: str,
    inputDirName: str
) -> Dict[str, pd.DataFrame]:
    """
    Calculate contact distances for each pair of residues in the given trajectory
    and save the results as CSV files.

    Args:
        traj (mdtraj.Trajectory): The trajectory containing the residues.
        residuePairs (Dict[str, Dict[str, Any]]): A dictionary containing the pairs of residues for which to compute the contacts.
        sysAnalDir (str): The directory where the analysis files will be saved.
        inputDirName (str): The name of the input directory.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary containing the DataFrames with the contact distances for each pair of residues.
    """
    print("---->\tCalculating contact distances!")
    contactDfs = {}
    for resTag in residuePairs:
        print(f"-------->{resTag}")
        # Get contacts and labels from residuePairs dict
        pairwiseContacts = residuePairs[resTag]["contacts"]
        if len(pairwiseContacts) == 0:
            continue
        plotLabels = residuePairs[resTag]["labels"]

        ## Calculate contact distances
        # Compute the contacts between the specified residues in the trajectory
        contactDistances, _ = md.compute_contacts(traj, pairwiseContacts)

        # Create a DataFrame with the calculated contact distances
        contactDf = pd.DataFrame(contactDistances, columns=plotLabels)

        # Create the file path for the CSV file
        outCsv = p.join(sysAnalDir, f"contacts_{resTag}_{inputDirName}.csv")

        # Save the DataFrame as a CSV file
        contactDf.to_csv(outCsv)

        # Add the DataFrame to the dictionary of DataFrames
        contactDfs.update({resTag: contactDf})

    return contactDfs

#############################################################################################
def find_pairwise_residue_contacts(
    traj: md.Trajectory,
    pdbDf: pd.DataFrame,
    keyResidues: Dict[str, Dict[str, Any]]
) -> Dict[str, Dict[str, Any]]:
    """
    Finds pairwise contacts between residues in a trajectory and saves the results as a dictionary.

    Args:
        traj (mdtraj.Trajectory): The trajectory containing the residues.
        pdbDf (pandas.DataFrame): The DataFrame containing the PDB information.
        keyResidues (Dict[str, Dict[str, Any]]): A dictionary containing the pairs of residues for which to compute the contacts.

    Returns:
        Dict[str, Dict[str, Any]]: A dictionary containing the pairwise contacts for each residue pair.
    """
    print("---->\t Finding pairwise contacts!")

    # Get the mapping of PDB residue IDs to trajectory residue IDs
    residueMapping: Dict[int, int] = get_residue_id_mapping(traj, pdbDf)

    # Initialize the atoms of interest
    interestingAtoms: Dict[str, List[str]] = {}
    _, interestingAtoms, _ = init_atoms_of_interest()

    # Initialize the dictionary to store the pairwise contacts
    residuePairs: Dict[str, Dict[str, Any]] = {}

    # Iterate over the key residues
    for resTag in keyResidues:
        print(f"--------> {resTag}")

        # Get the target residue ID, name, and chain ID
        targetResId: int = int(keyResidues[resTag]["RES_ID"])
        targetResName: str = keyResidues[resTag]["RES_NAME"]
        targetResChain: str = keyResidues[resTag]["CHAIN_ID"]

        # Get the residue name
        resName: str = keyResidues[resTag]["RES_NAME"]

        # Get the DataFrame containing the PDB information for the target residue
        resDf: pd.DataFrame = pdbDf[(pdbDf["RES_ID"] == targetResId) &
                                        (pdbDf["ATOM_NAME"].isin(interestingAtoms[resName])) &
                                        (pdbDf["RES_NAME"] == targetResName)]

        # Check if the target residue is present in the PDB file
        if len(resDf) == 0:
            print(f"--------> No {resTag} in pdb")
            continue

        # Get the query indices for the target residue
        queryIndecies: np.ndarray = resDf["ATOM_ID"].values -  1 ## zero-index

        # Compute the neighbors within a distance of 0.5 nm
        neighbors: List[List[int]] = md.compute_neighbors(traj, 0.5, queryIndecies)

        # Get the unique neighboring atoms
        uniqueNeighborAtoms: List[int] = list(set([item + 1 for sublist in neighbors for item in sublist])) ## one-index!

        # Exclude the query residue and its 2 adjacent residues
        exclueResidues: List[int] = range(targetResId -2, targetResId +3)#[resId - 2, resId, resId + 2]
        neighborDf: pd.DataFrame = pdbDf[(pdbDf["ATOM_ID"].isin(uniqueNeighborAtoms)) & (~pdbDf["RES_ID"].isin(exclueResidues))]

        # Get the neighboring residues and their corresponding chain IDs
        neighborResAndChains: pd.DataFrame = neighborDf[["RES_ID", "CHAIN_ID"]].drop_duplicates()
        neighborResidues: List[int] = neighborResAndChains["RES_ID"].tolist()
        neighborChains: List[str] = neighborResAndChains["CHAIN_ID"].tolist()

        # Translate PDB residue IDs to mdTraj residue IDs
        targetResId_traj: int = pdbResToTrajRes(residueMapping, targetResId, targetResChain)
        nearbyResIds_traj: List[int] = [pdbResToTrajRes(residueMapping, neighborResId, neighborChain)
                                        for neighborResId, neighborChain in zip(neighborResidues, neighborChains)]

        # Create a list of pairwise contacts
        pairwiseContacts: List[Tuple[int, int]] = [(targetResId_traj, nearbyIndex) for nearbyIndex in nearbyResIds_traj]

        # Create plot labels
        neighborDf: pd.DataFrame = neighborDf.copy()
        neighborDf.loc[:,"plotLabels"] = resTag + "--" + neighborDf["RES_NAME"] + ":" + neighborDf["RES_ID"].astype(str)
        plotLabels: List[str] = neighborDf["plotLabels"].unique().tolist()

        # Add the pairwise contacts to the dictionary
        residuePairs.update({resTag:{"contacts": pairwiseContacts,
                                     "labels": plotLabels}})

    return residuePairs
#############################################################################################
def find_hydrogen_bonds(traj: md.Trajectory,
                        pdbDf: pd.DataFrame,
                        keyResidues: Dict[str, Dict[str, str]]) -> Tuple[Dict[str, Dict[str, Union[List[Tuple[int, int]], List[str]]]],
                                                                        Dict[str, Dict[str, Union[List[Tuple[int, int]], List[str]]]]]:
    """
    Finds hydrogen bonds in a trajectory using a provided PDB DataFrame and a dictionary of key residues.

    Args:
        traj (md.Trajectory): The trajectory to search for hydrogen bonds in.
        pdbDf (pd.DataFrame): The DataFrame containing the PDB information.
        keyResidues (Dict[str, Dict[str, str]]): A dictionary of key residues.

    Returns:
        Tuple[Dict[str, Dict[str, Union[List[Tuple[int, int]], List[str]]]],
              Dict[str, Dict[str, Union[List[Tuple[int, int]], List[str]]]]]:
              A tuple containing two dictionaries. The first dictionary contains information about donor-acceptor pairs,
              while the second dictionary contains information about acceptor-donor pairs.
    """
    print("---->\t Finding hydrogen bonds!")
    _, hydrogenBondDonors, hydrogenBondAcceptors = init_atoms_of_interest()

    # Loop through each key residue
    donorAcceptorPairs: Dict[str, Dict[str, Union[List[Tuple[int, int]], List[str]]]] = {}
    acceptorDonorPairs: Dict[str, Dict[str, Union[List[Tuple[int, int]], List[str]]]] = {}
    for resTag in keyResidues:
        print(f"--------> {resTag}")
        # Get relevant information about the residue
        targetResId: int = int(keyResidues[resTag]["RES_ID"])
        targetResName: str = keyResidues[resTag]["RES_NAME"]
        resName: str = keyResidues[resTag]["RES_NAME"]

        # Get dataframes with residue of interest's hydrogen bond donors and acceptors
        hydrogenDonorDf: pd.DataFrame = pdbDf[(pdbDf["RES_ID"] == targetResId) &
                                              (pdbDf["ATOM_NAME"].isin(hydrogenBondDonors[resName])) &
                                              (pdbDf["RES_NAME"] == targetResName)]

        hydrogenAcceptorDf: pd.DataFrame = pdbDf[(pdbDf["RES_ID"] == targetResId) &
                                                 (pdbDf["ATOM_NAME"].isin(hydrogenBondAcceptors[resName])) &
                                                 (pdbDf["RES_NAME"] == targetResName)]

        # Find neighbors for HBA's and HBD's in target residue
        if len(hydrogenDonorDf) > 0:
            for idx, donorRow in hydrogenDonorDf.iterrows():
                # Compute neighboring atoms - get a unique list of them - produce a dataframe
                queryIndex: np.ndarray = np.array([donorRow["ATOM_ID"] -  1]) ## zero-index
                neighbors: List[List[int]] = md.compute_neighbors(traj, 0.5, queryIndex)
                uniqueNeighborAtoms: List[int] = list(set([item + 1 for sublist in neighbors for item in sublist])) ## one-index!
                neighborDf: pd.DataFrame = pdbDf[(pdbDf["ATOM_ID"].isin(uniqueNeighborAtoms))]

                # Drop non-hydrogen bond acceptor atoms
                idxToDrop: List[int] = []
                for idx, acceptorRow in neighborDf.iterrows():
                    if not acceptorRow["ATOM_NAME"] in hydrogenBondAcceptors[acceptorRow["RES_NAME"]]:
                        idxToDrop.append(idx)
                acceptorDf: pd.DataFrame = neighborDf.drop(index=idxToDrop)

                # Create plot labels
                targetAtomTag: str = "-".join([resTag, donorRow["ATOM_NAME"]])
                acceptorDf["plotLabels"] =  acceptorDf["CHAIN_ID"] +":"+ acceptorDf["RES_NAME"] +":"+ acceptorDf["RES_ID"].astype(str) +":"+ acceptorDf["ATOM_NAME"]
                acceptorDf["plotLabels"] = targetAtomTag + " -- H -- " + acceptorDf["plotLabels"]
                plotLabels: List[str] = acceptorDf["plotLabels"].unique().tolist()

                # Create outputs to feed to compute_distances, returns a list of [(X, Y), (X, Z)] contacts
                pairwiseContacts: List[Tuple[int, int]] = [(donorRow["ATOM_ID"] -  1, acceptorIndex - 1) for acceptorIndex in acceptorDf["ATOM_ID"].to_list()] ## zero-index!
                donorAcceptorPairs.update({targetAtomTag:{"contacts": pairwiseContacts,
                                                    "labels": plotLabels}})

        # Find neighbors for HBA's and HBD's in target residue
        if len(hydrogenAcceptorDf) > 0:
            for idx, acceptorRow in hydrogenAcceptorDf.iterrows():
                queryIndex: np.ndarray = np.array([acceptorRow["ATOM_ID"] -  1]) ## zero-index
                neighbors: List[List[int]] = md.compute_neighbors(traj, 0.5, queryIndex)
                uniqueNeighborAtoms: List[int] = list(set([item + 1 for sublist in neighbors for item in sublist])) ## one-index!
                neighborDf: pd.DataFrame = pdbDf[(pdbDf["ATOM_ID"].isin(uniqueNeighborAtoms))]

                # Drop non-hydrogen bond acceptor atoms
                idxToDrop: List[int] = []
                for idx, donorRow in neighborDf.iterrows():
                    if not donorRow["ATOM_NAME"] in hydrogenBondDonors[donorRow["RES_NAME"]]:
                        idxToDrop.append(idx)
                donorDf: pd.DataFrame = neighborDf.drop(index=idxToDrop)

                # Create plot labels
                targetAtomTag: str = "_".join([resTag, acceptorRow["ATOM_NAME"]])
                donorDf["plotLabels"] =  donorDf["CHAIN_ID"] +":"+ donorDf["RES_NAME"] +":"+ donorDf["RES_ID"].astype(str) +":"+ donorDf["ATOM_NAME"]
                donorDf["plotLabels"] = targetAtomTag + " -- H -- " + donorDf["plotLabels"]

                plotLabels: List[str] = donorDf["plotLabels"].unique().tolist()

                # Compute contacts, returns a list of [(X, Y), (X, Z)] contacts
                pairwiseContacts: List[Tuple[int, int]] = [(acceptorRow["ATOM_ID"] -  1, donorIndex - 1) for donorIndex in donorDf["ATOM_ID"].to_list()] ## zero-index!
                acceptorDonorPairs.update({targetAtomTag:{"contacts": pairwiseContacts,
                                                    "labels": plotLabels}})

    return donorAcceptorPairs, acceptorDonorPairs


#############################################################################################
def init_atoms_of_interest() -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Initialize the atoms of interest.

    Returns:
        interestingAtoms (Dict[str, List[str]]): A dictionary containing the atoms of interest for each amino acid.
        hydrogenBondDonors (Dict[str, List[str]]): A dictionary containing the hydrogen bond donors for each amino acid.
        hydrogenBondAcceptors (Dict[str, List[str]]): A dictionary containing the hydrogen bond acceptors for each amino acid.
    """
    # Define the atoms of interest for each amino acid
    interestingAtoms: Dict[str, List[str]] = {
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

    # Define the hydrogen bond donors for each amino acid
    hydrogenBondDonors: Dict[str, List[str]] = {
        "ALA": [],
        "LEU": ["N"],
        "GLY": ["N"],
        "ILE": ["N"],
        "VAL": ["N"],
        "PHE": ["N"],
        "TRP": ["N", "NE1"],
        "TYR": ["N", "OH"],
        "ASP": ["N"],
        "GLU": ["N"],
        "ARG": ["N", "NH1", "NH2"],
        "HIS": ["N", "ND1", "NE2"],
        "LYS": ["N", "NZ"],
        "SER": ["N", "OG"],
        "THR": ["N", "OG1"],
        "ASN": ["N", "ND2"],
        "GLN": ["N", "NE2"],
        "CYS": ["N", "SG"],
        "MET": ["N"],
        "PRO": ["N"]
    }

    # Define the hydrogen bond acceptors for each amino acid
    hydrogenBondAcceptors: Dict[str, List[str]] = {
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
