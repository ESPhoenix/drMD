import mdtraj as md
from sklearn.cluster import KMeans
import glob
import numpy as np
import os
from os import path as p
import instruments.drSelector as drSelector
from typing import Dict, Union, Any, List
#######################################################################
def rmsd_clustering_protocol(inDir: str, clusterInfo: Dict[str, Union[int, str]]) -> None:
    """
    Clusters a trajectory based on the RMSD values of a subset of atoms.

    Args:
        inDir (str): The input directory containing the trajectory files.
        clusterInfo (Dict[str, Union[int, str]]): A dictionary containing the number of clusters (int)
            and the selection of atoms to be used for clustustering (str).

    Returns:
        None
    """
    print("Clustering trajectory...")

    ## make outDir if needed
    outDir: str = p.join(inDir,"cluster_pdbs")
    os.makedirs(outDir,exist_ok=True)

    ## unpack clusterInfo
    nClusters: int = clusterInfo["nClusters"]
    clusterSelection: str = clusterInfo["selection"]

    ## find output files
    dcdFile: str = p.join(inDir, "trajectory.dcd")
    pdbFile: str = glob.glob(p.join(inDir,"*.pdb"))[0]

    # Check if the MetaDynamics output files exist
    if not p.isfile(dcdFile) or not p.isfile(pdbFile):
        print("MetaDynamics output files not found!")
        print("Better call the Doctor!")
        return

    # Get the atom indexes for the selected atoms
    clusterSelectionAtomIndexes: List[int] = drSelector.get_atom_indexes(clusterSelection, pdbFile)

    # Load trajectory
    traj: md.Trajectory = md.load(dcdFile, top=pdbFile)
    # Optionally, superimpose all frames to the first to remove translational and rotational motions
    traj.superpose(traj[0])

    # Convert trajectory to RMSD matrix
    rmsdMatrix: np.ndarray = convert_traj_to_rmsdMatrix(traj,clusterSelectionAtomIndexes)

    # Perform clustering and save the clusters to PDB files
    kmeans_clusters_to_pdb(rmsdMatrix, nClusters, outDir, traj)

#######################################################################
def convert_traj_to_rmsdMatrix(traj: md.Trajectory, atomIndexes: List[int]) -> np.ndarray:
    """
    Converts a trajectory to an RMSD matrix.

    Args:
        traj (md.Trajectory): The input trajectory.
        atomIndexes (List[int]): The atom indexes to be included in the subtrajectory.

    Returns:
        np.ndarray: The RMSD matrix.
    """
    # Create a subtrajectory only containing the selected atoms
    sectionTraj = traj.atom_slice(atomIndexes)

    # Compute the pairwise RMSD matrix for all frames
    nFrames = sectionTraj.n_frames  # Get the number of frames
    rmsdMatrix = np.empty((nFrames, nFrames))  # Create an empty RMSD matrix
    for i in range(nFrames):
        rmsdMatrix[i] = md.rmsd(sectionTraj, sectionTraj, frame=i)  # Compute the RMSD for each pair of frames

    return rmsdMatrix


#######################################################################
def kmeans_clusters_to_pdb(rmsdMatrix: np.ndarray, bestK: int, outDir: str, traj: md.Trajectory) -> None:
    """
    Saves the cluster centers of a k-means clustering to PDB files.

    Args:
        rmsdMatrix (np.ndarray): The RMSD matrix.
        bestK (int): The number of clusters.
        outDir (str): The output directory.
        traj (md.Trajectory): The input trajectory.
    """
    kmeans = KMeans(n_clusters=bestK)  # Create a KMeans object
    _ = kmeans.fit_predict(rmsdMatrix)  # Fit the KMeans object to the RMSD matrix

    # Save cluster centers
    representative_frames: List[int] = []  # List to store the indices of representative frames
    for i in range(bestK):
        # Find the frame closest to each cluster center
        cluster_center: int = np.argmin(np.linalg.norm(rmsdMatrix - kmeans.cluster_centers_[i], axis=1))
        representative_frames.append(cluster_center)

        # Save the representative frame as a PDB file
        traj[cluster_center].save_pdb(p.join(outDir, f"cluster_center_{i+1}.pdb"))

#######################################################################