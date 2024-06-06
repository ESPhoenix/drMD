import mdtraj as md
from sklearn.cluster import KMeans
import glob
import numpy as np
import os
from os import path as p
import instruments.drSelector as drSelector

#######################################################################
def rmsd_clustering_protocol(inDir, clusterInfo):
    print("Clustering trajectory...")
    ## make outDir if needed
    outDir = p.join(inDir,"cluster_pdbs")
    os.makedirs(outDir,exist_ok=True)

    ## unpack clusterInfo
    nClusters = clusterInfo["nClusters"]
    clusterSelection = clusterInfo["selection"]
    ## find output files
    dcdFile = p.join(inDir, "trajectory.dcd")
    pdbFile = glob.glob(p.join(inDir,"*.pdb"))[0]

    if not p.isfile(dcdFile) or not p.isfile(pdbFile):
        print("MetaDynamics output files not found!")
        print("Better call the Doctor!")
        return

    clusterSelectionAtomIndexes = drSelector.get_atom_indexes(clusterSelection, pdbFile)

    # Load trajectory
    traj = md.load(dcdFile, top=pdbFile)
    # Optionally, superimpose all frames to the first to remove translational and rotational motions
    traj.superpose(traj[0])
    ## convert traj into a matrix of rmsd values, using only the subset of atoms in clusterBy
    rmsdMatrix = convert_traj_to_rmsdMatrix(traj,clusterSelectionAtomIndexes)
    ## use best value 
    kmeans_clusters_to_pdb(rmsdMatrix, nClusters, outDir, traj)

#######################################################################
def convert_traj_to_rmsdMatrix(traj, atomIndexes):
    # Create a subtrajectory only containing FMN atoms
    sectionTraj = traj.atom_slice(atomIndexes)

    # Compute the pairwise RMSD matrix for all frames
    nFrames = sectionTraj.n_frames
    rmsdMatrix = np.empty((nFrames, nFrames))
    for i in range(nFrames):
        rmsdMatrix[i] = md.rmsd(sectionTraj, sectionTraj, frame=i)

    return rmsdMatrix

#######################################################################
def kmeans_clusters_to_pdb(rmsdMatrix, bestK, outDir, traj):
    kmeans = KMeans(n_clusters=bestK)
    _ = kmeans.fit_predict(rmsdMatrix)

    # Save cluster centers
    representative_frames = []
    for i in range(bestK):
        # Find the frame closest to each cluster center
        cluster_center = np.argmin(np.linalg.norm(rmsdMatrix - kmeans.cluster_centers_[i], axis=1))
        representative_frames.append(cluster_center)
        traj[cluster_center].save_pdb(p.join(outDir,f"cluster_center_{i+1}.pdb"))

#######################################################################