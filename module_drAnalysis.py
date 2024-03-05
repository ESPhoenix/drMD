import os
from os import path as p
import pandas as pd
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
from PIL import Image
import numpy as np
from scipy.spatial import distance

### drMD modules ###
from module_pdbUtils import pdb2df
#############################################################################################
def check_RMSD(traj):
    rmsd = md.rmsd(traj, traj, 0)
    rmsdDf = pd.DataFrame(rmsd, columns = ["RMSD"])
    return rmsdDf
#############################################################################################
def check_RMSF(traj):
    rmsf = md.rmsf(traj, traj, 0)
    rmsfDf = pd.DataFrame(rmsf, columns = ["RMSF"])
    return rmsfDf
#############################################################################################
def find_contacts(traj, pdbFile, keyResidues):
    contactsPerTarget = find_nearby_residues(keyResidues, pdbFile)
    for resTag in contactsPerTarget:
        print(f"-->{resTag}")
        pairwiseContacts = contactsPerTarget[resTag]
        ## compute contacts
        contactDistances, residuePairs = md.compute_contacts(traj, pairwiseContacts)
        contactIds = residuePairs[:, 1].tolist()
        contactIds = [x+1 for x in contactIds] # un-Zero index
        contactDf = pd.DataFrame(contactDistances, columns = contactIds)
        #print(contactDf)
        ## calculate radial distribution function for pairs
        radii, rdf = md.compute_rdf(traj, pairwiseContacts, bin_width = 0.05)
        rdfDf = pd.DataFrame({'Radii': radii, 'RDF': rdf})
        # peak detection with a cutoff of 500
        peak_distances = [rdfDf['Radii'][i] for i in range(len(rdfDf['RDF'])) if rdfDf['RDF'][i] > 500]
        for i, pair in enumerate(pairwiseContacts):
            pair_distances = contactDistances[i, :]
            for peak_distance in peak_distances:
                # Check if this pair often has a distance close to the peak distance
                close_to_peak = np.isclose(pair_distances, peak_distance, atol=0.01)
                if close_to_peak.mean() > 0.05:
                    print(f"Pair {pair} is associated with peak at distance {round(peak_distance,4)}")


        plot_rdf(rdfDf)

#############################################################################################
def plot_rdf(rdfDf):
    # Plot the RDF
    plt.figure(figsize=(10, 6))
    plt.plot(rdfDf['Radii'], rdfDf['RDF'], label='RDF')
    plt.xlabel('Radii')
    plt.ylabel('Radial Distribution Function')
    plt.title('Radial Distribution Function vs Radii')
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.close()
#############################################################################################
def find_nearby_residues(keyResidues, pdbFile, cutoff = 6.2):
    pdbDf = pdb2df(pdbFile)
    pdbDf = pdbDf[pdbDf["RES_NAME"] != "HOH"]
    pdbDf = pdbDf[pdbDf["ELEMENT"] != "H"]
    contactsPerTarget = {}
    for resTag in keyResidues:
        resId = keyResidues[resTag]["RES_ID"]
        resDf = pdbDf[pdbDf["RES_ID"] == resId]
        searchDf = pdbDf[pdbDf["RES_ID"] != resId].copy()
        resCoords = resDf[["X", "Y", "Z"]].values
        searchCoords = searchDf[["X", "Y", "Z"]].values
        distances = distance.cdist(resCoords, searchCoords, "euclidean")
        minDists = np.amin(distances, axis=0)
        searchDf.loc[:,"MIN_DIST"] = minDists
        searchDf = searchDf[searchDf["MIN_DIST"] <= cutoff] 
        potentialContacts = searchDf["RES_ID"].unique().tolist()
        pairwiseContacts = [(resId , nearbyResidue) for nearbyResidue in potentialContacts] # Zero Index
        contactsPerTarget.update({resTag: pairwiseContacts})
    return contactsPerTarget

#############################################################################################
def main():
    simDir = "/home/esp/scriptDevelopment/drMD/02_outputs/cvFAP_WT_PLM_FAD_3/05_NpT_equilibriation"
    pdbFile = p.join(simDir, "NpT_final_geom.pdb")
    dcdFile = p.join(simDir, "trajectory.dcd")

    traj = md.load_dcd(dcdFile, top = pdbFile)

    keyResidues = {"C372": {"CHAIN_ID": "A",
                            "RES_ID": 372,
                            "RES_NAME": "CYS"},
                    "R391": {"CHAIN_ID": "A",
                            "RES_ID": 391,
                            "RES_NAME": "ARG"}}

    find_contacts(traj, pdbFile, keyResidues)
    check_RMSD(traj)
    check_RMSF(traj)


#############################################################################################
if  __name__ == "__main__":
    main()