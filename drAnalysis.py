import os
from os import path as p
import pandas as pd
import mdtraj as md
import argpass
import yaml
### drMD modules ###
from instruments.pdbUtils import pdb2df
import instruments.drPlot as drPlot
import instruments.drDiagnosis as drDiagnosis



#####################################################################################################
def  compute_radial_dristrubution_function():
    return
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
    # plotting_protocol(analMenu, analDir, keyResidues)
    # return
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
    # plotting_protocol(analMenu, analDir, keyResidues)

#############################################################################################
def plotting_protocol(analMenu, analDir, keyResidues):
    ############ read analysis menu and do  plotting that has been ordered ############
    keyResiAnal = analMenu["keyResidueAnalysis"]
    wholeAnal = analMenu["wholeProteinAnalysis"]
    sysAnalDirs = [p.join(analDir,dir) for dir in os.listdir(analDir) if p.isdir(p.join(analDir,dir))]
    print(sysAnalDirs)
    referenceSystem = analMenu["referenceSystem"]
    ## for interesting residues
    if keyResiAnal["contactDistances"]:
        for sysAnalDir in sysAnalDirs:
            for resTag in keyResidues: 
                drPlot.plot_distance_hist(sysAnalDir, resTag)
    ## for whole protein properties
    ## FOR WHOLE PROTEIN PROPERTIES
    if wholeAnal["RMSD"]:
        for sysAnalDir in sysAnalDirs:

            drPlot.plot_RMSD(sysAnalDir)
    if wholeAnal["RMSF"]:
        for sysAnalDir in sysAnalDirs:
            drPlot.plot_RMSF(sysAnalDir)
        drDiagnosis.compute_delta_RMSF(analDir,referenceSystem)
        if wholeAnal["deltaRMSF"]:
            drPlot.plot_delta_RMSF(analDir, referenceSystem)


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

    ## FOR INTERESTING RESIDUES
    if any(job for job in keyResiAnal.values() if job):
        residuePairs = drDiagnosis.find_pairwise_contacts(traj, pdbDf, keyResidues)
        if keyResiAnal["contactDistances"]:
            contactDf = drDiagnosis.compute_contact_distances(traj, residuePairs, sysAnalDir, inputDirName)
            if keyResiAnal["radialDistribution"]:
                drDiagnosis.compute_radial_distribution(traj, residuePairs, contactDf, sysAnalDir, inputDirName)
    ## FOR WHOLE PROTEIN PROPERTIES
    if wholeAnal["RMSD"]:
        drDiagnosis.check_RMSD(traj,sysAnalDir, inputDirName)
    if wholeAnal["RMSF"]:
        drDiagnosis.check_RMSF(traj,pdbDf, sysAnalDir, inputDirName)
#############################################################################################
if  __name__ == "__main__":
    main()