import os
from os import path as p
import pandas as pd
import mdtraj as md
import argpass
import yaml
### drMD modules ###
from pdbUtils import pdbUtils
import instruments.drPlot as drPlot
import instruments.drDiagnosis as drDiagnosis

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
    # plotting_protocol(analMenu, analDir, keyResidues)

#############################################################################################
def plotting_protocol(analMenu, analDir, keyResidues):
    ############ read analysis menu and do  plotting that has been ordered ############
    keyResiAnal = analMenu["keyResidueAnalysis"]
    wholeAnal = analMenu["wholeProteinAnalysis"]
    sysAnalDirs = [p.join(analDir,dir) for dir in os.listdir(analDir) if p.isdir(p.join(analDir,dir))]
    referenceSystem = analMenu["referenceSystem"]
    ## for interesting residues
    ## plot histograms
    drPlot.histogram_plotting_manager(sysAnalDir)

    # if keyResiAnal["contactDistances"]:
    #     for sysAnalDir in sysAnalDirs:
    #         for resTag in keyResidues: 
    #             drPlot.histogram_plotting_manager(sysAnalDir, idTag = resTag, dataTag = "contacts")
    # if keyResiAnal["findHydrogenBonds"]:
    #     for sysAnalDir in sysAnalDirs:
    #      for resTag in keyResidues:
    #         drPlot.histogram_plotting_manager(sysAnalDir, idTag = resTag, dataTag = "acceptor")
    #         drPlot.histogram_plotting_manager(sysAnalDir, idTag = resTag, dataTag = "donor")

            
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
    pdbDf = pdbUtils.pdb2df(pdbFile)
    pdbDf = pdbDf[~pdbDf["RES_NAME"].isin(["Cl-","Na+","HOH"])].copy()

    ############ read analysis menu and do ordered analysis ############
    keyResiAnal = analMenu["keyResidueAnalysis"]
    wholeAnal = analMenu["wholeProteinAnalysis"]

    analysisData = {}

    ## FOR INTERESTING RESIDUES
    if any(job for job in keyResiAnal.values() if job):
        residuePairs = drDiagnosis.find_pairwise_residue_contacts(traj, pdbDf, keyResidues)
        if keyResiAnal["contactDistances"]:
            contactDf = drDiagnosis.compute_contact_distances(traj, residuePairs, sysAnalDir, inputDirName)
            if keyResiAnal["radialDistribution"]:
                rdfDf = drDiagnosis.compute_radial_distribution(traj, residuePairs, contactDf, sysAnalDir, inputDirName)
                analysisData.update({"rdf": rdfDf})
        if keyResiAnal["findHydrogenBonds"]:
            donorAcceptorPairs, acceptorDonorPairs = drDiagnosis.find_hydrogen_bonds(traj, pdbDf, keyResidues)
            donorAcceptorDistancesDfs = drDiagnosis.compute_atomic_distances(traj, donorAcceptorPairs, sysAnalDir, inputDirName, "donor")
            acceptorDonorDistancesDfs = drDiagnosis.compute_atomic_distances(traj, acceptorDonorPairs, sysAnalDir, inputDirName, "acceptor")
            analysisData.update({"donorAcceptor": donorAcceptorDistancesDfs,
                                 "acceptorDonor": acceptorDonorDistancesDfs})


    ## FOR WHOLE PROTEIN PROPERTIES
    if wholeAnal["RMSD"]:
        rmsdDf = drDiagnosis.check_RMSD(traj,sysAnalDir, inputDirName)
        analysisData.update({"rmsd": rmsdDf})

    if wholeAnal["RMSF"]:
        rmsfDf = drDiagnosis.check_RMSF(traj,pdbDf, sysAnalDir, inputDirName)
        analysisData.update({"rmsf": rmsfDf})

    return analysisData
#############################################################################################
if  __name__ == "__main__":
    main()