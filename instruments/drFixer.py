
import os
from os import path as p
from instruments import pdbUtils
import pandas as pd

##################################################################################################
def reset_chains_residues(goodPdb, badPdb):
    ## load pdb files as dataframes - separate out waters and ions
    goodDf = pdbUtils.pdb2df(goodPdb)
    badDf = pdbUtils.pdb2df(badPdb)


    hetAtomNames = ["HOH", "WAT", "TIP3",
                    "Na+", "Cl-", "Mg2+", "Ca2+"]
    
    goodDf = goodDf[goodDf["ATOM"] == "ATOM"]
    badHetAtoms = badDf[(badDf["ATOM"] == "HETATM") | (badDf["RES_NAME"].isin(hetAtomNames))]

    badDf = badDf.drop(badHetAtoms.index, errors='ignore')

    chainIds = goodDf["CHAIN_ID"].unique().tolist()
    ## create a list of lists containing residue numbers and associated chains
    goodResidues = []       ## stores unique residue ids for each chain
    goodChains = []     ## stores chain ids associated with the above
    for chainId in chainIds:
        goodChainDf = goodDf[goodDf["CHAIN_ID"] == chainId]
        goodChainResidues = list(set(goodChainDf["RES_ID"].to_list()))
        goodResidues.append(goodChainResidues)
        goodChainChains = [chainId for _ in range(len(goodChainResidues))]
        goodChains.append(goodChainChains)


    ## get list of residues in bad pdb
    badResidues = badDf["RES_ID"].to_list()
    ## create mappings that translate bad pdb to good pdb chains and residues
    residueMapping = {}
    chainMapping = {}
    previousBadResidue = badResidues[0]
    residueCount = 0
    chainCount = 0


    for badResidue in badResidues:
        if residueCount == len(goodResidues[chainCount]):
            if chainCount  + 1 < len(goodResidues):
                chainCount += 1
                residueCount = 0
                residueMapping.update({previousBadResidue:goodResidues[chainCount][residueCount]})
                chainMapping.update({previousBadResidue:goodChains[chainCount][residueCount]})
            else:
                residueMapping.update({previousBadResidue:goodResidues[chainCount][residueCount]})
                chainMapping.update({previousBadResidue:goodChains[chainCount][residueCount]})
                break

        elif badResidue != previousBadResidue:

            residueMapping.update({previousBadResidue:goodResidues[chainCount][residueCount]})
            chainMapping.update({previousBadResidue:goodChains[chainCount][residueCount]})
            residueCount += 1
            previousBadResidue  = badResidue

    residueMapping.update({badResidue:goodResidues[chainCount][residueCount]})
    chainMapping.update({badResidue:goodChains[chainCount][residueCount]})

    outputDf = badDf.copy()
    outputDf["CHAIN_ID"] = outputDf["RES_ID"].map(chainMapping)

    outputDf["RES_ID"] = outputDf["RES_ID"].map(residueMapping)

    outputDf = pd.concat([outputDf, badHetAtoms], axis = 0)

    pdbUtils.df2pdb(outputDf, badPdb)
    return badPdb

##################################################################################################