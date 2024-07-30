
import os
from os import path as p
from pdbUtils import pdbUtils
import pandas as pd

##################################################################################################
def  reset_atom_numbers(pdbFile: str) -> str:
    """
    Resets the atom numbers in a PDB file.

    Parameters:
        pdbFile (str): Path to the PDB file.

    Returns:
        str: Path to the modified PDB file.
    """

    pdbDf = pdbUtils.pdb2df(pdbFile)
    pdbDf["ATOM_ID"] = range(1, len(pdbDf) + 1)

    pdbUtils.df2pdb(pdbDf, pdbFile)

    return pdbFile
##################################################################################################
def reset_chains_residues(goodPdb: str, badPdb: str) -> str:
    """
    Resets the chains and residues in a PDB file to match another PDB file.

    Parameters:
        goodPdb (str): Path to the PDB file with the correct chains and residues.
        badPdb (str): Path to the PDB file with incorrect chains and residues.

    Returns:
        str: Path to the modified PDB file.
    """
    ## load pdb files as dataframes - separate out waters and ions
    # Load the good and bad PDB files as dataframes
    goodDf: pd.DataFrame = pdbUtils.pdb2df(goodPdb)
    badDf: pd.DataFrame = pdbUtils.pdb2df(badPdb)

    ## drop waters and ions from both good and bad dataframes - we don't need to re-do these!
    # Drop waters and ions from the good and bad dataframes
    hetAtomNames: List[str] = ["HOH", "WAT", "TIP3",
                    "Na+", "Cl-", "Mg2+", "Ca2+"]
    
    goodWaterAndIons: pd.DataFrame = goodDf[goodDf["RES_NAME"].isin(hetAtomNames)]
    goodDf: pd.DataFrame = goodDf.drop(goodWaterAndIons.index, errors='ignore')

    badWaterAndIons: pd.DataFrame = badDf[badDf["RES_NAME"].isin(hetAtomNames)]

    badDf: pd.DataFrame = badDf.drop(badWaterAndIons.index, errors='ignore')

    chainIds: list = goodDf["CHAIN_ID"].unique().tolist()
    ## create a list of lists containing residue numbers and associated chains
    # Create a list of lists containing unique residue IDs and associated chain IDs for each chain
    goodResidues: list = []       ## stores unique residue ids for each chain
    goodChains: list = []     ## stores chain ids associated with the above
    for chainId in chainIds:
        goodChainDf: pd.DataFrame = goodDf[goodDf["CHAIN_ID"] == chainId]
        goodChainResidues: list = list(set(goodChainDf["RES_ID"].to_list()))
        goodResidues.append(goodChainResidues)
        goodChainChains: list = [chainId for _ in range(len(goodChainResidues))]
        goodChains.append(goodChainChains)

    ## get list of residues in bad pdb
    # Get the list of residue IDs in the bad PDB file
    badResidues: list = badDf["RES_ID"].to_list()
    ## create mappings that translate bad pdb to good pdb chains and residues
    # Create mappings to translate residue IDs and chain IDs from the bad PDB file to the good PDB file
    residueMapping: dict = {}
    chainMapping: dict = {}
    previousBadResidue: int = badResidues[0]
    residueCount: int = 0
    chainCount: int = 0

    for badResidue in badResidues:
        if residueCount == len(goodResidues[chainCount]):
            if chainCount  + 1 < len(goodResidues):
                chainCount += 1
                residueCount: int = 0
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

    # Modify the chain and residue IDs in the bad PDB file to match the good PDB file and save the modified PDB file
    outputDf: pd.DataFrame = badDf.copy()
    outputDf["CHAIN_ID"] = outputDf["RES_ID"].map(chainMapping)

    outputDf["RES_ID"] = outputDf["RES_ID"].map(residueMapping)

    outputDf: pd.DataFrame = pd.concat([outputDf, badWaterAndIons], axis = 0)
    ## write output pdb
    pdbUtils.df2pdb(outputDf, badPdb)
    return badPdb

##################################################################################################