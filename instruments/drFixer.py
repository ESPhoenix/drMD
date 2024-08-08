
## BASIC LIBS
import os
from os import path as p
import pandas as pd
## CUSTOM LIBS
from pdbUtils import pdbUtils
## CLEAN CODE
from typing import List, Dict, Union, Any
try:
    from instruments.drCustomClasses import FilePath, DirectoryPath
except:
    from drCustomClasses import FilePath, DirectoryPath


##################################################################################################
def reset_atom_numbers(pdbFile: str) -> str:
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
    nonProtLigNames: List[str] = ["HOH", "WAT", "TIP3",
                    "Na+", "Cl-", "Mg2+", "Ca2+"]
    
    goodDf: pd.DataFrame = goodDf[~goodDf["RES_NAME"].isin(nonProtLigNames)]

    badNoProtLigDf: pd.DataFrame = badDf[badDf["RES_NAME"].isin(nonProtLigNames)]
    badDf: pd.DataFrame = badDf[~badDf["RES_NAME"].isin(nonProtLigNames)]

    # Create a dictionary to store the mappings of residue IDs and chain IDs from the bad PDB file to the good PDB file
    residueMapping: Dict[int, int] = {}
    chainMapping: Dict[int, str] = {}

    # Iterate through the unique chain IDs in the good PDB file
    for chainId in goodDf["CHAIN_ID"].unique():
        # Get the unique residue IDs for the current chain in the good PDB file
        goodChainResidues: pd.Series = goodDf[goodDf["CHAIN_ID"] == chainId]["RES_ID"].unique()
        # Get the corresponding residue IDs for the current chain in the bad PDB file
        badChainResidues: pd.Series = badDf[badDf["CHAIN_ID"] == chainId]["RES_ID"].unique()
        # Create mappings for the residue IDs and chain IDs
        residueMapping.update(dict(zip(badChainResidues, goodChainResidues)))
        chainMapping.update(dict(zip(badChainResidues, [chainId] * len(badChainResidues))))

    # Modify the chain and residue IDs in the bad PDB file to match the good PDB file
    badDf["CHAIN_ID"] = badDf["RES_ID"].map(chainMapping)
    badDf["RES_ID"] = badDf["RES_ID"].map(residueMapping)

    # Concatenate the modified bad PDB dataframe with the waters and ions from the original bad PDB dataframe
    outputDf: pd.DataFrame = pd.concat([badDf, badNoProtLigDf], axis = 0)
    # Write the modified PDB file
    pdbUtils.df2pdb(outputDf, badPdb)
    return badPdb

##################################################################################################