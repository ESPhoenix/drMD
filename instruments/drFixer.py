
## BASIC LIBS
import os
from os import path as p
import pandas as pd
## CUSTOM LIBS
from pdbUtils import pdbUtils
## CLEAN CODE
from typing import List, Dict, Union, Any, Optional
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



def reset_chains(refPdb: FilePath, inputPdb: FilePath, ligandNames = []) -> FilePath:

    refDf = pdbUtils.pdb2df(refPdb)
    inputDf = pdbUtils.pdb2df(inputPdb)

    refLigDf = refDf[refDf["RES_NAME"].isin(ligandNames)]

    inputLigDf = inputDf[inputDf["RES_NAME"].isin(ligandNames)]




    ligChainMap: dict = {}
    for chainId, chainDf in refLigDf.groupby("CHAIN_ID"):
        for resId, resDf in chainDf.groupby("RES_ID"):
            ligChainMap[resId] = chainId

    for chainId, chainDf in refLigDf.groupby("CHAIN_ID"):
        for resId, resDf in chainDf.groupby("RES_ID"):
            pass





##################################################################################################
def reset_chains_residues(goodPdb: str, badPdb: str, ligandNames = []) -> str:
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
    solventAndIonNames: List[str] = ["HOH", "WAT", "TIP3",
                    "Na+", "Cl-", "Mg2+", "Ca2+"]
    
    goodDf: pd.DataFrame = goodDf[~goodDf["RES_NAME"].isin(solventAndIonNames)]

    solventAndIonsDf: pd.DataFrame = badDf[badDf["RES_NAME"].isin(solventAndIonNames)]
    solventAndIonsDf["CHAIN_ID"] = " "
    badDf: pd.DataFrame = badDf[~badDf["RES_NAME"].isin(solventAndIonNames)]


    badDf["CHAIN_ID"] = goodDf["CHAIN_ID"]
    badDf["RES_ID"] = goodDf["RES_ID"]


    recombinedDf = pd.concat([badDf, solventAndIonsDf])

    fixedPdb = p.splitext(badPdb)[0] + "_fixed.pdb"
    pdbUtils.df2pdb(recombinedDf, fixedPdb)
    return fixedPdb

def old():
    pass

    # # Iterate through the unique chain IDs in the good PDB file
    # for chainId in goodDf["CHAIN_ID"].unique():
    #     # Get the unique residue IDs for the current chain in the good PDB file
    #     goodChainResidues: pd.Series = goodDf[goodDf["CHAIN_ID"] == chainId]["RES_ID"].unique()
    #     # Get the corresponding residue IDs for the current chain in the bad PDB file
    #     badDf['RES_ID'] = badDf['RES_ID'].fillna(-1)
    #     badChainResidues: pd.Series = badDf[badDf["CHAIN_ID"] == chainId]["RES_ID"].unique()
    #     # Create mappings for the residue IDs and chain IDs
    #     residueMapping.update(dict(zip(badChainResidues, goodChainResidues)))
    #     chainMapping.update(dict(zip(badChainResidues, [chainId] * len(badChainResidues))))
    # badDf["RES_ID"] = badDf["RES_ID"].map(residueMapping)


    # for goodChain, chainDf in goodDf.groupby("CHAIN_ID"):
    #     for resId, resDf in chainDf.groupby("RES_ID"):
    #         chainMapping[chainDf["RES_ID"]] = goodChain

    # # Modify the chain and residue IDs in the bad PDB file to match the good PDB file
    # badDf["CHAIN_ID"] = badDf["RES_ID"].map(chainMapping)
    # badDf["RES_ID"] = badDf["RES_ID"].map(residueMapping)


    # # Concatenate the modified bad PDB dataframe with the waters and ions from the original bad PDB dataframe
    # outputDf: pd.DataFrame = pd.concat([badDf, badNoProtLigDf], axis = 0)
    # # Write the modified PDB file
    # fixedPdb = p.splitext(badPdb)[0] + "_fixed.pdb"


##################################################################################################
def fix_atom_names(df): 
    # deal with unwanted apostrophies (prime)
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("'", "")
    # deal with numbers at the beginning of atom names
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].replace(r'^(\d+)(.+)$', r'\2\1', regex=True)
    # deal with "A" at the start of atom name
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].apply(lambda x: x.lstrip('A') if x.startswith('A') else x)

    ## ensure unique names
    count_series = df.groupby('ATOM_NAME').cumcount()
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'] + "_" +count_series.astype(str)
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("_0", "")
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("_", "")

    return df 

##################################################################################################

if __name__ == "__main__":
    goodPdb = "/home/esp/scriptDevelopment/drMD/03_test_outputs/A0A0D2XFD3_TPA_1/00_prep/WHOLE/A0A0D2XFD3_TPA_1.pdb"
    badPdb = "/home/esp/scriptDevelopment/drMD/03_test_outputs/A0A0D2XFD3_TPA_1/00_prep/WHOLE/A0A0D2XFD3_TPA_1_solvated.pdb"
    reset_chains_residues(goodPdb, badPdb)