from pdbUtils import pdbUtils
import pandas as pd
from typing import List
#######################################################################
def get_atom_indexes(selection: dict, pdbFile: str) -> List[int]:
    """
    This function takes a selection dictionary and a PDB file as input,
    and returns a list of atom indexes that correspond to the selection.

    Args:
        selection (dict): A dictionary containing the selection type and input.
        pdbFile (str): The path to the PDB file.

    Returns:
        list: A list of atom indexes that correspond to the selection.
    """
    # Read PDB file into a DataFrame
    pdbDf: pd.DataFrame = pdbUtils.pdb2df(pdbFile)

    # Initialize lists of amino acid residue names, backbone atom names,
    # solvent residue names, and ion residue names
    aminoAcidResNames, backboneAtomNames, solventResNames, ionResNames = init_name_lists()

    # Initialize list of atom indexes
    atomIndexes: List[int] = []
    print(selection)
    # Check selection type and find required atom indexes
    if selection["type"] == "backbone":
        # Find indexes for all backbone atoms
        backboneDf: pd.DataFrame = pdbDf[pdbDf["ATOM_NAME"].isin(backboneAtomNames)]
        atomIndexes = backboneDf.index.tolist()
    elif selection["type"] == "protein":
        # Find indexes for all protein atoms
        proteinDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        atomIndexes = proteinDf.index.tolist()
    elif selection["type"] == "water":
        # Find indexes for water molecules
        waterDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(solventResNames)]
        atomIndexes = waterDf.index.tolist()
    elif selection["type"] == "ions":
        # Find indexes for ions
        ionDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(ionResNames)]
        atomIndexes = ionDf.index.tolist()
    elif selection["type"] == "ligands":
        # Find indexes for all ligands / organics / oddball molecules
        ligandDf: pd.DataFrame = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcidResNames+solventResNames+ionResNames)]
        atomIndexes = ligandDf.index.tolist()
    elif selection["type"] == "residue":
        # Find indexes for whole residue selections
        selectionInput = selection["input"]
        for residue in selectionInput:
            residueDf: pd.DataFrame = pdbDf[(pdbDf["CHAIN_ID"] == residue[0])
                              & (pdbDf["RES_NAME"] == residue[1])
                              & (pdbDf["RES_ID"] == int(residue[2])) ]
            atomIndexes += residueDf.index.tolist()
    elif selection["type"] == "atom":
        # Find indexes for atom-by-atom selections
        selectionInput = selection["input"]
        for atom in selectionInput:
            atomDf: pd.DataFrame = pdbDf[(pdbDf["CHAIN_ID"] == atom[0])
                              & (pdbDf["RES_NAME"] == atom[1])
                              & (pdbDf["RES_ID"] == int(atom[2]))
                              & (pdbDf["ATOM_NAME"] == atom[3]) ]
            atomIndexes += atomDf.index.tolist()

    return atomIndexes

#######################################################################
def init_name_lists():
    aminoAcidResNames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN',
            'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
                'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR',
                'VAL']
    backboneAtomNames = ["N","CA","C","O"]
    solventResNames = ["HOH", "WAT"]
    ionResNames = ["NA", "K", "LI", "RB", "CS", "MG", "CA", "ZN", "CD", "HG", "MN"]
    return aminoAcidResNames, backboneAtomNames, solventResNames, ionResNames

#######################################################################