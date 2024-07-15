from pdbUtils import pdbUtils
import pandas as pd
from typing import List
import numpy as np
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
    # Check selection type and find required atom indexes
    if selection["keyword"] == "all":
        # Find indexes for all atoms
        atomIndexes = pdbDf.index.tolist()
    elif selection["keyword"] == "backbone":
        # Find indexes for all backbone atoms
        backboneDf: pd.DataFrame = pdbDf[pdbDf["ATOM_NAME"].isin(backboneAtomNames) &
                                          pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        atomIndexes = backboneDf.index.tolist()
    elif selection["keyword"] == "protein":
        # Find indexes for all protein atoms
        proteinDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        atomIndexes = proteinDf.index.tolist()
    elif selection["keyword"] == "water":
        # Find indexes for water molecules
        waterDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(solventResNames) &
                                       ~pdbDf["RES_NAME"].isin([aminoAcidResNames])]
        atomIndexes = waterDf.index.tolist()
    elif selection["keyword"] == "ions":
        # Find indexes for ions
        ionDf: pd.DataFrame = pdbDf[pdbDf["RES_NAME"].isin(ionResNames) &
                                    ~pdbDf["RES_NAME"].isin([aminoAcidResNames])]

        atomIndexes = ionDf.index.tolist()
    elif selection["keyword"] == "ligands":
        # Find indexes for all ligands / organics / oddball molecules
        ligandDf: pd.DataFrame = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcidResNames+solventResNames+ionResNames)]
        atomIndexes = ligandDf.index.tolist()
    elif selection["keyword"] == "custom":
        # Find indexes for whole residue selections
        customSelection = selection["customSelection"]
        ## check for "all" selections in selectionSytax
        for selection in customSelection:
            ## create selection conditions
            ## init empty list to store conditions
            selectionConditions: list  = []
            ## for each selection key...
            for selctionKey in ["CHAIN_ID", "RES_NAME", "RES_ID", "ATOM_NAME"]:
                ## use list inputs to create a .isin() condition
                if isinstance(selection[selctionKey], list):
                    selectionConditions.append(pdbDf[selctionKey].isin(range(selection[selctionKey])))
                ## use single input to create a == condition
                elif selection[selctionKey] != "all":
                    selectionConditions.append(pdbDf[selctionKey] == selection[selctionKey])
                ## "all" inputs create no condition
            selectionConditionsCombined = np.logical_and.reduce(selectionConditions)
            selectionDf = pdbDf[selectionConditionsCombined]
            ## add to atomIndexes list
            atomIndexes += selectionDf.index.tolist()

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