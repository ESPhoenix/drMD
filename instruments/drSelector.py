from pdbUtils import pdbUtils
#######################################################################
def get_atom_indexes(selection,  pdbFile):
    pdbDf = pdbUtils.pdb2df(pdbFile)

    aminoAcidResNames, backboneAtomNames, solventResNames, ionResNames = init_name_lists()

    ## reads input selection - finds required atom indexes
    atomIndexes = []
## look for keywords ##
    ## for all backbone atoms
    if selection["type"] == "backbone":
        backboneDf = pdbDf[pdbDf["ATOM_NAME"].isisn(backboneAtomNames)]
        atomIndexes = backboneDf.index.tolist()
    ## for all protein atoms
    elif selection["type"] == "protein":
        proteinDf = pdbDf[pdbDf["RES_NAME"].isin(aminoAcidResNames)]
        atomIndexes = proteinDf.index.tolist()
    ## for water molecules
    elif selection["type"] == "water":
        waterDf = pdbDf[pdbDf["RES_NAME"].isin(solventResNames)]
        atomIndexes = waterDf.index.tolist()
    ## for ions
    elif selection["type"] == "ions":
        ionDf = pdbDf[pdbDf["RES_NAME"].isin(ionResNames)]
        atomIndexes = ionDf.index.tolist()
    ## for all ligands / organics / oddball molecules
    elif selection["type"] == "ligands":
        ligandDf = pdbDf[~pdbDf["RES_NAME"].isin(aminoAcidResNames+solventResNames+ionResNames)]
        atomIndexes = ligandDf.index.tolist()

## if residue or atom style selections, get indexes for required atoms
    ## for whole residue selections
    elif  selection["type"] == "residue":
        selectionInput = selection["input"]
        for residue in selectionInput:
            print(residue)
            residueDf = pdbDf[(pdbDf["CHAIN_ID"] == residue[0])
                              & (pdbDf["RES_NAME"] == residue[1])
                              & (pdbDf["RES_ID"] == int(residue[2])) ]
            atomIndexes += residueDf.index.tolist()
    ## for atom-by-atom selections
    elif  selection["type"] == "atom":
        selectionInput = selection["input"]
        for atom in selectionInput:
            atomDf = pdbDf[(pdbDf["CHAIN_ID"] == atom[0])
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