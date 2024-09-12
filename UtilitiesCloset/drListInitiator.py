from typing import List


##################################################################################
def get_amino_acid_residue_names() -> List[str]:
    """
    Returns a list of the amino acid names.
    """
    return  {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN',
            'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
            'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'ASH', 'GLH', 'HIP', 'HIE', 'HID', 'CYX', 'CYM', 'LYN'} ## oddball protonations
##################################################################################
def get_ion_residue_names() -> List[str]:
    """ 
    Returns a list of the ion atom names.
    """
    return {"LI", "NA", "K", "RB", "CS", "TL", "CU", "AG", "NH4", "H3O", "F", "CL", "BR", "I",
        "BE2", "CU2", "NI2", "PT2", "ZN2", "CO2", "PD2", "AG2", "CR2", "FE2", 
        "MG2", "V2", "MN2", "HG2", "CD2", "YB2", "CA2", "SN2", "PB2", "EU2", 
        "SR2", "SM2", "BA2", "RA2", "AL3", "FE3", "CR3", "IN3", "TL3", "Y3", 
        "LA3", "CE3", "PR3", "ND3", "SM3", "EU3", "GD3", "TB3", "DY3", "ER3", 
        "TM3", "LU3", "HF4", "ZR4", "CE4", "U4", "PU4", "TH4"}

##################################################################################
def get_backbone_atom_names() -> List[str]:
    """
    Returns a list of the backbone atom names.
    """
    return {"N","CA","C","O"}
##################################################################################
def get_solvent_residue_names() -> List[str]:
    """
    Returns a list of the solvent residue names.
    """
    return {"HOH", "WAT"}
##################################################################################
