## OPENMM LIBS
import simtk.openmm as openmm
import simtk.unit as unit
import xml.etree.ElementTree as ET
from os import path as p
from pdbUtils import pdbUtils
import math
###########################################################################################
def constraints_handler(system, prmtop, inpcrd, sim, saveFile, pdbFile):

    if "restraints" in sim:
        restraintInfo = sim["restraints"]

        kNumber = 0        
        for restraint in restraintInfo:
            selection = restraint["selection"]
            if restraint["type"] == "position":
                system = create_position_restraint(system, prmtop, inpcrd, selection, kNumber, pdbFile)
            elif restraint["type"] == "distance":
                parameters = restraint["parameters"]
                system = create_distance_restraint(system, prmtop, selection, parameters, kNumber, pdbFile)
            elif restraint["type"] == "angle":
                parameters = restraint["parameters"]
                system = create_angle_restraint(system, prmtop, selection, parameters, kNumber, pdbFile)
            elif restraint["type"] == "torsion":
                parameters = restraint["parameters"]
                system = create_torsion_restraint(system, prmtop, selection, parameters, kNumber, pdbFile)
            kNumber += 1
    
    else:
        print("Running with no restraints...")
        if p.splitext(saveFile)[1] == ".chk":
            return system
        clear_all_restraints(saveFile)
    return system
###########################################################################################
def create_position_restraint(system, prmtop, inpcrd, selection, kNumber, pdbFile):
    ## create a position restraint
    positionRestraint = openmm.CustomExternalForce(f"k{str(kNumber)}*periodicdistance(x, y, z, x0, y0, z0)^2")
    positionRestraint.addGlobalParameter(f"k{str(kNumber)}", 1000.0 * unit.kilojoules_per_mole / unit.nanometer)
    positionRestraint.addPerParticleParameter("x0")
    positionRestraint.addPerParticleParameter("y0")
    positionRestraint.addPerParticleParameter("z0")
    system.addForce(positionRestraint)

    restraintAtomIndexes = get_restraint_atom_selection(selection, prmtop, pdbFile)
    for restraintAtomIndex in restraintAtomIndexes:
        positionRestraint.addParticle(restraintAtomIndex,inpcrd.getPositions()[restraintAtomIndex])
    return system

###########################################################################################
def create_distance_restraint(system, prmtop, selection, parameters, kNumber, pdbFile):
    ## Create a distance restraint
    k_force_constant = 1000.0 * unit.kilojoules_per_mole / unit.nanometer**2
    distanceRestraint = openmm.CustomBondForce(f"0.5 * k{str(kNumber)} * (r - r0)^2")
    ## add per bond parameters, k for force constant, r0 for desired distance
    distanceRestraint.addPerBondParameter(f"k{str(kNumber)}")
    distanceRestraint.addPerBondParameter("r0")
    ## add force to system
    system.addForce(distanceRestraint)

    # Get the indices of the two atoms to be restrained
    restraintAtomIndexes = get_restraint_atom_selection(selection, prmtop, pdbFile)
    if len(restraintAtomIndexes) != 2:
        raise ValueError("Expected exactly two atom indices for a distance restraint.")

    # Get target distance from the parameters dictionary, and convert from angstroms to nanometers
    kForceConstant = parameters["k"]
    targetDistance_nm = parameters["r0"] * unit.angstroms * 0.1

    # Add the atom pair and the calculated target distance in nanometers to the bond restraint
    distanceRestraint.addBond(restraintAtomIndexes[0], restraintAtomIndexes[1], [kForceConstant, targetDistance_nm])

    return system
###########################################################################################
def create_angle_restraint(system, prmtop, selection, parameters, kNumber, pdbFile):
    ## Create an angle restraint
    angleRestraint = openmm.CustomAngleForce(f"0.5 * k{str(kNumber)} * (theta - theta0)^2")

    ## add per angle parameters, k for force constant, theta0 for desired angle
    angleRestraint.addPerAngleParameter(f"k{str(kNumber)}")
    angleRestraint.addPerAngleParameter("theta0")

    ## add force to system
    system.addForce(angleRestraint)

    # Get the indices of the three atoms to be restrained
    restraintAngleAtoms = get_restraint_atom_selection(selection, prmtop, pdbFile)
    if len(restraintAngleAtoms) != 3:
        raise ValueError("Expected exactly three atom indices for an angle restraint.")

    # Get target angle from the parameters dictionary and convert from degrees to radians
    kForceConstant = parameters["k"] * unit.kilojoules_per_mole / unit.radians**2
    targetAngle_rad = parameters["theta0"] * unit.degrees * (math.pi / 180.0)

    # Add the atom triplet and the calculated target angle in radians to the angle restraint
    angleRestraint.addAngle(restraintAngleAtoms[0],
                             restraintAngleAtoms[1],
                               restraintAngleAtoms[2],
                                 [kForceConstant, targetAngle_rad])

    return system
###########################################################################################
def create_torsion_restraint(system, prmtop, selection, parameters, kNumber, pdbFile):
    ## Create a torsion restraint
    # Note: Ensure 'phi' is the variable used for torsion calculation in OpenMM
    torsionRestraint = openmm.CustomTorsionForce(f"0.5*k{str(kNumber)}*(1-cos(theta-theta0))")
    
    ## Add parameters: k for force constant, phi0 for desired torsion angle
    torsionRestraint.addPerTorsionParameter(f"k{str(kNumber)}")
    torsionRestraint.addPerTorsionParameter("theta0")
    
    ## Add force to system
    system.addForce(torsionRestraint)
    
    ## Get indices of the four atoms to be restrained
    restraintTorsionAtoms = get_restraint_atom_selection(selection, prmtop, pdbFile)
    if len(restraintTorsionAtoms) != 4:
        raise ValueError("Expected exactly four atom indices for a torsion restraint.")
    
    ## Get target torsion angle from parameters (in radians)
    kForceConstant = parameters.get("k", 1000) * unit.kilojoules_per_mole / unit.radians**2  # Default k value if not specified
    targetTorsion_rad = parameters["phi0"] * math.pi / 180.0
    
    ## Add the atom quartet and settings to the torsion restraint
    torsionRestraint.addTorsion(restraintTorsionAtoms[0], restraintTorsionAtoms[1],
                                restraintTorsionAtoms[2], restraintTorsionAtoms[3],
                                [kForceConstant, targetTorsion_rad])
    
    return system
###########################################################################################
def get_restraint_atom_selection(selection, prmtop, pdbFile):
    ## init some residue name lists for the preset options
    aminoAcids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN',
                   'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
                     'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR',
                       'VAL']
    solventResNames = ["HOH", "WAT"]
    ionResNames = ["CL","NA"]
    backboneAtomNames = ["N","HN","CA","O"]
    restraintAtomsIndexes = []
    ## deal with presets
    if selection == "protein":
        for atom in prmtop.topology.atoms():
            if atom.residue.name in aminoAcids:
                restraintAtomsIndexes.append(atom.index)

    elif selection == "water":
        for atom in prmtop.topology.atoms():
            if atom.residue.name in solventResNames:
                restraintAtomsIndexes.append(atom.index)

    elif selection == "ions":
        for atom in prmtop.topology.atoms():
            if atom.residue.name in ionResNames:
                restraintAtomsIndexes.append(atom.index)

    elif selection == "backbone":
        for atom in prmtop.topology.atoms():
            if atom.name not in backboneAtomNames:
                restraintAtomsIndexes.append(atom.index)

    elif selection == "ligands":
        for atom in prmtop.topology.atoms():
            if (atom.residue.name not in aminoAcids and
                atom.residue.name not in solventResNames and
                atom.residue.name not in ionResNames):
                restraintAtomsIndexes.append(atom.index)

    # deal with atom selection as a [CHAIN_ID, RES_NAME, RES_ID, ATOM_NAME]
    elif isinstance(selection,list):
        pdbDf = pdbUtils.pdb2df(pdbFile)
        for atomSelection in selection:
            atomDf = pdbDf[(pdbDf["CHAIN_ID"] == atomSelection[0])
                           & (pdbDf["RES_NAME"] == atomSelection[1])
                           & (pdbDf["RES_ID"] == int(atomSelection[2]))
                           & (pdbDf["ATOM_NAME"] == atomSelection[3])]
            atomIndex = atomDf.index.values[0]
            restraintAtomsIndexes.append(atomIndex)

    return restraintAtomsIndexes

###########################################################################################
def clear_all_restraints(saveXml):
    ## remove any leftover force constants (just k for now)
    tree = ET.parse(saveXml)
    root = tree.getroot()
    parametersElement = root.find("Parameters")
    # Safely remove the 'k' attribute if it exists

    paramsToPop = []
    if parametersElement is not None:
        for param in parametersElement.attrib:
            if param.startswith("k"):
                paramsToPop.append(param)
    for param in paramsToPop:
        parametersElement.attrib.pop(param)
    tree.write(saveXml)
    
# ###########################################################################################
# def heavy_atom_position_restraints(system, prmtop, inpcrd):
#     return restrain_all_atom_names_except_list(system, prmtop, inpcrd, ["H"])
# ###########################################################################################
# def constrain_all_atom_names_except_list(system, prmtop, unrestrictedAtomSymbolList):
#     # Sets the mass of all atoms to zero.
#     # Except atoms in unrestrictedAtomSymbolList and water

#     atomIndicesToFreeze = []
#     # Iterate over all atoms in the topology
#     for atom in prmtop.topology.atoms():
#         # Skip water, as massless objects cannot participate in openmm constraints (H is bound to Oxygen with constraint)
#         if atom.residue.name in ['HOH', 'WAT']:
#             continue
#         # Check if the atom is not in the unrestricted list
#         if atom.element.symbol not in unrestrictedAtomSymbolList:
#             atomIndicesToFreeze.append(atom.index)
#     # Set the mass of these atoms to 0.0 to freeze (constrain) them
#     for index in atomIndicesToFreeze:
#         system.setParticleMass(index, 0.0)

#     return system

###########################################################################################
