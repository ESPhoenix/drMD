## OPENMM LIBS
import simtk.openmm as openmm
import simtk.unit as unit
import xml.etree.ElementTree as ET
from os import path as p
###########################################################################################
def constraints_handler(system, prmtop,inpcrd, sim, saveFile):

    constraintOptions = ["relaxWaters"]
    clearRestraints = True
    if "relaxWaters" in sim:
        if sim["relaxWaters"]:
            print("Relaxing water positions...")
            clearRestraints = False
            system = restrain_protein(system, prmtop, inpcrd)
            return system, clearRestraints
    
    if clearRestraints:
        print("Running with no restraints...")
        if p.splitext(saveFile)[1] == ".chk":
            return system, clearRestraints
        clear_all_restraints(saveFile)
        return system, clearRestraints

###########################################################################################
def restrain_protein(system, prmtop, inpcrd):
    restraint = openmm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    system.addForce(restraint)
    restraint.addGlobalParameter("k", 1000.0 * unit.kilojoules_per_mole / unit.nanometer)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    solvent_residues = ['HOH', 'WAT']  
    for atom in prmtop.topology.atoms():
        # Check if the residue of the atom is not part of the solvent
        if atom.residue.name not in solvent_residues:
            restraint.addParticle(atom.index, inpcrd.getPositions()[atom.index])
    return system


###########################################################################################
def restrain_all_atom_names_except_list(system, prmtop, inpcrd, unrestrictedAtomSymbolList):
    # assert all names in unrestrictedAtomNamesList are in an allowed list.
    restraint = openmm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    system.addForce(restraint)
    restraint.addGlobalParameter("k", 1000.0 * unit.kilojoules_per_mole / unit.nanometer)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for atom in prmtop.topology.atoms():
        if atom.element.symbol not in unrestrictedAtomSymbolList:
            restraint.addParticle(atom.index, inpcrd.getPositions()[atom.index])
    return system

###########################################################################################
def clear_all_restraints(saveXml):
    ## remove any leftover force constants (just k for now)
    tree = ET.parse(saveXml)
    root = tree.getroot()
    parametersElement = root.find("Parameters")
    # Safely remove the 'k' attribute if it exists
    if parametersElement is not None and 'k' in parametersElement.attrib:
        parametersElement.attrib.pop('k')
    tree.write(saveXml)
    
###########################################################################################
def heavy_atom_position_restraints(system, prmtop, inpcrd):
    return restrain_all_atom_names_except_list(system, prmtop, inpcrd, ["H"])
###########################################################################################
def constrain_all_atom_names_except_list(system, prmtop, unrestrictedAtomSymbolList):
    # Sets the mass of all atoms to zero.
    # Except atoms in unrestrictedAtomSymbolList and water

    atomIndicesToFreeze = []
    # Iterate over all atoms in the topology
    for atom in prmtop.topology.atoms():
        # Skip water, as massless objects cannot participate in openmm constraints (H is bound to Oxygen with constraint)
        if atom.residue.name in ['HOH', 'WAT']:
            continue
        # Check if the atom is not in the unrestricted list
        if atom.element.symbol not in unrestrictedAtomSymbolList:
            atomIndicesToFreeze.append(atom.index)
    # Set the mass of these atoms to 0.0 to freeze (constrain) them
    for index in atomIndicesToFreeze:
        system.setParticleMass(index, 0.0)

    return system

###########################################################################################
