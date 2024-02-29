## OPENMM LIBS
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

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
def restrain_all_atom_names_except_list(system, prmtop, inpcrd, unrestrictedAtomSymbolList):
    # assert all names in unrestrictedAtomNamesList are in an allowed list.

    restraint = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    system.addForce(restraint)
    restraint.addGlobalParameter("k", 1000.0 * kilojoules_per_mole / nanometer)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for atom in prmtop.topology.atoms():
        if atom.element.symbol not in unrestrictedAtomSymbolList:
            restraint.addParticle(atom.index, inpcrd.getPositions()[atom.index])

    return system

###########################################################################################
def heavy_atom_position_restraints(system, prmtop, inpcrd):
    return restrain_all_atom_names_except_list(system, prmtop, inpcrd, ["H"])
###########################################################################################