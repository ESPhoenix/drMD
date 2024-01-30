## OPENMM LIBS
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


###########################################################################################
def restrain_all_atom_names_except_list(system, prmtop, inpcrd, unrestrictedAtomNamesList):
    # assert all names in unrestrictedAtomNamesList are in an allowed list.
    
    restraint = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    system.addForce(restraint)
    restraint.addGlobalParameter("k", 1000.0 * kilojoules_per_mole / nanometer)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for atom in prmtop.topology.atoms():
        if not atom.element.name in unrestrictedAtomNamesList:
            restraint.addParticle(atom.index, inpcrd.getPositions()[atom.index])

    return system

def heavy_atom_position_restraints(system, prmtop, inpcrd):
    return restrain_all_atom_names_except_list(
        system, prmtop, inpcrd, ["hydrogen"]
    )