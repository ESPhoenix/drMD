## BASIC LIBS
import xml.etree.ElementTree as ET
from os import path as p
from pdbUtils import pdbUtils
import math
## OPENMM LIBS
import simtk.openmm as openmm
import simtk.unit as unit
## CUSTOM drMD LIBS
import instruments.drSelector as drSelector
###########################################################################################
def constraints_handler(system, prmtop, inpcrd, sim, saveFile, pdbFile):

    if "restraints" in sim:
        restraintInfo = sim["restraints"]

        kNumber = 0        
        for restraint in restraintInfo:
            selection = restraint["selection"]
            if restraint["type"] == "position":
                system = create_position_restraint(system, inpcrd, selection, kNumber, pdbFile)
            elif restraint["type"] == "distance":
                parameters = restraint["parameters"]
                system = create_distance_restraint(system, selection, parameters, kNumber, pdbFile)
            elif restraint["type"] == "angle":
                parameters = restraint["parameters"]
                system = create_angle_restraint(system, selection, parameters, kNumber, pdbFile)
            elif restraint["type"] == "torsion":
                parameters = restraint["parameters"]
                system = create_torsion_restraint(system, selection, parameters, kNumber, pdbFile)
            kNumber += 1
    
    else:
        print("Running with no restraints...")
        if p.splitext(saveFile)[1] == ".chk":
            return system
        clear_all_restraints(saveFile)
    return system
###########################################################################################
def create_position_restraint(system, inpcrd, selection, kNumber, pdbFile):
    ## create a position restraint
    positionRestraint = openmm.CustomExternalForce(f"k{str(kNumber)}*periodicdistance(x, y, z, x0, y0, z0)^2")
    positionRestraint.addGlobalParameter(f"k{str(kNumber)}", 1000.0 * unit.kilojoules_per_mole / unit.nanometer)
    positionRestraint.addPerParticleParameter("x0")
    positionRestraint.addPerParticleParameter("y0")
    positionRestraint.addPerParticleParameter("z0")
    system.addForce(positionRestraint)

    restraintAtomIndexes = drSelector.get_atom_indexes(selection, pdbFile)
    for restraintAtomIndex in restraintAtomIndexes:
        positionRestraint.addParticle(restraintAtomIndex,inpcrd.getPositions()[restraintAtomIndex])
    return system

###########################################################################################
def create_distance_restraint(system, selection, parameters, kNumber, pdbFile):
    ## Create a distance restraint
    k_force_constant = 1000.0 * unit.kilojoules_per_mole / unit.nanometer**2
    distanceRestraint = openmm.CustomBondForce(f"0.5 * k{str(kNumber)} * (r - r0)^2")
    ## add per bond parameters, k for force constant, r0 for desired distance
    distanceRestraint.addPerBondParameter(f"k{str(kNumber)}")
    distanceRestraint.addPerBondParameter("r0")
    ## add force to system
    system.addForce(distanceRestraint)

    # Get the indices of the two atoms to be restrained
    restraintAtomIndexes = drSelector.get_atom_indexes(selection, pdbFile)
    if len(restraintAtomIndexes) != 2:
        raise ValueError("Expected exactly two atom indices for a distance restraint.")

    # Get target distance from the parameters dictionary, and convert from angstroms to nanometers
    kForceConstant = parameters["k"]
    targetDistance_nm = parameters["r0"] * unit.angstroms * 0.1

    # Add the atom pair and the calculated target distance in nanometers to the bond restraint
    distanceRestraint.addBond(restraintAtomIndexes[0], restraintAtomIndexes[1], [kForceConstant, targetDistance_nm])

    return system
###########################################################################################
def create_angle_restraint(system, selection, parameters, kNumber, pdbFile):
    ## Create an angle restraint
    angleRestraint = openmm.CustomAngleForce(f"0.5 * k{str(kNumber)} * (theta - theta0)^2")

    ## add per angle parameters, k for force constant, theta0 for desired angle
    angleRestraint.addPerAngleParameter(f"k{str(kNumber)}")
    angleRestraint.addPerAngleParameter("theta0")

    ## add force to system
    system.addForce(angleRestraint)

    # Get the indices of the three atoms to be restrained
    restraintAngleAtoms = drSelector.get_atom_indexes(selection, pdbFile)
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
def create_torsion_restraint(system, selection, parameters, kNumber, pdbFile):
    ## Create a torsion restraint
    # Note: Ensure 'phi' is the variable used for torsion calculation in OpenMM
    torsionRestraint = openmm.CustomTorsionForce(f"0.5*k{str(kNumber)}*(1-cos(theta-theta0))")
    
    ## Add parameters: k for force constant, phi0 for desired torsion angle
    torsionRestraint.addPerTorsionParameter(f"k{str(kNumber)}")
    torsionRestraint.addPerTorsionParameter("theta0")
    
    ## Add force to system
    system.addForce(torsionRestraint)
    
    ## Get indices of the four atoms to be restrained
    restraintTorsionAtoms = drSelector.get_atom_indexes(selection, pdbFile)
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
def clear_all_restraints(saveXml):
    ## remove any leftover force constants 
    tree = ET.parse(saveXml)
    root = tree.getroot()
    parametersElement = root.find("Parameters")

    # Safely remove any custom force constants
    paramsToPop = []
    if parametersElement is not None:
        for param in parametersElement.attrib:
            if param.startswith("k"):
                paramsToPop.append(param)
    for param in paramsToPop:
        parametersElement.attrib.pop(param)
    tree.write(saveXml)
    
###########################################################################################