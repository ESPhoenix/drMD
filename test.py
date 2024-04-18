from instruments import drFixer
import os
from os import path as p



targetPdb = "/home/esp/scriptDevelopment/drMD/01_inputs/Con1_wildtype.pdb"
inputPdb = "/home/esp/scriptDevelopment/drMD/03_outputs/Con1_wildtype/01_energy_minimisation/minimised_geom.pdb"

newDf = drFixer.openMM_pdb_fix(targetPdb, inputPdb)