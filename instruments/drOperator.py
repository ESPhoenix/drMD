## BASIC LIBS
import os
from os import path as p
## drMD UTILS
from instruments import drPrep
## drMD simulation
from instruments import drSim
from pdbUtils import pdbUtils
import instruments.drConfigInspector as drConfigInspector
#####################################################################################
def drMD_protocol(configYaml):
    config = drConfigInspector.read_config(configYaml)
    outDir = config["pathInfo"]["outputDir"]
    os.makedirs(outDir,exist_ok=True)
        
    mergedPdb, inputCoords, amberParams = drPrep.prep_protocol(config)

    drSim.run_simulation(config = config,
                   outDir = outDir,
                   inputCoords=inputCoords,
                   amberParams=amberParams,
                   pdbFile = mergedPdb)
#####################################################################################
