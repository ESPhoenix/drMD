import os
from os import path as p
from pdbUtils import pdbUtils

def inputs():
    prepDir = "/home/esp/scriptDevelopment/drMD/04_PET_proj_outputs/6eqe_PET-tetramer_1/00_prep/WHOLE"
    simDir = "/home/esp/scriptDevelopment/drMD/04_PET_proj_outputs/6eqe_PET-tetramer_1/031_production_md"
    trajPdb = p.join(simDir,"trajectory.pdb")

    return simDir, trajPdb

def main():
    simDir, trajPdb = inputs()
    trajDf = pdbUtils.pdb2df(trajPdb)
    print(trajDf)


if __name__ == "__main__":
    main()