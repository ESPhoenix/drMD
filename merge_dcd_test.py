## basic python libraries
import os
from os import path as p
from functools import wraps
import pandas as pd
## openMM libraries
import openmm
from openmm import app
from openmm import OpenMMException
import  simtk.unit  as unit
from openmm.app.dcdfile import DCDFile


import mdtraj as md
## drMD libraries
from instruments import drSim
from instruments import drSpash
## clean code
from typing import Tuple, Union, Dict, List
from os import PathLike
####################################################################
def merge_dcd_files(dcdFiles: list[Union[PathLike, str]],
                    pdbFile: Union[PathLike, str],
                    outputDcd: Union[PathLike, str]) -> Union[PathLike, str]:
    
    traj = md.load_dcd(dcdFiles[0], top = pdbFile)
    print("first frame")
    print(traj.n_frames)
    for file in dcdFiles[1:]:

        newTraj = md.load_dcd(file, top = pdbFile)
        print("new frame")

        print(newTraj.n_frames)

        traj = traj + newTraj
        print("concat frame")

        print(traj.n_frames)



    traj.save_dcd(outputDcd)


####################################################################
simDir = "/home/esp/scriptDevelopment/drMD/03_outputs/4a29_FMN_1/03_this_will_explode"

dcdFiles = [p.join(simDir,file) for file in os.listdir(simDir) if p.splitext(file)[1] == ".dcd" and "partial" in file]
dcdFiles = sorted(dcdFiles)
dcdFiles.append(p.join(simDir, "trajectory.dcd"))

merge_dcd_files(dcdFiles, pdbFile= "/home/esp/scriptDevelopment/drMD/03_outputs/4a29_FMN_1/03_this_will_explode/4a29_FMN_1.pdb", outputDcd= p.join(simDir, "merged.dcd") )