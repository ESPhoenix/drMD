import os
from os import path as p
import pandas as pd
import mdtraj as md


def main():
    dcdFile = "/scratch/collated_xml/Equlibriation_step/cvFAP_T169C_PLM_FAD_2_eq.dcd"
    pdbFile = "/scratch/collated_xml/topologies/cvFAP_T169C_PLM_FAD_2.pdb"
    check_vitals(dcdFile, pdbFile)

def  check_vitals(dcdFile,pdbFile):
    traj = md.load_dcd(dcdFile, top=pdbFile)
    print(traj.temperature)
    exit()    
    rmsd = md.rmsd(traj, traj)
    density = md.density(traj)

    print(traj.temperature)
main()