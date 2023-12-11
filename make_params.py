from audioop import add
import argpass
import os
from os import path as p
from subprocess import call

## REQUIRES AMBERTOOLS

def inputs():
    parser = argpass.ArgumentParser()
    parser.add_argument("-pdb")
    parser.add_argument("-charge")
    parser.add_argument("-addH")
    args = parser.parse_args()
    pdbFile = args.pdb
    charge = str(args.charge)
    name = p.splitext(pdbFile)[0]
    addH = args.addH
    return pdbFile, name, charge, addH

def main():
    pdbFile, name, charge, addH = inputs()
    
    frcmodFile = make_params(pdbFile, name, charge, addH)

def make_params(pdbFile,name,charge, addH):
    if addH == "y":
        ## add hydrogens 
        pdbH = f"{name}_h.pdb"
        addHcommand = ["reduce",{pdbFile},">",pdbH]
        call(addHcommand)
    else:
        pdbH = pdbFile
    ## convert to mol2 with antechamber
    mol2File = f"{name}.mol2"
    antechamberCommand = ["antechamber",
                           "-i", pdbH,"-fi","pdb",
                            "-o", mol2File,"-fo", "mol2",
                            "-c", "bcc", "-s", "2", "-nc", charge]
    call(antechamberCommand)
    ## generate parameters with parmchk2
    frcmodFile = f"{name}.frcmod"
    parmchk2Command = ["parmchk2",
                       "-i", mol2File,"-f", "mol2",
                       "-o", frcmodFile]
    call(parmchk2Command)
    return frcmodFile
main()