## basic python libraries
import os
from os import path as p
from functools import wraps

## openMM libraries
import openmm
from openmm import app
from openmm import OpenMMException
import  simtk.unit  as unit
## drMD libraries
from instruments import drSim
from instruments import drOperator

## clean code
from typing import Tuple, Union
from os import PathLike


#######################################################################
def triage_handler(fallback_function, max_retries=5):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            retries = 0
            while retries <= max_retries:
                try:
                    saveFile = func(*args, **kwargs)
                    if retries > 1:
                        print(f"Success after {retries} tries.")
                        merge_partial_outputs(kwargs["outDir"], kwargs["refPdb"])
                    return saveFile
                except (OpenMMException, ValueError) as e:
                    retries += 1
                    print(f"Attempting fallback triage, try {retries} of {max_retries}")


                    outDir = kwargs["outDir"]

                    triageDir, explodedCheckpoint = pre_triage_processing(outDir, kwargs["sim"], retries)

                    kwargs["outDir"] = triageDir
                    kwargs["saveFile"] = explodedCheckpoint
                    kwargs["triageTries"] = retries

                    saveFile = fallback_function(*args, **kwargs)


                    kwargs["outDir"] = outDir
                    del kwargs["triageTries"]
               
            print("Max retries reached. Stopping.")
            return None
        return wrapper
    return decorator
#######################################################################
def pre_triage_processing(outDir: Union[PathLike, str], sim, retries: int) -> Tuple:
    simDir = p.join(outDir, sim["stepName"])

    ## look for checkpoint file of exploded simulation
    explodedCheckpoint = p.join(simDir, "checkpoint.chk")
    if not p.isfile(explodedCheckpoint):
        raise FileNotFoundError(f"Checkpoint file not found at {explodedCheckpoint}")
    
    ## make a triage directory - this will be where each triage step is performed
    triageDir = p.join(simDir, "triage")
    os.makedirs(triageDir, exist_ok=True)

    ## rename output files in exploded directory
    rename_output_files(simDir, retries)

    ## find renamed checkpoint
    explodedCheckpoint = p.join(simDir, f"checkpoint_partial_{str(retries)}.chk")

    return triageDir, explodedCheckpoint
#######################################################################

def rename_output_files(outDir: Union[PathLike, str], retries: int) -> None:
    extensions = [".dcd", ".csv", ".chk"]
    for file in os.listdir(outDir):
        if "partial" in file:
            continue
        if p.splitext(file)[1] in extensions:
            os.rename(p.join(outDir, file),
                       p.join(outDir,f"{p.splitext(file)[0]}_partial_{str(retries)}{p.splitext(file)[1]}"))
#######################################################################



def merge_partial_outputs(simDir, pdbFile) -> None:
    ...
#######################################################################
def run_triage_simulation(prmtop: str, inpcrd: str, sim: dict, saveFile: str, outDir: str, platform: openmm.Platform, refPdb: str, triageTries: int) -> None:

    triageSimInfo =  {
        "stepName" : f"triage_step_{triageTries}",
        "simulationType": "NPT",
        "duration" : "10 ps",
        "timestep" : "0.5 ps",
        "temp" : 300,
        "logInterval" : "2 ps",
        "maxIterations" : 1000
        }
    triageSimInfo = drSim.process_sim_data(triageSimInfo)

    ## run triage simulation
    triageSaveFile = drSim.run_energy_minimisation(prmtop = prmtop,
                                                   inpcrd = inpcrd,
                                                   sim = triageSimInfo,
                                                    saveFile= saveFile,
                                                     outDir = outDir,
                                                     platform = platform,
                                                      refPdb= refPdb)

    # triageSaveFile = drSim.run_molecular_dynamics(prmtop = prmtop,
    #                                                inpcrd = inpcrd,
    #                                                sim = triageSimInfo,
    #                                                 saveFile= saveFile,
    #                                                  outDir = outDir,
    #                                                  platform = platform,
    #                                                   refPdb= refPdb)

    return triageSaveFile