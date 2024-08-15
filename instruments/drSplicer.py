import os
from os import path as p
from typing import Union, Dict
import pandas as pd

## OPEN MM LIBS
import openmm
from openmm import unit

from instruments.drCustomClasses import FilePath, DirectoryPath
from instruments import drLogger
#####################################################################################
def merge_partial_outputs(simDir: DirectoryPath, prmtop: FilePath, simInfo: Dict) -> None:
    """
    After firstAid protocols have been run, we need to merge the partial reports and trajectories

    Args:
        simDir (Union[PathLike, str]): The path to the simulation directory
        pdbFile (Union[PathLike, str]): The path to the pdb file

    Returns:
        None
    """
    drLogger.log_info(f"-->{' '*4}Merging partial outputs...", True)
    ## merge vitals reports
    vitalsDf = merge_partial_reports(simDir, "vitals_report", removePartials=True)
    vitalsDf = fix_merged_vitals(vitalsDf, simInfo)
    vitalsDf.to_csv(p.join(simDir, "vitals_report.csv"))
    ## merge progress reports
    merge_partial_reports(simDir, "progress_report", removePartials=True)
    ## merge trajectories
    merge_partial_trajectories(simDir, prmtop, removePartials=True)

#####################################################################################
def merge_partial_reports(simDir: DirectoryPath, matchString: str, removePartials: bool = False) -> None:
    """
    Merges partial reports into a single report

    Args:
        simDir (Union[PathLike, str]): The path to the simulation directory
        matchString (str): The string to match
        removePartials (bool, optional): Whether to remove partial reports. Defaults to False.  
    """
    ## rename the last report to be made
    lastReport = p.join(simDir, f"{matchString}.csv")
    if p.isfile(lastReport):
        os.rename(lastReport, p.join(simDir, f"{matchString}_partial_99.csv"))

    ## collect all partial reports into a list exept for last one written
    reports = []
    for file in os.listdir(simDir):
        if ("partial" in file) and (file.startswith(matchString)) and (p.splitext(file)[1] == ".csv") and (p.getsize(p.join(simDir, file)) > 0):
            reports.append(p.join(simDir, file))

    ## sort the list so that reports are in chronological order
    reports = sorted(reports)
    ## merge reports
    dfsToConcat = []
    for report in reports:
        dfsToConcat.append(pd.read_csv(report))

    ## concatonate | write back to csv
    ## remove partial reports to tidy up 
    [os.remove(report) for report in reports if removePartials]
    ## concat dataframes and write to file
    df = pd.concat(dfsToConcat, ignore_index=True)

    return df
#######################################################################
def merge_dcd_files(dcdFiles: list[FilePath],
                    pdbFile: FilePath,
                    outputDcd: FilePath) -> FilePath:
    
    traj = md.load_dcd(dcdFiles[0], top = pdbFile)
    for file in dcdFiles[1:]:
        newTraj = md.load_dcd(file, top = pdbFile)
        traj = traj + newTraj
    traj.save_dcd(outputDcd)
#######################################################################
def merge_partial_trajectories(simDir: DirectoryPath,
                                pdbFile: FilePath,
                                  removePartials: bool = False) -> None:
    """
    Merges partial trajectories into a single trajectory

    Args:
        simDir (Union[PathLike, str]): The path to the simulation directory
        pdbFile (Union[PathLike, str]): The path to the pdb file
        removePartials (bool, optional): Whether to remove partial trajectories. Defaults to False.  
    """
    ## rename last trajectory to be made
    lastTrajectory = p.join(simDir, "trajectory.dcd")
    if p.isfile(lastTrajectory):
        os.rename(lastTrajectory, p.join(simDir, "trajectory_partial_99.dcd"))

    ## collect all partial trajectories into a list exept for last one written
    ## collect all partial reports into a list exept for last one written
    trajectories = []
    for file in os.listdir(simDir):
        if ("partial" in file) and (file.startswith("trajectory")) and (p.splitext(file)[1] == ".dcd"):
            trajectories.append(p.join(simDir, file))

    ## sort the list so that trajectories are in chronological order
    trajectories = sorted(trajectories)


    ## merge trajectories
    merge_dcd_files(trajectories, pdbFile,  lastTrajectory)

    ## delete partial trajectories
    if removePartials:
        for trajectory in trajectories:
            if p.basename(trajectory) == "trajectory.dcd":
                continue
            os.remove(trajectory)

#######################################################################
def fix_merged_vitals(vitalsDf: pd.DataFrame, simInfo: Dict) -> pd.DataFrame:

    ## read stuff from simInfo
    logInterval_ps: int = round(simInfo["logInterval"])
    timeStep: openmm.Quantity = simInfo["timestep"]
    duration: openmm.Quantity = simInfo["duration"]

    ## convert to ints 
    duration_ps: int = round(duration.value_in_unit(unit.picoseconds))
    timeStep_ps: int = round(timeStep.value_in_unit(unit.picoseconds))
    timeRange_ps = range(logInterval_ps, duration_ps + logInterval_ps, logInterval_ps)

    stepsPerLog = int(logInterval_ps / timeStep_ps)

    stepsRange = [val * stepsPerLog for val in timeRange_ps]    

    vitalsDf['Time (ps)'] = timeRange_ps
    vitalsDf['#"Step"'] = stepsRange
    return vitalsDf
#######################################################################