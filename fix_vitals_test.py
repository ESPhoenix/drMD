from simtk import unit
import openmm
import pandas as pd
from typing import Dict, Union
import os
from os import path as p
from os import PathLike



#######################################################################
def merge_partial_outputs(simDir: Union[PathLike, str], prmtop: Union[PathLike, str], simInfo: Dict) -> None:
    """
    After firstAid protocols have been run, we need to merge the partial reports and trajectories

    Args:
        simDir (Union[PathLike, str]): The path to the simulation directory
        pdbFile (Union[PathLike, str]): The path to the pdb file

    Returns:
        None
    """
    print("-->\tMerging partial outputs...")
    ## merge vitals reports
    vitalsDf = merge_partial_reports(simDir, "vitals_report", removePartials=False)
    vitalsDf = fix_merged_vitals(vitalsDf, simInfo)
    vitalsDf.to_csv(p.join(simDir, "vitals_report.csv"))
    ## merge progress reports
    merge_partial_reports(simDir, "progress_report", removePartials=False)
    ## merge trajectories
    # merge_partial_trajectories(simDir, prmtop, removePartials=True)
#######################################################################

def fix_merged_vitals(vitalsDf: pd.DataFrame, simInfo: Dict) -> pd.DataFrame:

    ## read stuff from simInfo
    logInterval_ps: int = simInfo["logInterval"]
    timeStep: openmm.Quantity = simInfo["timestep"]
    duration: openmm.Quantity = simInfo["duration"]

    ## convert to ints 
    duration_ps: int = duration.value_in_unit(unit.picoseconds)
    timeStep_ps: int = timeStep.value_in_unit(unit.picoseconds)
    timeRange_ps = range(logInterval_ps, duration_ps + logInterval_ps, logInterval_ps)

    stepsPerLog = int(logInterval_ps / timeStep_ps)

    stepsRange = [val * stepsPerLog for val in timeRange_ps]    

    vitalsDf['Time (ps)'] = timeRange_ps
    vitalsDf['#"Step"'] = stepsRange
    return vitalsDf
#######################################################################
def merge_partial_reports(simDir: Union[PathLike, str], matchString: str, removePartials: bool = False) -> None:
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
        if ("partial" in file) and (file.startswith(matchString)) and (p.splitext(file)[1] == ".csv"): #and (p.getsize(p.join(simDir, file)) > 0):
            reports.append(p.join(simDir, file))

    print(reports)

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





simDir = "/home/esp/scriptDevelopment/drMD/03_outputs/4a29_FMN_1/03_this_will_explode"

simInfo = {
    "logInterval": 10 ,
    "duration": 300 * unit.picoseconds,
    "timestep": 2 * unit.femtoseconds}

merge_partial_outputs(simDir, None, simInfo)