## BASIC PYTHON LIBRARIES
import os
from os import path as p
import pandas as pd
import numpy as np
from functools import wraps
import textwrap
from shutil import move
import warnings

## PLOTTING LIBRARIES
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

## PDF LIBS
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML, CSS

## MDTRAJ AND MDANALYSIS LIBRARIES
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import rms

## drMD LIBRARIES
from ExaminationRoom import drLogger
from UtilitiesCloset import drSplicer, drSelector

## CLEAN CODE
from typing import Union, Dict, Tuple, List
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

## PDB // DATAFRAME UTILS
from pdbUtils import pdbUtils

## DISABLE WARNINGS
import logging
logging.getLogger('weasyprint').setLevel(logging.ERROR)
warnings.filterwarnings('ignore')

######################################################################
def check_vitals(simDir: DirectoryPath,
                  vitalsFiles: Dict[str, FilePath]) -> None:
    """
    Looks at the following properties of the simulation:
        - RMSD
        - Energies (Potential, Kinetic, Total)
        - Properties (Temperature, Volume, Density)
    Checks that each of these properties have converged
    
    
    Also gathers time data for the simulation
    Plots these data and creates a report PDF

    Args:
        simDir (DirectoryPath): The directory of the simulation
        vitalsFiles (Dict[str, FilePath]): A dictionary of the vitals files
    """
    drLogger.log_info(f"-->{' '*4}Checking vitals...")
    ## read OpenMM reporters into dataframes
    vitalsDf: pd.DataFrame = pd.read_csv(vitalsFiles["vitals"])
    progressDf: pd.DataFrame = pd.read_csv(vitalsFiles["progress"])

    ## skip this if the reporters are not present, or are less than 5 entries long
    if len(vitalsDf) < 5:
        return
    
    if len(progressDf) < 5:
        return

    ## use mdtraj to calculate RMSD for non water and ions, get that data into a dataframe
    rmsdDf: pd.DataFrame = calculate_rmsd_mda(trajectoryDcd = vitalsFiles["trajectory"],
                                                trajectoryPdb = vitalsFiles["pdb"])
    ## plot RMSD as a function of time
    rmsdPng: FilePath = plot_rmsd(rmsdDf, simDir, "Backbone_RMSD")
    ## get time data, plot a table
    timeDf: pd.DataFrame = extract_time_data(vitalsDf, progressDf)
    timePng: FilePath = plot_time_data(timeDf, simDir)
    
    ## plot energy as a function of time 
    energyList: list = ["Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)"]
    energyPlot: FilePath = plot_vitals(vitalsDf, simDir, energyList, "Energies")
    ## plot properties as a function of time 
    propertiesList: list = ["Temperature (K)", "Box Volume (nm^3)", "Density (g/mL)"]
    propertiesPlot: FilePath = plot_vitals(vitalsDf, simDir, propertiesList, "Properties")

    ## combine all feature names
    allFeaturesList: list = energyList + propertiesList + ["Backbone RMSD"]
    ## combine all features into one dataframe
    allFeaturesDf: pd.DataFrame = pd.concat([vitalsDf,rmsdDf],axis=1)
    ## check convergance for all features
    converganceDf: pd.DataFrame  = check_convergance(allFeaturesDf,allFeaturesList)
    ## plot convergance booleans as a table
    convergancePng: FilePath = plot_converged(simDir, converganceDf)
    ## create report PDF using all the plots created above
    create_vitals_pdf(simDir)
    ## tidy up reporters to avoid clutter (dont do this while testing!)
    if not __name__ == "__main__":
        tidy_up(simDir)
######################################################################
def tidy_up(simDir: DirectoryPath):
    """
    Tidies up the reporters in the simulation directory

    Args:
        simDir (DirectoryPath): The directory of the simulation
    """
    ## make a new directory to tidy up reporters and png files
    tidyDir = p.join(simDir, "00_reporters_and_plots")
    os.makedirs(tidyDir, exist_ok=True)
    ## move all reporters and png files to the new directory
    for file in os.listdir(simDir):
        if p.splitext(file)[1] in [".csv", ".png"]:
            move(p.join(simDir, file), p.join(tidyDir, file))

    plt.close("all")
######################################################################
def create_vitals_pdf(simDir: DirectoryPath):
    """
    Uses Jinja2 to create a vitals report PDF from 
    the PNG files generated in this script

    Args:   
        simDir (DirectoryPath): The directory of the simulation
    """
    ## get the instruments directory path
    instrumentsDir: DirectoryPath = p.dirname(__file__)
    env: Environment = Environment(loader=FileSystemLoader(instrumentsDir))
    template = env.get_template("vitals_template.html")
    
    # Render the template with any context variables you need
    context = {
        # Add your context variables here
    }
    rendered_html = template.render(context)
    
    # Generate the PDF with error handling
    outPdf: FilePath = p.join(simDir, "vitals_report.pdf")
    try:
        base_url = simDir # Set the base URL to the current working directory
        HTML(string=rendered_html, base_url=base_url).write_pdf(outPdf, stylesheets=[CSS(string='@page { margin: 0; }')])
    except Exception as e:
        drLogger.log_info(f"-->{' '*4}Error generating vitals report: {e}")
######################################################################
def check_up_handler():
    def decorator(simulationFunction):
        @wraps(simulationFunction)
        def wrapper(*args, **kwargs):
            saveFile: FilePath = simulationFunction(*args, **kwargs)

            vitalsFiles, simDir = find_vitals_files(simInfo=kwargs["sim"],
                                                     outDir=kwargs["outDir"], 
                                                     pdbFile=kwargs["refPdb"],
                                                     config=kwargs["config"])

            try:
                check_vitals(simDir = simDir,
                            vitalsFiles = vitalsFiles)
            except Exception as e:
                drLogger.log_info(f"-->{' '*4}Error running checkup: {e}", True)
            return saveFile
        return wrapper
    return decorator
######################################################################
def find_vitals_files(simInfo: Dict,
                       outDir: DirectoryPath,
                       pdbFile: FilePath,
                       config: Dict) -> Tuple[Dict[str, FilePath], DirectoryPath]:
    
    ## get the simulation directory
    simDir: DirectoryPath= p.join(outDir, simInfo["stepName"])

    ## check to see if multiple partial trajectories exist
    trajectoryDcds: list = [p.join(simDir, file) for file in os.listdir(simDir) if file.endswith(".dcd")]
    ## if more than one trajectory is found, merge partial outputs before running the health check
    if len(trajectoryDcds) > 1:
        drSplicer.merge_partial_outputs(simDir = simDir,
                                            pdbFile=pdbFile,
                                            config=config,
                                            simInfo = simInfo)

    ## find vitals reporter file
    vitalsReport: FilePath = p.join(simDir, "vitals_report.csv")
    if not p.isfile(vitalsReport):
        raise FileNotFoundError(f"->\tReporter file not found at {vitalsReport}")
    ## find progress reporter file
    progressReport: Union[PathLike, str] = p.join(simDir, "progress_report.csv")
    if not p.isfile(progressReport):
        raise FileNotFoundError(f"->\tReporter file not found at {progressReport}")
    ## find trajectory file
    trajectoryDcd: Union[PathLike, str] = p.join(simDir, "trajectory.dcd")
    if not p.isfile(trajectoryDcd):
        raise FileNotFoundError(f"->\Trajectory file not found at {trajectoryDcd}")
    
    ## find the pdb file
    pdbFile = p.join(simDir, "trajectory.pdb")
    if not p.isfile(pdbFile):
        raise FileNotFoundError(f"->\tPDB file not found at {simDir}")
    
    ## collect files into a dictionary
    vitalsFiles = {"vitals": vitalsReport, "progress": progressReport, "trajectory": trajectoryDcd, "pdb": pdbFile}
    
    return vitalsFiles, simDir

######################################################################
def cusum_test(series: pd.Series, threshold: float=0.05):
    """
    Performs the Cumulative Sum Test

    Args:
        series (pd.Series): The series to test
        threshold (float, optional): The threshold for the test. Defaults to 0.05.
    """
    ## remove any inf and -inf values
    seriesClean = series.replace([np.inf, -np.inf], np.nan).dropna()
    ## check to see if series is long enough to perform the test
    if len(seriesClean) < 2:
        return False
    ## perform the test
    cumSumSeries = np.cumsum(seriesClean - seriesClean.mean())

    return np.max(np.abs(cumSumSeries)) < threshold

######################################################################
def check_convergance(df: pd.DataFrame, columns: list, windowSize: int = 3) -> pd.DataFrame:
    """
    Checks for convergence in the dataframe

    Args:
        df (pd.DataFrame): The dataframe to check
        columns (list): The columns to check
        windowSize (int, optional): The window size for the rolling average. Defaults to 5.

    Returns:
        pd.DataFrame: The dataframe with the convergence status
    """

    # define the conversion tolerances TODO: make these a bit more scientific
    conversionTolerances = {
        "Potential Energy (kJ/mole)" : 3000, 
        "Kinetic Energy (kJ/mole)": 3000,
        "Total Energy (kJ/mole)": 3000,
        "Temperature (K)": 5,
        "Box Volume (nm^3)" : 1,
        "Density (g/mL)" : 1,
        "Backbone RMSD" : 3                     

    }


    ## create a dictionary to store the results
    convergedDict = {}
    ## loop through the columns in our dataframe

    for column in columns:
        ## skip if column length is 1 (we can't do any checks on that!)
        if len(column) == 1:
            continue
        ## get the rolling average
        runningAverage = df[column].rolling(window=windowSize).mean()
        ## check for convergence using cumsum test
        isConverged = cusum_test(runningAverage, threshold=conversionTolerances[column])
        ## update the dictionary
        entryLabel = column.split("(")[0].split()[0]
        convergedDict.update({entryLabel:isConverged})
    ## convert results to dataframe
    convergedData = [(key, value) for key, value in convergedDict.items()]
    convergedDf = pd.DataFrame(convergedData, columns=["Property", "Converged"])    

    return convergedDf


######################################################################
def plot_converged(simDir: DirectoryPath, convergedDf: pd.DataFrame) -> FilePath:
    """
    Plots a table of ticks and crosses for each feature we have tested for convergance

    Args:
        simDir (DirectoryPath): The simulation directory
        convergedDf (pd.DataFrame): The dataframe with the convergence status
    """

    ## create the plot
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.axis('off')
    table = ax.table(cellText=convergedDf.values, colLabels=convergedDf.columns, loc='center',
                     cellLoc='center', rowLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.5, 2)

    brightGreen = '#00FF00'  # Brighter green
    brightRed = '#FF0000'  # Brighter red
    darkGrey = '#1a1a1a'  # Much darker grey

    fig.patch.set_facecolor(darkGrey)
    ax.set_facecolor(darkGrey)

    # Replace True/False with ticks and crosses in the second column
    for i in range(1, len(convergedDf) + 1):
        if convergedDf.iloc[i - 1, 1] == True:
            table[(i, 1)].get_text().set_text('✔')
            table[(i, 1)].get_text().set_color(brightGreen)
        else:
            table[(i, 1)].get_text().set_text('✘')
            table[(i, 1)].get_text().set_color(brightRed)

    for key, cell in table.get_celld().items():
        cell.set_linewidth(0.5)
        cell.set_edgecolor(brightGreen)
        cell.set_facecolor(darkGrey)
        if key[1] != 1:  # Skip the second column for color setting
            cell.get_text().set_color(brightGreen)

    for i in range(len(convergedDf.columns)):
        table[(0, i)].set_facecolor(brightGreen)
        table[(0, i)].set_alpha(0.3)  # Set low alpha value
        table[(0, i)].get_text().set_color(brightGreen)
        table[(0, i)].get_text().set_fontweight('bold')


    # Save the plot as a PNG file
    savePng = p.join(simDir, f"convergance_checks.png")
    plt.savefig(savePng, bbox_inches='tight', facecolor=darkGrey)
    plt.close()
    return savePng



######################################################################
def convert_seconds(seconds: int) -> str:
    """
    Converts seconds as an int to HH:MM:SS format

    Args:
        seconds (int): The number of seconds
    """
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02}:{:02}:{:02}".format(int(hours), int(minutes), int(seconds))

######################################################################
def extract_time_data(vitalsDf: pd.DataFrame, progressDf: pd.DataFrame) -> pd.DataFrame:
    """
    Reads vitals and progress data to extract key time-related data
    
    Args:
        vitalsDf (pd.DataFrame): The vitals dataframe
        progressDf (pd.DataFrame): The progress dataframe

    Returns:
        timeDf (pd.DataFrame): The time dataframe
    """

    # get elapsed time
    timeElapsed: int = progressDf["Elapsed Time (s)"].tail(1).values[0]
    timeElapsed: str = convert_seconds(timeElapsed)
    # get average speed
    averageSpeed: str = str(round(progressDf["Speed (ns/day)"].mean(),2))
    # get nSteps
    nSteps: str = str(vitalsDf["#\"Step\""].tail(1).values[0])
    # get simulation duration
    duration: str = str(round(vitalsDf["Time (ps)"].tail(1).values[0],2))

    # create a DataFrame
    timeDf = pd.DataFrame({'Elapsed Time (HH:MM:SS)': [timeElapsed],
                            'Average Speed (ns/day)': [averageSpeed],
                            'Total Steps': [nSteps],
                            'Simulation Duration (ps)': [duration]})

    timeDf: pd.DataFrame = timeDf.transpose()
    timeDf.columns = ['Value']  # rename the column
    timeDf.reset_index(level=0, inplace=True)  # reset the index

    return timeDf

######################################################################
def plot_time_data(timeDf: pd.DataFrame, outDir: FilePath):
    """
    Plots time data into a table
    
    Args:
        timeDf (pd.DataFrame): The time dataframe
        outDir (FilePath): The output directory
    """

    ## set up plot size and convert cm to inches
    widthInch: float = 8.94 / 2.54
    heightInch: float = 6.09 / 2.54
    ## set up plot
    fig, ax = plt.subplots(figsize=(widthInch, heightInch))
    fig.patch.set_facecolor('#1a1a1a')  # Much darker grey
    brightGreen: str = '#00FF00'  # Brighter green
    ax.axis('off')
    
    # Wrap text in each cell to ensure it fits within the cell
    wrappedText: str = timeDf.map(lambda x: '\n'.join(textwrap.wrap(str(x), width=15)))

    # Create the table without column labels
    table = ax.table(cellText=wrappedText.values, cellLoc='center', loc='center', colLabels=None)
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.2)

    # Set a fixed height for all rows to be two lines deep
    for key, cell in table.get_celld().items():
        cell.set_height(0.2)  # Adjust this value as needed for two lines of text

    # Set table colors
    for key, cell in table.get_celld().items():
        cell.set_edgecolor(brightGreen)
        cell.set_facecolor('#1a1a1a')
        cell.get_text().set_color(brightGreen)
    plt.margins(0.01, 0.01)

    # Save the plot as a PNG image with reduced border
    savePng: FilePath = p.join(outDir, "time_info.png")
    # Adjust layout to remove extra space around the table
    plt.savefig(savePng, bbox_inches="tight", pad_inches=0.01, facecolor='#1a1a1a')
    plt.close()

    return savePng


######################################################################

def plot_vitals(vitalsDf: pd.DataFrame,
                 outDir: FilePath,
                   yData:pd.Series,
                     tag: str) -> None:
    """
    Plots vitals energy or properties into a 1x3 subplot of line graphs
    
    Args:
        vitalsDf (pd.DataFrame): The vitals dataframe
        outDir (FilePath): The output directory
        yData (pd.Series): The y data
        tag (str): The tag
    """
    # Convert cm to inches
    widthInch: float = 29.7 / 2.54
    heightInch: float = 10 / 2.54
    
    # Set up the figure and axis
    fig, axes = plt.subplots(nrows=1, ncols=len(yData), figsize=(widthInch, heightInch), constrained_layout=True)
    fig.patch.set_facecolor('#1a1a1a')  # Much darker grey
    brightGreen: str = '#00FF00'  # Brighter green
    fig.suptitle(f'Simulation {tag} vs Time', fontsize=12, y=1.05, color=brightGreen)  # Adjust y to move the title up
    
    # Ensure axes is iterable even if there's only one plot
    if len(yData) == 1:
        axes = [axes]
    
    for i, ax in enumerate(axes):
        lineColor = brightGreen
        ax.set_facecolor('#1a1a1a')  # Much darker grey
        ax.set_xlabel('Time (ps)', fontsize=10, color=brightGreen)
        ax.set_ylabel(yData[i], color=brightGreen, fontsize=10, labelpad=5)
        ax.plot(vitalsDf['Time (ps)'], vitalsDf[yData[i]], label=yData[i],
                linestyle='-', color=lineColor, linewidth=1)
        ax.tick_params(axis='y', labelcolor=brightGreen, labelsize=8)
        ax.tick_params(axis='x', labelcolor=brightGreen, labelsize=8)
        for label in ax.get_xticklabels():
            label.set_fontsize(8)
            label.set_color(brightGreen)
        for label in ax.get_yticklabels():
            label.set_fontsize(8)
            label.set_color(brightGreen)
        ax.grid(True, linestyle='--', alpha=0.7, color=brightGreen,
                linewidth=0.5, which='both')
        ax.minorticks_on()  # Enable minor ticks
        ax.grid(which='minor', linestyle=':', linewidth=0.5, color=brightGreen)
        legend = ax.legend(loc='best', fontsize=8, facecolor='#1a1a1a', edgecolor=brightGreen)
        for text in legend.get_texts():
            text.set_color(brightGreen)
    
    # Save the plot as a PNG image
    savePng = p.join(outDir, f"Vitals_{tag}.png")
    plt.savefig(savePng, bbox_inches="tight", facecolor='#1a1a1a')
    plt.close(fig)
    return savePng


######################################################################

def plot_rmsd(rmsdDf: pd.DataFrame,
               outDir: DirectoryPath,
                   tag: str) -> FilePath:
    """
    Plots rmsd trace as a line graph

    Args:
        rmsdDf (pd.DataFrame): The rmsd dataframe
        outDir (DirectoryPath): The output directory
        tag (str): The tag

    Returns:
        FilePath: The path to the plot  
    """

    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(3.54, 3.54))
    fig.patch.set_facecolor('#1a1a1a')  # Much darker grey
    brightGreen: str = '#00FF00'  # Brighter green
    brightRed: str = '#FF0000'  # Brighter red
    fig.suptitle(f'Simulation {tag} vs Time', fontsize=12, y=0.98, color=brightGreen)
    
    ax.set_facecolor('#1a1a1a')  # Much darker grey
    ax.set_xlabel('Timestep (ps)', fontsize=10, color=brightGreen)
    ax.set_ylabel('Backbone RMSD', color=brightGreen, fontsize=10, labelpad=15)

    ax.plot(rmsdDf['Timestep (ps)'], rmsdDf["Backbone RMSD"],
            linestyle='-', color=brightRed, linewidth=1)
    ax.tick_params(axis='y', labelcolor=brightGreen, labelsize=8)
    ax.tick_params(axis='x', labelcolor=brightGreen, labelsize=8)
    for label in ax.get_xticklabels():
        label.set_fontsize(8)
        label.set_color(brightGreen)
    for label in ax.get_yticklabels():
        label.set_fontsize(8)
        label.set_color(brightGreen)
    ax.grid(True, linestyle='--', alpha=0.7, color=brightGreen,
            linewidth=0.5, which='both')
    ax.minorticks_on()  # Enable minor ticks
    ax.grid(which='minor', linestyle=':', linewidth=0.5, color=brightGreen)
    
    
    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    # save
    savePng = p.join(outDir, f"Vitals_{tag}.png")
    plt.savefig(savePng, bbox_inches="tight", facecolor='#1a1a1a')
    plt.close()
    return savePng


#########################################################################################################
def calculate_rmsd_mda(trajectoryDcd: FilePath, trajectoryPdb: FilePath) -> pd.DataFrame:

    universe = mda.Universe(trajectoryPdb, trajectoryDcd)

    rmsdCalculation = mda.analysis.rms.RMSD(universe, select="backbone")

    rmsdCalculation.run()


    rmsdDf = pd.DataFrame(rmsdCalculation.rmsd)
    rmsdDf.columns = ["Frame Index", "Timestep (ps)", "Backbone RMSD"]

    rmsdDf.drop(columns=["Frame Index"], inplace=True)
    rmsdDf['Backbone RMSD'] = rmsdDf['Backbone RMSD'].round(2)

    return rmsdDf



######################################################################
def calculate_rmsd(trajectoryDcd: FilePath, pdbFile: FilePath, trajectorySelections: List[Dict], outDir: DirectoryPath) -> pd.DataFrame:
    """
    Calculate RMSD of a trajectory against its first frame, excluding water and ions.

    Args:
    traj_file (str): Path to the trajectory file.
    top_file (str): Path to the topology file (e.g., PDB file).
    atom_indices (list, optional): List of atom indices to consider for RMSD calculation.

    Returns:
    pd.DataFrame: DataFrame containing frame indices and corresponding RMSD values.
    """

    dcdAtomSelection: List = []
    for selection in trajectorySelections:
        dcdAtomSelection.extend(drSelector.get_atom_indexes(selection["selection"], pdbFile))


    pdbDf = pdbUtils.pdb2df(pdbFile)
    dcdDf = pdbDf.iloc[dcdAtomSelection]

    subsetPdb = p.join(outDir, "trajectory.pdb")
    pdbUtils.df2pdb(dcdDf, subsetPdb)

    # Load the trajectory with the topology file
    traj = md.load(trajectoryDcd, top=subsetPdb)
    # Select atoms that are not water or ions
    non_water_ions = traj.topology.select('not (resname HOH or resname WAT or resname NA or resname CL)')

    # Use the first frame as the reference
    ref = traj[0]

    # Align the trajectory to the first frame using the selected atoms
    traj.superpose(ref, atom_indices=non_water_ions)

    # Calculate RMSD using the selected atoms
    rmsdValues = md.rmsd(traj, ref, atom_indices=non_water_ions)

    # Create a DataFrame
    rmsdDf = pd.DataFrame({'Frame': range(len(rmsdValues)), 'RMSD': rmsdValues})
    ## clean up temporary pdb file
    return rmsdDf

######################################################################
if __name__ == "__main__":


    simDir = "/home/esp/scriptDevelopment/drMD/04_PET_proj_outputs/6eqe_PET-tetramer_1/021_NVT_warmup"


    for file in os.listdir(simDir):
        if p.splitext(file)[1] in [".pdf", ".png"]:
            os.remove(p.join(simDir, file))


    vitalsCsv = p.join(simDir,"vitals_report.csv")
    progressCsv = p.join(simDir,"progress_report.csv")
    trajectoryDcd = p.join(simDir,"trajectory.dcd")
    pdbFile = p.join(simDir,"trajectory.pdb")

    vitals = {"vitals": vitalsCsv, "progress": progressCsv, "trajectory": trajectoryDcd, "pdb": pdbFile}

    check_vitals(simDir,
                  vitals,
                    [{"selection":{"keyword": "protein"}}, {"selection":{"keyword": "ligand"}}])