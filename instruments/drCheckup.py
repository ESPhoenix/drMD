## BASIC LIBS
import os
from os import path as p
import pandas as pd
import numpy as np
from functools import wraps
import textwrap

## PLOTTING LIBS
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

## PDF LIBS
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML, CSS

## Molecular Dynamics LIBS
import mdtraj as md

## CLEAN CODE
from typing import Union, Dict, Tuple
from os import PathLike
try:
    from instruments.drCustomClasses import FilePath, DirectoryPath
except:
    from drCustomClasses import FilePath, DirectoryPath


######################################################################
def check_vitals(simDir: Dict, vitalsFiles: Dict[str, FilePath]) -> None:
    """
    Looks at the following properties of the simulation:
        - RMSD
        - Energies (Potential, Kinetic, Total)
        - Properties (Temperature, Volume, Density)
    Checks that each of these properties have converged
    
    
    Also gathers time data for the simulation
    Plots these data and creates a report PDF

    Args:
        simDir (Dict): The directory of the simulation
        vitalsFiles (Dict[str, FilePath]): A dictionary of the vitals files
    """

    ## read OpenMM reporters into dataframes
    vitalsDf = pd.read_csv(vitalsFiles["vitals"])
    progressDf = pd.read_csv(vitalsFiles["progress"])

    ## skip this if the reporters are not present, or are less than 5 entries long
    if len(vitalsDf) < 5:
        return
    
    if len(progressDf) < 5:
        return
    
    ## use mdtraj to calculate RMSD for non water and ions, get that data into a dataframe
    rmsdDf = calculate_rmsd(trajectoryDcd = vitalsFiles["trajectory"], pdbFile = vitalsFiles["pdb"])

    ## plot RMSD as a function of time
    rmsdPng = plot_rmsd(rmsdDf, simDir, ["RMSD"], "RMSD")

    ## get time data, plot a table
    timeDf = extract_time_data(vitalsDf, progressDf)
    timePng = plot_time_data(timeDf, simDir)
    
    ## plot energy as a function of time 
    energyList = ["Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)"]
    energyPlot = plot_vitals(vitalsDf, simDir, energyList, "Energies")
    ## plot properties as a function of time 
    propertiesList = ["Temperature (K)", "Box Volume (nm^3)", "Density (g/mL)"]
    propertiesPlot = plot_vitals(vitalsDf, simDir, propertiesList, "Properties")

    ## combine all feature names
    allFeaturesList = energyList + propertiesList + ["RMSD"]
    ## combine all features into one dataframe
    allFeaturesDf = pd.concat([vitalsDf,rmsdDf],axis=1)
    ## check convergance for all features
    converganceDf  = check_convergance(allFeaturesDf,allFeaturesList,simDir,"Convergance Checks")
    ## plot convergance booleans as a table
    convergancePng = plot_converged(simDir, converganceDf)

    ## create report PDF using all the plots created above
    create_vitals_pdf(simDir)
######################################################################
def create_vitals_pdf(simDir):
    instrumentsDir = p.dirname(__file__)
    env = Environment(loader=FileSystemLoader(instrumentsDir))
    template = env.get_template("vitals_template.html")
    
    # Render the template with any context variables you need
    context = {
        # Add your context variables here
    }
    rendered_html = template.render(context)
    
    # Generate the PDF with error handling
    outPdf = p.join(simDir, "vitals_report.pdf")
    try:
        base_url = simDir # Set the base URL to the current working directory
        HTML(string=rendered_html, base_url=base_url).write_pdf(outPdf, stylesheets=[CSS(string='@page { margin: 0; }')])
    except Exception as e:
        print(f"Error generating PDF: {e}")
######################################################################
def check_up_handler():
    def decorator(simulationFunction):
        @wraps(simulationFunction)
        def wrapper(*args, **kwargs):
            saveFile: Union[PathLike, str] = simulationFunction(*args, **kwargs)

            vitalsFiles, simDir = find_vitals_files(kwargs["sim"], kwargs["outDir"])

            check_vitals(simDir = simDir,
                          vitalsFiles = vitalsFiles)

            return saveFile
        return wrapper
    return decorator
######################################################################
def find_vitals_files(simInfo: Dict, outDir: Union[PathLike, str]):
    simDir: Union[PathLike, str] = p.join(outDir, simInfo["stepName"])

    vitalsReport: Union[PathLike, str] = p.join(simDir, "vitals_report.csv")
    if not p.isfile(vitalsReport):
        raise FileNotFoundError(f"->\tReporter file not found at {vitalsReport}")

    progressReport: Union[PathLike, str] = p.join(simDir, "progress_report.csv")
    if not p.isfile(progressReport):
        raise FileNotFoundError(f"->\tReporter file not found at {progressReport}")


    trajectoryDcd: Union[PathLike, str] = p.join(simDir, "trajectory.dcd")
    if not p.isfile(trajectoryDcd):
        raise FileNotFoundError(f"->\Trajectory file not found at {trajectoryDcd}")
    

    trajectoryDcd: Union[PathLike, str] = p.join(simDir, "trajectory.dcd")
    if not p.isfile(trajectoryDcd):
        raise FileNotFoundError(f"->\Trajectory file not found at {trajectoryDcd}")
    
    pdbFile = False
    for file in os.listdir(simDir):
        if file.endswith(".pdb"):
            pdbFile = p.join(simDir, file)

    if not pdbFile:
        raise FileNotFoundError(f"->\tPDB file not found at {simDir}")
    

    vitalsFiles = {"vitals": vitalsReport, "progress": progressReport, "trajectory": trajectoryDcd, "pdb": pdbFile}
    
    return vitalsFiles, simDir


######################################################################
def check_convergance(df, columns, simDir, tag, windowSize = 5):
    convergedDict = {}
    for column in columns:
        if len(column) == 1:
            continue
        dataRange = df[column].max() - df[column].min()
        converganceTolerance = dataRange * 0.05 
        runningAverage = df[column].rolling(window=windowSize).mean()
        lastTwoAverages = runningAverage.tail(2).values
        isConverged = np.abs(lastTwoAverages[0] - lastTwoAverages[1]) < converganceTolerance
        entryLabel = column.split("(")[0].split()[0]
        convergedDict.update({entryLabel:isConverged})

    convergedData = [(key, value) for key, value in convergedDict.items()]
    convergedDf = pd.DataFrame(convergedData, columns=["Property", "Converged"])    

    return convergedDf
######################################################################
def plot_converged(simDir, convergedDf):
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
def convert_seconds(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02}:{:02}:{:02}".format(int(hours), int(minutes), int(seconds))

######################################################################
def extract_time_data(vitalsDf, progressDf):
    # get elapsed time
    timeElapsed = progressDf["Elapsed Time (s)"].tail(1).values[0]
    timeElapsed = convert_seconds(timeElapsed)
    # get average speed
    averageSpeed = str(round(progressDf["Speed (ns/day)"].mean(),2))
    # get nSteps
    nSteps = str(vitalsDf["#\"Step\""].tail(1).values[0])
    # get simulation duration
    duration = str(round(vitalsDf["Time (ps)"].tail(1).values[0],2))

    # create a DataFrame
    df = pd.DataFrame({'Elapsed Time (HH:MM:SS)': [timeElapsed],
        'Average Speed (ns/day)': [averageSpeed],
        'Total Steps': [nSteps],
        'Simulation Duration (ps)': [duration]})

    df = df.transpose()
    df.columns = ['Value']  # rename the column
    df.reset_index(level=0, inplace=True)  # reset the index

    return df

######################################################################
def plot_time_data(timeDf, outDir):
    # Convert cm to inches
    width_inch = 8.97 / 2.54
    height_inch = 6.09 / 2.54
    
    fig, ax = plt.subplots(figsize=(width_inch, height_inch))
    fig.patch.set_facecolor('#1a1a1a')  # Much darker grey
    brightGreen = '#00FF00'  # Brighter green
    ax.axis('off')
    
    # Wrap text in each cell to ensure it fits within the cell
    wrapped_text = timeDf.map(lambda x: '\n'.join(textwrap.wrap(str(x), width=15)))

    # Create the table without column labels
    table = ax.table(cellText=wrapped_text.values, cellLoc='center', loc='center', colLabels=None)
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

    # Save the plot as a PNG image with reduced border
    savePng = p.join(outDir, "time_info.png")
    plt.savefig(savePng, bbox_inches="tight", pad_inches=0.05, facecolor='#1a1a1a')
    plt.close()
    return savePng


######################################################################

def plot_vitals(vitalsDf, outDir, yData, tag):
    # Convert cm to inches
    width_inch = 29.7 / 2.54
    height_inch = 10 / 2.54
    
    # Set up the figure and axis
    fig, axes = plt.subplots(nrows=1, ncols=len(yData), figsize=(width_inch, height_inch), constrained_layout=True)
    fig.patch.set_facecolor('#1a1a1a')  # Much darker grey
    brightGreen = '#00FF00'  # Brighter green
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

def plot_rmsd(vitalsDf, outDir, yData, tag):
    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(3.54, 3.54))
    fig.patch.set_facecolor('#1a1a1a')  # Much darker grey
    brightGreen = '#00FF00'  # Brighter green
    brightRed = '#FF0000'  # Brighter red
    fig.suptitle(f'Simulation {tag} vs Time', fontsize=12, y=0.98, color=brightGreen)
    
    lineColor = brightRed
    ax.set_facecolor('#1a1a1a')  # Much darker grey
    ax.set_xlabel('Frame', fontsize=10, color=brightGreen)
    ax.set_ylabel(yData[0], color=brightGreen, fontsize=10, labelpad=15)
    ax.plot(vitalsDf['Frame'], vitalsDf[yData[0]], label=yData[0],
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
    
    # Dynamically place the legend to avoid overlapping with the trace
    legend = ax.legend(loc='best', fontsize=8, facecolor='#1a1a1a', edgecolor=brightGreen)
    for text in legend.get_texts():
        text.set_color(brightGreen)
    
    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    # save
    savePng = p.join(outDir, f"Vitals_{tag}.png")
    plt.savefig(savePng, bbox_inches="tight", facecolor='#1a1a1a')
    fig.clf()
    return savePng
######################################################################
def calculate_rmsd(trajectoryDcd, pdbFile):
    """
    Calculate RMSD of a trajectory against its first frame, excluding water and ions.

    Args:
    traj_file (str): Path to the trajectory file.
    top_file (str): Path to the topology file (e.g., PDB file).
    atom_indices (list, optional): List of atom indices to consider for RMSD calculation.

    Returns:
    pd.DataFrame: DataFrame containing frame indices and corresponding RMSD values.
    """
    # Load the trajectory with the topology file
    traj = md.load(trajectoryDcd, top=pdbFile)

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

    return rmsdDf

######################################################################
if __name__ == "__main__":


    simDir = "/home/esp/scriptDevelopment/drMD/03_outputs/A/02_NVT_pre-equilibraition"
    for file in os.listdir(simDir):
        if p.splitext(file)[1] in [".pdf", ".png"]:
            os.remove(p.join(simDir, file))



    vitalsCsv = p.join(simDir,"vitals_report.csv")
    progressCsv = p.join(simDir,"progress_report.csv")
    trajectoryDcd = p.join(simDir,"trajectory.dcd")
    pdbFile = p.join(simDir,"A.pdb")

    vitals = {"vitals": vitalsCsv, "progress": progressCsv, "trajectory": trajectoryDcd, "pdb": pdbFile}
    check_vitals(simDir, vitals)