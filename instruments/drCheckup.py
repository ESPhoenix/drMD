from os import path as p
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from functools import wraps

from matplotlib.colors import to_rgba
from PIL import Image
import numpy as np
from fpdf import FPDF


from typing import Union, Dict, Tuple
from os import PathLike
######################################################################
def check_up_handler():
    def decorator(simulationFunction):
        @wraps(simulationFunction)
        def wrapper(*args, **kwargs):
            saveFile: Union[PathLike, str] = simulationFunction(*args, **kwargs)

            vitalsReport, progressReport, simDir = find_reporters(kwargs["sim"], kwargs["outDir"])

            check_vitals(simDir = simDir, vitalsCsv = vitalsReport, progressCsv = progressReport)

            return saveFile
        return wrapper
    return decorator
######################################################################
def find_reporters(simInfo: Dict, outDir: Union[PathLike, str]) -> Tuple[Union[PathLike, str], Union[PathLike, str]]:
    simDir: Union[PathLike, str] = p.join(outDir, simInfo["stepName"])

    vitalsReport: Union[PathLike, str] = p.join(simDir, "vitals_report.csv")
    if not p.isfile(vitalsReport):
        raise FileNotFoundError(f"->\tReporter file not found at {vitalsReport}")

    progressReport: Union[PathLike, str] = p.join(simDir, "progress_report.csv")
    if not p.isfile(progressReport):
        raise FileNotFoundError(f"->\tReporter file not found at {progressReport}")
    
    return vitalsReport, progressReport, simDir


######################################################################
def create_vitals_pdf(simDir, basicPng, energyConvTable, energyPlot, propertiesConvTable, propertiesPlot):
    """
    Creates a PDF document containing simulation vitals using the given images and plots.
    
    Args:
        simDir (str): The directory where the PDF will be saved.
        basicPng (str): The path to the basic information image file.
        energyConvTable (str): The path to the energy convergence table image file.
        energyPlot (str): The path to the energy plot image file.
        propertiesConvTable (str): The path to the properties convergence table image file.
        propertiesPlot (str): The path to the properties plot image file.
        
    Returns:
        None
        
    This function creates a PDF document using the FPDF library with the given images and plots. The PDF document contains the following sections:
    
    - Title: "Simulation Vitals"
    - Basic Information: The basic information image file is added to the PDF with a specified width and height.
    - Energy Convergence Table: The energy convergence table image file is added to the PDF with a specified width and height.
    - Properties Convergence Table: The properties convergence table image file is added to the PDF with a specified width and height.
    - Energy Plot: The energy plot image file is added to the PDF with a specified width and height.
    - Properties Plot: The properties plot image file is added to the PDF with a specified width and height.
    
    The PDF document is saved in the specified simulation directory with the filename "Simulation_Vitals.pdf".
    """
    pdf = FPDF()
    # Add a page
    pdf.add_page()
    # Set font
    pdf.set_font("Arial", size = 36)
    ##### TITLE ####
    pdf.cell(200, 10, txt = "Simulation Vitals", ln = True, align = 'L')
    ##### BASIC INFO #####
    # Open the image file
    with Image.open(basicPng) as bI:
        bI_width, bI_height = bI.size
    # Calculate the desired width and height in mm
    desired_width = pdf.w * 0.6 # 90% of the page width
    desired_height = pdf.h  # 2nd to 5th of the page height
    # Calculate scale factors
    scale_width = desired_width / bI_width
    scale_height = desired_height / bI_height
    # Use the smaller scale factor to ensure the image fits both dimensions
    scale = min(scale_width, scale_height)
    # Calculate the scaled dimensions
    bI_width *= scale
    bI_height *= scale
    # Add the image to the pdf
    pdf.image(basicPng, 30, 25, bI_width, bI_height)

    nowHeight = 15 + bI_height
    # ##### ENERGY CONV #####
    with Image.open(energyConvTable) as eC:
        eC_width, eC_height = eC.size
    desired_width  = pdf.w * 0.3
    desired_height  = pdf.h
    scale_width = desired_width / eC_width
    scale_height = desired_height / eC_height
    scale = min(scale_width, scale_height)
    # Calculate the scaled dimensions
    eC_width *= scale
    eC_height *= scale
    pdf.image(energyConvTable, 30 , nowHeight, eC_width, eC_height)

    ##### PROPERTIES CONV #####
    with Image.open(propertiesConvTable) as pC:
        pC_width, pC_height = pC.size
    desired_width  = pdf.w * 0.3
    desired_height  = pdf.h
    scale_width = desired_width / pC_width
    scale_height = desired_height / pC_height
    scale = min(scale_width, scale_height)
    # Calculate the scaled dimensions
    pC_width *= scale
    pC_height *= scale
    pdf.image(propertiesConvTable, 35 + eC_width, nowHeight, pC_width, pC_height)

    nowHeight = nowHeight + 2 + eC_height
   ##### ENERGY PLOT #####
    with Image.open(energyPlot) as eP:
        eP_width, eP_height = eP.size
    desired_width  = pdf.w * 0.75
    desired_height  = pdf.h
    scale_width = desired_width / eP_width
    scale_height = desired_height / eP_height
    scale = min(scale_width, scale_height)
    # Calculate the scaled dimensions
    eP_width *= scale
    eP_height *= scale
    pdf.image(energyPlot, 15, nowHeight, eP_width, eP_height)

    nowHeight = nowHeight + 2 + eP_height
    ##### PROPERTIES PLOT #####
    with Image.open(propertiesPlot) as pP:
        pP_width, pP_height = pP.size
    desired_width  = pdf.w * 0.75
    desired_height  = pdf.h
    scale_width = desired_width / pP_width
    scale_height = desired_height / pP_height
    scale = min(scale_width, scale_height)
    # Calculate the scaled dimensions
    pP_width *= scale
    pP_height *= scale
    pdf.image(propertiesPlot, 15, nowHeight, pP_width, pP_height)


    ##### SAVE PDF #####
    pdf.output(p.join(simDir,"Simulation_Vitals.pdf"), "F")



######################################################################
def check_convergance(df, columns, simDir, tag, windowSize = 5,):
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
    convergedDf = pd.DataFrame(convergedData, columns=["Property", "Converged?"])    

    fig, ax = plt.subplots(figsize=(2, 2))
    ax.axis('off')
    table = ax.table(cellText=convergedDf.values, colLabels=convergedDf.columns, loc='center',
                     cellLoc='center', rowLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.5, 2)

    for key, cell in table.get_celld().items():
        cell.set_linewidth(0.5)

    for i in range(len(convergedDf.columns)):
        table[(0, i)].set_facecolor("#40466e")
        table[(0, i)].get_text().set_color('white')

    for key, cell in table.get_celld().items():
        cell.set_linewidth(0.5)


    # Save the plot as a PNG file
    savePng = p.join(simDir, f"{tag}.png")
    plt.savefig(savePng, bbox_inches='tight')
    plt.close()
    return savePng

######################################################################
def convert_seconds(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02}:{:02}:{:02}".format(int(hours), int(minutes), int(seconds))
######################################################################
def check_vitals(simDir, vitalsCsv, progressCsv):

    vitalsDf = pd.read_csv(vitalsCsv)
    progressDf = pd.read_csv(progressCsv)

    basicPng = extract_basic_data(vitalsDf, progressDf, simDir)
    energyList = ["Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)"]
    energyConvTable = check_convergance(vitalsDf,energyList,simDir,"energy_convergance")
    energyPlot = plot_vitals(vitalsDf, simDir, energyList, "Energies", colorOffset=0)

    propertiesList = ["Temperature (K)", "Box Volume (nm^3)", "Density (g/mL)"]
    propertiesConvTable = check_convergance(vitalsDf, propertiesList, simDir, "properties_convergance")
    propertiesPlot = plot_vitals(vitalsDf, simDir, propertiesList, "Properties", colorOffset=3)


    create_vitals_pdf(simDir, basicPng, energyConvTable, energyPlot, propertiesConvTable, propertiesPlot)
######################################################################
def extract_basic_data(vitalsDf, progressDf, outDir):
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
    df = pd.DataFrame({'Elapsed Time': [timeElapsed],
        'Average Speed': [averageSpeed],
        'Total Steps': [nSteps],
        'Simulation Duration (ps)': [duration]})

    df = df.transpose()
    df.columns = ['Value']  # rename the column
    df.reset_index(level=0, inplace=True)  # reset the index

    fig, ax = plt.subplots(figsize=(5, 2))
    ax.axis('off')
    table = ax.table(cellText=df.values, cellLoc = 'center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.2)

    # save the plot as a png image
    savePng = p.join(outDir,"basic_info.png")
    plt.savefig(savePng,bbox_inches="tight")
    plt.close()
    return savePng


######################################################################
def plot_vitals(vitalsDf, outDir, yData, tag, colorOffset):
    # Set up the figure and axis
    fig, ax1 = plt.subplots(figsize=(12, 6))
    vitalsColors = init_colors()
    # Plot First line on ax1
    lineColor = vitalsColors[0]
    ax1.set_xlabel('Time (ps)', fontsize=18)
    ax1.set_ylabel(yData[0], color=lineColor, fontsize=18)
    ax1.plot(vitalsDf['Time (ps)'], vitalsDf[yData[0]], label=yData[0], linestyle='-', color=lineColor, linewidth=2)
    ax1.tick_params(axis='y', labelcolor=lineColor)
    for label in ax1.get_xticklabels():
        label.set_fontsize(18)
    for label in ax1.get_yticklabels():
        label.set_fontsize(18)
    # get labels
    lines, labels = ax1.get_legend_handles_labels()
    for i in range(1,len(yData)):
        data = yData[i]
        lineColor = vitalsColors[i+colorOffset]
    # Create a new  y-axis 
        ax2 = ax1.twinx()
         # Adjusting the position of the new y-axis to avoid overlap
        offset = 120 * (i-1)  # Adjusting this value changes how far away each subsequent axis is
        ax2.spines['right'].set_position(('outward', offset))
        ax2.set_ylabel(data, color=lineColor, fontsize=18)
        ax2.plot(vitalsDf['Time (ps)'], vitalsDf[data], label=data, linestyle='-', color=lineColor, linewidth=2)
        ax2.tick_params(axis='y', labelcolor=lineColor)
        for label in ax2.get_xticklabels():
            label.set_fontsize(18)
        for label in ax2.get_yticklabels():
            label.set_fontsize(18)
        # get labels and append
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines += lines2
        labels += labels2
    # Customize the appearance to resemble an ECG readout
    ax1.set_title(f'Simulation {tag} vs Time', loc = "right", fontsize=25)
    ax1.grid(True, linestyle='--', alpha=0.7, color=vitalsColors[1],linewidth=0.5, which='both')
    # legend
    ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.35, 1.35), fontsize=18)
    # save
    savePng = p.join(outDir, f"Vitals_{tag}.png")
    plt.savefig(savePng, bbox_inches="tight")
    fig.clf()
    return savePng
######################################################################
def init_colors():
    black = to_rgba((0, 0, 0, 1.0))  # Black
    darkRed = to_rgba((0.545, 0, 0, 1.0))  # Dark Red
    brightRed = to_rgba((1.0, 0, 0, 1.0))  # Bright Red
    darkBlue = to_rgba((0, 0, 0.545, 1.0))  # Dark Blue
    brightBlue = to_rgba((0, 0, 1.0, 1.0))  # Bright Blue
    darkGreen = to_rgba((0, 0.5, 0, 1.0))  # Dark Green
    brightGreen = to_rgba((0, 1.0, 0, 1.0))  # Bright Green
    purple = to_rgba((0.5, 0, 0.5, 1.0))  # Purple
    orange = to_rgba((1.0, 0.5, 0, 1.0))  # Orange

    vitalsColors = [black, darkRed, brightRed, darkBlue, brightBlue, darkGreen, brightGreen, purple, orange]
    return vitalsColors
######################################################################
if __name__ == "__main__":
    simDir = "/home/esp/scriptDevelopment/drMD/02_outputs/cvFAP_WT_PLM_FAD_3/05_NpT_equilibriation"
    vitalsCsv = p.join(simDir,"vitals_report.csv")
    progressCsv = p.join(simDir,"progress_report.csv")
    check_vitals(simDir, vitalsCsv, progressCsv)