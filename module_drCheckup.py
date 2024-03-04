import os
from os import path as p
import pandas as pd
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
from PIL import Image
import numpy as np

######################################################################
def convert_seconds(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02}:{:02}:{:02}".format(int(hours), int(minutes), int(seconds))

######################################################################
def check_vitals(simDir, vitalsCsv, progressCsv):

    basicPng = extract_basic_data(vitalsCsv, progressCsv, simDir)
    energyList = ["Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)"]
    energyPng = plot_vitals(vitalsCsv, simDir, energyList, "Energies", colorOffset=0)

    propertiesList = ["Temperature (K)", "Box Volume (nm^3)", "Density (g/mL)"]
    propertiesPng = plot_vitals(vitalsCsv, simDir, propertiesList, "Properties", colorOffset=3)

    pngFiles = [basicPng, energyPng, propertiesPng]

    create_vitals_pdf(simDir, pngFiles)
######################################################################
def extract_basic_data(vitalsCsv, progressCsv, outDir):
    vitalsDf = pd.read_csv(vitalsCsv)
    progressDf = pd.read_csv(progressCsv)
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

    fig, ax = plt.subplots(figsize=(10, 1))
    ax.axis('off')
    table = ax.table(cellText=df.values, cellLoc = 'center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)

    # save the plot as a png image
    savePng = p.join(outDir,"basic_info.png")
    plt.savefig(savePng,bbox_inches="tight")

    return savePng
######################################################################
def create_vitals_pdf(outDir, pngFiles):
    titleTag = p.basename(outDir)
    savePng = p.join(outDir, "Simulation_Vitals.png")

    basicTable = mpimg.imread(pngFiles[0])
    h, w, _ = basicTable.shape  # get dimensions of image

    vitalsPng = Image.new("RGB", (2480, 3508))
    vitalsPng.paste(Image.fromarray((basicTable * 255).astype(np.uint8)),(100, 100, 100+w, 100+h))
    vitalsPng.save(savePng)
######################################################################
def plot_vitals(vitalsCsv, outDir, yData, tag, colorOffset):
    offRed = to_rgba((204/255, 102/255, 102/255), alpha=1)
    vitalsDf = pd.read_csv(vitalsCsv)
    # Set up the figure and axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    vitalsColors = init_colors()
    # Plot First line on ax1
    lineColor = vitalsColors[0]
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel(yData[0], color=lineColor)
    ax1.plot(vitalsDf['Time (ps)'], vitalsDf[yData[0]], label=yData[0], linestyle='-', color=lineColor, linewidth=2)
    ax1.tick_params(axis='y', labelcolor=lineColor)
    # get labels
    lines, labels = ax1.get_legend_handles_labels()
    for i in range(1,len(yData)):
        data = yData[i]
        lineColor = vitalsColors[i+colorOffset]
    # Create a new  y-axis 
        ax2 = ax1.twinx()
         # Adjusting the position of the new y-axis to avoid overlap
        offset = 80 * (i-1)  # Adjusting this value changes how far away each subsequent axis is
        ax2.spines['right'].set_position(('outward', offset))
        ax2.set_ylabel(data, color=lineColor)
        ax2.plot(vitalsDf['Time (ps)'], vitalsDf[data], label=data, linestyle='-', color=lineColor, linewidth=2)
        ax2.tick_params(axis='y', labelcolor=lineColor)

        # get labels and append
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines += lines2
        labels += labels2
    # Customize the appearance to resemble an ECG readout
    ax1.set_title(f'Simulation {tag} vs Time')
    ax1.grid(True, linestyle='--', alpha=0.7, color=vitalsColors[1],linewidth=0.5, which='both')

    # legend
    ax1.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.35, 1.2))

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
# simDir = "/home/esp/scriptDevelopment/drMD/02_outputs/cvFAP_WT_PLM_FAD_3/03_NpT_pre-equilibriation"
# vitalsCsv = p.join(simDir,"vitals_report.csv")
# progressCsv = p.join(simDir,"progress_report.csv")
# check_vitals(simDir, vitalsCsv, progressCsv)