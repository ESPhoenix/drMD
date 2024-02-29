import os
from os import path as p
import pandas as pd
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
##########################################################################################################################
def main():
    vitalsCsv = "/home/eugene/drMD/02_outputs/cvFAP_WT_PLM_FAD_3/02_NVT_pre-equilibraition/vitals_report.csv"
    outDir = "/home/eugene/drMD/02_outputs/cvFAP_WT_PLM_FAD_3/02_NVT_pre-equilibraition"

    energyList = ["Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)"]
    plot_vitals(vitalsCsv, outDir, energyList, "Energies")
    propertiesList = ["Temperature (K)", "Box Volume (nm^3)", "Density (g/mL)"]
    plot_vitals(vitalsCsv, outDir, propertiesList, "Properties")
##########################################################################################################################
def plot_vitals(vitalsCsv, outDir, yData, tag):
    offRed = to_rgba((204/255, 102/255, 102/255), alpha=1)

    vitalsDf = pd.read_csv(vitalsCsv)
    # Set up the figure and axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    vitalsColors = init_colors()

    # Plot First line on ax1
    color = vitalsColors[0]
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel(yData[0], color=color)
    ax1.plot(vitalsDf['Time (ps)'], vitalsDf[yData[0]], label=yData[0], linestyle='-', color=color, linewidth=2)
    ax1.tick_params(axis='y', labelcolor=color)

    # get labels
    lines, labels = ax1.get_legend_handles_labels()

    for i in range(1,len(yData)):
        data = yData[i]
        lineColor = vitalsColors[i]
    # Create a new  y-axis 
        ax2 = ax1.twinx()
        color = lineColor

         # Adjusting the position of the new y-axis to avoid overlap
        offset = 80 * (i-1)  # Adjusting this value changes how far away each subsequent axis is
        ax2.spines['right'].set_position(('outward', offset))
        ax2.set_ylabel(data, color=color)
        ax2.plot(vitalsDf['Time (ps)'], vitalsDf[data], label=data, linestyle='-', color=color, linewidth=2)
        ax2.tick_params(axis='y', labelcolor=color)

        # get labels and append
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines += lines2
        labels += labels2

    # Customize the appearance to resemble an ECG readout
    ax1.set_title(f'Simulation {tag} vs Time')
    ax1.grid(True, linestyle='--', alpha=0.7, color=offRed,linewidth=0.5, which='both')

    # legend
    ax1.legend(lines, labels, loc='upper left')
    # save
    plt.savefig(p.join(outDir, f"Vitals_{tag}.png"), bbox_inches="tight")
    fig.clf()
##########################################################################################################################
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


##########################################################################################################################
if __name__ == "__main__":
    main()