import os
from os import path as p
import pandas as pd
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

def main():
    vitalsCsv = "/home/eugene/drMD/02_outputs/cvFAP_WT_PLM_FAD_3/02_NVT_pre-equilibraition/vitals_report.csv"
    outDir = "/home/eugene/drMD/02_outputs/cvFAP_WT_PLM_FAD_3/02_NVT_pre-equilibraition"
    plot_vitals(vitalsCsv, outDir)

def plot_vitals(vitalsCsv, outDir):
    offRed = to_rgba((204/255, 102/255, 102/255), alpha=1)

    vitalsDf = pd.read_csv(vitalsCsv)
    # Set up the figure and axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot Potential Energy on the first y-axis
    color = 'tab:blue'
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Potential Energy (kJ/mole)', color=color)
    ax1.plot(vitalsDf['Time (ps)'], vitalsDf['Potential Energy (kJ/mole)'], label='Potential Energy', linestyle='-', color=color, linewidth=2)
    ax1.tick_params(axis='y', labelcolor=color)

    # Create a second y-axis for Kinetic Energy
    ax2 = ax1.twinx()
    color = 'tab:green'
    ax2.set_ylabel('Kinetic Energy (kJ/mole)', color=color)
    ax2.plot(vitalsDf['Time (ps)'], vitalsDf['Kinetic Energy (kJ/mole)'], label='Kinetic Energy', linestyle='-', color=color, linewidth=2)
    ax2.tick_params(axis='y', labelcolor=color)

    # Create a third y-axis for Total Energy
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))
    color = 'tab:red'
    ax3.set_ylabel('Total Energy (kJ/mole)', color=color)
    ax3.plot(vitalsDf['Time (ps)'], vitalsDf['Total Energy (kJ/mole)'], label='Total Energy', linestyle='-', color=color, linewidth=2)
    ax3.tick_params(axis='y', labelcolor=color)

    # Customize the appearance to resemble an ECG readout
    ax1.set_title('Energy vs Time - ECG Readout Style')
    ax1.grid(True, linestyle='--', alpha=0.7, color=offRed,linewidth=0.5, which='both')

    # Combine legends for all lines
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines3, labels3 = ax3.get_legend_handles_labels()
    lines += lines2 + lines3
    labels += labels2 + labels3
    ax1.legend(lines, labels, loc='upper left')
    # save
    plt.savefig(p.join(outDir, "vitals_energy_ecg.png"), bbox_inches="tight")
if __name__ == "__main__":
    main()