import matplotlib.pyplot as plt
import os
from os import path as p
import pandas as pd
import numpy as np
import math
import yaml
from adjustText  import adjust_text
#########################################################################################################
def plot_delta_RMSF(analDir, referenceSystem):
    deltaDf = pd.read_csv(p.join(analDir, "delta_RMSF.csv"),index_col="RES_ID")
    peaksYaml = p.join(analDir,"delta_RMSF_peaks.yaml")
    with open(peaksYaml,"r") as file:
        deltaPeaks = yaml.safe_load(file)

    # Create a figure and axes
    nSubplots = len(deltaDf.columns)
    fig, axs = plt.subplots(nSubplots, figsize=(8, 6*nSubplots))

    # Plot each column on a different axis and annotate the peaks
    for ax, col in zip(axs, deltaDf.columns):
        texts = []  # to store all the text objects
        ax.plot(deltaDf[col], label=col)
        if col in deltaPeaks:
            for peak in deltaPeaks[col]:
                if deltaDf[col].iloc[peak] > 0:
                    va = 'bottom'
                else:
                    va = 'top'
                texts.append(ax.text(peak, deltaDf[col].iloc[peak], str(peak), color='red', fontsize=8, ha='left', va=va))
                adjust_text(texts, ax = ax, x = deltaDf.index, y = deltaDf[col])
        ax.set_xlabel('Residue Number')
        ax.set_ylabel('RMSF Values')
        if referenceSystem:
            ax.set_title(f'Delta RMSF of {col} compared to {referenceSystem}')
        else:
            ax.set_title(f'Delta RMSF of {col}')
        ax.legend()

    rmsfPng = p.join(analDir,"delta_RMSF.png")
    plt.savefig(rmsfPng,bbox_inches="tight")

#########################################################################################################
def plot_distance_hist(sysAnalDir, resTag):
    print(sysAnalDir)
    print(resTag)
    ## load contact dfs from sysAnalDir | concat into one big df
    dfsToConcat = []
    for file in os.listdir(sysAnalDir):
        if not p.splitext(file)[1] == ".csv":
            continue
        if  file.startswith("contacts") and file.split("_")[1] == resTag:
            runDf = pd.read_csv(p.join(sysAnalDir,file))
    
            dfsToConcat.append(runDf)
    if len(dfsToConcat) == 0:
        return
    df = pd.concat(dfsToConcat, axis = 0, ignore_index=True)
    df.drop(columns = ["Unnamed: 0"], inplace = True)
    # Generate a list of unique colors
    colors = plt.cm.viridis(np.linspace(0, 1, len(df.columns)))

    # Calculate the number of rows needed for the subplots
    num_rows = math.ceil(len(df.columns) / 4)

    # Create a subplot for each column in a 4xN grid
    fig, axs = plt.subplots(num_rows, 4, figsize=(20, 5*num_rows))

    # Increase the font size of the suptitle and adjust its y position
    plt.suptitle(f"Interaction distances between {resTag} and nearby residues (in Ang)", fontsize=32, y=1.02)

    for i, column in enumerate(df.columns):
        row = i // 4
        col = i % 4
        # Convert distances from nm to Angstrom
        distances_in_angstrom = df[column] * 10
        axs[row, col].hist(distances_in_angstrom, bins=30, alpha=0.5, label=column, color=colors[i])
        axs[row, col].legend()
        # Set x-axis limits
        axs[row, col].set_xlim([2, 10])
        axs[row, col].set_ylim([0,100])

    # Remove empty subplots
    [fig.delaxes(ax) for ax in axs.flatten() if not ax.has_data()]

    plt.tight_layout()
    plt.savefig(p.join(sysAnalDir, f"distances_{resTag}.png"), bbox_inches="tight")
    plt.close()
#############################################################################################
def plot_RMSF(sysAnalDir):
    print(sysAnalDir)
    ## load contact dfs from sysAnalDir | concat into one big df
    dfsToConcat = []
    for file in os.listdir(sysAnalDir):
        if not p.splitext(file)[1] == ".csv":
            continue
        if  file.startswith("RMSF"):
            runDf = pd.read_csv(p.join(sysAnalDir,file),index_col="RES_ID")
            dfsToConcat.append(runDf)
    if len(dfsToConcat) == 0:
        return
    df = pd.concat(dfsToConcat, axis = 1, ignore_index=False) 
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot each column with a different color
    for col in df:
        plt.plot(df[col], label=col)

    # Calculate the mean for each row
    mean_values = df.mean(axis=1)

    # Plot the mean in bold black
    plt.plot(mean_values, color='black', linewidth=2, label='Mean')

    # Customize the plot
    plt.xlabel('Residue Number')
    plt.ylabel('RMSF Values')
    plt.title('RMSF Per Residue')
    plt.legend()
    rmsfPng = p.join(sysAnalDir,"RMSF.png")
    plt.savefig(rmsfPng,bbox_inches="tight")

#############################################################################################
def plot_RMSD(sysAnalDir):
    ## load contact dfs from sysAnalDir | concat into one big df
    dfsToConcat = []
    for file in os.listdir(sysAnalDir):
        if not p.splitext(file)[1] == ".csv":
            continue
        if  file.startswith("RMSD"):
            runDf = pd.read_csv(p.join(sysAnalDir,file),index_col="Unnamed: 0")
            dfsToConcat.append(runDf)
    if len(dfsToConcat) == 0:
        return
    df = pd.concat(dfsToConcat, axis = 1, ignore_index=False) 
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot each column with a different color
    for col in df:
        plt.plot(df[col], label=col)

    # Calculate the mean for each row
    mean_values = df.mean(axis=1)

    # Plot the mean in bold black
    plt.plot(mean_values, color='black', linewidth=2, label='Mean')

    # Customize the plot
    plt.xlabel('Data Points')
    plt.ylabel('RMSD Values')
    plt.title('RMSD Values Over Time')
    plt.legend()
    rmsdPng = p.join(sysAnalDir,"RMSD.png")
    plt.savefig(rmsdPng,bbox_inches="tight")

#############################################################################################
def plot_rdf(rdfDf, outDir, resTag):
    # Plot the RDF
    plt.figure(figsize=(10, 6))
    plt.plot(rdfDf['Radii'], rdfDf['RDF'], label='RDF')
    plt.xlabel('Radii')
    plt.ylabel('Radial Distribution Function')
    plt.title('Radial Distribution Function vs Radii')
    plt.legend()
    plt.grid(True)
    plt.savefig(p.join(outDir,f"RDF_{resTag}.png"))
    plt.close()