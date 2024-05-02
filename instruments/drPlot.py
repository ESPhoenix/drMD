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
def histogram_plotting_manager(sysAnalDir):
    ## read through analysis output files | gather info from filenames | package into a dataframe
    dfsToConcat = []
    for file in os.listdir(sysAnalDir):
        fileData = p.splitext(file)
        if not fileData[1] == ".csv":
            continue
        fileData = fileData[0].split("_")
        dataTag = fileData[0]                   ## get type of data
        resAtomTag = fileData[1]   ## identifier of residue or atom of interest
        protName = "_".join(fileData[2:])                 ## protein name (contains repeats information)
        csvFile = p.join(sysAnalDir,file)

        if not dataTag in ["contacts", "donor", "acceptor"]:
            continue

        row = pd.DataFrame([[dataTag, resAtomTag, protName, csvFile]], columns=["dataTag", "resAtomTag", "protName", "csvPath"])
        dfsToConcat.append(row)
    
    managerDf = pd.concat(dfsToConcat, axis=0, ignore_index=True)

    uniqueDataTags = managerDf["dataTag"].unique().tolist()
    for dataTag in uniqueDataTags:
        groupedDf = managerDf[managerDf["dataTag"] == dataTag]

        plot_distance_hist(groupedDf, dataTag, sysAnalDir)



#########################################################################################################
def plot_distance_hist(groupedDf, dataTag, sysAnalDir):
    print(dataTag)
    ## load csv files and extract columnNames
    csvFiles = groupedDf["csvPath"].to_list()
    dfs = []
    columnNames = []
    for csvFile in csvFiles:
        df = pd.read_csv(csvFile, index_col="Unnamed: 0")
        dfs.append(df)
        columnNames += df.columns.to_list()
    columnNames = list(set(columnNames))

    # Generate a list of unique colors
    colors = plt.cm.viridis(np.linspace(0, 1, len(columnNames)))

    # Calculate the number of rows needed for the subplots
    num_rows = math.ceil(len(columnNames) / 4)

    # Create a subplot for each column in a 4xN grid
    fig, axs = plt.subplots(num_rows, 4, figsize=(20, 5*num_rows))
    axs = axs.flatten()

    # Increase the font size of the suptitle and adjust its y position
    i = 0
    for colName in columnNames:
        for df in dfs:
            if not colName in df.columns:
                continue
            # Convert distances from nm to Angstrom
            distances_in_angstrom = df[colName] * 10
            print(distances_in_angstrom)
            plotColor = colors[i]
            axs[i].hist(distances_in_angstrom, bins=30, alpha=0.5, label=colName, color=plotColor)
            axs[i].legend()
            # Set x-axis limits
            axs[i].set_xlim([2, 10])
            axs[i].set_ylim([0,100])
        axs[i].set_title(colName)

        i += 1
        # plt.suptitle(f"{colName}", fontsize=32, y=1.02)

    # Remove empty subplots
    [fig.delaxes(ax) for ax in axs.flatten() if not ax.has_data()]

    plt.tight_layout()
    plt.savefig(p.join(sysAnalDir, f"{dataTag}.png"), bbox_inches="tight")
    plt.close()
#############################################################################################
def plot_RMSF(sysAnalDir):
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