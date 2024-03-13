import matplotlib.pyplot as plt
import os
from os import path as p
import pandas as pd
import numpy as np
import math
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
            print(file)
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
#############################################################################################
def plot_rmsf(df, outDir):
        # Plot the RMSF
    plt.figure(figsize=(8, 6))
    plt.plot(df["RMSF"], label="RMSF")
    plt.xlabel("Residue Index")
    plt.ylabel("RMSF (Å)")
    plt.title("Root Mean Square Fluctuation (RMSF) vs. Residue Index")
    plt.grid(True)
    plt.legend()
    plt.savefig(p.join(outDir,"RMSF.png"))
#############################################################################################
def plot_rmsd(df, outDir):
    # Plot the RMSD
    plt.figure(figsize=(8, 6))
    plt.plot(df["RMSD"], label="RMSD")
    plt.xlabel("Frame Index")
    plt.ylabel("RMSD (Å)")
    plt.title("Root Mean Square Deviation (RMSD) vs. Frame Index")
    plt.grid(True)
    plt.legend()
    plt.savefig(p.join(outDir,"RMSD.png"))

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