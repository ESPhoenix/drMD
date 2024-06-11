import matplotlib.pyplot as plt
import os
from os import path as p
import pandas as pd
import numpy as np
import math
import yaml
from adjustText  import adjust_text
#########################################################################################################
def plot_delta_RMSF(analDir: str, referenceSystem: str) -> None:
    """
    Plots delta RMSF values for each column in the delta_RMSF.csv file.

    Args:
        analDir (str): Directory containing delta_RMSF.csv and delta_RMSF_peaks.yaml
        referenceSystem (str): Name of the reference system to compare against
    """
    # Read delta_RMSF.csv into a dataframe
    deltaDf: pd.DataFrame = pd.read_csv(p.join(analDir, "delta_RMSF.csv"),index_col="RES_ID")
    # Read delta_RMSF_peaks.yaml and extract peaks data
    peaksYaml: str = p.join(analDir,"delta_RMSF_peaks.yaml")
    with open(peaksYaml,"r") as file:
        deltaPeaks: dict = yaml.safe_load(file)

    # Create a figure and axes for each column
    nSubplots: int = len(deltaDf.columns)
    fig, axs = plt.subplots(nSubplots, figsize=(8, 6*nSubplots))

    # Plot each column on a different axis and annotate peaks
    for ax, col in zip(axs, deltaDf.columns):
        texts = []  # to store all the text objects
        ax.plot(deltaDf[col], label=col)
        if col in deltaPeaks:
            for peak in deltaPeaks[col]:
                if deltaDf[col].iloc[peak] > 0:
                    va: str = 'bottom'
                else:
                    va: str = 'top'
                texts.append(ax.text(peak, deltaDf[col].iloc[peak], str(peak), color='red', fontsize=8, ha='left', va=va))
                adjust_text(texts, ax = ax, x = deltaDf.index, y = deltaDf[col])
        ax.set_xlabel('Residue Number')
        ax.set_ylabel('RMSF Values')
        # Set title with reference system if specified
        if referenceSystem:
            ax.set_title(f'Delta RMSF of {col} compared to {referenceSystem}')
        else:
            ax.set_title(f'Delta RMSF of {col}')
        ax.legend()

    # Save the figure as delta_RMSF.png
    rmsfPng: str = p.join(analDir,"delta_RMSF.png")
    plt.savefig(rmsfPng,bbox_inches="tight")

#########################################################################################################
def histogram_plotting_manager(sysAnalDir: str) -> None:
    """
    Read through analysis output files, gather info from filenames, and package into a dataframe.
    
    Args:
        sysAnalDir (str): Directory containing analysis output files.
    """
    # Initialize list to store dataframes to concatenate
    dfsToConcat = []
    
    # Iterate through files in sysAnalDir
    for file in os.listdir(sysAnalDir):
        # Split file extension from filename
        fileData: tuple = p.splitext(file)
        
        # Skip if file is not a CSV
        if not fileData[1] == ".csv":
            continue
        
        # Split filename by underscores and extract relevant information
        fileData: list = fileData[0].split("_")
        dataTag: str = fileData[0]                   # Get type of data
        resAtomTag: str = fileData[1]   # Identifier of residue or atom of interest
        protName: str = "_".join(fileData[2:])                 # Protein name (contains repeats information)
        csvFile: str = p.join(sysAnalDir,file)   # Full path to CSV file
        
        # Skip if dataTag is not "contacts", "donor", or "acceptor"
        if not dataTag in ["contacts", "donor", "acceptor"]:
            continue
        
        # Create a dataframe row with relevant information and add to list
        row: pd.DataFrame = pd.DataFrame([[dataTag, resAtomTag, protName, csvFile]], columns=["dataTag", "resAtomTag", "protName", "csvPath"])
        dfsToConcat.append(row)
    
    # Concatenate dataframes in dfsToConcat into a single dataframe
    managerDf: pd.DataFrame = pd.concat(dfsToConcat, axis=0, ignore_index=True)
    
    # Get unique dataTags from managerDf
    uniqueDataTags: list = managerDf["dataTag"].unique().tolist()
    
    # Iterate through each unique dataTag
    for dataTag in uniqueDataTags:
        # Filter managerDf to get rows corresponding to dataTag
        groupedDf: pd.DataFrame = managerDf[managerDf["dataTag"] == dataTag]
        
        # Plot histogram of distances for each dataTag
        plot_distance_hist(groupedDf, dataTag, sysAnalDir)



#########################################################################################################
def plot_distance_hist(groupedDf: pd.DataFrame, dataTag: str, sysAnalDir: str) -> None:
    """
    Plots the histogram of distances for each column in the given dataframes.

    Parameters:
        groupedDf (pd.DataFrame): Dataframe containing information about the CSV files.
        dataTag (str): The type of data being plotted.
        sysAnalDir (str): Directory containing the CSV files.
    """
    
    ## load csv files and extract columnNames
    # Get a list of CSV file paths
    csvFiles: list = groupedDf["csvPath"].to_list()
    
    # Create empty lists to store dataframes and column names
    dfs: list = []
    columnNames: list = []
    
    # Iterate through each CSV file path
    for csvFile in csvFiles:
        # Read CSV file into a dataframe
        df: pd.DataFrame = pd.read_csv(csvFile, index_col="Unnamed: 0")
        
        # Append dataframe to dfs list and column names to columnNames list
        dfs.append(df)
        columnNames += df.columns.to_list()
    
    # Get unique column names
    columnNames: list = list(set(columnNames))

    # Generate a list of unique colors
    colors: list = plt.cm.viridis(np.linspace(0, 1, len(columnNames)))

    # Calculate the number of rows needed for the subplots
    num_rows: int = math.ceil(len(columnNames) / 4)

    # Create a subplot for each column in a 4xN grid
    fig, axs = plt.subplots(num_rows, 4, figsize=(20, 5*num_rows))
    axs = axs.flatten()

    # Iterate through each column name
    i: int = 0
    for colName in columnNames:
        for df in dfs:
            # Skip if column name is not in dataframe
            if not colName in df.columns:
                continue
            
            # Convert distances from nm to Angstrom
            distances_in_angstrom = df[colName] * 10
            
            # Get the color for the current column
            plotColor: str = colors[i]
            
            # Plot histogram of distances with specified color and label
            axs[i].hist(distances_in_angstrom, bins=30, alpha=0.5, label=colName, color=plotColor)
            axs[i].legend()
            
            # Set x-axis limits
            axs[i].set_xlim([2, 10])
            axs[i].set_ylim([0,100])
        
        # Set title for the current column
        axs[i].set_title(colName)

        i += 1
    
    # Remove empty subplots
    [fig.delaxes(ax) for ax in axs.flatten() if not ax.has_data()]

    # Adjust figure layout and save figure
    plt.tight_layout()
    plt.savefig(p.join(sysAnalDir, f"{dataTag}.png"), bbox_inches="tight")
    plt.close()
#############################################################################################
def plot_RMSF(sysAnalDir: str) -> None:
    """
    Plots the Root Mean Square Fluctuations (RMSF) per residue for all CSV files
    in the given directory.

    Args:
        sysAnalDir (str): Directory containing CSV files with RMSF data.
    """
    ## load contact dfs from sysAnalDir | concat into one big df
    # Initialize an empty list to store dataframes
    dfsToConcat: list = []
    
    # Loop through each file in the directory
    for file in os.listdir(sysAnalDir):
        # Skip if the file does not have a .csv extension
        if not p.splitext(file)[1] == ".csv":
            continue
        
        # Check if the file starts with "RMSF"
        if  file.startswith("RMSF"):
            # Read the CSV file into a dataframe and set the index column
            runDf: pd.DataFrame = pd.read_csv(p.join(sysAnalDir,file), index_col="RES_ID")
            
            # Append the dataframe to the list
            dfsToConcat.append(runDf)
    
    # If there are no dataframes to concat, return
    if len(dfsToConcat) == 0:
        return
    
    # Concatenate all dataframes into one big dataframe
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
    
    # Save the figure as RMSF.png
    rmsfPng: str = p.join(sysAnalDir,"RMSF.png")
    plt.savefig(rmsfPng,bbox_inches="tight")

#############################################################################################
def plot_RMSD(sysAnalDir: str) -> None:
    """
    Plots the Root Mean Square Displacement (RMSD) per time point for all CSV files
    in the given directory.

    Args:
        sysAnalDir (str): Directory containing CSV files with RMSD data.
    """
    ## load contact dfs from sysAnalDir | concat into one big df
    dfsToConcat: list = []  # Initialize an empty list to store dataframes
    
    # Loop through each file in the directory
    for file in os.listdir(sysAnalDir):
        # Skip if the file does not have a .csv extension
        if not p.splitext(file)[1] == ".csv":
            continue
        
        # Check if the file starts with "RMSD"
        if  file.startswith("RMSD"):
            # Read the CSV file into a dataframe and set the index column
            runDf: pd.DataFrame = pd.read_csv(p.join(sysAnalDir,file), index_col="Unnamed: 0")
            
            # Append the dataframe to the list
            dfsToConcat.append(runDf)
    
    # If there are no dataframes to concat, return
    if len(dfsToConcat) == 0:
        return
    
    # Concatenate all dataframes into one big dataframe
    df: pd.DataFrame = pd.concat(dfsToConcat, axis = 1, ignore_index=False) 
    
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot each column with a different color
    for col in df:
        plt.plot(df[col], label=col)

    # Calculate the mean for each row
    mean_values: pd.Series = df.mean(axis=1)

    # Plot the mean in bold black
    plt.plot(mean_values, color='black', linewidth=2, label='Mean')

    # Customize the plot
    plt.xlabel('Data Points')
    plt.ylabel('RMSD Values')
    plt.title('RMSD Values Over Time')
    plt.legend()
    
    # Save the figure as RMSD.png
    rmsdPng: str = p.join(sysAnalDir,"RMSD.png")
    plt.savefig(rmsdPng,bbox_inches="tight")

#############################################################################################
def plot_rdf(rdfDf: pd.DataFrame, outDir: str, resTag: str) -> None:
    """
    Plots the Radial Distribution Function (RDF) given a DataFrame containing
    the RDF data and saves the plot as a PNG file.

    Parameters:
        rdfDf (pd.DataFrame): DataFrame containing the RDF data.
        outDir (str): Output directory where the plot will be saved.
        resTag (str): Tag used in the output file name.

    Returns:
        None
    """
    # Plot the RDF
    plt.figure(figsize=(10, 6))  # Set figure size
    plt.plot(rdfDf['Radii'], rdfDf['RDF'], label='RDF')  # Plot the RDF
    plt.xlabel('Radii')  # Set x-axis label
    plt.ylabel('Radial Distribution Function')  # Set y-axis label
    plt.title('Radial Distribution Function vs Radii')  # Set plot title
    plt.legend()  # Add legend
    plt.grid(True)  # Enable grid
    
    # Save the plot as a PNG file
    output_file = p.join(outDir, f"RDF_{resTag}.png")
    plt.savefig(output_file, bbox_inches="tight")
    
    plt.close()  # Close the plot to free up memory
