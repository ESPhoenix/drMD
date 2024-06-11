# drMD
Automated workflow for running molecular dynamics simulations with Amber and Openmm

# Installation
## Clone this repository
```bash
git clone https://github.com/ESPhoenix/drMD.git
```
## Create and activate conda environment
```bash
conda create -n drMD python=3.10
```
```bash
conda activate drMD
```
## Install AmberTools (needs to be before OpenMM)
```bash
conda install -c conda-forge ambertools=23
``` 
## Install OpenMM
```bash
conda install -c omnia openmm
``` 
## Install OpenBabel
```bash
conda install -c conda-forge openbabel
```

## Install other python libraries
```bash
pip install -r requirements.txt
```

# Usage

Now that you have sucessfully set up the dependancies for drMD, you are nearly ready to run some biomolecular simulations!

To the drMD script takes a single config.yaml file as a command line argument, using the "--config" flag:
```bash
python drMD.py --config config.yaml
```
This config file contains all of the user inputs drMD needs to run a series of biomolecular simulations.
The following section will detail the correct formatting of this config.yaml file

## Config syntax
The config.yaml file is in the YAML format *(https://en.wikipedia.org/wiki/YAML)* 
Inputs are grouped by theme and are stored as nested dictionaries and lists according as needed

### pathInfo
The pathInfo entry in the config file is a dictionary containing two parameters:
- inputDir:     This is the absoloute path towards a directory containing PDB files that will be used as starting points for your simulations.

                >Note that if you whish to perform multiple repeats of your simulation protocols upon one PDB file, simply duplicate it within 
                >this directory with suitable changes to the filenames

- outputDir:    This is the absoloute path towards a directory that you want your drMD outputs to be written to.

                >Note that this file does not need to exist at the point of running drMD, the script will create outputDir if it does not already exist

                >Note that within outputDir, a directory will be created for each PDB file contained in inputDir, 
                >in this document, these subdirectories will be refered to as runDirs*

Example pathInfo:
```yaml
pathInfo:
  inputDir: "/home/esp/scriptDevelopment/drMD/01_inputs"
  outputDir: "/home/esp/scriptDevelopment/drMD/03_outputs"
```
        
### generalInfo
The generalInfo entry in the config file is a dictionary containing three parameters:
- platform:         This is the platform that will be used to run simulations in OpenMM. Accepted arguments are "CUDA", "OpenCL", and "CPU"

- paralellCPU:      This is the number (int) of simulations that will be run in paralell, if using "CUDA" or "OpenCL" options in the platform 
                    parameter, it is recommended to set this to 1. 

- subprocessCpus:  This is the number (int) of cpu cores that will be allocated to each simulation. It is recommended to set this to 1. 

Example generalInfo:
```yaml
generalInfo:
  parallelCPU: 1
  platform: "CUDA"
  subprocessCpus: 1
```

### ligandInfo
The ligandInfo entry in the config file is optional and must be used if your PDB files have organic ligands or cofactors.
These small molecules will not have parameters in the AMBER forcefield, drMD will run an automated protocol to generate these parameters for you.
To do this, you will need to tell drMD some things about each ligand you whish tp simulate.
ligandInfo contains two parameters:
- nLigands:         This is the number (int) of ligands in your PDB files

- ligands:          This is a list of dictionaries containing information about each ligand. 
                    You can add as many ligands as you wish, with one dictionary per ligand
                    Each ligand dictionary contains the following parameters:

    - ligandName:   This is the three letter name of the ligand, this will be used to match the residue names in your PDB files

    - protons:      This is a boolean (TRUE/FALSE) to tell drMD whether you have protons on your ligand. 
                    If set to FALSE, drMD will run an automated protonation protocol to add protons to your ligand
                    *Note that the automatic protonation protocol only works reliably for simple organic ligands.
                    If your ligand is complex it is recommended to manually add protons in your input PDB file prior to running drMD*

    - charge:       This is the charge of the ligand (int)

    - toppar:       This is a boolean (TRUE/FALSE) to tell drMD whether you have an frcmod file for your ligand already made.
                    If you already have one, it must be located in the 01_ligand_parameters directory within your outputDir

    - mol2:         This is a boolean (TRUE/FALSE) to tell drMD whether you have a mol2 file for your ligand already made.
                    If you already have one, it must be located in the 01_ligand_parameters directory within your outputDir
Example ligandInfo:
```yaml
ligandInfo:
  nLigands: 1
  ligands:
    - ligandName: "FMN"
      protons: True
      charge: -1
      toppar: False
      mol2: False
```

### simulationInfo
This is the real meat and potatoes of drMD. 
The simulationInfo entry in the config file is a list of dictionaries containing information about each simulation.
Each simulation detailed in simulationInfo will be run in sequence, with the output of the previous simulation being the starting point for the next simulation.
#### Manditory parameters
Each simulation dictionary contains the following parameters:
- stepName:         This is the name of the step that will be used to create a subdirectory in the runDir, we reccomend prefixing these names with
                    numbers to make them order nicely

- type:             This is the type of simulation that will be run. Accepted arguments are:

    - "EM":         This will run a steepest-decent Energy Minimisation step. 
                    *Note that it is reccomended to run one of these steps before any other simulation steps*
    - "NVT":        This will run an NVT molecular dynamics simulation
    - "NPT":        This will run an NPT molecular dynamics simulation
    - "META":       This will run a Metadynamics simulation

- temp:             This is the temperature (int) of the simulation in Kelvin 

Depending on the type of simulation specified, additional parameters will be required.
#### Energy Minimisation Pararameters
For Energy Minimisation steps, the following additional parameters are required:
- maxIterations:    This is the maximum number (int) of iterations that will be run in the Energy Minimisation step
                    If this parameter is set to -1, the step will run until the energy converges

Example Energy Minimisation syntax:
```yaml
simulationInfo:
  - stepName: "01_energy_minimisation"
    type: "EM"
    temp: 300
    maxIterations: -1
```
This will run a energy minimisation until the energy converges

#### Generic Simulation Parameters
For "normal" simulations using NVT or NPT ensembles, as well as for Metadynamics simulations, the following additional parameters are required:
- duration:         This is the duration of the simulation step, as a string "int unit" eg. "1000 ps"

- timestep:         This is the timestep of the simulation, as a string "int unit" eg. "2 fs"

- logInterval:      This is the frequency that the simulation will write to file using built-in OpemMM reporters. As a string "int unit" eg. "100 ps"


Example NVT simulation syntax:
```yaml
simulationInfo:
  - stepName: "02_NVT_pre-equilibraition"
    type: "NVT"
    duration: "100 ps"
    timestep: "2 fs"
    temp: 300
    logInterval: "10 ps"
```
This will run a 100 ps NVT molecular dynamics simulation with a timestep of 2 fs, a temp of 300 and a logInterval of 10 ps

#### Adding Restraints with drMD
If you whish to perform simulations (not Energy Minimisations) with restraints, the following additional parameters are required:
- restraints:       This is a list of dictionaries containing information about each restraints. 
                    You can add as many restraints as you wish, with one dictionary per restraints
                    Each restraints dictionary contains the following parameters:

    - type:         This is the type of restraints that will be added. Accepted arguments are: "distance", "angle", "dihedral", "position"

    - parameters:   This is a dictionary containing the parameters for the restraints.

    All restraints need to have the following parameter:
        -k :        This is the force constant of the restraint (int)

    Distance restraints require the following parameter:
        -r0:        This is the distance in Angstroms that the restraint should be applied at (int/float)

    Angle restrainst require the following parameter:
        -theta0:    This is the angle in degrees that the angle should be constrained at (int/float)

    Dihedral restraints require the following parameter:
        -phi0:      This is the angle in degrees that the dihedral should be constrained at (int/float)

    - selection:    This is a dictionary containing information on the selection of atoms that the restraints will be applied to.
                    It contains the following parameters:
        - type:     This can be set the following presets: "backbone", "protein", "ligand", "water", or "ions".
                    If a preset is used, atoms will be selected automatically for you and no other parameters need to be set.
                    Other options for more granular selections are "residue" or "atom".

    If "residue" has been selected and additional parameter must be used:
        - input:    This is a list of lists containing the information needed to select resiudes. This must use the following format:
                    [[chainId, residueName, residueNumber], [chainId, residueName, residueNumber], ...]      

    If "atom" has been selected and additional parameter must be used:
        - input:    This is a list of lists containing the information needed to select atoms. This must use the following format:
                    [[chainId, residueName, residueNumber, atomName], [chainId, residueName, residueNumber, atomName], ...]      

Example restraints syntax:
```yaml
    restraints:
    - type: "position"
      selection:
        type: "protein"
        parameters:
          k: 1000
    - type: "distance"
      selection: 
        type: "atoms"  
        parameters:
          k: 1000
          r0: 3
        input:
          - ["A", "ALA", 1, "CA"]
          - ["A", "ALA", 2, "CA"]
```
This example will add a position restraints to the protein atoms and 
a 3 Angstrom distance restraint between the CA atoms of residues 1 and 2 of the protein

#### Running Metadynamics with drMD
To run a metadynamics simulation, first set the simulation type to "META".
To do this, you will need to specify the following parameters:
- metaDynamicsParameters:   This is a dictionary containing the parameters for the Metadynamics simulation,
                            requried parameters within this dictionary are:
    - height:               This is the height (int) parameter used in the Metadynamics simulation

    - biasFactor:           This is the bias factor (int) parameter used in the Metadynamics simulation

You will also need to specify at least one biasVariable for the simulation to sample.
You can specify as many biasVariables as you wish, with one dictionary per biasVariable
This is within the following biases list:
- biases:    This is a list of dictionaries containing information about each biasVariable.

Within each biasVariable dictionary, you need to specify the following parameters:  
    - biasVar:      This is the type of biasVariable that will be added. Accepted arguments are: "RMSD", "torsion"
                    TODO: add "distance" and "angle" options
    - selection:    This is a dictionary containing information on the selection of atoms that the biasVariable will be applied to.
                    The selection syntax is identical to that used for the restraints.
                    *Note that for distance bias variables, the selection type must be "atoms" with two atoms selected*
                    *Note that for distance bias variables, the selection type must be "atoms" with three atoms selected*
                    *Note that for torsion bias variables, the selection type must be "atoms" with four atoms selected*

Example MetaDynamics syntax:
```yaml
    metaDynamicsParams:
      height: 2
      biasFactor: 5
    biases: 
      - biasVar: "RMSD"
        selection: 
          type: "backbone"
    - biasVar: "torsion"
        selection: 
          type: "atom"
          input:
            - ["A", "ALA", 1, "CA"]
            - ["A", "ALA", 2, "CA"]
            - ["A", "ALA", 3, "CA"]
            - ["A", "ALA", 4, "CA"]
```
This example will add a RMSD bias to the backbone of the protein
and a torsion bias between the CA atoms of residues 1, 2, 3, and 4 of the protein

#### Clustering your outputs
Molecular dyamics simulations can generate very large output files that can become rather unweildy and difficulat to anaylse.
drMD offers a functionality to perform k-means clustering on your output trajectories. 
This will return a managable number of individual PDB files for you to analyse or use in the next step of your pipeline. 

To perform this clustering, you may include the follwoing optional parameters in any (not EM) of you simulation dictionaries:
- clusterTrajectory:       This is a dictionary containing the parameters for the trajectory clustering.

Within this dictionary, you need to specify the following parameters:
    - clusterBool:          This is a boolean (TRUE/FALSE) to tell drMD whether you want to perform clustering
    - nClusters:            This is the number (int) of clusters PDB files that will be generated
    - selection:            The clustering is performed using RMSD, this can be calculated on any selection of atoms
                            The selection syntax is identical to that used for the restraints, described above.


### cleanUpInfo
Once all of your simulationsare complete, this entry contains optionals for processing your output files
The cleanUpInfo entry in the config file is a dictionary containing four parameters:
- getEndpointPdbs:        This is a boolean (TRUE/FALSE) to tell drMD whether you want to get the pdb files at the end of each simulation
- removeWaters:           This is a boolean (TRUE/FALSE) to tell drMD whether you want to remove the water molecules from your endpoint PDB files
- removeIons:             This is a boolean (TRUE/FALSE) to tell drMD whether you want to remove the ions from your endpoint PDB files
- keepFileNames:          This is a boolean (TRUE/FALSE) to tell drMD whether you want to keep the original filenames in your endpoint PDB files

This section is only relevant if getEndpointPdbs is set to TRUE
TODO: maybe move this section to within each sim dictionary?? Is this redundant with clustering?

