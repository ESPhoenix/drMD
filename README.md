# drMD
Automated workflow for running molecular dynamics simulations with Amber and Openmm
# Installation
1. Clone this repository
```bash
git clone https://github.com/ESPhoenix/drMD.git
```
2. Create and activate conda environment
```bash
conda create -n drMD python=3.10
```
```bash
conda activate drMD
```
3. Install AmberTools (needs to be before OpenMM) with conda
```bash
conda install -c conda-forge ambertools=23
``` 
4. Install OpenMM with conda
```bash
conda install -c omnia openmm
``` 
5. Install OpenBabel with conda
```bash
conda install -c conda-forge openbabel
```
6. Install other python libraries with pip
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
Inputs are grouped by theme and are stored as nested dictionaries and lists
The next few sections will detail the correct formatting of the config.yaml file

### pathInfo
The pathInfo entry in the config file is a dictionary containing two parameters:
- inputDir:     This is the absoloute path towards a directory containing PDB files that will be used as starting points for your simulations.
            
  > :medical_symbol:
  > **To Perform Replicate** simulations, simply create copies of your starting PDB files in the inputDir, with each copy
  >named with a unique number. For example, your inputDir could contain my_protein_1.pdb, my_protein_2.pdb, etc.


- outputDir:    This is the absoloute path towards a directory that you want your drMD outputs to be written to.

  > The outputDir will be created if it does not already exist at the point of running drMD

  > Within outputDir, a directory will be created for each PDB file contained in inputDir, in this document, thesesubdirectories will be refered to as **runDirs**

Example pathInfo:
```yaml
pathInfo:
  inputDir: "/home/esp/scriptDevelopment/drMD/01_inputs"
  outputDir: "/home/esp/scriptDevelopment/drMD/02_outputs"
```
        
### hardwareInfo
This config entry tells drMD about your computer hardware and how you want to use it to run your simulations
The hardwareInfo entry in the config file is a dictionary containing three parameters:
- platform: This is the platform that will be used to run simulations in OpenMM. 
  > Accepted arguments for **platform** are "CUDA", "OpenCL", and "CPU"
- paralellCPU:      This is the number (int) of simulations that will be run in paralell
- subprocessCpus:  This is the number (int) of cpu cores that will be allocated to each simulation. It is recommended to set this to 1. 
  > The total CPU usage will be parallelCPU * subprocessCpus, so make sure you have enough CPus when you set these parameters

Example hardwareInfo:
```yaml
hardwareInfo:
  parallelCPU: 1
  platform: "CUDA"
  subprocessCpus: 1
```

### ligandInfo
The ligandInfo entry in the config file is optional and must be used if your PDB files have organic ligands or cofactors.
These small molecules will not have parameters in the AMBER forcefield, drMD will run an automated protocol to generate these parameters for you.
To do this, you will need to tell drMD some things about each ligand you whish tp simulate.


>[!NOTE]
>The ligandInfo entry is optional. drMD will automatically detect ligands in your PDB files. It will also detect parameter 
>files in your input directory

ligandInfo is a list of dictionaries that contain the following parameters:

- ligandName:   This is the three letter name of the ligand, this will be used to match the residue names in your PDB files

- protons:      This is a boolean (TRUE/FALSE) to tell drMD whether you have protons on your ligand. 
                If set to FALSE, drMD will run an automated protonation protocol to add protons to your ligand
  > [!CAUTION]
  > The automatic protonation protocol only works reliably for simple organic ligands.

  > For more complex ligands, we recommended that you manually add protons in your input PDB file prior to running drMD

- charge:       This is the charge of the ligand (int)

- toppar:       This is a boolean (TRUE/FALSE) to tell drMD whether you have an frcmod file for your ligand already made.
                If you already have one, it must be located in the 01_ligand_parameters directory within your outputDir

- mol2:         This is a boolean (TRUE/FALSE) to tell drMD whether you have a mol2 file for your ligand already made.
                If you already have one, it must be located in the 01_ligand_parameters directory within your outputDir

Example ligandInfo:
```yaml
ligandInfo:
  - ligandName: "FMN"
    protons: True
    charge: -1
    toppar: False
    mol2: False
  - ligandName: "TPA"
    protons: True
    charge: -2
    toppar: False
    mol2: False
```
---

### simulationInfo
This is the real meat and potatoes of the drMD config file. 

The simulationInfo entry in the config file is a list of dictionaries containing information about each simulation.

Each simulation detailed in simulationInfo will be run in sequence, with the output of the previous simulation being the starting point for the next simulation.
#### Selecting Simulation Type
Each simulation dictionary contains the following parameters:
- stepName:         This is the name of the step that will be used to create a subdirectory in the runDir, we reccomend prefixing these names with numbers to make them order nicely

- simulationType: This is the type of simulation that will be run. Accepted arguments are:

    - "EM":         This will run a steepest-decent Energy Minimisation step. 
    > 💡 **TIP:**
    > We reccomended that you run one of these steps before any other simulation steps
    - "NVT":        This will run an NVT (constant volume) molecular dynamics simulation
    - "NPT":        This will run an NPT (constant pressure) molecular dynamics simulation
    > 💡 **TIP:**
    > For the majority of protein simulations, the NPT ensemble is used for production MD simulations, while the NVT ensemble is oonly used in equilibration steps
    - "META":       This will run a Metadynamics simulation

#### Selecting simulation temperature 
For most simulations, a constant temperature is used. In this case the following parameter is required:

- temperature: This is the temperature (int) of the simulation in Kelvin 

If you wish to change the temperature throughout the simulation, the following parameter is required:

- temperatureRange: This is a list of integers (again, in Kelvin) that will be used to change the temperature throughout the simulation. 

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
For "normal" MD simulations using NVT or NpT ensembles, as well as for Metadynamics simulations, the following additional parameters are required:
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
- restraintInfo:       This is a list of dictionaries containing information about each restraints. 
                      You can add as many restraints as you wish, with one dictionary per restraints
                      Each restraints dictionary contains the following parameters:

- restraintType: This is the type of restraints that will be added. Accepted arguments are: "distance", "angle", "dihedral", "position"

- parameters:   This is a dictionary containing the parameters for the restraints.

 All restraints need to have the following parameters:

- k :        This is the force constant of the restraint (int)

- selection: This is the selection of atoms that the restraints will be applied to. 

>The selection method is shared between multiple different inputs in the drMD config file. This is described in more detail in the next section

 Distance restraints require the following parameter:
 - r0:        This is the distance in Angstroms that the restraint should be applied at (int/float)

Angle restrainst require the following parameter:
- theta0:    This is the angle in degrees that the angle should be constrained at (int/float)

Dihedral restraints require the following parameter:
- phi0:      This is the angle in degrees that the dihedral should be constrained at (int/float)

All restraints require the selection parameter. This tells drMD what atoms to apply the restraint to
    - selection:    This is a dictionary containing information on the selection of atoms that the restraints will be applied to.
 

Example restraints syntax:
```yaml
    restraints:
    - restraintType: "position"
      selection:
        keyword: "protein"
        parameters:
          k: 1000
    - type: "distance"
      selection: 
        keyword: "custom"  
        selectionSyntax:
          - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 1, ATOM_NAME: "CA"}
          - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 2, ATOM_NAME: "CA"}
      parameters:
        k: 1000
        r0: 3

```

This example will add a position restraints to the protein atoms and 
a 3 Angstrom distance restraint between the CA atoms of residues 1 and 2 of the protein

#### drMD Selection syntax

When creating restraints, metadynamics bias variables or running post-simulation clustering, you will need to specify the selection of atoms that the restraints will be applied to. To do this, you will need to supply a "selection" dictionary. This dictionary must contain the following parameter:

- keyword: This is the keyword that will be used to specify the selection. Accepted arguments are:
  - **protein** : This will select all protein atoms in your system
  - **water** : This will select all water molecules in your system
  - **ions**: This will select all ions in your system
  - **ligand** : This will select all non-protein, non-water and non-ion atoms in your system
  - **custom** : This will select all atoms that match the selectionSyntax

If you have used the **custom** keyword, you will need to use an additional parameter the selection dictionary:

- customSelection: This is a list of dictionaries containing details of the atoms that will be selected.
Each dictionary in the list must contain the following parameters:
  - CHAIN_ID: (str) This is the chain ID of the atom to be selected
  - RES_NAME: (str) This is the three-letter residue name of the atom to be selected
  - RES_ID: This (int) is the residue ID of the atom to be selected
  - ATOM_NAME: (str) This is the atom name of the atom to be selected

> The customSelection syntax accepts the wildcards by using "all" in the place of a dictionary value. 

Example customSelection syntax:

```yaml
selection:
  keyword: "custom"
  customSelection:
    - {CHAIN_ID: "A", RES_NAME: "all", RES_ID: "all", ATOM_NAME: "CA"}
    - {CHAIN_ID: "B", RES_NAME: "FMN", RES_ID: "all", ATOM_NAME: "all"}
```
In the above example, a selection containing all CA atoms in chain A as well as all atoms in the residue FMN in chain B

#### Running Metadynamics with drMD
To run a metadynamics simulation, first set the simulation type to "META".
To do this, you will need to specify the following parameters:
- metaDynamicsInfo:   This is a dictionary containing the parameters for the Metadynamics simulation,
    *Note The requried parameters within this dictionary are:
    - height:               This is the height (int) parameter used in the Metadynamics simulation

    - biasFactor:           This is the bias factor (int) parameter used in the Metadynamics simulation

    - biases:               This is a list of dictionaries containing information about each biasVariable.

> You will also need to specify at least one biasVariable for the simulation to sample.
> You can specify as many biasVariables as you wish, with one dictionary per biasVariable

Each dictionary in biases must contain the following parameters:

- biasVar:      This is the type of biasVariable that will be added. Accepted arguments are "RMSD", "torsion"
              TODO: add "distance" and "angle" options
- selection:    This is a dictionary containing information on the selection of atoms that the biasVariable will be applied to.
              The selection syntax is identical to that used for the restraints. (described above)

> Note that for distance bias variables, the selection type must be "atoms" with two atoms selected

> Note that for angle bias variables, the selection type must be "atoms" with three atoms selected

> Note that for torsion bias variables, the selection type must be "atoms" with four atoms selected

Example MetaDynamics syntax:
```yaml
    metaDynamicsInfo:
      height: 2
      biasFactor: 5
      biases: 
        - biasVar: "RMSD"
          selection: 
            keyword: "backbone"
        - biasVar: "torsion"
          selection: 
            keyword: "custom"
          selectionSyntax:
            - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 1, ATOM_NAME: "CA"}
            - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 2, ATOM_NAME: "CA"}
            - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 3, ATOM_NAME: "CA"}
            - {CHAIN_ID: "A", RES_NAME: "ALA", RES_ID: 4, ATOM_NAME: "CA"}

```
This example will add a RMSD bias to the backbone of the protein
and a torsion bias between the CA atoms of residues 1, 2, 3, and 4 of the protein

#### Post-simulation processing
After all of your simulations have been run, drMD contains some simple utilities for organising your output files and deleting any unwanted files.

If you want to do any post-processing, you will need to provide the following parameter in your config file:
- postSimulationInfo : This is a dictionary containing the parameters for the post-simulation processing

If you wish to collect PDB files that represent the last frame of each simulation, you may include the following parameter in postSimulationInfo:
- endpointInfo: This is a dictionary containing the parameters the following parameters:

  - stepNames : This is a list of strings containing the names of the steps in the simulation, these should match the stepNames that you have used in your simulationInfo dictionary (described above). Endpoint PDB files will be gathered for these steps

  - removeAtoms: This is a list of dictionaries containing the selections of atoms to be removed from the PDB files. The format for these selections is the same as that used for the restraints, metadynamics, etc. (described above)

Moleculatr Dyamics simulations can generate very large output files that can become rather unweildy and difficulat to anaylse. One way to quickly see the most important parts of your simulation is to perform clustering on your simulation trajectories. To do this with drMD, include the following parameter in your config file:

- clusterTrajectory: This is a dictionary containing the parameters following parameters:

  - stepNames: This is a list of strings containing the names of the steps in the simulation, these should match the stepNames that you have used in your simulationInfo dictionary (described above). Clustering will be performed on trajectories of these steps
  - nClusters: This is the number (int) of clusters PDB files that will be generated
  - clusterBy: This is a list of selections of atoms (described above) to cluster on. 
> By selecting only the parts of our system that you are interested in with the clusterBy parameter, you can generated cluster PDB files where these atoms are most separated by RMSD.
