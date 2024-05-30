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
## Run example pdbs
```bash
python path/to/drMD/batch_drMD.py --config path/to/drMD/batch_config.yaml
```
#### To run your own pdb files, modify `config.yaml` file 