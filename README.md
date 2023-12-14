# drMD
Automated workflow for running molecular dynamics simulations with Amber and Openmm

# Installation

##Create and activate conda environment
```bash
conda create -n drMD python=3.10
conda activate drMD
```
##Install AmberTools (needs to be before OpenMM)
```bash
conda install -c conda-forge ambertools=23
``` 
##Install OpenMM
```bash
conda install -c omnia openmm
``` 
##Install other python libraries
```bash
pip install argpass
pip install pyyaml
pip install pandas
pip install propka
```