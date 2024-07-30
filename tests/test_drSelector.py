## import basic libraries
import os
from os import path as p
import sys
import pandas as pd
## import testing libraries
import pytest

## import openMM modules
import openmm
from openmm import app
import  simtk.unit  as unit


## add instruments to path
sys.path.insert(0, p.abspath(p.join(p.dirname(__file__), '..')))
## import drMD modules
from instruments import drSelector

########

# Mock pdbUtils module
class MockPdbUtils:
    @staticmethod
    def pdb2df(pdbFile):
        # Create a mock DataFrame representing a PDB file
        data = {
            "ATOM_NAME": ["N", "CA", "C", "O", "H", "O", "H", "NA", "CL"],
            "RES_NAME": ["ALA", "ALA", "ALA", "ALA", "HOH", "HOH", "HOH", "NA", "CL"],
            "CHAIN_ID": ["A", "A", "A", "A", "A", "A", "A", "A", "A"],
            "RES_ID": [1, 1, 1, 1, 2, 2, 2, 3, 4]
        }
        return pd.DataFrame(data)

# Mock init_name_lists function
def mock_init_name_lists():
    aminoAcidResNames = ["ALA", "GLY", "SER"]
    backboneAtomNames = ["N", "CA", "C", "O"]
    solventResNames = ["HOH"]
    ionResNames = ["NA", "CL"]
    return aminoAcidResNames, backboneAtomNames, solventResNames, ionResNames

########

@pytest.fixture
def mock_pdb_utils(monkeypatch):
    monkeypatch.setattr(drSelector.pdbUtils, "pdb2df", MockPdbUtils.pdb2df)
    monkeypatch.setattr(drSelector, "init_name_lists", mock_init_name_lists)

########

def test_get_atom_indexes_all(mock_pdb_utils):
    selection = {"keyword": "all"}
    pdbFile = "mock.pdb"
    atom_indexes = drSelector.get_atom_indexes(selection, pdbFile)
    assert atom_indexes == [0, 1, 2, 3, 4, 5, 6, 7, 8]

########

def test_get_atom_indexes_backbone(mock_pdb_utils):
    selection = {"keyword": "backbone"}
    pdbFile = "mock.pdb"
    atom_indexes = drSelector.get_atom_indexes(selection, pdbFile)
    assert atom_indexes == [0, 1, 2, 3]

########

def test_get_atom_indexes_protein(mock_pdb_utils):
    selection = {"keyword": "protein"}
    pdbFile = "mock.pdb"
    atom_indexes = drSelector.get_atom_indexes(selection, pdbFile)
    assert atom_indexes == [0, 1, 2, 3]

########

def test_get_atom_indexes_water(mock_pdb_utils):
    selection = {"keyword": "water"}
    pdbFile = "mock.pdb"
    atom_indexes = drSelector.get_atom_indexes(selection, pdbFile)
    assert atom_indexes == [4, 5, 6]

########

def test_get_atom_indexes_ions(mock_pdb_utils):
    selection = {"keyword": "ions"}
    pdbFile = "mock.pdb"
    atom_indexes = drSelector.get_atom_indexes(selection, pdbFile)
    assert atom_indexes == [7, 8]

########

def test_get_atom_indexes_ligands(mock_pdb_utils):
    selection = {"keyword": "ligands"}
    pdbFile = "mock.pdb"
    atom_indexes = drSelector.get_atom_indexes(selection, pdbFile)
    assert atom_indexes == []

########

def test_get_atom_indexes_custom(mock_pdb_utils):
    selection = {
        "keyword": "custom",
        "customSelection": [
            {"CHAIN_ID": "A", "RES_NAME": "ALA", "RES_ID": 1, "ATOM_NAME": "N"},
            {"CHAIN_ID": "A", "RES_NAME": "HOH", "RES_ID": 2, "ATOM_NAME": "O"}
        ]
    }
    pdbFile = "mock.pdb"
    atom_indexes = drSelector.get_atom_indexes(selection, pdbFile)
    assert atom_indexes == [0, 5]

########
def test_get_atom_indexes_custom_wildcards(mock_pdb_utils):
    selection = {
        "keyword": "custom",
        "customSelection": [
            {"CHAIN_ID": "A", "RES_NAME": "ALA", "RES_ID": 1, "ATOM_NAME": "all"}
        ]
    }
    pdbFile = "mock.pdb"
    atom_indexes = drSelector.get_atom_indexes(selection, pdbFile)
    assert atom_indexes == [0, 1, 2, 3]
########

if __name__ == "__main__":
    pytest.main()