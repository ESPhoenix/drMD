
## import basic libraries
import os
from os import path as p
import sys
## import testing libraries
import pytest


## import openMM modules
import openmm
from openmm import app
import  simtk.unit  as unit

## add instruments to path
sys.path.insert(0, p.abspath(p.join(os.path.dirname(__file__), '..')))
## import drMD modules
from instruments import drSim
from instruments import drFirstAid
from instruments import drCheckup
from instruments import drRestraints
from instruments import drFixer
from instruments import drClusterizer
########
## INIT SOME MOCK CLASSES
class MockAmberPrmtopFile:
    def createSystem(self, nonbondedMethod, nonbondedCutoff, constraints):
        # Create a mock system with minimal setup
        system = openmm.System()
        
        # Add a mock force to the system
        force = openmm.CustomExternalForce("0")
        system.addForce(force)
        
        return system
## MAKE SOME MOCK CLASSES
class MockSimulation:
    def __init__(self, topology, system, integrator, platform):
        self.topology = topology
        self.system = system
        self.integrator = integrator
        self.platform = platform
        self.context = MockContext()
        self.reporters = []
        self.currentStep = 0
    
    def step(self, steps):
        self.currentStep += steps
    
    def saveState(self, file):
        self.saved_state = file
    
    def loadCheckpoint(self, file):
        self.loaded_checkpoint = file
    
    def loadState(self, file):
        self.loaded_state = file

class MockContext:
    def setTime(self, time):
        self.time = time
    
    def setStepCount(self, step_count):
        self.step_count = step_count
    
    def getState(self, getPositions=False, getEnergy=False):
        return MockState()

class MockState:
    def getPositions(self):
        return []

