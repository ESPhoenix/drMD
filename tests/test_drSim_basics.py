
## import basic libraries
import os
from os import path as p
import sys
## import testing libraries
import pytest

## add instruments to path
sys.path.insert(0, p.abspath(p.join(os.path.dirname(__file__), '..')))
## import openMM modules
import openmm
from openmm import app
import  simtk.unit  as unit


## import drMD modules
from instruments import drSim
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

class MockSimulation:
    def __init__(self):
        self.context = MockContext()
    
    def loadCheckpoint(self, file):
        self.loaded_checkpoint = file
    
    def loadState(self, file):
        self.loaded_state = file

class MockContext:
    def setTime(self, time):
        self.time = time
    
    def setStepCount(self, step_count):
        self.step_count = step_count



########

def test_init_system():
    # Create a mock prmtop object
    mock_prmtop = MockAmberPrmtopFile()
    
    # Call the function
    system = drSim.init_system(mock_prmtop)
    
    # Check the type of the returned object
    assert isinstance(system, openmm.System)
    
    # Check if the system has forces
    assert system.getNumForces() > 0

########

def test_process_sim_data():
    # Create a mock simulation dictionary
    sim = {
        "timestep": "2 fs",
        "duration": "1000 ps",
        "logInterval": "10 ps",
        "temp": 300
    }
    
    # Call the function
    processed_sim = drSim.process_sim_data(sim)
    
    # Check if the dictionary is marked as processed
    assert processed_sim["processed"] == True
    
    # Check the processed timestep
    expected_timestep = 2 * unit.femtoseconds
    assert processed_sim["timestep"] == expected_timestep
    
    # Check the processed duration
    expected_duration = 1000 * unit.picoseconds
    assert processed_sim["duration"] == expected_duration
    
    # Check the number of steps
    expected_nSteps = int(expected_duration / expected_timestep)
    assert processed_sim["nSteps"] == expected_nSteps
    
    # Check the number of log steps
    expected_logInterval = 10 * unit.picoseconds
    expected_logIntervalInSteps = int(round(expected_logInterval / expected_timestep))
    assert processed_sim["logInterval"] == expected_logIntervalInSteps
    
    # Check the number of log steps (every 500 steps)
    expected_nLogSteps = round(expected_nSteps / 500)
    assert processed_sim["nLogSteps"] == expected_nLogSteps
    
    # Check the processed temperature
    expected_temp = 300 * unit.kelvin
    assert processed_sim["temp"] == expected_temp


########
def test_init_reporters(tmp_path):
    # Create a mock simulation directory
    sim_dir = tmp_path / "simulation"
    sim_dir.mkdir()
    
    # Define the number of steps and report interval
    n_steps = 1000
    report_interval = 100
    
    # Call the function
    reporters = drSim.init_reporters(str(sim_dir), n_steps, report_interval)
    
    # Check if the dictionary contains the expected keys
    expected_keys = {"vitals", "progress", "trajectory", "checkpoint"}
    assert set(reporters.keys()) == expected_keys
    
    # Check the vitals reporter
    vitals_reporter, vitals_file = reporters["vitals"]
    assert isinstance(vitals_reporter, app.StateDataReporter)
    assert vitals_file == p.join(sim_dir, "vitals_report.csv")
    
    # Check the progress reporter
    progress_reporter, progress_file = reporters["progress"]
    assert isinstance(progress_reporter, app.StateDataReporter)
    assert progress_file == p.join(sim_dir, "progress_report.csv")
    
    # Check the trajectory reporter
    trajectory_reporter, trajectory_file = reporters["trajectory"]
    assert isinstance(trajectory_reporter, app.DCDReporter)
    assert trajectory_file == p.join(sim_dir, "trajectory.dcd")
    
    # Check the checkpoint reporter
    checkpoint_reporter, checkpoint_file = reporters["checkpoint"]
    assert isinstance(checkpoint_reporter, app.CheckpointReporter)
    assert checkpoint_file == p.join(sim_dir, "checkpoint.chk")


########
def test_load_simulation_state_chk(tmp_path):
    # Create a mock simulation object
    simulation = MockSimulation()
    
    # Create a mock checkpoint file
    chk_file = tmp_path / "state.chk"
    chk_file.write_text("mock checkpoint data")
    
    # Call the function
    modified_simulation = drSim.load_simulation_state(simulation, str(chk_file))
    
    # Check if the checkpoint was loaded
    assert modified_simulation.loaded_checkpoint == str(chk_file)
    
    # Check if the time and step count were reset
    assert modified_simulation.context.time == 0.0
    assert modified_simulation.context.step_count == 0

########

def test_load_simulation_state_xml(tmp_path):
    # Create a mock simulation object
    simulation = MockSimulation()
    
    # Create a mock XML file
    xml_file = tmp_path / "state.xml"
    xml_file.write_text("mock XML data")
    
    # Call the function
    modified_simulation = drSim.load_simulation_state(simulation, str(xml_file))
    
    # Check if the state was loaded
    assert modified_simulation.loaded_state == str(xml_file)
    
    # Check if the time and step count were reset
    assert modified_simulation.context.time == 0.0
    assert modified_simulation.context.step_count == 0

########
if __name__ == "__main__":
    pytest.main()