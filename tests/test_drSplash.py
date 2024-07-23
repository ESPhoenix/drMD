
## import basic libraries
import os
from os import path as p
import sys
## import testing libraries
import pytest

## add instruments to path
sys.path.insert(0, p.abspath(p.join(os.path.dirname(__file__), '..')))
## import drMD modules
from instruments import drSplash

#############################################################################################
def test_print_drMD_logo(capsys):
    drSplash.print_drMD_logo()
    captured = capsys.readouterr()
    assert "Molecular Dynamics: Just what the Doctor Ordered!" in captured.out
#############################################################################################
def test_print_config_error_no_error(capsys):
    with pytest.raises(SystemExit):
        drSplash.print_config_error()
    captured = capsys.readouterr()

    assert "⚕⚕⚕⚕⚕⚕⚕" in captured.out
#############################################################################################
def test_print_config_error_no_error(capsys):
    with pytest.raises(SystemExit):
        drSplash.print_config_error()
    captured = capsys.readouterr()
    assert "⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕" in captured.out
#############################################################################################
def test_print_config_error_with_error(capsys):
    error_message = "Sample error message"
    with pytest.raises(SystemExit):
        drSplash.print_config_error(error_message)
    captured = capsys.readouterr()
    assert "⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕" in captured.out
    assert f"-->\t{error_message}" in captured.out
#############################################################################################