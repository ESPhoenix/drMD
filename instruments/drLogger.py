## BASIC LIBS
import logging
import sys
import time
import pandas as pd
from functools import wraps
import threading
import shutil
from os import path as p
## CLEAN CODE
from typing import List, Dict, Union, Any, Tuple
try:
    from instruments.drCustomClasses import FilePath, DirectoryPath
except:
    from drCustomClasses import FilePath, DirectoryPath

class OverwriteStreamHandler(logging.StreamHandler):
    def emit(self, record):
        try:

            terminalWidth = shutil.get_terminal_size().columns
            msg = self.format(record)
            stream = self.stream
            stream.write('\r' + msg + " "*(terminalWidth-len(msg)))
            stream.flush()
        except Exception:
            self.handleError(record)

def setup_logging(logFile):
    """
    Set up logging configuration to write to the specified log file.
    """
    # Get the root logger
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)
    

    # Create the file handler
    fileHandler = logging.FileHandler(logFile)
    fileFormatter = logging.Formatter('%(message)s')
    fileHandler.setFormatter(fileFormatter)
    
    # Remove any existing handlers
    for handler in rootLogger.handlers[:]:
        rootLogger.removeHandler(handler)
    
    # Add the file handler to the root logger
    rootLogger.addHandler(fileHandler)

def log_info(message, printToTerminal=False, persist=False):
    """
    Log an info message to the log file and optionally print to the terminal.
    """
    # Log to the file
    logging.info(message)
    
    # Optionally print to the terminal
    if printToTerminal:
        if persist:
            streamHandler = logging.StreamHandler(sys.stdout)
        else:
            streamHandler = OverwriteStreamHandler(sys.stdout)
        streamFormatter = logging.Formatter('%(message)s')
        streamHandler.setFormatter(streamFormatter)
        streamHandler.emit(logging.LogRecord(name='',
                                              level=logging.INFO,
                                                pathname='',
                                                  lineno=0, 
                                                  msg=message,
                                                    args=(),
                                                      exc_info=None))

def close_logging():
    """
    Close and remove all handlers from all loggers.
    """
    # Get the root logger
    root_logger = logging.getLogger()
    
    # Remove and close all handlers associated with the root logger
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()
    
    # Also remove and close handlers from all other loggers
    for logger_name in logging.root.manager.loggerDict:
        logger = logging.getLogger(logger_name)
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            handler.close()


def read_simulation_progress(progressReporterCsv: FilePath) -> Tuple[str, str]:
    try:
        progressDf: pd.DataFrame = pd.read_csv(progressReporterCsv)
        progressPercent: str = str(progressDf['#"Progress (%)"'].iloc[-1])
        timeRemaining: str = str(progressDf["Time Remaining"].iloc[-1])
        averageSpeed: str = str(round(progressDf["Speed (ns/day)"].mean()))
        return progressPercent, timeRemaining, averageSpeed
    except (FileNotFoundError, pd.errors.EmptyDataError) as e:
        return "N/A", "N/A", "N/A"
    

def monitor_progress_decorator(checkInterval=20):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            startTime = time.time()
            monitoring = threading.Event()

            def monitor_progress():
                stepName = kwargs["sim"]["stepName"]
                protName = kwargs["config"]["proteinInfo"]["proteinName"]
                outDir = kwargs["outDir"]
                simDir = p.join(outDir, stepName)
                progressReporterCsv = p.join(simDir, "progress_report.csv")
                while not monitoring.is_set():
                    progressPercent, timeRemaining, averageSpeed = read_simulation_progress(
                        progressReporterCsv)
                    log_info(f"-->{' '*4}Running {stepName} Step for: "
                             f"{protName} | Progress: {progressPercent} | "
                             f"Time Remaining: {timeRemaining} | "
                             f"Average Speed: {averageSpeed} ns/day ", True)
                    time.sleep(checkInterval)

            monitorThread = threading.Thread(target=monitor_progress)
            monitorThread.start()
            try:
                result = func(*args, **kwargs)
            except SystemExit as e:
                if e.code == 1:
                    monitoring.set()
                    monitorThread.join()
                    exit(1)
                else:
                    raise
            monitoring.set()
            monitorThread.join()
            return result
        return wrapper
    return decorator

          
    




if __name__ == '__main__':
    logFile = 'exampleLogFile.log'
    setup_logging(logFile)

    # Log messages using the custom log_info function
    log_info('This is the first message written to file only.')
    time.sleep(1)  
    log_info('This is the second message written to file and terminal.', printToTerminal=True)
    time.sleep(1) 
    log_info('This is the third message written to file and terminal.', printToTerminal=True)
    time.sleep(1)  
    log_info('This is the final message written to file and terminal.\n', printToTerminal=True)  # Add a newline at the end to finalize