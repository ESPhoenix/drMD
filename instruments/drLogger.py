## BASIC LIBS
import logging
import sys
import time

## CLEAN CODE
from typing import List, Dict, Union, Any

# from instruments.drCustomClasses import FilePath, DirectoryPath


class OverwriteStreamHandler(logging.StreamHandler):
    def emit(self, record):
        try:
            msg = self.format(record)
            stream = self.stream
            stream.write('\r' + msg)
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
        streamHandler.emit(logging.LogRecord(name='', level=logging.INFO, pathname='', lineno=0, msg=message, args=(), exc_info=None))

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