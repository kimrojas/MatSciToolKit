from abc import ABC, abstractmethod
from matscitoolkit.analysis_workflow.logger import get_log


def test_logger():
    log = get_log()
    # Example logging messages
    log.debug('This is a debug message')
    log.info('This is an info message')
    log.warning('This is a warning message')
    log.error('This is an error message')
    log.critical('This is a critical message')


