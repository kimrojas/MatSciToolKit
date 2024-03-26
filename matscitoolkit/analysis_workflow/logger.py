import logging


def get_log(name="log", logfile=None):
    log = logging.getLogger(name)

    # Configure a formatter
    formatter = logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s", "%Y-%m-%d %H:%M:%S")

    # Create a file handler and configure it
    if logfile is not None:
        file_handler = logging.FileHandler(logfile)
        file_handler.setLevel(logging.DEBUG)  # Set the logging level
        file_handler.setFormatter(formatter)  # Set the formatter
        log.addHandler(file_handler)  # Add handler to the logger

    # Create a stream handler and configure it
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)  # Set the logging level
    stream_handler.setFormatter(formatter)  # Set the formatter
    log.addHandler(stream_handler)  # Add handler to the logger

    # Set the logging level for the logger
    log.setLevel(logging.DEBUG)

    return log
