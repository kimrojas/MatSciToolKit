import logging


class logger:
    def __init__(self, name="log", logfile=None, debug=True):
        self.log = logging.getLogger(name)

        # Configure a formatter
        formatter = logging.Formatter("%(asctime)s [%(levelname)s]-- %(message)s", "%Y-%m-%d %H:%M:%S")

        # Create a file handler and configure it
        if logfile is not None:
            file_handler = logging.FileHandler(logfile)
            file_handler.setLevel(logging.DEBUG)  # Set the logging level
            file_handler.setFormatter(formatter)  # Set the formatter
            self.log.addHandler(file_handler)  # Add handler to the logger

        # Create a stream handler and configure it
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.DEBUG)  # Set the logging level
        stream_handler.setFormatter(formatter)  # Set the formatter
        self.log.addHandler(stream_handler)  # Add handler to the logger

        # Set the logging level for the logger
        if debug:
            self.log.setLevel(logging.DEBUG)
        else:
            self.log.setLevel(logging.WARNING)

    def debug(self, *message):
        self.log.debug(message)

    def info(self, *message):
        self.log.info(*message)

    def warning(self, *message):
        self.log.warning(*message)

    def error(self, *message):
        self.log.error(*message)

    def critical(self, *message):
        self.log.critical(*message)

    def assert_(self, condition, message):
        try:
            assert condition, message
        except AssertionError:
            self.log.error(message)
            raise AssertionError
        
    def close(self):
        for handler in self.log.handlers[:]:
            self.log.removeHandler(handler)
