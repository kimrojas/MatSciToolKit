import logging

def logger(name="log"):
    log = logging.getLogger(name)
    log.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler('example.log'),  # Output to file
                            logging.StreamHandler()  # Output to terminal
                        ])
    
    return log
