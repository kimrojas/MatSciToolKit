from pprint import pprint

# Accepts string and writes to logfile
def print_log(string: str, logfile: str = "MatciSciToolkit.log"):
    print(string)
    with open(logfile, "a") as f:
        f.write(string + "\n")


