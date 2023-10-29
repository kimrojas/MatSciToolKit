import numpy as np

def get_high_sympoints(filepath):
    """Gets the high symmetry points from a bands.x output file"""
    with open(filepath, "r") as f:
        lines = f.read().splitlines()

    data_lines = [line for line in lines if "high-symmetry point" in line]

    data_xcoord = []
    for line in data_lines:
        x = line.split()[-1]
        data_xcoord.append(x)

    data = np.array(data_xcoord, dtype=float)

    return data

def get_fermi(filepath):
    """Gets the Fermi energy from a pwscf output file"""
    efermi = None
    
    with open(filepath, "r") as f:
        lines = f.read().splitlines()
    
    for line in lines:
        if "Fermi energy" in line:
            efermi = float(line.split()[-2])
            print(f"Efermi found: {efermi}")
            break
        
    return efermi
        