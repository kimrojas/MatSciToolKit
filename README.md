# MatSciToolKit
Set of utilities for processing input and output data for Mat Sci softwares


## Installation

1. Create an environment
```
conda install --channel conda-forge --name mtkenv python=3.9
conda activate mtkenv
```

2. Install the MatSciToolKit package
```
pip install git+https://github.com/kimrojas/MatSciToolKit.git
```

3. Install ASE package with specific version
```
pip install git+https://gitlab.com/ase/ase.git@f1b37b76dda641bcdd7dc3f41a5aa243659f4a99
```


## Test

1. Have a required Quantum Espresso executable
```
which pw.x
```

2. Set the environment variable
```
export OMP_NUM_THREADS=1
```


3. Run the test
```
pytest -v --parallel=4
```


