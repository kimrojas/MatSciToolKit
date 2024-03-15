from setuptools import setup, find_packages
import re

# Version Number:
version_file = "matscitoolkit/__version__.py"
with open(version_file) as f:
    lines = f.readlines()

for line in lines:
    if "__version_info__" in line:
        result = re.findall("\d+", line)
        result = [int(x) for x in result]
        version = "{}.{}.{}".format(*result)
        break


setup(
    name="matscitoolkit",
    python_requires=">3.9",
    version=version,
    packages=find_packages(),
    install_requires=[
        "numpy==1.22",
        "matplotlib",
        "pandas",
        "tqdm",
        "scipy",
        "pytest",
        "multiprocess",
    ],
)
