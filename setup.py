from setuptools import setup, find_packages

setup(
    name="matscitoolkit",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "pandas",
        "tqdm",
        "scipy",
        "agox",
        "pytest"
    ]
)