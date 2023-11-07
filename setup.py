from setuptools import setup, find_packages

setup(
    name="matscitoolkit",
    python_requires='>3.9',
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy==1.22",
        "matplotlib",
        "pandas",
        "tqdm",
        "scipy",
        "pytest",
        "multiprocess",
    ]
)