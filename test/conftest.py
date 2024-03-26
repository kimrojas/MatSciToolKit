import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--parallel", type=int, action="store", default=1, help="max parallelization (default: 1)"
    )