
# clear; py -m pytest --runmatlab --tb=short

import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--runmatlab", action="store_true", default=False, help="run matlab tests"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runmatlab"):
        # --runmatlab given in cli: do not skip slow tests
        return
    skip_matlab = pytest.mark.skip(reason="need --runmatlab option to run")
    for item in items:
        if "matlab" in item.keywords:
            item.add_marker(skip_matlab)


def pytest_configure(config):
    """
    Allows plugins and conftest files to perform initial configuration.
    This hook is called for every plugin and initial conftest
    file after command line options have been parsed.
    """
    
    config.addinivalue_line("markers", "matlab: mark test as matlab to run")

    if config.getoption("--runmatlab"):
        print("\nImporting MATLAB engine into environment... ", end="")
        import matlab.engine
        print("done", end="")
        print("\nStarting MATLAB instance... ", end="")
        pytest.matlabEng = matlab.engine.start_matlab()
        print("done", end="")
        print("\nChanging MATLAB path... ", end="")
        pytest.matlabEng.cd(r'matlab/unit', nargout=0)
        pytest.matlabEng.addpath('../unitpy_test_helpers', nargout=0)
        print("done", end="")
        print("\n")
    

def pytest_unconfigure(config):
    """
    called before test process is exited. 
    """

    if config.getoption("--runmatlab"): 
        pytest.matlabEng.quit()


""" 
06-05-2025 :
    Got pytest to work with matlab.engine yada yada
    Will start implementing tests for the UNIT start*.m functions next
"""