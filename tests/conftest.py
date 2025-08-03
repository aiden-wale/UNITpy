
# py -m pytest -v --tb=short --runmatlab

import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--runmatlab", action="store_true", default=False, help="run matlab tests"
    )
#endfunction


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runmatlab"):
        # --runmatlab given in cli: do not skip matlab tests
        return
    #endif

    skip_matlab = pytest.mark.skip(reason="need --runmatlab option to run")
    for item in items:
        if "matlab" in item.keywords:
            item.add_marker(skip_matlab)
        #endif
    #endfor
#endfunction


def pytest_configure(config):    
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
    #endif
#endfunction
    

def pytest_unconfigure(config):
    if config.getoption("--runmatlab"): 
        pytest.matlabEng.quit()
    #endif
#endfunction

