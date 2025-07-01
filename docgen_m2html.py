import sys

def docgen_m2html_run(path_m2html, path_unit):
    print("\nImporting MATLAB engine into environment... ", end="")
    import matlab.engine
    print("done", end="")
    print("\nStarting MATLAB instance... ", end="")
    matlabEng = matlab.engine.start_matlab()
    print("done", end="")
    print("\nUpdating MATLAB path... ", end="")
    matlabEng.addpath(path_m2html, nargout=0)
    matlabEng.addpath(path_unit, nargout=0)
    print("done", end="")
    print("\n")
    
    print("\nm2html: Generating documentation for "+path_unit+"... ", end="")
    matlabEng.workspace['path_unit']    = path_unit
    matlabEng.workspace['path_doc']     = "doc"
    matlabEng.eval("m2html('mfiles', path_unit, 'htmldir', path_doc);", nargout=0)
    print("done", end="")
    print("\n")

    matlabEng.quit()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\nMust be two arguments, example use: \n\tpy -m docgen_m2html /path/to/unit /path/to/m2html\n\n")
        raise Exception("WRONG NUMBER OF ARGUMENTS")
    path_m2html = str(sys.argv[1])
    path_unit   = str(sys.argv[2])
    docgen_m2html_run(path_unit, path_m2html)

