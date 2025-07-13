import sys
#import subprocess

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

    docDirName = "doc"
    
    print("\nm2html: Generating documentation for "+path_unit+"...\n ", end="")
    matlabEng.workspace['path_unit']    = path_unit
    matlabEng.workspace['path_doc']     = docDirName
    matlabEng.eval("m2html('mfiles', path_unit, 'htmldir', path_doc, 'graph', 'on');", nargout=0)
    print("\nm2html: Finished generating documentation", end="")
    # print("\ngraphviz: Creating graph.pdf from dot file", end="")
    # print("\n")
    # bashcmd = "dot -Tpdf "+ docDirName + "/matlab/unit/graph.dot -o " + docDirName + "/matlab/unit/graph.pdf"
    # process = subprocess.Popen(bashcmd.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()
    # if error:
    #     raise Exception("ERROR OCCURED DURING GRAPH GENERATION")
    # print(output.decode())
    print("\n")

    matlabEng.quit()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\nMust be two arguments, example use: \n\tpy -m docgen_m2html /path/to/unit /path/to/m2html\n\n")
        raise Exception("WRONG NUMBER OF ARGUMENTS")
    path_m2html = str(sys.argv[1])
    path_unit   = str(sys.argv[2])
    docgen_m2html_run(path_unit, path_m2html)

