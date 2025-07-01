To install the package, navigate to the UNITpy directory and do
```
py -m pip install .
```
or
```
py -m pip install -e .
```
to install the package in an editable state.

To run tests do
```
py -m pytest -v --tb=short
```
or
```
py -m pytest --runmatlab -v --tb=short
```
to also run matlab tests (which will be slow on first run after fresh boot).

if you wish to generate documentation for UNIT (MATLAB) you will need m2html:
https://github.com/gllmflndn/m2html
Then navigate to the `/UNITpy/` parent directory and do
```
py -m docgen_m2html /matlab/unit /path/to/m2html
```
in a terminal or
```matlab
m2html('mfiles', "matlab/unit", 'htmldir', "doc");
```
in MATLAB

