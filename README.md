# cif2json

Make an exess input file (json) from a cif file. A sphere of molecules are created with a given distance from the central fragment. If ions, the molecules can be paired together into [cation+anion] neutral fragments. This has been adapted from code by Anh LP Nguyen.

## Dependences
- qcelemental
- numpy
- qcelemental
- QCP (https://github.com/zoeseeger/qcp-python-app) and change sys.path.append(path_to_qcp)
