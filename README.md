# cif2json

Make an exess input file (json) from a cif file. A sphere of molecules are created with a given distance from the central fragment. If ions, the molecules can be paired together into [cation+anion] neutral fragments. This has been adapted from code by Anh LP Nguyen.


## Dependences
- qcelemental
- numpy
- PyCifRW
- QCP (https://github.com/zoeseeger/qcp-python-app) and change sys.path.append(path_to_qcp)

## Usage

For a sphere of 100Å where the cif file is in the same directory:
`cif2json.py 100`

or if no pairing of molecules in desired:

`cif2json.py 100 file.cif none`

Molecules are formed using QCP. To add molecules to QCP edit ~/myMolecules.qcp which will be created after first use.