# cif2json

Make an exess input file (json) from a cif file. A sphere of molecules are created with a given distance from the central fragment. If ions, the molecules can be paired together into [cation+anion] neutral fragments. This has been adapted from code by Anh LP Nguyen.

## Dependences
- qcelemental
- numpy
- qcelemental
- QCP (https://github.com/zoeseeger/qcp-python-app) and change sys.path.append(path_to_qcp)

## Usage

`cif2json.py inputfile.inp`

### Input file
```
// [C1mpyr][NTf2]
Title        : c1mpyr ntf2
Cif          : ../P11NTF2.cif
Ionic        : Y
Cation       : C 6 H 14 N 1
Anion        : C 2 F 6 N 1 O 4 S 2
Rgiven       : Y
Rsphere      : 22.0
Rdim         : 0.0
Rtrim        : 0.0
```
