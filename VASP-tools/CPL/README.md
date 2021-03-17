# Chemical Potential Landscape Plotter

- An object oriented program for computing the chemical potential landscape (CPL) of a target compound
- See Journal of Applied Physics 117, 035702 (2015); https://doi.org/10.1063/1.4906065 for details

## Setup

- CPL_OOP.py for the current implementation
- General plotting parameters are given in the code and may be modified
- A test is provided for CuBiW2O8
  - J. Mater. Chem. A,2021, 9,1643–1654; https://dx.doi.org/10.1039/d0ta07653h
  - This test shows an example of plotting a quarternary metal oxide by slicing its 3D CPL at a value represented by a free-parameter (Cu in the example setup)
  - See test files for input format

## Defect Formation Analysis

- The original implementation (CPL.py) allows for plotting of defect formation energies at different chemical potential conditions. The current OOP implementation of this functionality is in progress.

## Note

- I am not the author of widgets.py, though I did modify it slightly to work for the purpose of CPL plotting

## Questions

- All questions may be directed to edan.bainglass@gmail.com
