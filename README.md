```
________     ___________      __________________                                     
____  _/________  __/_(_)________(_)_⚡ ____/_/ /________ ______________ _____________
 __  /__  __ \_  /_ _  /_  __ \_  /_/ /   __/ __ \  ՛__ `/_/ ___/_  __ `/  _ \_  ___/
__/ / _  / / /  __/_  /_  / / /  / / /___ _/ / / // /_/ /_/ / ___  /_/ //  __/(__  ) 
/___/ /_/ /_//_/   /_/ /_/ /_//_/  \____/ /_/ /_/ \__,_/ /_/   __\__, / \___//____/  
                                                                /____/               

A program to generate partial charges for periodic systems.
Copyright (C) 2015  Andrea Gabrieli and Marco Sant
```
################################################################################
 
    This software is distributed under the GNU General Public License.
 
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 
    InfiniCharges, including its sources and pointers to the authors
    can be found at http://www.physchem.uniss.it/cpc/

    Contact Address:
    agabrieli@uniss.it
    msant@uniss.it

################################################################################

    InfiniCharges is a computer program for generating reliable
    partial charges for molecular simulations in periodic systems.
    It relies on the DM-REPEAT method where the stability of the
    resulting charges, over a large set of fitting regions, is obtained
    through the simultaneous fit of multiple electrostatic potential (ESP)
    configurations together with the total dipole fluctuations (TDF).

    This program performs the following kinds of fits:
    M-REPEAT (also standard REPEAT)
    DM-REPEAT (also D-REPEAT)
    PARABOLIC RESTRAINED M-REPEAT
    "RAPPE-GODDARD LIKE" RESTRAINED M-REPEAT

################################################################################

The InfiniCharges distribution includes the following files and directories:

- README.md, this file
- LICENSE, the GNU General Public License (GPL)
- InfiniCharges.py, main program
- settings.py, input processing tools
- dir_atautils, contains the fortran library
- rappe_goddard.py, values needed for the rappe-gooddard like restraints
- examples, simple test problems

For a full reference guide look at:
UserManual.pdf

################################################################################
Installation

1) Prerequisites, the following programs should be installed:
- GNU make
- Python 2.7 (already available in all recent Linux distributions)
- Numpy (and F2Py if not automatically installed with numpy)
- a recent Fortran 90 compiler (e.g., gfortran, ifort)

For example in Ubuntu Linux:
> sudo apt-get install make gfortran python-numpy

2) Build the atautils.so library:
> cd dir_atautils
> make

################################################################################
Running the Examples

- enter the wanted directory (“esp”, for example)
> cd examples/esp

- and run the program
> nohup ../../InfiniCharges.py > output &

- check the progress
> tail -f output

################################################################################
Running the Program

- create a working directory and enter it
> mkdir run1
> cd run1

- copy the input file InfiniInput.py 
> cp ../examples/sample_InfiniInput.py InfiniInput.py 

- modify the input as needed
- create a directory named “cubes” and copy here the cube files
> mkdir cubes
> cp path_to_your_cubes/cubefile.cube cubes

- alternatively, to avoid wasting space, you can create a link to the cube files
> ln -s path_to_your_cubes/cubefile.cube cubes

- only in case of “ESP+TDF” fit, create a directory named “tdf” and put here
  the TDF signal file (either a “.tdf” or a CP2K output file) and the trajectory
  in .xyz format
> mkdir tdf
> cp path_to_your_tdf/signal.tdf tdf
> cp path_to_your_tdf/trajectory.xyz tdf

- now it is possible to run
> nohup ../InfiniCharges.py > output &

- check the progress
> tail -f output  


