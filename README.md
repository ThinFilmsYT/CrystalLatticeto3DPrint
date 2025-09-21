# CrystalLatticeto3DPrint
This repository is meant to be an open source and approachable way to 3D print a model of any crystal lattice. The project is designed to be a resource for teachers, hobbyists, and nerds to print and then investigate a 3D structure in a uniquely enriching way.

Spefically, this software takes a CIF (Crystallographic Information File) file and generate a STL file that can be printed and then assembled.

## Step 1: Download CIF Online
There are many online resources where CIF files can be downloaded.

[MaterialsProject](https://materialsproject.org/) (free after registration) (recommended)

[Crystallography Open Database](http://crystallography.net/cod/) (free, no registration)

[Inorganic Crystal Structure Database](https://icsd.fiz-karlsruhe.de/) (not free, but universities often pay for it) (recommended if university affiliated)


## Step 2: Download OpenSCAD
The free [OpenSCAD software](https://openscad.org/) can be controlled with python to generate 3D printable files. After downloading, install the OpenSCAD. For windows, the program will automatically look to find OpenSCAD at the path "C:/Program Files/OpenSCAD/openscad.exe", which is the default installation path. (Note, this should be improved; I am thinking to do a check for OpenSCAD at that path and then have a pop up come up if it can't find it. Maybe add a file that contains OpenSCAD path for later use.)

## Step 3: Download LatticetoSTL.zip
