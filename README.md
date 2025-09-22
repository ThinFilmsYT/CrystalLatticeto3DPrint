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
Download LatticetoSTL.zip from this repository; this version is currently only for windows. After Download extract it into a new folder (Note, all files in the folder must remain there for the .exe to work)

## Step 4: Run LatticetoSTL.exe
The software is now ready to go! Here will give a brief tutorial on using the software as well as clarifying details about the functionality of the program, as well as some general tips and tricks.

### Generating a Lattice

Use the Browse button to find the CIF file that you have downloaded, then generate the lattice with the generate button.

Note, the log is at the bottom of the window and is used to keep track of what has happened and is happening within the software.

After a lattice has been generated, a model of the lattice will appear on the right of the window. This can be rotated for 360 viewing. There is also a show bonds button to display the calculated bonds (see Bond Range for more details). The lattice can be regenerated with the regnerate button to reflect parameter changes.

### Adjsutable Parameters
nx, ny, and nz: These adjust number of a, b, and c lattice parameters long that will be printed: i.e. nx*ny*nz = number of unit cells.

Curve Resolution: Number of fragments that are used when rendering final STL shapes. This is effectively the "smoothness of the spheres". I would not reccomend going much below 10, but increasing too high creates extremely large render times. (Note: this has no effect on the display resolution. This only changes rendering/STL smoothness.)

Buildplate Length/Width: This gives the max dimensions of the eventual STL file(s). Note, the default unit is mm.

Radius: This is the maximum radius of the printed nuclei. (Note, for lattices composed of multiple sized atoms, the largest atom radius will be set to radius, and the other atoms will be proportioned appropiately. Note, the default unit is mm.

Nucleus Scale: This is a scaling factor of the nucleus. As is common with many ball and stick type atomic models, the nuclear radii can be displayed smaller than reality for clarity. This adjustment allows the scaling of the nucleus from 0 to 1 (0% and 100% of the Radius).

### Outputting Files
Use the Browse button next to the Output Path bar to select where the output files will be places. Then export will begin the STL creation procedure. The program will first make one or more OpenSCAD files. Because traditinoal filiment 3D printing struggles to produce clean spheres, each atom is cut into two hemispheres. The hemispheres are male and female pairs eventually fit together with a friction fit and no additional pieces. The bonds are simply cylinders that fit into cut outs on the hemispheres. An algotithm uses simulated annealing in order to cut the hemispsheres such that the bonds are as far from the hemisphere split as possible. These bonds are laid out after the hemispheres. For many atoms or large radius prints, multiple OpenSCAD files will be generated. In order to keep track of which atom is which within a the lattice, a txt file(s) containing atom data (AtomInfo#.txt) is exported as well. A new txt file will be generated for each OpenSCAD file containing atoms. Similarly BondInfo#.txt will be generated cooresponding to each OpenSCAD files containing bonds.

After the OpenSCAD files have been generated, OpenSCAD is used to convert each file to an STL. This process can be slow, but time is variable and strongly dependent on number of atoms and Curve Resolution. A very rudimentary time estimate is provided in the log to give a guess of how long the conversion will take, but expect time results to only be good to about an order of magnitude.

While exporting, the program will appear frozen and will not update until exporting is completed.

### Advanced Mode
Advanced mode displays many extra adjustable parameters for those with more unsatified with the defaults. Although, the defaults are desogned to be well-suited for most cases.

