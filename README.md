# CrystalLatticeto3DPrint
This repository is meant to be an open source and approachable way to 3D print a model of any crystal lattice. The project is designed to be a resource for teachers, hobbyists, and nerds to print and then investigate a 3D structure in a uniquely enriching way.

Spefically, this software takes a CIF (Crystallographic Information File) file and generate a STL file that can be printed and then assembled.

## Step 1: Download CIF Online
There are many online resources where CIF files can be downloaded.

[MaterialsProject](https://materialsproject.org/) (free after registration) (recommended)

[Crystallography Open Database](http://crystallography.net/cod/) (free, no registration)

[Inorganic Crystal Structure Database](https://icsd.fiz-karlsruhe.de/) (not free, but universities often pay for it) (recommended if university affiliated)


## Step 2: Download OpenSCAD
The free [OpenSCAD software](https://openscad.org/) can be controlled with python to generate 3D printable files. After downloading, install the OpenSCAD. For windows, the program will automatically look to find OpenSCAD at the path "C:/Program Files/OpenSCAD/openscad.exe", which is the default installation path. If the program cannot find OpenSCAD, it will prompt the user to find the OpenSCAD.exe file.

## Step 3: Download LatticetoSTL.zip
Download LatticetoSTL.zip from this repository; this version is currently only for windows. After Download extract it into a new folder (Note, all files in the folder must remain there for the .exe to work)

## Step 4: Run LatticetoSTL.exe
The software is now ready to go! Here will give a brief tutorial on using the software as well as clarifying details about the functionality of the program, as well as some general tips and tricks.

### Generating a Lattice

Use the Browse button to find the CIF file that you have downloaded, then generate the lattice with the generate button.

Note, the log is at the bottom of the window and is used to keep track of what has happened and is happening within the software.

After a lattice has been generated, a model of the lattice will appear on the right of the window. This can be rotated for 360 viewing. There is also a show bonds button to display the calculated bonds (see Bond Range for more details). The lattice can be regenerated with the regnerate button to reflect parameter changes.

There is also a Separate atom types and bonds button. When toggled on, a separate STL file(s) will be generated for each atom in the lattice as well as an additional file(s) for the bonds. This allows for simple color changing of the bonds and atom types while printing.

### Adjsutable Parameters
nx, ny, and nz: These adjust number of a, b, and c lattice parameters long that will be printed: i.e. nx*ny*nz = number of unit cells.

Curve Resolution: Number of fragments that are used when rendering final STL shapes. This is effectively the "smoothness of the spheres". I would not reccomend going much below 10, but increasing too high creates extremely large render times. (Note: this has no effect on the display resolution. This only changes rendering/STL smoothness.)

Buildplate Length/Width: This gives the max dimensions of the eventual STL file(s). Note, the default unit is mm.

Radius: This is the maximum radius of the printed nuclei. (For lattices composed of multiple sized atoms, the largest atom radius will be set to radius, and the other atoms will be proportioned appropiately.) The default unit is mm. 

Nucleus Scale: This is a scaling factor of the nucleus. As is common with many ball and stick type atomic models, the nuclear radii can be displayed smaller than reality for clarity. This adjustment allows the scaling of the nucleus from 0 to 1 (0% and 100% of the Radius). This reduces the printed size of the radius of each atoms for printing: a Radius of 30 and Nucleus Scale of 0.5 prints with a radius of 15.

### Outputting Files
Use the Browse button next to the Output Path bar to select where the output files will be places. Then export will begin the STL creation procedure. The program will first make an OpenSCAD file(s). Because traditional filiment 3D printing struggles to produce clean spheres, each atom is cut into two hemispheres. The hemispheres are male and female pairs eventually fit together with a friction fit and no additional pieces. The bonds are simply cylinders that fit into cut outs on the hemispheres. An algotithm uses simulated annealing in order to cut the hemispsheres such that the bonds are as far from the hemisphere split as possible. These bonds are laid out after the hemispheres. For many atoms or large radius prints, multiple OpenSCAD files will be generated. In order to keep track of which atom is which within a the lattice, a txt file(s) containing atom data (AtomInfo#.txt) is exported as well. A new txt file will be generated for each OpenSCAD file containing atoms. Similarly BondInfo#.txt will be generated cooresponding to each OpenSCAD files containing bonds.

After the OpenSCAD files have been generated, OpenSCAD is used to convert each file to an STL. This process can be slow, but time is variable and strongly dependent on number of atoms and Curve Resolution. A very rudimentary time estimate is provided in the log to give a guess of how long the conversion will take, but expect time results to only be good to about an order of magnitude.

While exporting, the program will appear frozen and will not update until exporting is completed. The log will periodically update as checkpoints are reached.

### Advanced Mode
Advanced mode displays many extra adjustable parameters for those with more unsatified with the defaults. Although, the defaults are designed to be well-suited for most cases.

Bond Range: This is a scaling factor that adjusts which atoms are connected by a bond. In the program, bonds are found by first finding the distance to all the other atoms in the lattice. A bond is always created to the nearest neighbor to each atom. Then, any other atoms within a radius of d_nn * (1 + Bond Range) will also be included in the bonds. This can be useful for complex lattices where one could want to display different bond lengths between several atoms.

Tolerance: Tolerance only effects the 3D printing part of the program. In order to get a nice friction fit, it is useful to be able to make the radius of the male parts slightly smaller than the radius of the female parts. This parameter reduces diameter and length of male parts by Tolerance for good fitting.

Cut Depth Ratio: This controls the depth of each bond hole in the hemisphere. Cut Depth Ratio is the ratio of the cut depth to Radius. If too large, this can interfere with other bonds or female hemisphere interior.

Cut Radius Ratio: This controls the radius of each bond hole in the hemisphere. Cut Depth Radius is the ratio of the cut radius to Radius. If too large, this can interfere with other bonds or female hemisphere interior. It is important to consider tolerance especially with this dimension, as the tolerance is extremely important to the fitting of the bonds within the hemisphere cut holes.

Trapezoid Side Ratio: This and the next two parameters specify the geometry of the trapezoid that connects the two hemispheres. The sphere of the nucleus is split into two hemispheres: one male and one female. The male hemisphere has an isosceles prism extruding from it bottom, while the female has a cavity that allows the male piece to insert. An isosceles trapeziod is chosen because it breaks rotational symmetry whereas a shape like a sqaure would allow the bonds to be rotated four different ways. (Nearly right angles are the easiest to remove supports from as well.) Specifically, Trapezoid Side Ratio is the ratio of the largest side of the trapeziod (base) to the radius of each atom.

Trapeziod Depth Ratio: Trapezoid Depth Ratio is the ratio of the the trapeziodal prismic cut/extrusion depth to the radius of each atom.

Trapezoid Oblong Ratio: Trapezoid Side Ratio is the ratio of the smallest size of the trapeziod (top) to the the largest size of the trapeziod (base). 1 will make a square

Overhang Angle: To avoid generating supports in the female part, the trapeziodal prismic cut is topped with a pyramid. The base of this pyramid is the trapezoid and the other faces are determined by Overhang Angle. This controls the angle between the base of the pyramid and each upright face. Most extrusion 3D printers can comfortably print features with overhang angles of 45 degrees without supports.

## General Tips and Tricks
It is generally a good idea to try a test atom and bond. This can be useful to adjust tolerances and sizing. One common problem is [elephant footing](https://help.prusa3d.com/article/elephant-foot-compensation_114487). This is a common 3D printing issue caused by the nozzle's Z position being too close to the bed. This can cause the male part of the hemisphere or any vertically printed bonds to be widened at the bottom. However, this can be useful to provide nice tight, reinsertable bonds or hemispheres. Testing is useful to get a feeling for how this will work before you print an entire lattice.

For cleanest nucleui, it is useful to separate the atoms and bonds. Then, the hemispheres can be printed with no brim, but with supports. When printing vertical bonds, it is useful to use a brim, especially for longer bonds.
