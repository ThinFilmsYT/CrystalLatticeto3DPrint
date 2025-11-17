# CrystalLatticeto3DPrint
This repository is meant to be an open source and approachable way to 3D print a model of any crystal lattice. The project is designed to be a resource for teachers, hobbyists, and nerds to print and then investigate a 3D structure in a uniquely enriching way.

Specifically, this repo contains a free, downloadable software that takes a CIF (Crystallographic Information File) file and generates a STL file that can be 3D printed and then assembled.

## Step 1: Download CIF Online
There are many online resources where CIF files can be downloaded.

[MaterialsProject](https://materialsproject.org/) (free after registration) (recommended)

[Crystallography Open Database](http://crystallography.net/cod/) (free, no registration)

[Inorganic Crystal Structure Database](https://icsd.fiz-karlsruhe.de/) (not free, but universities often pay for it) (recommended if university affiliated)


## Step 2: Download OpenSCAD
The free [OpenSCAD software](https://openscad.org/) can be controlled with python to generate 3D printable files. After downloading, install the OpenSCAD. For windows, the program will automatically try to find OpenSCAD at the path "C:/Program Files/OpenSCAD/openscad.exe", which is the default installation path. For mac, the program will look for OpenSCAD at /Applications/ folder. If the program cannot find OpenSCAD, it will prompt the user to find the .exe or .app file.

## Step 3: Download LatticetoSTL
For windows, download LatticetoSTL.zip from this repository. After downloading, extract it into a new folder (Note, all files in the folder must remain in the same directory as the .exe for the program to work)

For mac, download the .dmg file, and drag the LatticetoSTL into the Applications folder.

## Step 4: Run LatticetoSTL
The software is now ready to go! Here will give a brief tutorial on using the software as well as clarifying details about the functionality of the program, as well as some general tips and tricks. 

Important note, the software can be slow on start up. It needs to open some unfortunately bulky python packages but after opening should be much faster.

### Generating a Lattice

After opening the software, use the Browse button to find the CIF file that you have downloaded, then generate the lattice with the generate button.

Note, the log is at the bottom of the window and is used to keep track of what has happened and is happening within the software.

After a lattice has been generated, a model of the lattice will appear on the right of the window. This can be rotated for full 3D viewing. There is also a "show bonds" button to display the bonds that the program will print (see Bond Range for more details about bond adjustments). The lattice can be regenerated with the regenerate button to reflect parameter changes or inputting a different CIF..

There is also a "Separate atom types and bonds" button. When toggled on, a separate STL file(s) will be generated for each atom in the lattice as well as an additional file(s) for the bonds. This allows for easy color changing of the bonds and atom types for eventual 3D printing.

### Adjustable Parameters
nx, ny, and nz: These adjust number of a, b, and c lattice parameters long that will be printed: i.e. nx*ny*nz = number of unit cells.

Curve Resolution: Number of fragments that are used when rendering final STL shapes. This is effectively the "smoothness of the spheres". I would not recommend going much below 10, but increasing too high creates extremely large render times. (Note: this has no effect on the display resolution. This only changes rendering/STL smoothness.)

Buildplate Length/Width: This gives the max dimensions of the eventual STL file(s). The default unit is mm.

Radius: This is the maximum radius of the printed nuclei. (For lattices composed of multiple sized atoms, the largest atom radius will be set to Radius, and the other atoms will be proportioned appropriately.) This is the only parameter that determines the size of the final printed lattice, but is not the only parameter that determines an atom's printed size (see Nucleus Scale). The default unit is mm. 

Nucleus Scale: This is a scaling factor of the nucleus. As is common with many ball and stick type atomic models, the nuclear radii can be displayed smaller than reality for clarity. This adjustment allows the scaling of the nucleus from 0 to 1 (0% and 100% of Radius). This reduces the printed size of each atom for printing: a Radius of 30 and Nucleus Scale of 0.5 prints with a radius of 15. 

Bond Cross Section: This is a togglable selection between Circle and D. Circle prints a vertical cylinder that is hollow. The ratio of the ID to the OD is 0.7. Circle is useful if rotating the atoms about their bonds or if simplicity is desired. Oppositely, if no bond rotation is desired, D ensures correct relative positions of all atoms/bonds. The D shafts are printed horizontally with the flat D face down. The depth of the D cut is controlled by Overhang Angle. It is important to note that the atoms are cut with the D shape as well; all of the flat faces of the D cuts for each atom points toward a vector that is normal to the crystal lattice x-y plane and that passes through the center of the atom (e.g. for orthogonal crystals, the flats will point towards the z axis that passes through the center of the atom it is on). If the bond points exactly normal to the x-y plane, the flat face will face towards the x axis passing through the atom instead. 

### Outputting Files
Use the Browse button next to the Output Path bar to select where the output files will be placed. Then "export" will begin the STL creation procedure. The program will first make an OpenSCAD file(s). Because traditional filament 3D printing struggles to produce clean spheres, each atom is cut into two hemispheres. The hemispheres are male and female pairs that eventually fit together with a friction fit, making a clean sphere with no additional pieces. The bonds are rods that fit into cuts on the hemispheres. To find a plane to split the sphere into hemispheres, simulated annealing is used. The plane is optimized such that the bonds are as far from the cutting plane as possible. 

During printing, hemispheres are laid out in pairs first and then the bonds second. For many atoms or large radius prints, multiple OpenSCAD files will be generated. In order to keep track of which atom is which within a the lattice, a txt file(s) containing atom data (AtomInfo#.txt) is created for reference. A new txt file will be generated for each OpenSCAD file containing atoms. Similarly BondInfo#.txt will be generated corresponding to each OpenSCAD file containing bonds.

After the OpenSCAD files have been generated, OpenSCAD is used to convert each file to an STL. This process can be slow, but time is variable and strongly dependent on the number of atoms, bonds and Curve Resolution. A very rudimentary time estimate is provided in the log to give a guess of how long the conversion will take, but expect time results to only be good to about an order of magnitude.

While exporting, the program will appear frozen and will not update until exporting is completed. The log will periodically update as checkpoints are reached.

### Advanced Mode
Advanced mode displays many extra adjustable parameters for those who are unsatisfied with the defaults. Although, the defaults are designed to be well-suited for most cases.

Bond Range: This is a scaling factor that adjusts which atoms are connected by a bond. In the program, bonds are found by first finding the distance to all the other atoms in the lattice. A bond is always created to the nearest neighbor to each atom. Then, any other atoms within a radius of d_nn * (1 + Bond Range) will also be included in the bonds. This can be useful for complex lattices where one could want to display different bond lengths between several types of atoms.

Tolerance: To get a nice friction fit, it is useful to be able to make the dimensions of the male parts slightly smaller than the dimensions of the female parts. This parameter reduces diameter of the bonds and height and width of the of male trapezoidal extrusion by length, Tolerance, for good friction fits. Tolerance only affects the 3D printing part of the program and should be adjusted for your printer.The default unit is mm. 

Cut Depth Ratio: This controls the depth of each bond hole in the hemisphere and the corresponding male length. Cut Depth Ratio is the ratio of the cut depth to the size of the atom (for larger atomic radii, the depth will be proportionally deeper). If too large, this can interfere with other bonds or the female hemisphere’s interior.

Cut Radius Ratio: This controls the radius of each bond hole in the hemisphere. Cut Depth Radius is the ratio of the cut radius to Radius. If too large, this can interfere with other bonds or the female hemisphere’s interior. It is important to consider Tolerance especially with this dimension, as the tolerance is extremely important to the fitting of the bonds within the hemisphere cut holes.

Trapezoid Side Ratio: This and the next two parameters specify the geometry of the trapezoid that connects the two hemispheres. The sphere of the nucleus is split into two hemispheres: one male and one female. The male hemisphere has an isosceles trapezoidal prism extruding from its bottom, while the female has a cavity that allows the male piece to insert. An isosceles trapezoid is chosen because it breaks rotational symmetry whereas a shape like a square would allow the bonds to be rotated four different ways. (Nearly right angles are the easy to remove supports from as well.) Specifically, Trapezoid Side Ratio is the ratio of the largest side of the trapezoid (base) to the radius of each atom.

Trapezoid Depth Ratio: Trapezoid Depth Ratio is the ratio of the trapezoidal prismic cut/extrusion depth to the radius of each atom. Like Cut Depth Ratio, this similarly scales with atom size.

Trapezoid Oblong Ratio: Trapezoid Side Ratio is the ratio of the smallest size of the trapezoid (top) to the largest size of the trapezoid (base). 1 will make a square, 0 will probably make an error lol.

Overhang Angle: To avoid the need for supports inside the female part, the trapezoidal-prism-cut is topped with a pyramid. The base of this pyramid is the trapezoid and the slope of the other faces are determined by Overhang Angle. This controls the angle between the base of the pyramid and each upright face. Most extrusion 3D printers can comfortably print features with overhang angles of 45 degrees without supports. Units are in degrees

## General Tips and Tricks
It is generally a good idea to try a test atom and bond before printing an entire lattice. This can be useful to adjust tolerances and sizing. One common problem is [elephant footing](https://help.prusa3d.com/article/elephant-foot-compensation_114487). This is a common 3D printing issue caused by the nozzle's Z position being too close to the bed during the first layer or two. This can cause the male part of the hemisphere or any vertically printed bonds to be widened at the bottom. However, this widening can be useful to provide nice tight, reinsertable hemispheres (but sadly only affects one side of the bonds). Testing is useful to get a feeling for how this will work before you print an entire lattice.

It is quite easy to print too small. Printed atomic radii of about 15mm is a good place to start. Going much smaller than that can cause problems depending on 3D printer resolution, and going too much bigger takes a lot of time and filament.

For cleanest nuclei with circle bonds, it is useful to separate the atoms and bonds. Then, the hemispheres can be printed with no brim, but with supports. When printing vertical bonds, it is useful to use a brim, especially for the tall vertically printed bonds.

In order to use the provided txt docs to their full potential, it is often useful to label the atoms. I like to number the flat faces of the hemispheres with a permanent marker before removing the supports from the male parts. This is especially useful for the D shaft bond mode. The atoms are numbered with increasing numbers down the column (negative y direction). Once a column is full, a new column is started to the right (x direction). OpenSCAD has built-in axes to verify atom numbering. Bond numbering follows the same scheme.


