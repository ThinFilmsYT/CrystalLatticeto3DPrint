from solid import *
# from solid.utils import *
import subprocess
import os
import importlib.metadata as importlib_metadata  # Modern alternative for metadata
import sys
import numpy as np
import math



# === Configuration ===
# Set this to the OpenSCAD command or full path to the executable on your system.
# If 'openscad' is in your PATH, no change is needed. Otherwise, use the full path:
openscad_cmd = "C:/Program Files/OpenSCAD/openscad.exe"

# Ensure solidpython is installed
try:
    version = importlib_metadata.version("solidpython")
    # print(f"Using solidpython version {version}")
except importlib_metadata.PackageNotFoundError:
    print("solidpython is not installed. Install with 'pip install solidpython'.")
    sys.exit(1)


def hemisphere(r):
    # Cube spans x:[-r,r], y:[-r,r], z:[r,3r]
    cutter = translate([0, 0, r/2])(cube([2*r, 2*r, r], center=True))#starting with corner at 
    return intersection()(sphere(r), cutter) #moves hemisphere center to origin

def cut_cylinder_on_hemisphere(radius, polar, azimuth, cut_depth, cyl_radius):

    # Direction vector from center of sphere
    x = np.sin(polar) * np.cos(azimuth)
    y = np.sin(polar) * np.sin(azimuth)
    z = np.cos(polar)

    # Position cylinder so its base starts just outside the hemisphere and cuts inward
    center_offset = radius - cut_depth / 2
    position = [center_offset * x, center_offset * y, center_offset * z]
    # print(position)
    # Create a long cylinder, then align and position it
    cyl = cylinder(h=cut_depth, r=cyl_radius, center=True)

    # Rotate from Z to target vector using rotate(a, v), where `a` is angle, `v` is axis
    # But OpenSCAD/SolidPython uses Euler angles, so we need to rotate properly
    # First rotate around Y, then Z
    # Z-rotate by azimuth, then Y-rotate by polar angle from Z

    return translate(position)(rotate(a=[0,polar*180/np.pi,azimuth*180/np.pi])(cyl))

def trapozoidal_prism(side, ratio, depth):
    """Return an trapozoid in the XY plane extruded downward in Z direction."""
    trap = polygon([
        [-side, -side],
        [-side, side],
        [side, ratio*side],
        [side, -ratio*side]        
    ])
    return linear_extrude(height=depth)(trap)

def pyramid_from_points(base_pts, apex):
    """Create a pyramid from 4 coplanar base points and 1 apex point using triangulated base."""
    if len(base_pts) != 4:
        raise ValueError("Base must have exactly 4 points")
    
    points = base_pts + [apex]
    faces = [
        [0, 1, 2], [0, 2, 3],  # base as two triangles
        [0, 4, 1],
        [1, 4, 2],
        [2, 4, 3],
        [3, 4, 0]
    ] #very subtle detail here is that the points must be ordered clockwise around the face from outside looking in
    # print(points)
    return polyhedron(points=points, faces=faces,convexity=2)

def make_trap_pyramid(trap_side,trap_depth, trap_oblong_ratio, overhang_angle):
    trapozoidal_base_points = [
        [-trap_side, -trap_side,trap_depth],
        [-trap_side, trap_side,trap_depth],
        [trap_side, trap_oblong_ratio*trap_side,trap_depth],
        [trap_side, -trap_oblong_ratio*trap_side,trap_depth]        
    ]
    trapozoidal_pyramid_height = float(trap_side / np.tan(np.radians(overhang_angle)))
    trapozoidal_pyramid_apex_point = [0,0,trapozoidal_pyramid_height+trap_depth]
    return pyramid_from_points(trapozoidal_base_points,trapozoidal_pyramid_apex_point)

def reorient_angles(normal_v, angles):
    """
    Rotate spherical coordinates so that the polar axis aligns with vector v.
    
    Parameters:
    -----------
    normal_v : tuple or list of 3 floats
        Vector (vx, vy, vz) defining the new polar axis.
    angles : list of tuples
        List of (azimuth, polar) angles in radians.
        
    Returns:
    --------
    list of tuples
        List of (azimuth', polar') angles in the rotated frame.
    """
    normal_v = np.array(normal_v, dtype=float)
    normal_v /= np.linalg.norm(normal_v)  # normalize
    
    # Find rotation axis and angle to align v with z-axis
    z_axis = np.array([0, 0, 1.0])
    axis = np.cross(normal_v, z_axis)
    norm_axis = np.linalg.norm(axis)
    if norm_axis < 1e-12:  # v already along z
        R = np.eye(3)
    else:
        axis /= norm_axis
        angle = np.arccos(np.clip(np.dot(normal_v, z_axis), -1.0, 1.0))
        
        # Rodrigues' rotation formula
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

    list_lt = []  # polar < pi/2
    list_ge = []  # polar >= pi/2

    for bond_leng, pol, az in angles:
        # Convert to Cartesian
        x = np.sin(pol) * np.cos(az)
        y = np.sin(pol) * np.sin(az)
        z = np.cos(pol)
        vec = np.array([x, y, z])
        
        # Apply rotation
        vec_rot = R @ vec
        
        # Back to spherical
        r = np.linalg.norm(vec_rot)
        pol_new = np.arccos(np.clip(vec_rot[2]/r, -1.0, 1.0))
        az_new = np.arctan2(vec_rot[1], vec_rot[0])

        # Classify
        if pol_new < np.pi/2:
            list_lt.append((pol_new, az_new))
        else:
            list_ge.append((pol_new, az_new))
    list_flipped = []
    for pol, az in list_ge:#flip hemisphere upright for better 3d printing
        if az>=np.pi:
            list_flipped.append((np.pi-pol,az-np.pi))
        else:
            list_flipped.append((np.pi-pol,az+np.pi))

    return list_lt, list_flipped

def bond_cylinder(bond_x, bond_y, bond_z,bond_length, cut_radius, tolerance,ID_ratio = 0.8):
    return translate([bond_x, bond_y, bond_z])(cylinder(h=bond_length, r=cut_radius-tolerance/2, center=True) - cylinder(h=bond_length, r=cut_radius*ID_ratio, center=True))

def makeSTL(STL_path,atoms,cut_planes,plate_dim = (210,250),parameters = [],log = None):
    failed = True
    
    if len(parameters) == 10:
        try:
            res = str(int(parameters[0]))
            radius = float(parameters[1])
            nucleus_scale = float(parameters[2])
            tolerance = float(parameters[3])
            cut_depth_ratio = float(parameters[4])
            cut_radius_ratio = float(parameters[5])
            trap_side_ratio = float(parameters[6])
            trap_depth_ratio = float(parameters[7])
            trap_oblong_ratio = float(parameters[8])
            overhang_angle = float(parameters[9])
            failed = False
        except:
            print('Inputted parameters failed to convert to floats')
            if log:
                log.log_message('Inputted parameters failed to convert to floats', error=True)
    if failed:
        # Parameters
        res = '50' #number of fragments for each 3D printed piece. basically 3D printed circle resolution
        radius = 10 #max radius size in mm
        nucleus_scale = 1 # scale from 0 to 1 to decrease nucleus size in order to see structure more clearly
        tolerance = 0.1 #add tolerance to male parts
        cut_depth_ratio = 1/4       # Depth of the cylindrical cut for bonds
        cut_radius_ratio = 1/6       # Radius of the cutting cylinder

        trap_side_ratio = 1/3 #ratio of half long trap leg to radius
        trap_depth_ratio = 1/5 # ratio of trap prism cut depth to radius
        trap_oblong_ratio = 0.9 #ratio of height of short trap leg to long trap leg

        overhang_angle = 45 #degrees
    
    cut_depth = cut_depth_ratio*radius
    cut_radius = cut_radius_ratio*radius

    #add a check to see if bond cuts could intersect interior features
    #Max bond cut depth is just r_cut <= r_sphere * [1 - (trap_depth_ratio + trap_side_ratio*tan(overhang_angle))]
    #this is just the top of the pyramid
    #max bond cut depth should also be r_cut <= r_sphere * [1 - sqrt(trap_depth_ratio^2 + trap_side_ratio^2)]
    #this is the corner of the pyramid base
    max_cut_ratio = max(1 - (trap_depth_ratio + trap_side_ratio*np.tan(overhang_angle)),1 - np.sqrt(trap_depth_ratio**2 + trap_side_ratio**2))
    if cut_depth/radius> max_cut_ratio:
        print('Warning: Bond cuts might intersect interior features.')
        if log:
            log.log_message('Warning: Bond cuts might intersect interior features.', error=False)
        
    elif cut_depth>radius:
        print('Error: Cut depth way too big')
        if log:
            log.log_message('Error: Cut depth way too big', error=True)
        return


        ''' #first pass, bad layout
        #find how many atoms fit on build plate
        max_bonds = 0 #max number of bonds for a single atom
        max_radius = 0 # in units of Anstroms
        for atom in atoms:
            if len(atom.bonds)>max_bonds:
                max_bonds = len(atom.bonds)
            if atom.radius > max_radius:
                max_radius = atom.radius
        width = (max_bonds+1)//2 * 5 * cut_radius + 5 * radius #floor((max_bonds+1)/2) to get max bonds on either side of hemispheres
        row_max = plate_dim[0]//width
        col_max = plate_dim[1]//(radius*2.5)
        
        if (atoms_per_plate:=row_max*col_max)>0:
            plate_total = math.ceil(len(atoms)/atoms_per_plate)
        else:
            print(f'Error: plate too small for radius = {radius},{width}')

        print(f'Writing {plate_total} .stl files.')
        atom_idx = 0 # counter must keep counting regardless of plate number
        for p in range(plate_total): #allow for multiple plates adn thus multiple stls
            p_models = [] #stores all the models for a single plate
            while atom_idx<atoms_per_plate*(p+1) and atom_idx < len(atoms):
                atom = atoms[atom_idx] #get atom
                atom_x = (atom_idx//col_max) * width #these are the x,y coords of the center between the hemispheres
                atom_y = (atom_idx%col_max) * radius*2.5

                #get cut angles
                cut_angles1,cut_angles2 = reorient_angles(cut_planes[atom_idx],atom.bonds)
                atom_print_r = atom.radius/max_radius*radius*nucleus_scale #convert atom radius in anstroms into units of radius (mm)
                cuts1 = [cut_cylinder_on_hemisphere(atom_print_r, pol, az, cut_depth, cut_radius) for pol, az in cut_angles1]
                cuts2 = [cut_cylinder_on_hemisphere(atom_print_r, pol, az, cut_depth, cut_radius) for pol, az in cut_angles2]

                trap_side = atom_print_r * trap_side_ratio
                trap_depth = atom_print_r * trap_depth_ratio

                #make female hemispher
                model1 = translate([radius*1.25 + atom_x, atom_y, -trap_depth])(difference()(hemisphere(radius*atom.radius/max_radius*nucleus_scale), *cuts1) - trapozoidal_prism(trap_side,trap_oblong_ratio,trap_depth) - make_trap_pyramid(trap_side,trap_depth, trap_oblong_ratio, overhang_angle)) #female
                #make male hemisphere
                model2 = translate([-radius*1.25 + atom_x, atom_y, 0])(difference()(hemisphere(radius*atom.radius/max_radius*nucleus_scale), *cuts2) + translate([0, 0, -trap_depth])(trapozoidal_prism(trap_side-tolerance/2,trap_oblong_ratio,trap_depth-tolerance))) #male
                
                # make bond models
                bond_models = []
                for j, bond in enumerate(atom.bonds):
                    bond_length = (bond[0] - (atom.bond_atom_sizes[j][0] + atom.bond_atom_sizes[j][1])*nucleus_scale)*atom_print_r/max_radius + 2*cut_depth
                    if j == 0: #first bond at x = 0
                        bond_x = atom_x
                    elif j%2 == 1: #for odd, put on right side
                        bond_x = 2.5*cut_radius*(j-1)/2 + 2.5 * radius + atom_x
                    else:# for even put on left
                        bond_x = -radius*2.5 - 2.5*cut_radius*(j-2)/2 + atom_x
                    bond_z = bond_length/2-trap_depth
                    bond_models.append(bond_cylinder(bond_x,atom_y,bond_z,bond_length,cut_radius,tolerance))

                #finalize models
                model3 = union()(*bond_models)
                p_models.append(model1 + model2 + model3)
                #next atom
                atom_idx+=1

            #combone all models on buildplate and make STL
            model = union()(*p_models)
            os.makedirs(STL_path, exist_ok=True)
            scad_file = os.path.join(STL_path, f"hemispheres{p+1}.scad")
            stl_file = os.path.join(STL_path, f"hemispheres{p+1}.stl")
            scad_render_to_file(model, scad_file, file_header="$fn=50;")

            print(f"Wrote SCAD: {scad_file}")
            # Export STL using OpenSCAD
            subprocess.run([openscad_cmd, "-o", stl_file, scad_file], check=True)
            print(f"Exported STL: {stl_file}")'''
    #layout algorithm. not perfect, but fairly space effective if the nuclear radii are similar
    #for big prints, this will break up the print into different stl files if 3D printer build plate fills
    # will also make info txt files that describes each set of atom. each atom build plate will get a txt describing each atom
    # bonds will also make a txt describing the layout, but this index starts at the first build plate that contains bonds   
    max_radius = 0 # in units of Anstroms
    for atom in atoms:
        if atom.radius > max_radius:
            max_radius = atom.radius

    row_max = plate_dim[0]//(radius*5)
    col_max = plate_dim[1]//(radius*2.5)
    if (atoms_per_plate:=row_max*col_max)>0:
        plate_total = math.ceil(len(atoms)/atoms_per_plate)
    else:
        print(f'Error: plate too small for radius = {radius}')
        if log:
            log.log_message(f'Error: plate too small for radius = {radius}', error=True)
        return
        

    
    
    #phase one: add only hemisphere pairs to the buildplate
    atom_idx = 0 # counter must keep counting regardless of plate number
    plates = []
    for p in range(plate_total): #allow for multiple plates adn thus multiple stls
        output_info = open(os.path.join(STL_path,f"AtomInfo{p+1}.txt"), "w")       
        output_info.write("Atom Index,\tAtom,\tx,\ty,\tz,\tBond Count,\tColumn\n")
        plates.append([]) #stores all the models for a single plate
        while atom_idx<atoms_per_plate*(p+1) and atom_idx < len(atoms):
            atom = atoms[atom_idx] #get atom
            output_info.write(f"{atom_idx+1},\t{atom.element}{atom.charge},\t{atom.x},\t{atom.y},\t{atom.z},\t{len(atom.bonds)},\t{(atom_idx//col_max)+1}\n")
            atom_x = ((atom_idx%atoms_per_plate)//col_max) * radius*5 + 2.5*radius#these are the x,y coords of the center between the hemispheres
            atom_y = ((atom_idx%atoms_per_plate)%col_max) * radius*2.5 + 1.25 * radius

            #get cut angles
            cut_angles1,cut_angles2 = reorient_angles(cut_planes[atom_idx],atom.bonds)
            atom_print_r = atom.radius/max_radius*radius*nucleus_scale #convert atom radius in anstroms into units of radius (mm)
            cuts1 = [cut_cylinder_on_hemisphere(atom_print_r, pol, az, cut_depth, cut_radius) for pol, az in cut_angles1]
            cuts2 = [cut_cylinder_on_hemisphere(atom_print_r, pol, az, cut_depth, cut_radius) for pol, az in cut_angles2]

            trap_side = atom_print_r * trap_side_ratio
            trap_depth = atom_print_r * trap_depth_ratio

            #make female hemispher
            model1 = translate([radius*1.25 + atom_x, -atom_y, -trap_depth])(difference()(hemisphere(atom_print_r), *cuts1) - trapozoidal_prism(trap_side,trap_oblong_ratio,trap_depth) - make_trap_pyramid(trap_side,trap_depth, trap_oblong_ratio, overhang_angle)) #female
            #make male hemisphere
            model2 = translate([-radius*1.25 + atom_x, -atom_y, 0])(difference()(hemisphere(atom_print_r), *cuts2) + translate([0, 0, -trap_depth])(trapozoidal_prism(trap_side-tolerance/2,trap_oblong_ratio,trap_depth-tolerance))) #male
            

            #add hemisphere models
            plates[p].append(model1 + model2)
            #next atom
            atom_idx+=1
        output_info.close()
        
            

    #phase 2: add bond models
    bond_idx = 0
    bond_data = []
    for atom in atoms:#collect bond info
        for atom_bond_idx in range(len(atom.bonds)):
            curr_atom_idx = atom.bond_atom_idxs[atom_bond_idx][0]
            other_atom_idx = atom.bond_atom_idxs[atom_bond_idx][1]
            if curr_atom_idx<other_atom_idx:
                atom_print_r = atom.radius/max_radius*radius*nucleus_scale #convert atom radius in anstroms into units of radius (mm)
                bond_length = (atom.bonds[atom_bond_idx][0] - (atom.bond_atom_sizes[atom_bond_idx][0] + atom.bond_atom_sizes[atom_bond_idx][1])*nucleus_scale)*atom_print_r/max_radius + 2*cut_depth
                bond_data.append((bond_idx,
                                    str(atom.element)+str(atom.charge),
                                    str(atoms[other_atom_idx].element)+str(atoms[other_atom_idx].charge),
                                    curr_atom_idx,
                                    other_atom_idx,
                                    atom.bonds[atom_bond_idx][0], #bond length in angstroms
                                    bond_length
                                                                            ))
                bond_idx +=1
    bond_idx = 0

    #first complete column of bonds with
    atoms_on_unfinished_plate = len(atoms)%atoms_per_plate
    remaining_length = plate_dim[1] - radius*2.5*(atoms_on_unfinished_plate%col_max)
    rem_cols = remaining_length // (cut_radius*2.5)
    rem_rows = radius *5 // (cut_radius*2.5)
    rem_bonds = rem_cols*rem_rows
    new_x_O = (atoms_on_unfinished_plate//col_max) * radius*5
    new_y_O = (atoms_on_unfinished_plate%col_max) * radius*2.5
    
    while bond_idx<min(len(bond_data),rem_bonds):
        if bond_idx == 0: #Start new bond info file. There will be one file per plate, and all these bonds are nessicaroly on the same plate
            output_info = open(os.path.join(STL_path,f"BondInfo1.txt"), "w")
            output_info.write("Bond Index,\tAtom 1,\tAtom 2,\tAtom 1 Index,\tAtom 2 Index,\tBond Length (Angstrom),\tPrint Bond Length,\tColumn\n")
        bond_x = (bond_idx//rem_cols) * cut_radius*2 + 1.25*cut_radius + new_x_O#these are the x,y coords of the center between the hemispheres
        bond_y = (bond_idx%rem_cols) * cut_radius*2.5 + 1.25 * cut_radius + new_y_O#note, the "new origin" is moved to the location of the len(atoms) + 1th location
        bond_z = bond_data[bond_idx][6]/2-trap_depth
        plates[-1].append(bond_cylinder(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,tolerance))
        output_info.write(f"{bond_data[bond_idx][0]},\t{bond_data[bond_idx][1]},\t{bond_data[bond_idx][2]},\t{bond_data[bond_idx][3]},\t{bond_data[bond_idx][4]},\t{bond_data[bond_idx][5]},\t{bond_data[bond_idx][6]},\t{bond_idx//rem_cols+1}\n")
        bond_idx += 1


        # now full length columns for the remaining bonds on the same plate
    remaining_width = plate_dim[0] - radius*5*(atoms_on_unfinished_plate//col_max)
    fin_rows = remaining_width // (cut_radius*2.5)
    fin_cols = radius *5 // (cut_radius*2.5)
    fin_bonds = fin_cols*fin_rows

    new_x_O2 = (atoms_on_unfinished_plate//col_max+1) * radius*5
    while bond_idx-rem_bonds<min(len(bond_data)-rem_bonds,fin_bonds): # autoskips if bonds are already made
        if bond_idx == 0: #Start new bond info file. There will be one file per plate
            output_info = open(os.path.join(STL_path,f"BondInfo1.txt"), "w")
            output_info.write("Bond Index,\tAtom 1,\tAtom 2,\tAtom 1 Index,\tAtom 2 Index,\tBond Length (Angstrom),\tPrint Bond Length,\tColumn\n")
        bond_x = (bond_idx//fin_cols) * cut_radius*2 + 1.25*cut_radius + new_x_O2#these are the x,y coords of the center between the hemispheres
        bond_y = (bond_idx%fin_cols) * cut_radius*2.5 + 1.25 * cut_radius
        bond_z = bond_data[bond_idx][6]/2-trap_depth
        plates[-1].append(bond_cylinder(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,tolerance))
        output_info.write(f"{bond_data[bond_idx][0]},\t{bond_data[bond_idx][1]},\t{bond_data[bond_idx][2]},\t{bond_data[bond_idx][3]},\t{bond_data[bond_idx][4]},\t{round(bond_data[bond_idx][5],4)},\t{round(bond_data[bond_idx][6],4)},\t{bond_idx//rem_cols+1}\n")
        bond_idx += 1


    # now start new plate(s) is that is required
    bond_row_max = plate_dim[0]//(cut_radius*2.5)
    bond_col_max = plate_dim[1]//(cut_radius*2.5)
    bonds_per_plate = bond_row_max*bond_col_max
    bond_plate_total = math.ceil((len(bond_data)-bond_idx)/bonds_per_plate) #number of extra plates for just bonds
    if bond_idx == 0:#handle if no bond files have been created
        filenum = 0
    else:
        filenum = 1
        output_info.close()#close the previous bond file if it exists
    for p in range(bond_plate_total):
        plates.append([])
        output_info = open(os.path.join(STL_path,f"BondInfo{filenum + p + 1}.txt"), "w")
        output_info.write("Bond Index,\tAtom 1,\tAtom 2,\tAtom 1 Index,\tAtom 2 Index,\tBond Length (Angstrom),\tPrint Bond Length,\tColumn\n")
        while bond_idx-rem_bonds-fin_bonds<min(len(bond_data)-rem_bonds-fin_bonds,bonds_per_plate):
            bond_x = (bond_idx//bond_col_max) * cut_radius*2 + 1.25*cut_radius #these are the x,y coords of the center between the hemispheres
            bond_y = (bond_idx%bond_col_max) * cut_radius*2.5 + 1.25 * cut_radius
            bond_z = bond_data[bond_idx][6]/2-trap_depth
            plates[-1].append(bond_cylinder(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,tolerance))
            output_info.write(f"{bond_data[bond_idx][0]},\t{bond_data[bond_idx][1]},\t{bond_data[bond_idx][2]},\t{bond_data[bond_idx][3]},\t{bond_data[bond_idx][4]},\t{bond_data[bond_idx][5]},\t{bond_data[bond_idx][6]},\t{bond_idx//rem_cols+1}\n")
            bond_idx += 1   
        output_info.close()
    #phase 3: combone all models on buildplate and make STL

    if log:
        if len(plates) == 1:
            log.log_message(f'Generating 1 buildplate', error=False)
        else:
            log.log_message(f'Generating {len(plates)} buildplates', error=False)
    
    for p_num, p_models in enumerate(plates):
        model = union()(*p_models)
        os.makedirs(STL_path, exist_ok=True)
        scad_file = os.path.join(STL_path, f"BuildPlate{p_num+1}.scad")
        stl_file = os.path.join(STL_path, f"BuildPlate{p_num+1}.stl")
        scad_render_to_file(model, scad_file, file_header=f"$fn={res};",include_orig_code=False)# if you include orig code, it tries to call a file that doesnt exist and creates an error

        print(f"Wrote SCAD: {scad_file}")
        if log:
            log.log_message(f'Wrote OpenSCAD file {p_num+1}: {scad_file}', error=False)
            t = log._render_time_estimate()
            log.log_message(f'Rendering STL file {p_num+1}', error=False)
            if math.ceil(t)>3:
                log.log_message(f'Warning. Long rendering time, roughly {math.ceil(t)} minutes', error=True)
            elif t>1:
                log.log_message(f'Rendering time, roughly {math.ceil(t)} minutes', error=False)
            else:
                log.log_message(f'Rendering time, less than a minute', error=False)

        # Export STL using OpenSCAD
        a = subprocess.run([openscad_cmd, "-o", stl_file, scad_file], check=True,capture_output=True,text=True)
        if log:
            log.log_message(f'Wrote STL file {p_num+1}: {stl_file}', error=False)
            log.log_message(f'{a.stderr.splitlines()[4]}\n{a.stderr.splitlines()[7]}', error=False)

        print(f"Exported STL: {stl_file}")
        if log:
            log.log_message(f'STL OpenSCAD file {p_num+1}: {stl_file}', error=False)
    output_info.close()
    if log:
        log.log_message(f'Export Successful', error=False)