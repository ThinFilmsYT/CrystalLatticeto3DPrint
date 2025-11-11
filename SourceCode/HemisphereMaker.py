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

    return translate([position[0],position[1],position[2]])(rotate(a=[0,polar*180/np.pi,azimuth*180/np.pi])(cyl))

def cut_dshaft(radius, theta_v, phi_v, cut_depth, cyl_radius, overhang_angle, theta_z, phi_z, theta_x, phi_x):
    """
    Creates a D-shaped prism aligned with a rotated local coordinate system (x', y', z').
    The prism is rotated so that its x' axis points in direction (theta_v, phi_v),
    and its local xy-plane is minimally tilted so its normal passes through +z.

    Parameters:
        length: length of prism along x' axis
        radius: radius of the circular part of the D
        flat_cut: distance from center to flat side (0 < flat_cut < radius)
        theta_xp, phi_xp: polar and azimuthal angles of x' in global coords
        theta_yp, phi_yp: polar and azimuthal angles of y' in global coords
        theta_v, phi_v: target direction (in global coords) that x' will align to
    """
    flat_cut = cyl_radius*(1-np.cos(overhang_angle))
    prism = d_prism(cyl_radius,flat_cut,cut_depth)

    # --- compute gamma ---
    # D shaft flat edge is represented by (1,0,0) vector.
    # Like every other point will be rotated gamma about z, theta about y and phi about z
    # this puts them at vector position: 
    # vf_x = np.cos(gamma) * np.cos(phi) * np.cos(theta) - np.sin(gamma) * np.sin(phi)
    # vf_y = np.sin(gamma) * np.cos(phi) + np.cos(gamma) * np.sin(phi) * np.cos(theta)
    # vf_z = -np.cos(gamma) * np.sin(theta)

    # The z' unit vector is 
    # z'_x = np.sin(theta_p) * np.cos(phi_p)
    # z'_y = np.sin(theta_p) * np.sin(phi_p)
    # z'_z = np.cos(theta_p)

    # the distance between these vectors is d=sqrt(1-(Acos(gamma)+Bsin(gamma))^2) with A,B defined below.
    # minimize this distance using
    # 2 (B cos(x) - A sin(x)) (A cos(x) + B sin(x)) = 0
    # with 4 solutions in the x range of 0-2pi (gamma 1-4 below)



    # compute Δφ
    delta_phi = phi_v - phi_z
    
    # define coefficients
    A = np.sin(theta_z) * np.cos(theta_v) * np.cos(delta_phi) - np.sin(theta_v) * np.cos(theta_z)
    B = -np.sin(theta_z) * np.sin(delta_phi)
    
    # print('A',A,'B',B)

    # gamma1 = 2*(np.arctan2(B-np.sqrt(A**2+B**2),A)) % (2*np.pi)
    # gamma2 = 2*(np.arctan2(B+np.sqrt(A**2+B**2),A)) % (2*np.pi)

    gamma3 = 2*(np.arctan2(-A-np.sqrt(A**2+B**2),B)) % (2*np.pi)
    gamma4 = 2*(np.arctan2(-A+np.sqrt(A**2+B**2),B)) % (2*np.pi)
    # print(gamma1,gamma2,gamma3, gamma4)
    # print(np.sqrt(1 - (A * np.cos(gamma1) + B * np.sin(gamma1))**2),np.sqrt(1 - (A * np.cos(gamma2) + B * np.sin(gamma2))**2), np.sqrt(1 - (A * np.cos(gamma3) + B * np.sin(gamma3))**2),np.sqrt(1 - (A * np.cos(gamma4) + B * np.sin(gamma4))**2))
    dot3 = A * np.cos(gamma3) + B * np.sin(gamma3)
    dot4 = A * np.cos(gamma4) + B * np.sin(gamma4)
    # print(dot3,dot4)#,theta_v,theta_z, phi_v, phi_z)
    if dot3< dot4: #check which z position of the rotated (1,0,0) vector is higher
        gamma = gamma3
    else:
        gamma = gamma4
    # to make sure exactly top and bottom facing bonds globally match, always make them point to x axis
    bond_dot_zp = np.sin(theta_v)*np.sin(theta_z)*np.cos(phi_v-phi_z) + np.cos(theta_v)*np.cos(theta_z)
    if abs(bond_dot_zp) >= 1.0 - 1e-6: # if bond is aligned or antialigned to z prime direction, align to global x instead
        delta_phi2 = phi_v - phi_x
    
        # define coefficients
        A2 = np.sin(theta_x) * np.cos(theta_v) * np.cos(delta_phi2) - np.sin(theta_v) * np.cos(theta_x)
        B2 = -np.sin(theta_x) * np.sin(delta_phi2)
        gamma3 = 2*(np.arctan2(-A2-np.sqrt(A2**2+B2**2),B2)) % (2*np.pi)
        gamma4 = 2*(np.arctan2(-A2+np.sqrt(A2**2+B2**2),B2)) % (2*np.pi)
        # print(gamma1,gamma2,gamma3, gamma4)
        # print(np.sqrt(1 - (A * np.cos(gamma1) + B * np.sin(gamma1))**2),np.sqrt(1 - (A * np.cos(gamma2) + B * np.sin(gamma2))**2), np.sqrt(1 - (A * np.cos(gamma3) + B * np.sin(gamma3))**2),np.sqrt(1 - (A * np.cos(gamma4) + B * np.sin(gamma4))**2))
        dot3 = A2 * np.cos(gamma3) + B2 * np.sin(gamma3)
        dot4 = A2 * np.cos(gamma4) + B2 * np.sin(gamma4)
        # print(dot3,dot4)#,theta_v,theta_z, phi_v, phi_z)
        if dot3< dot4: #check which z position of the rotated (1,0,0) vector is higher
            gamma = gamma3
        else:
            gamma = gamma4

    x = np.sin(theta_v) * np.cos(phi_v)
    y = np.sin(theta_v) * np.sin(phi_v)
    z = np.cos(theta_v)

    # Position cylinder so its base starts just outside the hemisphere and cuts inward
    center_offset = radius - cut_depth
    position = [center_offset * x, center_offset * y, center_offset * z]
    # print(position)
    # Apply rotation to orient the prism
    return translate([position[0],position[1],position[2]])(rotate(a=[0, math.degrees(theta_v), math.degrees(phi_v)])(rotate(a=[0, 0, math.degrees(gamma)])(prism)))

def see_z(theta_z, phi_z):
    return rotate(a=[0, math.degrees(theta_z), math.degrees(phi_z)])(d_prism(.1, 0, 40))

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

def reorient_angles(normal_v, angles,get_orig_coords = False):
    """
    Rotate spherical coordinates so that the polar axis aligns with vector v.
    
    Parameters:
    -----------
    normal_v : tuple or list of 3 floats
        Vector (vx, vy, vz) defining the new polar axis.
    angles : list of tuples
        List of (polar, phi_v) angles in radians.
        
    Returns:
    --------
    list of tuples
        List of (polar', azimuth') angles in the rotated frame.
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
    orig_lt = []
    orig_ge = []

    for bond_leng, pol, az in angles:
        # print('before',pol,az)
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
        az_new = np.arctan2(vec_rot[1], vec_rot[0]) % (2*np.pi)
        # print('after',pol_new,az_new)

        x_rot = R @ np.array([1,0,0])
        r_x = np.linalg.norm(x_rot)
        pol_x = np.arccos(np.clip(x_rot[2]/r_x, -1.0, 1.0))
        az_x = np.arctan2(x_rot[1], x_rot[0]) % (2*np.pi)

        z_rot = R @ np.array([0,0,1])
        r_z = np.linalg.norm(z_rot)
        pol_z = np.arccos(np.clip(z_rot[2]/r_z, -1.0, 1.0))
        az_z = np.arctan2(z_rot[1], z_rot[0]) % (2*np.pi)

        # Classify
        if pol_new < np.pi/2:
            list_lt.append((pol_new, az_new))
            orig_lt.append((pol_z,az_z,pol_x,az_x))
        else:
            list_ge.append((np.pi-pol_new,2*np.pi-az_new))#flip hemisphere upright for better 3d printing
            orig_ge.append((np.pi-pol_z,2*np.pi-az_z,np.pi-pol_x,2*np.pi-az_x))
            
    if get_orig_coords:
        return list_lt, list_ge, orig_lt, orig_ge
    else:
        return list_lt, list_ge

def bond_cylinder(bond_x, bond_y, bond_z,bond_length, cut_radius, tolerance,ID_ratio = 0.7):
    return translate([bond_x, bond_y, bond_z])(cylinder(h=bond_length, r=cut_radius-tolerance/2, center=True) - cylinder(h=bond_length, r=cut_radius*ID_ratio, center=True))

def bond_dshaft(bond_x, bond_y, bond_z,bond_length, cut_radius, overhang_angle, tolerance):
    r = cut_radius - tolerance/2
    f = r*(1-np.cos(overhang_angle)) #this is how deeply the circle is cut into. the overhang determines how deep the cut needs to be so the printer can make it
    return (translate([bond_x - bond_length/2, bond_y , bond_z + cut_radius - f])(rotate([0,-90,180])(d_prism(r,f,bond_length))))

def d_prism(radius=1.0, flat_offset=0.0, length=10.0, circle_segments=128):
    """
    Create a local primitive prism that:
      - extends from x=0 to x=L (local coordinates)
      - cross-section lies in the local y-z plane and is a D-shaft:
          circle of radius `Radius` centered at (y=0,z=0) cut by chord
          whose plane is at z = +flat_offset (flat faces +z direction).
    NOTE: flat_offset must satisfy -Radius < flat_offset < Radius (reasonable).
    """

    R = float(radius)
    d = float(flat_offset)
    if d >= R or d <= -R:
        raise ValueError("flat_offset must be strictly between -R and R for a valid D.")

    # Build a circle in 2D XY plane (this 2D plane will encode (u,v) = (y,z))
    # We'll create the portion of the circle with u <= -d (this maps to z' >= d after final transforms).
    # To reason: after extrude & rotate, mapping makes this flat face point in +z' direction.
    circ = circle(r=R, segments=circle_segments)

    # Create a big rectangle whose max x is X_max = -d (so it keeps u <= -d).
    # rectangle lower-left at x = -2R, y = -2R, width = ( -d - (-2R) ) = 2R - d
    rect_w = 2*R - d
    rect_h = 2*R
    rect = translate([d/2, 0])(square([rect_w, rect_h], center=True))

    # Intersect circle with that rectangle -> keeps x <= -d region of circle (flat at x = -d).
    dshape_2d = intersection()(circ, rect)

    return linear_extrude(height=length)(dshape_2d)

def sph_to_cart(theta, phi):
        return [
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(theta)
        ]

def cartesian_to_spherical(x, y, z):
    """
    Convert Cartesian coordinates (x, y, z) to spherical coordinates (r, theta, phi).
    
    Returns:
        r: radial distance
        theta: polar angle (from +z axis)
        phi: azimuthal angle (from +x axis in x-y plane)
    """
    r = math.sqrt(x**2 + y**2 + z**2)
    if r == 0:
        theta = 0.0
        phi = 0.0
    else:
        theta = math.acos(z / r)  # polar angle
        phi = math.atan2(y, x)    # azimuthal angle
        if phi < 0:
            phi += 2 * math.pi    # ensure phi is in [0, 2pi)
    return r, theta, phi

def makeSTL(STL_path,atoms,cut_planes,plate_dim = (210,250),parameters = [],log = None, sep = False, openscad_cmd = "C:/Program Files/OpenSCAD/openscad.exe", bondtype = 'd'):
    try:
        failed = True
        # bondtype options: 'circle', 'd'
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
        
        # cut_depth = cut_depth_ratio*radius
        cut_radius = cut_radius_ratio*radius

        max_radius = 0 # in units of Anstroms
        min_radius = np.inf
        for atom in atoms:
            if atom.radius > max_radius:
                max_radius = atom.radius
            if atom.radius < min_radius:
                min_radius = atom.radius

        #add a check to see if bond cuts could intersect interior features
        #Max bond cut depth is just r_cut <= r_sphere * [1 - (trap_depth_ratio + trap_side_ratio*tan(overhang_angle))]
        #this is just the top of the pyramid
        #max bond cut depth should also be r_cut <= r_sphere * [1 - sqrt(trap_depth_ratio^2 + trap_side_ratio^2)]
        #this is the corner of the pyramid base
        max_cut_ratio = min(1 - (trap_depth_ratio + trap_side_ratio*np.tan(np.radians(overhang_angle))),1 - np.sqrt(trap_depth_ratio**2 + trap_side_ratio**2))
        if cut_depth_ratio> max_cut_ratio:
            print('Warning: Bond cuts might intersect interior features.')
            if log:
                log.log_message(f'Warning: Bond cuts might intersect interior features. Try using Cut Depth Ratio of {max_cut_ratio} or less', error=True)
            # return
        elif cut_radius_ratio*max_radius > 0.75 * min_radius:
            print('Error: Cut radius too big.')
            if log:
                log.log_message('Error: Cut radius too big.', error=True)
            # return

        time_total = 0
        if log:
            t = log._render_time_estimate()
            if math.ceil(t)>3:
                log.log_message(f'Warning. Long rendering time, roughly {math.ceil(t)} minutes', error=True)
            elif t>1:
                log.log_message(f'Rendering time, roughly {math.ceil(t)} minutes', error=False)
            else:
                log.log_message(f'Rendering time, less than a minute', error=False)
            
        

        bond_idx = 0
        bond_data = []
        for atom in atoms:#collect bond info
            for atom_bond_idx in range(len(atom.bonds)):
                curr_atom_idx = atom.bond_atom_idxs[atom_bond_idx][0]
                other_atom_idx = atom.bond_atom_idxs[atom_bond_idx][1]
                if curr_atom_idx<other_atom_idx:
                    atom_print_r = atom.radius/max_radius*radius*nucleus_scale #convert atom radius in anstroms into units of radius (mm)
                    cut_depth = cut_depth_ratio * atom_print_r
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
        
        if sep: # separate build plate for each type of atom, and separate build plate for the bonds
            #collect atom type info
            atom_kinds = []
            atom_lists = []
            atom_cuts = []
            for i, atom in enumerate(atoms):
                if atom.element in atom_kinds:
                    atom_lists[atom_kinds.index(atom.element)].append(atom)
                    atom_cuts[atom_kinds.index(atom.element)].append(cut_planes[i])
                else:
                    atom_kinds.append(atom.element)
                    atom_lists.append([atom])
                    atom_cuts.append([cut_planes[i]])

            #loop through atom types, making a new plate each time
            plates = []
            atom_plates = 0
            f_ind = 0
            for i, atom_list in enumerate(atom_lists):
                atom_idx = 0 # counter must keep counting regardless of plate number
                atom_print_r = atom_list[0].radius/max_radius*radius*nucleus_scale
                cut_depth = cut_depth_ratio * atom_print_r

                row_max = plate_dim[0]//(atom_print_r*5)
                col_max = plate_dim[1]//(atom_print_r*2.5)
                if (atoms_per_plate:=row_max*col_max)>0:
                    plate_total = math.ceil(len(atom_list)/atoms_per_plate)
                else:
                    print(f'Error: plate too small for radius = {radius}')
                    if log:
                        log.log_message(f'Error: plate too small for radius = {radius}', error=True)
                    return
                for p in range(plate_total): #allow for multiple plates and thus multiple stls
                    f_ind += 1
                    output_info = open(os.path.join(STL_path,f"AtomInfo{f_ind}.txt"), "w")       
                    output_info.write("Atom Index,\tAtom,\tx,\ty,\tz,\tBond Count,\tColumn\n")
                    plates.append([]) #stores all the models for a single plate
                    atom_plates += 1
                    while atom_idx<atoms_per_plate*(p+1) and atom_idx < len(atom_list):
                        atom = atom_list[atom_idx] #get atom
                        output_info.write(f"{atom_idx+1},\t{atom.element}{atom.charge},\t{atom.x},\t{atom.y},\t{atom.z},\t{len(atom.bonds)},\t{(atom_idx//col_max)+1}\n")
                        atom_x = ((atom_idx%atoms_per_plate)//col_max) * atom_print_r*5 + 2.5*atom_print_r#these are the x,y coords of the center between the hemispheres
                        atom_y = ((atom_idx%atoms_per_plate)%col_max) * atom_print_r*2.5 + 1.25 * atom_print_r

                        

                        #get cut angles
                        if bondtype == 'circle':
                            cut_angles1,cut_angles2 = reorient_angles(atom_cuts[i][atom_idx],atom.bonds,False)
                            cuts1 = [cut_cylinder_on_hemisphere(atom_print_r, pol, az, cut_depth, cut_radius) for pol, az in cut_angles1]
                            cuts2 = [cut_cylinder_on_hemisphere(atom_print_r, pol, az, cut_depth, cut_radius) for pol, az in cut_angles2]
                        elif bondtype == 'd':
                            cut_angles1,cut_angles2,orig_angles1,orig_angles2 = reorient_angles(atom_cuts[i][atom_idx],atom.bonds,True)
                            cuts1 = [cut_dshaft(atom_print_r, cut_angles1[i][0], cut_angles1[i][1], cut_depth, cut_radius, overhang_angle, orig_angles1[i][0], orig_angles1[i][1],orig_angles1[i][2],orig_angles1[i][3]) for i in range(len(cut_angles1))]
                            cuts2 = [cut_dshaft(atom_print_r, cut_angles2[i][0], cut_angles2[i][1], cut_depth, cut_radius, overhang_angle, orig_angles2[i][0], orig_angles2[i][1],orig_angles2[i][2],orig_angles2[i][3]) for i in range(len(cut_angles2))]
                        trap_side = atom_print_r * trap_side_ratio
                        trap_depth = atom_print_r * trap_depth_ratio

                        #make female hemispher
                        model1 = translate([atom_print_r*1.25 + atom_x, -atom_y, 0])(difference()(hemisphere(atom_print_r), *cuts1) - trapozoidal_prism(trap_side,trap_oblong_ratio,trap_depth) - make_trap_pyramid(trap_side,trap_depth, trap_oblong_ratio, overhang_angle)) #female
                        #make male hemisphere
                        model2 = translate([-atom_print_r*1.25 + atom_x, -atom_y, trap_depth-tolerance])(difference()(hemisphere(atom_print_r), *cuts2) + translate([0, 0, -(trap_depth-tolerance)])(trapozoidal_prism(trap_side-tolerance/2,trap_oblong_ratio,trap_depth-tolerance))) #male
                        plates[-1].append(model1 + model2)
                        atom_idx+=1
                    output_info.close()
            
            # now start new plate(s) for just bonds
            if bondtype == 'circle':
                bond_row_max = plate_dim[0]//(cut_radius*2.5)
                bond_col_max = plate_dim[1]//(cut_radius*2.5)
                bonds_per_plate = bond_row_max*bond_col_max
                bond_plate_total = math.ceil(len(bond_data)/bonds_per_plate) #number of extra plates for just bonds
                filenum = 0
                bond_idx = 0
                for p in range(bond_plate_total):
                    plates.append([])
                    output_info = open(os.path.join(STL_path,f"BondInfo{filenum + atom_plates + 1}.txt"), "w")
                    output_info.write("Bond Index,\tAtom 1,\tAtom 2,\tAtom 1 Index,\tAtom 2 Index,\tBond Length (Angstrom),\tPrint Bond Length,\tColumn\n")
                    while bond_idx-p*bonds_per_plate<min(len(bond_data)-p*bonds_per_plate,bonds_per_plate):
                        bond_x = ((bond_idx-p*bonds_per_plate)//bond_col_max) * cut_radius*2.5 + 1.25*cut_radius #these are the x,y coords of the center between the hemispheres
                        bond_y = ((bond_idx-p*bonds_per_plate)%bond_col_max) * cut_radius*2.5 + 1.25 * cut_radius
                        bond_z = bond_data[bond_idx][6]/2
                        plates[-1].append(bond_cylinder(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,tolerance))
                        output_info.write(f"{bond_data[bond_idx][0]},\t{bond_data[bond_idx][1]},\t{bond_data[bond_idx][2]},\t{bond_data[bond_idx][3]},\t{bond_data[bond_idx][4]},\t{bond_data[bond_idx][5]},\t{bond_data[bond_idx][6]},\t{bond_idx//bond_col_max+1}\n")
                        bond_idx += 1   
                    output_info.close()
                    
            elif bondtype == 'd':
                max_length = 0
                for b in bond_data:
                    if b[6] > max_length:
                        max_length = b[6]
                bond_row_max = plate_dim[0]//(max_length*1.25)
                bond_col_max = plate_dim[1]//(cut_radius*2.5)
                bonds_per_plate = bond_row_max*bond_col_max
                bond_plate_total = math.ceil(len(bond_data)/bonds_per_plate) #number of extra plates for just bonds
                filenum = 0
                bond_idx = 0
                for p in range(bond_plate_total):
                    plates.append([])
                    output_info = open(os.path.join(STL_path,f"BondInfo{filenum + atom_plates + 1}.txt"), "w")
                    output_info.write("Bond Index,\tAtom 1,\tAtom 2,\tAtom 1 Index,\tAtom 2 Index,\tBond Length (Angstrom),\tPrint Bond Length,\tColumn\n")
                    while bond_idx-p*bonds_per_plate<min(len(bond_data)-p*bonds_per_plate,bonds_per_plate):
                        bond_x = ((bond_idx-p*bonds_per_plate)//bond_col_max) * max_length*1.25 + 0.625 * max_length #these are the x,y coords of the center between the hemispheres
                        bond_y = ((bond_idx-p*bonds_per_plate)%bond_col_max) * cut_radius*2.5 + 1.25*cut_radius
                        bond_z = cut_radius*np.cos(overhang_angle)
                        plates[-1].append(bond_dshaft(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,overhang_angle, tolerance))
                        output_info.write(f"{bond_data[bond_idx][0]},\t{bond_data[bond_idx][1]},\t{bond_data[bond_idx][2]},\t{bond_data[bond_idx][3]},\t{bond_data[bond_idx][4]},\t{bond_data[bond_idx][5]},\t{bond_data[bond_idx][6]},\t{bond_idx//bond_col_max+1}\n")
                        bond_idx += 1   
                    output_info.close()


        else:
            #layout algorithm. not perfect, but fairly space effective if the nuclear radii are similar
            #for big prints, this will break up the print into different stl files if 3D printer build plate fills
            # will also make info txt files that describes each set of atom. each atom build plate will get a txt describing each atom
            # bonds will also make a txt describing the layout, but this index starts at the first build plate that contains bonds   
            row_max = plate_dim[0]//(radius*nucleus_scale*5)
            col_max = plate_dim[1]//(radius*nucleus_scale*2.5)
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
                    atom_x = ((atom_idx%atoms_per_plate)//col_max) * radius*nucleus_scale*5 + 2.5*radius*nucleus_scale#these are the x,y coords of the center between the hemispheres
                    atom_y = ((atom_idx%atoms_per_plate)%col_max) * radius*nucleus_scale*2.5 + 1.25 * radius*nucleus_scale

                    

                    #get cut angles
                    atom_print_r = atom.radius/max_radius*radius*nucleus_scale #convert atom radius in anstroms into units of radius (mm)
                    cut_depth = cut_depth_ratio * atom_print_r
                    if bondtype == 'circle':
                        cut_angles1,cut_angles2 = reorient_angles(cut_planes[atom_idx],atom.bonds,False)
                        cuts1 = [cut_cylinder_on_hemisphere(atom_print_r, pol, az, cut_depth, cut_radius) for pol, az in cut_angles1]
                        cuts2 = [cut_cylinder_on_hemisphere(atom_print_r, pol, az, cut_depth, cut_radius) for pol, az in cut_angles2]
                    elif bondtype == 'd':
                        cut_angles1,cut_angles2,orig_angles1,orig_angles2 = reorient_angles(cut_planes[atom_idx],atom.bonds,True)
                        cuts1 = [cut_dshaft(atom_print_r, cut_angles1[i][0], cut_angles1[i][1], cut_depth, cut_radius, overhang_angle, orig_angles1[i][0], orig_angles1[i][1],orig_angles1[i][2],orig_angles1[i][3]) for i in range(len(cut_angles1))]
                        cuts2 = [cut_dshaft(atom_print_r, cut_angles2[i][0], cut_angles2[i][1], cut_depth, cut_radius, overhang_angle, orig_angles2[i][0], orig_angles2[i][1],orig_angles2[i][2],orig_angles2[i][3]) for i in range(len(cut_angles2))]
                        
                    trap_side = atom_print_r * trap_side_ratio
                    trap_depth = atom_print_r * trap_depth_ratio

                    #make female hemispher
                    model1 = translate([radius*nucleus_scale*1.25 + atom_x, -atom_y, 0])(difference()(hemisphere(atom_print_r), *cuts1) - trapozoidal_prism(trap_side,trap_oblong_ratio,trap_depth) - make_trap_pyramid(trap_side,trap_depth, trap_oblong_ratio, overhang_angle)) #female
                    #make male hemisphere
                    model2 = translate([-radius*nucleus_scale*1.25 + atom_x, -atom_y, trap_depth-tolerance])(difference()(hemisphere(atom_print_r), *cuts2) + translate([0, 0, -(trap_depth-tolerance)])(trapozoidal_prism(trap_side-tolerance/2,trap_oblong_ratio,trap_depth-tolerance))) #male
                    #add hemisphere models
                    plates[p].append(model1 + model2)
                    #next atom
                    atom_idx+=1
                output_info.close()
                
                    

            #phase 2: add bond models

            #first complete column of bonds with
            atoms_on_unfinished_plate = len(atoms)%atoms_per_plate
            remaining_length = plate_dim[1] - radius*nucleus_scale*2.5*(atoms_on_unfinished_plate%col_max)
            if bondtype == 'circle':
                rem_cols = remaining_length // (cut_radius*2.5)
                rem_rows = radius*nucleus_scale *5 // (cut_radius*2.5)
                
            elif bondtype == 'd':
                max_length = 0
                for b in bond_data:
                    if b[6] > max_length:
                        max_length = b[6]
                rem_cols = remaining_length // (cut_radius*2.5)
                rem_rows = radius*nucleus_scale * 5 // (max_length*1.25)

            rem_bonds = rem_cols*rem_rows
            new_x_O = (atoms_on_unfinished_plate//col_max) * radius*nucleus_scale*5
            new_y_O = (atoms_on_unfinished_plate%col_max) * radius*nucleus_scale*2.5

            while bond_idx<min(len(bond_data),rem_bonds):
                if bond_idx == 0: #Start new bond info file. There will be one file per plate, and all these bonds are nessicaroly on the same plate
                    output_info = open(os.path.join(STL_path,f"BondInfo1.txt"), "w")
                    output_info.write("Bond Index,\tAtom 1,\tAtom 2,\tAtom 1 Index,\tAtom 2 Index,\tBond Length (Angstrom),\tPrint Bond Length,\tColumn\n")
                if bondtype == 'circle':
                    bond_x = (bond_idx//rem_cols) * cut_radius*2.5 + 1.25*cut_radius + new_x_O#these are the x,y coords of the center between the hemispheres
                    bond_y = (bond_idx%rem_cols) * cut_radius*2.5 + 1.25 * cut_radius + new_y_O#note, the "new origin" is moved to the location of the len(atoms) + 1th location
                    bond_z = bond_data[bond_idx][6]/2
                    plates[-1].append(bond_cylinder(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,tolerance))
                elif bondtype == 'd':
                    bond_x = (bond_idx//rem_cols) * max_length*1.25 + 0.625 * max_length + new_x_O#these are the x,y coords of the center between the hemispheres
                    bond_y = (bond_idx%rem_cols) * cut_radius*2.5 + 1.25*cut_radius + new_y_O
                    bond_z = cut_radius*np.cos(overhang_angle)
                    plates[-1].append(bond_dshaft(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,overhang_angle, tolerance))
                
                output_info.write(f"{bond_data[bond_idx][0]},\t{bond_data[bond_idx][1]},\t{bond_data[bond_idx][2]},\t{bond_data[bond_idx][3]},\t{bond_data[bond_idx][4]},\t{bond_data[bond_idx][5]},\t{bond_data[bond_idx][6]},\t{bond_idx//rem_cols+1}\n")
                bond_idx += 1


            #now fill with full length columns for the remaining bonds on the same plate
            remaining_width = plate_dim[0] - radius*nucleus_scale*5*math.ceil(atoms_on_unfinished_plate/col_max)
            fin_cols = plate_dim[1] // (cut_radius*2.5)
            if bondtype == 'circle':
                fin_rows = remaining_width // (cut_radius*2.5)
            elif bondtype == 'd':
                fin_rows = remaining_width // (max_length*1.25)
            fin_bonds = fin_cols*fin_rows

            new_x_O2 = (atoms_on_unfinished_plate//col_max+1) * radius*nucleus_scale*5
            while bond_idx-rem_bonds<min(len(bond_data)-rem_bonds,fin_bonds): # autoskips if bonds are already made
                if bond_idx == 0: #Start new bond info file. There will be one file per plate
                    output_info = open(os.path.join(STL_path,f"BondInfo1.txt"), "w")
                    output_info.write("Bond Index,\tAtom 1,\tAtom 2,\tAtom 1 Index,\tAtom 2 Index,\tBond Length (Angstrom),\tPrint Bond Length,\tColumn\n")
                if bondtype == 'circle':
                    bond_x = ((bond_idx-rem_bonds)//fin_cols) * cut_radius*2.5 + 1.25*cut_radius + new_x_O2#these are the x,y coords of the center between the hemispheres
                    bond_y = ((bond_idx-rem_bonds)%fin_cols) * cut_radius*2.5 + 1.25 * cut_radius
                    bond_z = bond_data[bond_idx][6]/2
                    plates[-1].append(bond_cylinder(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,tolerance))
                elif bondtype == 'd':
                    bond_x = ((bond_idx-rem_bonds)//fin_cols) * max_length*1.25 + 0.625 * max_length + new_x_O2#these are the x,y coords of the center between the hemispheres
                    bond_y = ((bond_idx-rem_bonds)%fin_cols) * cut_radius*2.5 + 1.25*cut_radius
                    bond_z = cut_radius*np.cos(overhang_angle)
                    plates[-1].append(bond_dshaft(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,overhang_angle, tolerance))
                output_info.write(f"{bond_data[bond_idx][0]},\t{bond_data[bond_idx][1]},\t{bond_data[bond_idx][2]},\t{bond_data[bond_idx][3]},\t{bond_data[bond_idx][4]},\t{round(bond_data[bond_idx][5],4)},\t{round(bond_data[bond_idx][6],4)},\t{bond_idx//fin_cols+1}\n")
                bond_idx += 1


            # now start new plate(s) that is required
            if bondtype == 'circle':
                bond_row_max = plate_dim[0]//(cut_radius*2.5)
            elif bondtype == 'd':
                bond_row_max = plate_dim[0]//(max_length*1.25)
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
                while bond_idx-rem_bonds-fin_bonds-p*bonds_per_plate<min(len(bond_data)-rem_bonds-fin_bonds-p*bonds_per_plate,bonds_per_plate):
                    if bondtype == 'circle':
                        bond_x = ((bond_idx-rem_bonds-fin_bonds-p*bonds_per_plate)//bond_col_max) * cut_radius*2.5 + 1.25*cut_radius #these are the x,y coords of the center between the hemispheres
                        bond_y = ((bond_idx-rem_bonds-fin_bonds-p*bonds_per_plate)%bond_col_max) * cut_radius*2.5 + 1.25 * cut_radius
                        bond_z = bond_data[bond_idx][6]/2
                        plates[-1].append(bond_cylinder(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,tolerance))
                    elif bondtype == 'd':
                        bond_x = ((bond_idx-rem_bonds-fin_bonds-p*bonds_per_plate)//bond_col_max) * max_length*1.25 + 0.625 * max_length#these are the x,y coords of the center between the hemispheres
                        bond_y = ((bond_idx-rem_bonds-fin_bonds-p*bonds_per_plate)%bond_col_max) * cut_radius*2.5 + 1.25*cut_radius
                        bond_z = cut_radius*np.cos(overhang_angle)
                        plates[-1].append(bond_dshaft(bond_x,-bond_y,bond_z,bond_data[bond_idx][6],cut_radius,overhang_angle, tolerance))
                    output_info.write(f"{bond_data[bond_idx][0]},\t{bond_data[bond_idx][1]},\t{bond_data[bond_idx][2]},\t{bond_data[bond_idx][3]},\t{bond_data[bond_idx][4]},\t{bond_data[bond_idx][5]},\t{bond_data[bond_idx][6]},\t{bond_idx//bond_col_max+1}\n")
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
                log.log_message(f'Rendering STL file {p_num+1}', error=False)

            # Export STL using OpenSCAD
            a = subprocess.run([openscad_cmd, "-o", stl_file, scad_file], check=True,capture_output=True,text=True)
            # print(a.stderr)
            # print(a.stdout)
            if log:
                log.log_message(f'Wrote STL file {p_num+1}: {stl_file}', error=False)
                log.log_message(f'{a.stderr.splitlines()[4]}\n{a.stderr.splitlines()[7]}', error=False)
                # print(a.stderr)
            _,hours,mins,secs = a.stderr.splitlines()[4].split(':')
            time_total += int(hours)*3600 + int(mins) * 60 + float(secs)
            print(f"Exported STL: {stl_file}")
            if log:
                log.log_message(f'STL OpenSCAD file {p_num+1}: {stl_file}', error=False)
            output_info.close()
            
        if log:
            log.log_message(f'Export Successful', error=False)
        log._register_render_time(res,time_total)

    
    except PermissionError:
        log.log_message(f"Error: Permission to write files denied. Try different output path.", error=True)
        return
    except Exception as e:
        log.log_message(f"Error in CAD file making: {e}", error=True)
        return