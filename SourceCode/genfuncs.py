import json
import numpy as np
import math

class Atom:
    def __init__(self,element,charge,x,y,z,occ=1):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.occupation = occ
        self.radius = self.get_radii_data(element,charge)
        self.charge = charge
        self.bonds= [] # [(bond_length,polar,azimuth),()]
        self.bond_atom_sizes = [] # list of tuples (atom 1 radius, atom 2 radius)
        self.bond_atom_idxs = [] # list of tuples (atom 1 index, atom 2 index)
        self.lattice_index = None
        self.coord_number = None

    def fractional_position(self):
        return (self.x, self.y, self.z)
    
    def read_ele_data(self):
        with open('ElementalRadiiData.json') as f:
            out = f.read()
        return json.loads(out)

    def get_radii_data(self, ele, charge, info = 'r_crystal', coord = '',spin_state = ''):
        charge = str(charge)
        if type(coord) is int:
            coord = int_to_roman(coord)
        data = self.read_ele_data()
        if charge == '0': #for covalent crystals
            return data[ele][charge][info] #simple. ignoring coordination number
        else: #for ions
            if coord != '':#given coordination number
                coord_info = data[ele][charge][coord]
                if len(coord_info) == 4: #no spin state data
                    return coord_info[info]
                elif len(coord_info) == 2: #has data for High and Low spin states
                    if spin_state != '':
                        return coord_info[spin_state][info]
                    else:
                        print('Warning: Returning average of two spin states')
                        return (coord_info['Low Spin'][info] + coord_info['High Spin'][info])/2
                else:
                    raise NameError("json file not in expected format. Try redownloading")
            else:#no coordination number
                charge_info = data[ele][charge]
                summ = 0
                unique_nums = 0
                for coord in charge_info:
                    coord_info = charge_info[coord]
                    if len(coord_info) == 4: #no spin state data
                        try:
                            summ += float(coord_info[info])
                            unique_nums += 1
                        except:
                            raise NameError("Need to include coordination number if not looking for radii data")
                    elif len(coord_info) == 2: #has data for High and Low spin states
                        if spin_state != '':
                            try:
                                summ += float(coord_info[spin_state][info])
                                unique_nums += 1
                            except:
                                raise NameError("Need to include coordination number if not looking for radii data")
                        else:
                            try:
                                summ += float(coord_info['Low Spin'][info]) + float(coord_info['High Spin'][info])
                                unique_nums += 2
                            except:
                                raise NameError("Need to include coordination number if not looking for radii data")
                            
                    else:
                        raise NameError("json file not in expected format. Try redownloading")
                if unique_nums>1:
                    print(f'Warning: Returning average of {unique_nums} states for {ele}{charge}')
                # print(summ/unique_nums)
                return summ/unique_nums

    def __repr__(self):
        return f"{self.element}: x={self.x}, y={self.y}, z={self.z}, radius={self.radius}, charge={self.charge}"

class Lattice:
    def __init__(self, cif_path,nx,ny,nz,bond_range = 0.1):
        self.cif_path = cif_path
        self.nx,self.ny,self.nz = nx,ny,nz
        self.bond_range = bond_range # error bars around bond range. for an atom with bond length 1.0 and 1.09 both are included for bond range = .1
        raw_atoms, params = self.read_cif(cif_path)
        self.raw_atoms = raw_atoms # raw atoms for just a unit cell of the lattice
        self.a = params.get('_cell_length_a', 0.0) #lengths in angstroms
        self.b = params.get('_cell_length_b', 0.0)
        self.c = params.get('_cell_length_c', 0.0)
        self.alpha = params.get('_cell_angle_alpha', 90.0) #all angles in degrees, because that is cif syntax
        self.beta  = params.get('_cell_angle_beta', 90.0)
        self.gamma = params.get('_cell_angle_gamma', 90.0)
        #initialize atom and bonds
        self.atoms = [] # list of Atom objects
        self.bonds = [] # list of tuples (r, theta, phi)

        atoms = [] # another list of Atom objects
        seen = set()
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    for raw_atom in raw_atoms: #raw atom is of the form: (label, xi, yi, zi, occ) where label is a string like 'Ru4' or 'O2-'
                        ele_name, charge = self.get_ele_charge(raw_atom[0])
                        fx, fy, fz = raw_atom[1] + ix, raw_atom[2] + iy, raw_atom[3] + iz
                        if not( (fx,fy,fz) in seen):
                            atoms.append(Atom(ele_name, charge,fx,fy,fz,raw_atom[4]))
                            seen.add((fx,fy,fz))
        
        self.add_atoms(atoms) #add Atom objects 
        self.bonds = self.find_bonds(self.atoms,self.atoms) #find bonds
        # self.refine_atomic_radii()

    def read_cif(self,path):
        with open(path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        params = {
            '_cell_length_a': None,
            '_cell_length_b': None,
            '_cell_length_c': None,
            '_cell_angle_alpha': None,
            '_cell_angle_beta': None,
            '_cell_angle_gamma': None
        }

        symmetry_ops = []
        atom_blocks = []
        current_headers = []
        current_rows = []
        mode = None

        for line in lines:
            if line.startswith('loop_'):
                # Commit previous loop if it's an atom block
                if mode == 'atom' and current_headers and current_rows:
                    atom_blocks.append((current_headers, current_rows))
                current_headers = []
                current_rows = []
                mode = None
                continue

            # Find Lattice constants
            if any(line.startswith(k) for k in params):
                for k in params:
                    if line.startswith(k):
                        val = line.split()[1]
                        if '(' in val:
                            val = val[:val.find('(')]
                        params[k] = float(val)
                continue

            # Initiate Symmetry mode
            if line.startswith('_symmetry_equiv_pos_as_xyz') or line.startswith('_space_group_symop_operation_xyz'):
                mode = 'sym'
                continue

            # Initaite Atom site mode
            if line.startswith('_atom_site_'):
                current_headers.append(line)
                mode = 'atom'
                continue

            # Collect symmetry operations
            if mode == 'sym' and not line.startswith('_'):
                symmetry_ops.append(line.replace("'", "").replace('"', '').strip())

            # Collect atom site data
            elif mode == 'atom' and current_headers and not line.startswith('_'):
                if not(line.startswith('#')): #don't include commented lines
                    current_rows.append(line.split())

        # Capture final block if needed
        if mode == 'atom' and current_headers and current_rows:
            atom_blocks.append((current_headers, current_rows))

        # Choose the best atom block: one with fract_x/y/z fields and consistent data
        best_atoms = []
        for headers, rows in atom_blocks:
            try:
                has_x = any(h.endswith('fract_x') for h in headers)
                has_y = any(h.endswith('fract_y') for h in headers)
                has_z = any(h.endswith('fract_z') for h in headers)
                if has_x and has_y and has_z and len(rows) > 0:
                    best_atoms = rows
                    break  # take first valid one
            except:
                continue
        
        best_atoms = self.extract_atoms_with_symmetry(headers, best_atoms, symmetry_ops)
        return best_atoms, params

    def add_atoms(self, atoms):
        i = len(self.atoms)
        for atom in atoms:
            atom.lattice_index = i
            self.atoms.append(atom)
            i += 1

    def get_atomic_distance(self, atom1, atom2):
        P = atom1.fractional_position()
        Q = atom2.fractional_position()
        # return np.sqrt((**(c1[0] - c2[0]))**2+(lattice.b*(c1[1] - c2[1])**2)+(lattice.c*(c1[2] - c2[2])**2))
        """
        Distance between two points P and Q in a 3D skewed basis (a,b,c).

        Points P and Q are given as coordinates (pa, pb, pc) in the (a,b,c) basis.
        """

        # Differences in basis coordinates
        da = P[0] - Q[0]
        db = P[1] - Q[1]
        dc = P[2] - Q[2]

        # Dot products between basis vectors

        a_dot_b = self.a * self.b * np.cos(np.radians(self.gamma))
        a_dot_c = self.a * self.c * np.cos(np.radians(self.beta))
        b_dot_c = self.b * self.c * np.cos(np.radians(self.alpha))

        # Expand distance^2 = vᵀ G v, where v = (da,db,dc)
        dist_sq = (
            (da * self.a)**2 +
            (db * self.b)**2 +
            (dc * self.c)**2 +
            2 * da * db * a_dot_b +
            2 * da * dc * a_dot_c +
            2 * db * dc * b_dot_c
        )

        return np.sqrt(dist_sq)

    def find_bonds(self, core_atoms, bondable_atoms):
        bonds = [] # this is a list of tuples with the atomic indices at make up a bond
        for i, atom1 in enumerate(core_atoms):
            a1_dists = {}
            for j, atom2 in enumerate(bondable_atoms):
                if i!=j:
                    a1_dists[f'{min(i,j)},{max(i,j)}'] = self.get_atomic_distance(atom1,atom2)
            # print(a1_dists)
            bond_min = min(a1_dists.values())
            bond_max = (1+self.bond_range)*bond_min
            for k in a1_dists:
                if a1_dists[k] <= bond_max:
                    inds = k.split(',')
                    i,j = int(inds[0]),int(inds[1])
                    if not((i,j) in bonds):
                        bonds.append((i,j)) # bonds of the form (i,j)

        #now bonds contains all the bond info, time to add the info to each atom
        for i, atom in enumerate(core_atoms):
            for b in bonds:
                if i in b:
                    if max(b)<len(core_atoms): #only take bonds of core atoms
                        if b[0] == i:
                            other_atm_idx = b[1]
                        else:
                            other_atm_idx = b[0]
                        atom.bonds.append(self.bond_coords(atom,self.atoms[other_atm_idx]))
                        atom.bond_atom_sizes.append((atom.radius,self.atoms[other_atm_idx].radius))
                        atom.bond_atom_idxs.append((i,other_atm_idx))
        return bonds
    
    def read_ele_data():
        with open('ElementalRadiiData.json') as f:
            out = f.read()
        return json.loads(out)
    
    def get_atom_cart_coords(self,atom):
        u,v,w = atom.fractional_position()
        
        
        ca = np.cos(np.radians(self.alpha))
        cb = np.cos(np.radians(self.beta))
        cg = np.cos(np.radians(self.gamma))
        sg = np.sin(np.radians(self.gamma))

        if abs(sg) < 1e-12:
            raise ValueError("gamma is 0 or 180 deg (sin(gamma)=0); formula singular.")

        # lattice vectors in Cartesian (columns of matrix)
        a_x, a_y, a_z = self.a, 0.0, 0.0
        b_x = self.b * cg
        b_y = self.b * sg
        b_z = 0.0

        # c vector components (standard crystallography formula)
        c_x = self.c * cb
        c_y = self.c * (ca - cb * cg) / sg
        # compute z-component safely
        arg = 1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg
        if arg < -1e-12:
            raise ValueError("Invalid lattice angles (negative argument inside sqrt).")
        arg = max(arg, 0.0)
        c_z = self.c * math.sqrt(arg) / sg

        # fractional -> Cartesian: r = u*a + v*b + w*c
        x = u * a_x + v * b_x + w * c_x
        y = u * a_y + v * b_y + w * c_y
        z = u * a_z + v * b_z + w * c_z

        return x, y, z
    
    def get_ele_charge(self, ele_label):
        for i,s in enumerate(ele_label):
            try:
                int(s)
                num_idx = i
                break
            except ValueError:
                pass
        ele = ele_label[:num_idx]
        if ele_label.endswith('-'):
            return ele, -1*int(ele_label[num_idx:-1])
        elif ele_label.endswith('+'):
            return ele, int(ele_label[num_idx:-1])
        else:
            raise ValueError(f'Unexpected format of atom: {ele_label}')
   
    def bond_coords(self, atom1,atom2):
        """
        Convert vector Q-P from skewed basis (a,b,c) to spherical coords in an
        orthonormal basis where:
            a is along +z axis
            b is in xz plane
        Output angles in radians.
        Returns (r, theta, phi).
        """
        P = atom1.fractional_position()
        Q = atom2.fractional_position()
        # Step 1: basis vectors in orthonormal Cartesian coords
        a_vec = (0.0, 0.0, self.a)

        b_vec = (self.b * np.sin(np.radians(self.gamma)), 0.0, self.b * np.cos(np.radians(self.gamma)))

        # Find c_vec coordinates:
        # Let c_vec = (cx, cy, cz)
        cz = self.c * np.cos(np.radians(self.beta)) # from a · c
        cx = (self.b * self.c * np.cos(np.radians(self.alpha)) - b_vec[2] * cz) / b_vec[0]
        cy_sq = self.c**2 - cx**2 - cz**2
        if cy_sq < 0:
            cy_sq = 0  # numerical tolerance
        cy = np.sqrt(cy_sq)

        c_vec = (cx, cy, cz)

        # Step 2: relative coords in skew basis
        dcoords = (
            Q[0] - P[0],
            Q[1] - P[1],
            Q[2] - P[2]
        )

        # Step 3: convert to Cartesian
        dx = dcoords[0]*a_vec[0] + dcoords[1]*b_vec[0] + dcoords[2]*c_vec[0]
        dy = dcoords[0]*a_vec[1] + dcoords[1]*b_vec[1] + dcoords[2]*c_vec[1]
        dz = dcoords[0]*a_vec[2] + dcoords[1]*b_vec[2] + dcoords[2]*c_vec[2]

        # Step 4: spherical conversion
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        theta = np.acos(dz / r) if r != 0 else 0.0  # polar angle
        phi = math.atan2(dy, dx)                      # azimuth

        return r, theta, phi
    
    def extract_atoms_with_symmetry(self, headers, atoms, sym_ops):
        def eval_sym_op(op, x, y, z):
            coords = {'x': x, 'y': y, 'z': z}
            op = op.replace(",", "")
            ops = op.lower().split()[1:]
            result = []
            for term in ops:
                # print(term)
                val = eval(term, {}, coords)
                result.append(val%1)
            return result
        def reflect_about_cell_boundaries(atom_list, decimal_rounding=4):
            """
            Reflects atoms across unit cell boundaries (x, y, z = 0 or 1)
            to ensure periodic coverage, including face, edge, and corner reflections.

            Parameters:
                atom_list: List of tuples (label, x, y, z, occ) with fractional coords.
                decimal_rounding: Rounding precision for comparisons.

            Returns:
                Sorted list of unique atoms, including mirrored copies.
            """
            seen = set()
            for label, x, y, z, occ in atom_list:
                x = round(x, decimal_rounding)
                y = round(y, decimal_rounding)
                z = round(z, decimal_rounding)
                seen.add((label, x, y, z, occ))

            result = set(seen)

            for label, x, y, z, occ in seen:
                x_opts = [x]
                y_opts = [y]
                z_opts = [z]

                if x == 0.0:
                    x_opts.append(round(1.0, decimal_rounding))
                elif x == 1.0:
                    x_opts.append(round(0.0, decimal_rounding))
                if y == 0.0:
                    y_opts.append(round(1.0, decimal_rounding))
                elif y == 1.0:
                    y_opts.append(round(0.0, decimal_rounding))
                if z == 0.0:
                    z_opts.append(round(1.0, decimal_rounding))
                elif z == 1.0:
                    z_opts.append(round(0.0, decimal_rounding))

                # Manually generate all combinations (up to 8 total)
                for xi in x_opts:
                    for yi in y_opts:
                        for zi in z_opts:
                            if (xi, yi, zi) != (x, y, z):
                                mirrored = (label, xi, yi, zi, occ)
                                if mirrored not in seen:
                                    result.add(mirrored)

            return sorted(result)
        
        label_idx = 0
        #extract indices
        for i,h in enumerate(headers):
            if h.endswith('fract_x'):
                x_idx = i
            elif h.endswith('fract_y'):
                y_idx = i
            elif h.endswith('fract_z'):
                z_idx = i
            elif h.endswith('occupancy'):
                occ_idx = i
            elif h.endswith('type_symbol'):
                label_idx = i
        full_atoms = []
        seen = set()
        for atom in atoms:
            label = atom[label_idx]
            fx = get_num(atom[x_idx])%1
            fy = get_num(atom[y_idx])%1
            fz = get_num(atom[z_idx])%1
            occ = get_num(atom[occ_idx])
            for op in sym_ops:
                x_new, y_new, z_new = eval_sym_op(op, fx, fy, fz)
                key = (label, round(x_new, 4), round(y_new, 4), round(z_new, 4))
                if key not in seen:
                    seen.add(key)
                    full_atoms.append((label, x_new, y_new, z_new, occ))
        atoms_with_corners = reflect_about_cell_boundaries(full_atoms)
        return atoms_with_corners

    def refine_atomic_radii(self):
        bondable_atoms = []
        seen = set()
        for ix in range(self.nx): # first, add all core atoms with correct indeces
            for iy in range(self.ny):
                for iz in range(self.nz):
                    for raw_atom in self.raw_atoms: #raw atom is of the form: (label, xi, yi, zi, occ) where label is a string like 'Ru4' or 'O2-'
                        ele_name, charge = self.get_ele_charge(raw_atom[0])
                        fx, fy, fz = raw_atom[1] + ix, raw_atom[2] + iy, raw_atom[3] + iz
                        if not( (fx,fy,fz) in seen):
                            bondable_atoms.append(Atom(ele_name, charge,fx,fy,fz,raw_atom[4]))
                            seen.add((fx,fy,fz))

        for ix in range(-1,self.nx+1): #now add peripheral atoms with large indeces
            for iy in range(-1,self.ny+1):
                for iz in range(-1,self.nz+1):
                    for raw_atom in self.raw_atoms: #raw atom is of the form: (label, xi, yi, zi, occ) where label is a string like 'Ru4' or 'O2-'
                        ele_name, charge = self.get_ele_charge(raw_atom[0])
                        fx, fy, fz = raw_atom[1] + ix, raw_atom[2] + iy, raw_atom[3] + iz
                        if not( (fx,fy,fz) in seen):
                            bondable_atoms.append(Atom(ele_name, charge,fx,fy,fz,raw_atom[4]))
                            seen.add((fx,fy,fz))
        additional_bonds = self.find_bonds(self.atoms,bondable_atoms)

        for i,atom in enumerate(self.atoms):# find coord number aka number of bonds, and try to find a more accurate radii
            atom.coord_number = 0
            for bond in additional_bonds:
                if i in bond:
                    atom.coord_number += 1
            atom.radius = atom.get_radii_data(atom.element,atom.charge,coord=atom.coord_number)
            print(f'Notice: Refined radii for {atom.element}{atom.charge}')

    def __repr__(self):
        return f"Lattice a={self.a}, b={self.b}, c={self.c}, α={self.alpha}, β={self.beta}, γ={self.gamma}, {len(self.atoms)} atoms"

    def all_nodes_connected(self):
        # Build adjacency list
        adj = {}
        for a, b in self.bonds:
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)

        # Pick an arbitrary starting node
        start = next(iter(adj))
        visited = set()
        stack = [start]

        # Depth-first search
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                stack.extend(adj[node] - visited)

        # Check if all nodes were visited
        return len(visited) == len(adj)
    

ELEMENT_COLORS = {
    "H": "#FFFFFF",  # white
    "He": "#D9FFFF",  # pale cyan
    "Li": "#CC80FF",  # violet
    "Be": "#C2FF00",  # green
    "B": "#FFB5B5",  # salmon
    "C": "#909090",  # dark gray
    "N": "#3050F8",  # blue
    "O": "#FF0D0D",  # red
    "F": "#90E050",  # green
    "Ne": "#B3E3F5",  # light blue
    "Na": "#AB5CF2",  # purple
    "Mg": "#8AFF00",  # lime green
    "Al": "#BFA6A6",  # gray
    "Si": "#F0C8A0",  # beige
    "P": "#FF8000",  # orange
    "S": "#FFFF30",  # yellow
    "Cl": "#1FF01F",  # bright green
    "Ar": "#80D1E3",  # cyan
    "K": "#8F40D4",  # purple
    "Ca": "#3DFF00",  # green
    "Sc": "#E6E6E6",  # light gray
    "Ti": "#BFC2C7",  # gray
    "V": "#A6A6AB",  # gray
    "Cr": "#8A99C7",  # light blue
    "Mn": "#9C7AC7",  # purple
    "Fe": "#E06633",  # brownish orange
    "Co": "#F090A0",  # pink
    "Ni": "#50D050",  # green
    "Cu": "#C88033",  # brown
    "Zn": "#7D80B0",  # lavender
    "Ga": "#C28F8F",  # gray
    "Ge": "#668F8F",  # teal
    "As": "#BD80E3",  # purple
    "Se": "#FFA100",  # orange
    "Br": "#A62929",  # dark red
    "Kr": "#5CB8D1",  # cyan
    "Rb": "#702EB0",  # purple
    "Sr": "#00FF00",  # bright green
    "Y": "#94FFFF",  # pale cyan
    "Zr": "#94E0E0",  # pale cyan
    "Nb": "#73C2C9",  # teal
    "Mo": "#54B5B5",  # teal
    "Tc": "#3B9E9E",  # teal
    "Ru": "#248F8F",  # teal
    "Rh": "#0A7D8C",  # teal
    "Pd": "#006985",  # dark teal
    "Ag": "#C0C0C0",  # silver
    "Cd": "#FFD98F",  # peach
    "In": "#A67573",  # brownish
    "Sn": "#668080",  # gray
    "Sb": "#9E63B5",  # purple
    "Te": "#D47A00",  # orange
    "I": "#940094",  # dark purple
    "Xe": "#429EB0",  # teal
    "Cs": "#57178F",  # purple
    "Ba": "#00C900",  # green
    "La": "#70D4FF",  # pale blue
    "Ce": "#FFFFC7",  # pale yellow
    "Pr": "#D9FFC7",  # pale green
    "Nd": "#C7FFC7",  # pale green
    "Pm": "#A3FFC7",  # pale green
    "Sm": "#8FFFC7",  # pale green
    "Eu": "#61FFC7",  # pale green
    "Gd": "#45FFC7",  # pale green
    "Tb": "#30FFC7",  # pale green
    "Dy": "#1FFFC7",  # pale green
    "Ho": "#00FF9C",  # green
    "Er": "#00E675",  # green
    "Tm": "#00D452",  # green
    "Yb": "#00BF38",  # green
    "Lu": "#00AB24",  # green
    "Hf": "#4DC2FF",  # light blue
    "Ta": "#4DA6FF",  # light blue
    "W": "#2194D6",  # blue
    "Re": "#267DAB",  # blue
    "Os": "#266696",  # blue
    "Ir": "#175487",  # dark blue
    "Pt": "#D0D0E0",  # light gray
    "Au": "#FFD123",  # gold
    "Hg": "#B8B8D0",  # light gray
    "Tl": "#A6544D",  # brown
    "Pb": "#575961",  # dark gray
    "Bi": "#9E4FB5",  # purple
    "Po": "#AB5C00",  # brown
    "At": "#754F45",  # brown
    "Rn": "#428296",  # blue
    "Fr": "#420066",  # purple
    "Ra": "#007D00",  # green
    "Ac": "#70ABFA",  # light blue
    "Th": "#00BAFF",  # cyan
    "Pa": "#00A1FF",  # cyan
    "U": "#008FFF",  # blue
    "Np": "#0080FF",  # blue
    "Pu": "#006BFF",  # blue
    "Am": "#545CF2",  # blue
    "Cm": "#785CE3",  # purple
    "Bk": "#8A4FE3",  # purple
    "Cf": "#A136D4",  # purple
}

num_chars = {'1','2','3','4','5','6','7','8','9','0','.','-'}

def get_num(num_string):
    output = ''
    for s in num_string:
        if s in num_chars:
            output += s
        else:
            return float(output)
    return float(output)

def int_to_roman(num: int) -> str:
    if not (1 <= num <= 3999):
        raise ValueError("Argument must be between 1 and 3999")

    val = [
        1000, 900, 500, 400,
        100,  90,  50,  40,
        10,   9,   5,   4,
        1
    ]
    syms = [
        "M", "CM", "D", "CD",
        "C", "XC", "L", "XL",
        "X", "IX", "V", "IV",
        "I"
    ]

    roman = []
    for i, v in enumerate(val):
        count, num = divmod(num, v)
        roman.append(syms[i] * count)

    return "".join(roman)
