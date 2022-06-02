#!/usr/bin/env python3
import glob
import sys
import math
import json
import copy
import time
import CifFile
import numpy as np
import qcelemental as qcel
from itertools import permutations
sys.path.append("/Users/zoes/apps/qcp-python-app/qcp")
sys.path.append("/g/data/k96/apps/qcp/qcp")
from system import systemData


class Timer:

    def __init__(self, to_print):
        self.start = time.time()
        self.t1 = time.time()
        self.print = to_print

    def stop(self, descr=None):

        if not self.print: return

        now = time.time()
        if descr:
            print(round(now - self.t1, 2), round(now - self.start, 2), descr)
        else:
            print(round(now - self.t1, 2), round(now - self.start, 2))
        self.t1 = time.time()


### CRYSTAL STRUCTURE --------------------------------------------


def find_cif(inp_dict):
    try:
        cif_file = inp_dict["Cif"]
    except KeyError:
        print("ERROR: cannot find 'Cif' in the input file.")

    return(cif_file)


def cell_lengths(r, cif_data):

    print("Radius: {} A".format(r))

    Nx = math.ceil(r/cif_data["_cell_length_a"])+1
    Ny = math.ceil(r/cif_data["_cell_length_b"])+1
    Nz = math.ceil(r/cif_data["_cell_length_c"])+1

    return Nx, Ny, Nz


def factors_convert_fract2cartes(cif_data):
    """
    Edge vectors (a, b, c) in fractional coordinate –> (x, y, z) in Cartesian coordinate

    cos(alpha) = b*c/(|b||c|)
    cos(beta) = a*c/(|a||c|)
    cos(gamma) = a*b/(|a||b|)

    a = (a, 0, 0)
    b = (bcos(gamma), bsin(gamma), 0)
    c = (cx, cy, cz)

    x = La*u + Lb*cos(gamma)*v + Lc*cos(beta)*w
    y = Lb*sin(gamma)*v + Lc*((cos(alpha)cos(gamma) - cos(alpha))/sin(gamma))*w
    z = Lc * (sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2*cos_a*cos_b*cos_g)/sin_g)*w
    """

    # Lengths of the unit cell
    La = cif_data["_cell_length_a"]
    Lb = cif_data["_cell_length_b"]
    Lc = cif_data["_cell_length_c"]

    # Angles in the unit cell
    alpha = math.radians(cif_data["_cell_angle_alpha"])
    beta = math.radians(cif_data["_cell_angle_beta"])
    gamma = math.radians(cif_data["_cell_angle_gamma"])

    cos_a = math.cos(alpha)
    sin_a = math.sin(alpha)

    cos_b = math.cos(beta)
    sin_b = math.sin(beta)

    cos_g = math.cos(gamma)
    sin_g = math.sin(gamma)

    ax = La
    # ay = az = 0

    bx = Lb * cos_g
    by = Lb * sin_g
    # bz = 0

    cx = Lc * cos_b
    cy = Lc * (cos_a - cos_g*cos_b)/sin_g
    cz = Lc * math.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2*cos_a*cos_b*cos_g)/sin_g

    # Use the volume to check that we calculated the vectors correctly
    V = ax * by * cz

    if abs(V - cif_data["_cell_volume"]) > 0.1:
        print("WARNING: Volume calculated with the real vectors is not the same as the volume in CIF file.")

    return({"ax": ax, "ay": 0, "az": 0, "bx": bx, "by": by, "bz": 0, "cx": cx, "cy": cy, "cz": cz})


def convert_fract2carte_atom(u, v, w, factors_dict, print_=False):
    ax = factors_dict["ax"]
    bx = factors_dict["bx"]
    cx = factors_dict["cx"]
    by = factors_dict["by"]
    cy = factors_dict["cy"]
    cz = factors_dict["cz"]

    x = ax*u + bx*v + cx*w
    y = by*v + cy*w
    z = cz*w

    if print_:
        print('convert_fract2carte_atom')
        print("u, v, w", u, v, w)
        print("ax, bx, cx, by, cy, cz", ax, bx, cx, by, cy, cz)
        print('x, y, z', x, y, z)

    return x, y, z


def convert_carte2fract_atom(x, y, z, factors_dict, print_=None):
    ax = factors_dict["ax"]
    bx = factors_dict["bx"]
    cx = factors_dict["cx"]
    by = factors_dict["by"]
    cy = factors_dict["cy"]
    cz = factors_dict["cz"]

    w = z/cz
    v = (y - cy * w) / by
    u = (x - bx*v - cx*w) / ax

    if print_:
        print('convert_carte2fract_atom')
        print('x, y, z', x, y, z)
        print("ax, bx, cx, by, cy, cz", ax, bx, cx, by, cy, cz)
        print("u, v, w", u, v, w)

    return(u, v, w)


def read_cif(cif_name):
    cif = CifFile.ReadCif(cif_name)

    for data in cif:
        cif_dblock = data
        break

    cif_data = {}

    # Extract CIF data and remove the square brackets in the numbers


    cif_data["_chemical_name"] = cif_dblock["_chemical_name_systematic"]
    if cif_dblock["_chemical_name_systematic"] == "?":
        cif_data["_chemical_name"] = cif_dblock["_chemical_name_common"]

    cif_data["_chemical_formula_moiety"] = cif_dblock["_chemical_formula_moiety"]
    cif_data["_cell_length_a"] = float(cif_dblock["_cell_length_a"].replace("(", "").replace(")", ""))
    cif_data["_cell_length_b"] = float(cif_dblock["_cell_length_b"].replace("(", "").replace(")", ""))
    cif_data["_cell_length_c"] = float(cif_dblock["_cell_length_c"].replace("(", "").replace(")", ""))
    cif_data["_cell_angle_alpha"] = float(cif_dblock["_cell_angle_alpha"].replace("(", "").replace(")", ""))
    cif_data["_cell_angle_beta"] = float(cif_dblock["_cell_angle_beta"].replace("(", "").replace(")", ""))
    cif_data["_cell_angle_gamma"] = float(cif_dblock["_cell_angle_gamma"].replace("(", "").replace(")", ""))
    cif_data["_cell_volume"] = float(cif_dblock["_cell_volume"].replace("(", "").replace(")", ""))

    # Extract the symmetry operations that define the space group
    '''
    In some cases, it might be called "_space_group_symop_operation_xyz".
    In the CIF file, the symmetry-equivalent position in the xyz format look like:
    ```
    loop_
    _symmetry_equiv_pos_as_xyz
        'x,y,z'
        'y,x,2/3-z'
        '-y,x-y,2/3+z'
        '-x,-x+y,1/3-z'
        '-x+y,-x,1/3+z'
        'x-y,-y,-z'
    ```
    Except for the space group P1, these data will be repeated in a loop.
    '''
    cif_data["_symmetry_equiv_pos_as_xyz"] = []

    try:
        sym_op = cif_dblock["_symmetry_equiv_pos_as_xyz"]
    except KeyError:
        try:
            sym_op = cif_dblock["_space_group_symop_operation_xyz"]
        except KeyError:
            print("\n ERROR: Cif file does not have an item: either \"_symmetry_equiv_pos_as_xyz\" or \"_space_group_symop_operation_xyz\".")
            sys.exit()

    for xyz_op in sym_op:
        cif_data["_symmetry_equiv_pos_as_xyz"].append(xyz_op)

    # Get the fractional coordinates u, v, w (x, y, z) of the atoms
    cif_data["_atom_site_label"] = cif_dblock["_atom_site_label"]
    cif_data["_atom_site_type_symbol"] = cif_dblock["_atom_site_type_symbol"]

    cif_data["_atom_site_fract_x"] = []
    for u in cif_dblock["_atom_site_fract_x"]:
        cif_data["_atom_site_fract_x"].append(float(u.split("(")[0]))

    cif_data["_atom_site_fract_y"] = []
    for v in cif_dblock["_atom_site_fract_y"]:
        cif_data["_atom_site_fract_y"].append(float(v.split("(")[0]))

    cif_data["_atom_site_fract_z"] = []
    for w in cif_dblock["_atom_site_fract_z"]:
        cif_data["_atom_site_fract_z"].append(float(w.split("(")[0]))

    return(cif_data)


def asym_unit(cif_data):

    # Atom labels
    atom_labels = cif_data["_atom_site_type_symbol"]
    # Atom coordinates
    atom_u = cif_data["_atom_site_fract_x"]
    atom_v = cif_data["_atom_site_fract_y"]
    atom_w = cif_data["_atom_site_fract_z"]

    asym_unit = []
    asym_unit = [(atom_labels[i], atom_u[i], atom_v[i], atom_w[i]) for i in range(len(atom_labels))]

    # Move atoms into a unit cell
    asym_unit = [(atom[0], atom[1]%1.0, atom[2]%1.0, atom[3]%1.0) for atom in asym_unit]

    return(asym_unit)


def unit_cell(atoms, cif_data):
    '''
    Use symmetry operations to create the unit cell

    The CIF file consists of a few atom positions and several "symmetry operations" that indicate the other atom positions within the unit cell.
    Using these symmetry operations, create copies of the atoms until no new copies can be made.

    For each atom, apply each symmetry operation to create a new atom.
    '''

    # Symmetry operation
    sym_op = cif_data["_symmetry_equiv_pos_as_xyz"]

    imax = len(atoms)
    i = 0

    atoms_uc = []

    while i < imax:
        label, x, y, z = atoms[i]
        # Keep x, y, z as they are! Cause they will be inserted in the eval(op) later.

        for op in sym_op:
            # eval will convert the string into a 3-tuple using the current values for x, y, z
            u, v, w = eval(op)

            # Move new atom into the unit cell
            u = u % 1.0
            v = v % 1.0
            w = w % 1.0

            # Check if the new position is actually new, or already exists
            # Two atoms are on top of each other if they are less than "eps" away.
            eps = 0.01

            new_atom = True

            for atom in atoms:
                if (abs(atom[1] - u) < eps) and (abs(atom[2] - v) < eps) and (abs(atom[3] - w) < eps):
                    new_atom = False

                    # Check that this is the same atom type.
                    if atom[0] != label:
                        print("\nERROR: Invalid CIF file: atom of type %s overlaps with atom of type %s" % (atom[0], label))

            if (new_atom):
                atoms.append((label, u, v, w))

        i = i + 1
        imax = len(atoms)

    atoms_uc = atoms

    return(atoms_uc)


def supercell(atoms_uc, Nx, Ny, Nz):
    atoms_sc = []

    for atom in atoms_uc:

        label, u, v, w = atom

        for i in range(Nx):
            uu = i + u

            for j in range(Ny):
                vv = j + v

                for k in range(Nz):
                    ww = k + w
                    atoms_sc.append((label, uu, vv, ww))

    return(atoms_sc)


def convert_supercell_tocartes(atoms_sc, cif_data, Nx, Ny, Nz, print_=None):
    factors_fract2carte_dict = factors_convert_fract2cartes(cif_data)

    # Check if we have a rectangular box
    rect_box = False

    bx = factors_fract2carte_dict["bx"]
    cx = factors_fract2carte_dict["cx"]
    cy = factors_fract2carte_dict["cy"]

    eps = 0.1

    if (bx < eps) and (cx < eps) and (cy < eps):
        rect_box = True

    # Calculate the box size
    Lx = Nx * cif_data['_cell_length_a']
    Ly = Ny * cif_data['_cell_length_b']
    Lz = Nz * cif_data['_cell_length_c']

    atoms_rsc = []

    for val, atom in enumerate(atoms_sc):
        label, xf, yf, zf = atom
        (xn1, yn1, zn1) = convert_fract2carte_atom(xf, yf, zf, factors_fract2carte_dict, print_)

        # if rect_box:
        if False:
            # only changes values below zero
            xn2 = (xn1 + Lx) % Lx
            yn2 = (yn1 + Ly) % Ly
            zn2 = (zn1 + Lz) % Lz
        else:
            xn2 = xn1
            yn2 = yn1
            zn2 = zn1

        atoms_rsc.append((label, xn2, yn2, zn2))

        if print_:
            print("convert_supercell_tocartes")
            print("label, xf, yf, zf", label, xf, yf, zf)
            print("label, xn, yn, zn", label, xn1, yn1, zn1)
            print("label, xn, yn, zn", label, xn2, yn2, zn2)

    return(atoms_rsc)


def convert_tofracts(atoms_sc, cif_data):
    factors_fract2carte_dict = factors_convert_fract2cartes(cif_data)

    atoms_rsc = []

    for label, xf, yf, zf in atoms_sc:
        (xn, yn, zn) = convert_carte2fract_atom(xf, yf, zf, factors_fract2carte_dict)
        atoms_rsc.append((label, xn, yn, zn))

    return(atoms_rsc)


def finalise_supercell(atoms_rsc):
    '''
    Clean the duplicates in the supercell coordinates
    Translate the supercell to the origin
    '''

    # Make sure there is no duplicate rows in the super cell coordinates
    clean_sc, uniq_idx_1 = np.unique(atoms_rsc, return_index=True, axis=0)

    sc_coord = []
    test_sc_coord = [] # testing if there are coordinates that are similar
    sc_elem = []
    #test_sc_elem = []

    for atom in clean_sc:
        x = round(float(atom[1]), 2)
        y = round(float(atom[2]), 2)
        z = round(float(atom[3]), 2)
        sc_coord.append([float(atom[1]), float(atom[2]), float(atom[3])])
        test_sc_coord.append([x, y, z])
        sc_elem.append(atom[0])

    sc_elem = np.array(sc_elem)

    test_sc_coord = np.array(test_sc_coord)
    clean_sc, uniq_idx_2 = np.unique(test_sc_coord, return_index=True, axis=0)

    fin_sc_coord = []
    fin_sc_elem = []

    for i in uniq_idx_2:
        atom = sc_coord[i]
        fin_sc_coord.append(sc_coord[i])
        fin_sc_elem.append(sc_elem[i])

    # Find the origin of the super cell
    center_sc = (np.max(fin_sc_coord, axis=0) - np.min(fin_sc_coord, axis=0))/2

    # Translate the supercell to the origin
    fin_sc_coord -= center_sc

    # Glue the coordinate and elements together
    sc_rcoord = [(fin_sc_elem[i], fin_sc_coord[i][0], fin_sc_coord[i][1], fin_sc_coord[i][2]) for i in range(len(fin_sc_elem))]

    return(sc_rcoord)


def fract_min_max_from_cart(coords, cif_data):
    """Convert carts to fracts and get min and max in x, y and z."""

    coords = convert_tofracts(coords, cif_data)

    ua, va, wa = [], [], []
    for sym, u, v, w in coords:
        ua.append(u)
        va.append(v)
        wa.append(w)

    (minu, maxu), (minv, maxv), (minw, maxw) = min_max(ua), min_max(va), min_max(wa)

    # SHOULD BE JUST UNDER 3
    if not 2.8 < maxu - minu < 3:
        print("WARNING: maxu-minu not as expected", maxu, minu)
    elif not 2.8 < maxv - minv < 3:
        print("WARNING: maxv-minv not as expected", maxv, minv)
    elif not 2.8 < maxw - minw < 3:
        print("WARNING: maxw-minw not as expected", maxw, minw)
    return minu, maxu, minv, maxv, minw, maxw


### COORDS/FRAGS/ATOMS --------------------------------------------

angstrom2bohr = 1.88973
hartree2kjmol = 2625.4996
osvdz2srs = 1.752
osvtz2srs = 1.64


def dist_threshold(atom_list, tolerance_value):
    rcov = {}

    for i in atom_list:
        rcov[i] = {}
        for j in atom_list:
            r = (qcel.covalentradii.get(i, units="angstrom") + qcel.covalentradii.get(j, units="angstrom")) * tolerance_value

            #r = round((covalent_radius_dict[atom_i] + covalent_radius_dict[atom_j]) * tolerance_value, 5)
            rcov[i][j] = r

    return rcov


def closest_distance(x, y, z, atmList):
    """Get closest distance between atmList and x, y, z."""

    min_d = 100000
    for atm in atmList:
        dist = distance([x, y, z], [atm['x'], atm['y'], atm['z']])
        min_d = min(dist, min_d)
    return min_d


def distance(i, j):
    dist = math.sqrt((i[0] - j[0])**2 + (i[1] - j[1])**2 + (i[2] - j[2])**2)
    return dist


def centroid(mol_coord):
    """
    mol_coord is a list of tuples. Each tuple holds an atomic coordinate.
    ('O', -0.8310492693341551, -8.864856732599998, 5.019346296047775)
    centoid
    """

    x, y, z, i = 0, 0, 0, 0
    for atom in mol_coord:
        i += 1
        x += atom[1]
        y += atom[2]
        z += atom[3]
    x = x/i
    y = y/i
    z = z/i

    return (x, y, z)


def mol_centroid_in_central_unit_cell(fragments, coords, cif_data, minu, minv, minw, atoms_uc):
    """Return list of atoms whose centroid is in central unit cell."""

    factors_fract2carte_dict = factors_convert_fract2cartes(cif_data)

    new_atoms = []
    for frag in fragments:
        atoms = [coords[i] for i in frag]
        cx, cy, cz = centroid(atoms)
        u, v, w = convert_carte2fract_atom(cx, cy, cz, factors_fract2carte_dict)
        if (minu+1 <= u < minu+2) and (minv+1 <= v < minv+2) and (minw+1 <= w < minw+2):
            print(u, v, w)
            new_atoms.extend(atoms)

    if len(new_atoms) != len(atoms_uc):
        sys.exit(f"Not the same number of atoms in the unit cell ({len(atoms_uc)}) as whole atom unit cell ({len(new_atoms)}). exiting ...")

    return new_atoms


def min_max(list_):
    """Min and max of 1D array."""

    return min(list_), max(list_)


def add_centroid_frags(fragList, atmList):
    """Add centroid to each fragment."""

    for frag in fragList:
        atoms = [atmList[i] for i in frag['ids']]
        frag['cx'], frag['cy'], frag['cz'] = centroid_atmList(atoms)

    return fragList


def centroid_atmList(atmList):
    """Get centroid for list of atoms."""

    i, x, y, z = 0, 0, 0, 0
    for atm in atmList:
        x += atm["x"]
        y += atm["y"]
        z += atm["z"]
        i += 1
    return x/i, y/i, z/i


def coords_midpoint(atom_list):
    """Midpoint of all points of xyz."""

    listx, listy, listz = [], [], []
    for atm in atom_list:
        listx.append(atm["x"])
        listy.append(atm["y"])
        listz.append(atm["z"])
    return midpoint(listx), midpoint(listy), midpoint(listz)


def midpoint(list_):
    """Return midpoint between a list of values in 1D."""

    return np.min(list_) + (np.max(list_) - np.min(list_))/2


def add_two_frags_together(fragList, atm_list, frag1_id, frag2_id):
    """Combine two fragments in fragList."""

    new_id = min(frag1_id, frag2_id)
    other_id = max(frag1_id, frag2_id)
    new_fragList = fragList[:new_id] # copy up to the combined one

    new_frag = { # combined frag
        'ids': fragList[frag1_id]['ids'] + fragList[frag2_id]['ids'],
        'syms': fragList[frag1_id]['syms'] + fragList[frag2_id]['syms'],
        'grp': new_id,
        'chrg': fragList[frag1_id]['chrg'] + fragList[frag2_id]['chrg'],
        'mult': fragList[frag1_id]['mult'] + fragList[frag2_id]['mult'] - 1,
        'name': fragList[new_id]['name'],
    }

    new_frag = add_centroid_frags([new_frag], atm_list)

    new_fragList.extend(new_frag) # add new frag

    # add up to removed frag
    new_fragList.extend(fragList[new_id+1:other_id])

    # change rest of values
    for i in range(other_id+1,len(fragList)):
        fragList[i]['grp'] = i-1
        fragList[i]['name'] = f"frag{i-1}"
        new_fragList.append(fragList[i])

    for i in range(len(new_fragList)):
        if i != new_fragList[i]["grp"]:
            print(i, "does not")

    return new_fragList, new_id


def combination_smallest_distance(fragList, combinations):
    """Return the list of anion-cation pairs that has the smallest distance."""

    comb_use = None
    min_dist = 1000
    for comb in combinations:

        tot_dist = 0

        # FOR EACH CATION, ANION PAIR
        for cat, an in comb:

            tot_dist += distance(
                [fragList[cat]['cx'], fragList[cat]['cy'], fragList[cat]['cz']],
                [fragList[an]['cx'], fragList[an]['cy'], fragList[an]['cz']]
            )

        if tot_dist < min_dist:
            min_dist = tot_dist
            comb_use = comb

    return comb_use, min_dist


def pair_ions_lowest_dist(fragList, atmList):
    """Get pairing of molecules that has lowest total distance."""

    # cation/anion lists
    cations, anions = [], []
    for i in range(len(fragList)):
        if fragList[i]['chrg'] == 1:
            cations.append(i)
        elif fragList[i]['chrg'] == -1:
            anions.append(i)
        else:
            sys.exit("Only written for singly charged species. exiting ...")

    anions = list(permutations(anions)) # perms of anions
    cations = [cations] * len(anions)   # make list of lists of cations

    # make combinations
    combinations = []
    for an_list, cat_list in zip(anions, cations):
        comb = []
        for an, cat in zip(an_list, cat_list):
            comb.append([cat, an])
        combinations.append(comb)

    # pair
    comb, min_dist = combination_smallest_distance(fragList, combinations)

    # sort combinations largest val to smallest so can combine frags safely
    comb_sorted = []
    starting_frags = len(fragList)
    for i in range(starting_frags-1, -1, -1):
        for _ in comb:
            if i in _:
                comb_sorted.append(_)
                comb.remove(_)
                break

    # combine frags
    for index1, index2 in comb_sorted:
        lines = []
        fragList, newid = add_two_frags_together(fragList, atmList, index1, index2)
        # for id in fragList[newid]['ids']:
        #     atm = atmList[id]
        #     lines.append(f"{atm['sym']} {atm['x']} {atm['y']} {atm['z']}\n")
        # write_xyz_zoe(f"{index1}-{index2}.xyz", lines)

    return fragList


def central_frag_with_charge(frag_list, atmList, midpointx, midpointy, midpointz, charge=0):
    """Returns the frag_id/grp of the central fragment with charge=charge by finding the average
    distance to the midpoint for each fragment."""

    min_dist = 10000
    min_ion  = None
    for frag in frag_list:
        if frag['chrg'] == charge or charge == "any":
            dist = 0
            for id in frag['ids']:
                # DIST
                dist += distance([midpointx, midpointy, midpointz],
                                 [atmList[id]['x'], atmList[id]['y'], atmList[id]['z']])
            # AVERAGE
            dist = dist / len(frag['ids'])
            # IF SMALLEST DIST
            if dist < min_dist:
                min_dist  = dist
                min_ion   = frag['grp']
    return min_ion


def pair_mols_by_type(fragList_uc, atmList_uc, mx, my, mz, pair_mols):
    """Pair ions by type and return new fragList and central fragment ID."""

    print(f"{len(fragList_uc)} molecules in whole unit cell")

    # Pair all ions
    if pair_mols == "all":

        # Add centroid
        fragList_uc = add_centroid_frags(fragList_uc, atmList_uc)

        # Pair molecules by lowest total pairing distance
        fragList_uc = pair_ions_lowest_dist(fragList_uc, atmList_uc)

        # central frag
        center_frag_id = central_frag_with_charge(fragList_uc, atmList_uc, mx, my, mz, 0)

    # Only pair central ion pair
    elif pair_mols == "central":

        # Find central anion and cations
        center_cat = central_frag_with_charge(fragList_uc, atmList_uc, mx, my, mz, 1)
        center_an = central_frag_with_charge(fragList_uc, atmList_uc, mx, my, mz, -1)
        fragList_uc, center_frag_id = add_two_frags_together(fragList_uc, atmList_uc, center_cat, center_an)

    # Do not pair any ions
    elif pair_mols == "none":
        center_frag_id = central_frag_with_charge(fragList_uc, atmList_uc, mx, my, mz, "any")

    print(f"{len(fragList_uc)} fragments in whole unit cell")
    return fragList_uc, center_frag_id


### FIND FRAGMENTS --------------------------------------------


def pairing_atoms(coords):
    """Put atoms in boxes and loop through boxes to get atom pairs."""

    uniq_atoms = list(set([i[0] for i in coords]))
    rcov_dict  = dist_threshold(uniq_atoms, 1.2)

    # BOX LENGTH
    box_length = 4

    # GET MIN AND MAX X, Y, Z
    group_dict = {}  # STORE GROUP NUMBER
    minx, miny, minz = 1000, 1000, 1000
    maxx, maxy, maxz = -1000, -1000, -1000
    for val, (sym, x, y, z) in enumerate(coords):
        if x < minx:
            minx = x
        if x > maxx:
            maxx = x
        if y < miny:
            miny = y
        if y > maxy:
            maxy = y
        if z < minz:
            minz = z
        if z > maxz:
            maxz = z

        group_dict[val] = None

    # NUMBER OF BOXES
    Nbx = math.ceil((maxx - minx) / box_length)
    Nby = math.ceil((maxy - miny) / box_length)
    Nbz = math.ceil((maxz - minz) / box_length)

    # BOXES CONTAINER
    Box_list = {}
    for bx in range(Nbx+1):
        for by in range(Nby+1):
            for bz in range(Nbz+1):
                Box_list[bx + Nbx*by + Nbx*Nby*bz] = []
    # print("Number of boxes = ", len(Box_list))

    # ASSIGN ATOMS TO BOXES
    for id, atom in enumerate(coords):
        bx = math.floor((atom[1] - minx) / box_length)
        by = math.floor((atom[2] - miny) / box_length)
        bz = math.floor((atom[3] - minz) / box_length)
        Box_list[bx + Nbx*by + Nbx*Nby*bz].append([id, atom])

    # LOOP OVER BOXES
    group = 0
    for bx in range(Nbx):
        for by in range(Nby):
            for bz in range(Nbz):
                box_id = bx + Nbx*by + Nbx*Nby*bz # BOX ID

                # FOR ATOMS IN BOX
                for id1, atom1 in Box_list[box_id]:

                    con_vwd = False # CONNECTED TO OTHER ATOMS

                    # LOOP OVER SURROUNDING BOXES
                    for sbx in range(max(0, bx-1), min(bx+2, Nbx+1)):
                        for sby in range(max(0, by-1), min(by+2, Nby+1)):
                            for sbz in range(max(0, bz-1), min(bz+2, Nbz+1)):

                                sbox_id = sbx + Nbx*sby + Nbx*Nby*sbz

                                for id2, atom2 in Box_list[sbox_id]:

                                    # DONT DO FOR SAME ATOM
                                    if not id2 == id1:

                                        # DO THE THING THAT NEEDS DOING FOR EACH PAIR
                                        group, group_dict, con_vwd = \
                                            update_groups(id1, id2, coords, group, group_dict, con_vwd, rcov_dict)

                    # IF NELLY NO MATES
                    if not con_vwd:
                        group_dict[id1] = group
                        group += 1

    # CONVERT TO LIST OF LISTS
    groups_list = listFromGroups(group, group_dict)

    return groups_list


def update_groups(atm1_id, atm2_id, coords, ngroups, group_dict, con_vwd, rcov_dict):
    """If two atoms are within vdW's radii update group_dict."""

    # CHECK IF DIST < VDW RADII
    if distance(coords[atm1_id][1:], coords[atm2_id][1:]) < rcov_dict[coords[atm1_id][0]][coords[atm2_id][0]]:

        # FOUND A MATE
        con_vwd = True

        # IF NEITHER ATOM1 NOR ATOM2 PART OF A GROUP ADD THEM TO A NEW GROUP
        if group_dict[atm1_id] is None and group_dict[atm2_id] is None:
            group_dict[atm1_id], group_dict[atm2_id] = ngroups, ngroups
            ngroups += 1

        # IF BOTH HAVE BEEN ASSIGNED TO A DIFF GROUP
        elif not group_dict[atm1_id] is None and not group_dict[atm2_id] is None:
            if group_dict[atm1_id] != group_dict[atm2_id]:
                grp_chng = group_dict[atm2_id]

                # CHANGE ALL IN SECOND GROUP TO FIRST
                for key, value in group_dict.items():
                    if value == grp_chng:
                        group_dict[key] = group_dict[atm1_id]

        # IF ATOM2 NOT ASSIGNED
        elif not group_dict[atm1_id] is None and group_dict[atm2_id] is None:
            group_dict[atm2_id] = group_dict[atm1_id]

        # IF ATOM1 NOT ASSIGNED
        elif not group_dict[atm2_id] is None and group_dict[atm1_id] is None:
            group_dict[atm1_id] = group_dict[atm2_id]

    return ngroups, group_dict, con_vwd


def listFromGroups(groups_, group_dict_):
    """Convert my dict of groups to list of lists."""

    # INITIATE EMPTY LIST
    list_ = [[] for i in range(groups_)]

    # ADD INDEX TO LIST GROUP NUMBER
    for atm_indx, group_no in group_dict_.items():
        list_[group_no].append(atm_indx)

    # REMOVE EMPTY LISTS
    list_ = [ele for ele in list_ if ele != []]

    return list_


### MAKE SPHERE --------------------------------------------


def make_sphere_from_whole_unit_cell(fragList_uc, atmList_uc, mx, my, mz, Nx, Ny, Nz, cif_data, r_thres, dist_cutoff):
    """Loop over number of boxes that will fit in radius with half on either side of initial unit cell and duplicate
    frag."""

    fragList = copy.deepcopy(fragList_uc)
    atmList  = copy.deepcopy(atmList_uc)
    n_frags = len(fragList)
    n_atoms = len(atmList)
    # print(cif_data["_cell_length_a"], cif_data["_cell_length_b"], cif_data["_cell_length_c"])

    # print(convert_fract2carte_atom(0, 0, 0, factors_convert_fract2cartes(cif_data)))
    # print(convert_fract2carte_atom(1, 0, 0, factors_convert_fract2cartes(cif_data)))
    # print(convert_fract2carte_atom(0, 1, 0, factors_convert_fract2cartes(cif_data)))
    # print(convert_fract2carte_atom(0, 0, 1, factors_convert_fract2cartes(cif_data)))
    # print(convert_fract2carte_atom(1, 1, 1, factors_convert_fract2cartes(cif_data)))
    for nx in range(-Nx, Nx):
        for ny in range(-Ny, Nx):
            for nz in range(-Nz, Nz):
                if nx == 0 and ny == 0 and nz == 0: continue
                x, y, z = convert_fract2carte_atom(nx, ny, nz, factors_convert_fract2cartes(cif_data))
                # loop over frags
                for frag in fragList_uc:
                    n_atoms_copy = n_atoms # reset number of atoms in case last frag was not added
                    new_atoms = [] # store atoms of new_frag
                    new_atom_ids = [] # store ids as lists for frag data
                    for id in frag['ids']:
                        # new atom
                        new_atoms.append({
                            'id'  : n_atoms_copy,
                            'sym' : atmList_uc[id]['sym'],
                            'x'   : atmList_uc[id]['x'] + x,
                            'y'   : atmList_uc[id]['y'] + y,
                            'z'   : atmList_uc[id]['z'] + z,
                            'nu'  : atmList_uc[id]['nu'],
                            'mas' : atmList_uc[id]['mas'],
                            'vdw' : atmList_uc[id]['vdw'],
                        })
                        new_atom_ids.append(n_atoms_copy)
                        n_atoms_copy += 1

                    # new frag
                    new_frag = {
                        'ids': new_atom_ids,
                        'syms': frag['syms'],
                        'grp': n_frags,
                        'chrg': frag['chrg'],
                        'mult': frag['mult'],
                        'name': f'frag{n_frags}',
                    }

                    # centroid
                    cx, cy, cz = centroid_atmList(new_atoms)

                    # if in sphere
                    if dist_cutoff == 'centroid':
                        dist = distance([mx, my, mz], [cx, cy, cz])

                    elif dist_cutoff == 'smallest':
                        dist = closest_distance(mx, my, mz, new_atoms)

                    if dist < r_thres: # distance of midpoint to frag
                        atmList.extend(new_atoms) # add atoms to atmList
                        n_atoms = n_atoms_copy # update number of atoms
                        fragList.append(new_frag) # add to fragList
                        n_frags += 1 # update number of frags


    print(f"{len(fragList)} fragments in sphere")

    write_xyz_atmList("sphere.xyz", atmList)

    return atmList, fragList


### JSON --------------------------------------------


def exess_mbe_template(frag_ids, frag_charges, symbols, geometry, method="RIMP2", nfrag_stop=None, basis="cc-pVDZ", auxbasis="cc-pVDZ-RIFIT", number_checkpoints=3, ref_mon=0):
    """Json many body energy exess template."""

    # FRAGS
    mons = len(frag_charges)
    # WARNING: this is not the correct number of nfrags if the trimers and tetramers have cutoffs
    total_frags = int(mons + mons-1 + (mons-1)*(mons-2)/2 + (mons-1)*(mons-2)*(mons-3)/6)

    if not nfrag_stop:
        nfrag_stop = total_frags

    # CHECKPOINTING
    ncheck = number_checkpoints + 1
    ncheck = int((nfrag_stop+ncheck)/ncheck)

    dict_ = {
        "driver"    : "energy",
        "model"     : {
            "method"        : method,
            "basis"         : basis,
            "aux_basis"     : auxbasis,
            "fragmentation" : True
        },
        "keywords"  : {
            "scf"           : {
                "niter"             : 100,
                "ndiis"             : 10,
                "dele"              : 1E-8,
                "rmsd"              : 1E-8,
                "dynamic_threshold" : 10,
                "debug"             : False,
            },
            "frag": {
                "method"                : "MBE",
                "level"                 : 4,
                "ngpus_per_group"       : 4,
                "lattice_energy_calc"   : True,
                "reference_monomer"     : ref_mon,
                "dimer_cutoff"          : 100*angstrom2bohr,
                "dimer_mp2_cutoff"      : 25*angstrom2bohr,
                "trimer_cutoff"         : 35*angstrom2bohr,
                "trimer_mp2_cutoff"     : 20*angstrom2bohr,
                "tetramer_cutoff"       : 20*angstrom2bohr,
                "tetramer_mp2_cutoff"   : 10*angstrom2bohr
            },
            "FMO": {
                "fmo_type": "CPF",
                "mulliken_approx": True,
                "esp_cutoff": 100000,
                "esp_maxit": 50
            },
            "check_rst": {
                "checkpoint": True,
                "restart": False,
                "nfrag_check": min(ncheck, total_frags),
                "nfrag_stop": min(nfrag_stop, total_frags)
            }
        },
        "molecule"  : {
            "fragments"     : {
                "nfrag"             : len(frag_charges),
                "fragid"            : frag_ids,
                "fragment_charges"  : frag_charges,
                "broken_bonds"      : [],
            },
            "symbols"       : symbols,
            "geometry"      : geometry,
        },
    }

    if number_checkpoints == 0:
        del dict_["keywords"]["check_rst"]

    return dict_


def make_json_from_frag_ids(frag_indexs, fragList, atmList, nfrag_stop=None, basis="cc-pVDZ", auxbasis="cc-pVDZ-RIFIT", number_checkpoints=1, ref_mon=0):

    symbols      = []
    frag_ids     = []
    frag_charges = []
    geometry     = []
    xyz_lines    = []
    num          = 0
    # FOR EACH FRAGMENT
    for index in frag_indexs:
        num += 1
        frag_charges.append(fragList[index]['chrg'])
        # FOR EACH ATOM OF THAT FRAG
        for id in fragList[index]['ids']:
            symbols.append(atmList[id]['sym'])
            frag_ids.append(num)
            geometry.extend([atmList[id]['x'], atmList[id]['y'], atmList[id]['z']])
            xyz_lines.append(f"{atmList[id]['sym']} {atmList[id]['x']} {atmList[id]['y']} {atmList[id]['z']}")
    # TO JSON
    json_dict = exess_mbe_template(frag_ids, frag_charges, symbols, geometry, "RIMP2", nfrag_stop, basis, auxbasis, number_checkpoints, ref_mon=ref_mon)
    json_lines = format_json_input_file(json_dict)

    return json_lines, xyz_lines


def format_json_input_file(dict_):
    """Put 5 items per line and 3 coords per line."""

    # GET JSON LINES
    lines = json.dumps(dict_, indent=4)
    # COMPACT LISTS
    newlines = []
    list_ = False
    geometry_ = False
    list_lines =  []
    for line in lines.split('\n'):

        if "]" in line and not '[]' in line:
            list_ = False

            # LIST OF STRINGS - 5 PER LINE
            if not geometry_:
                newline = ""
                for i in range(len(list_lines)):
                    if i % 5 == 0:
                        newline += list_lines[i]
                    elif i % 5 == 4:
                        newline += " " + list_lines[i].strip()
                        newlines.append(newline)
                        newline = ""
                    else:
                        newline += " " + list_lines[i].strip()
                newlines.append(newline)
                newline = ""

            # LIST OF NUMBERS THREE PER LINE
            else:
                newline = ""
                for i in range(len(list_lines)):
                    if i % 3 == 0:
                        newline += list_lines[i]
                    elif i % 3 == 2:
                        newline += " " + list_lines[i].strip()
                        newlines.append(newline)
                        newline = ""
                    else:
                        newline += " " + list_lines[i].strip()
                newlines.append(newline)
                newline = ""

            list_lines = []
            geometry_ = False

        if ": [" in line and not '[]' in line:
            newlines.append(line)
            list_ = True
            if "geometry" in line:
                geometry_ = True

        elif list_:
            list_lines.append(line)

        else:
            newlines.append(line)

    return newlines


### WRITE TO FILE --------------------------------------------


def write_xyz(output_fname, atom_list):
    """
    Write a coordinate file in the standard xyz format.

    output_fname := the name of the xyz file
        E.g. "aspirin.xyz"

    atom_list := the list of atoms and their coordinates
    It is a list of tuple.
    E.g.
    [('O', 0.62355, 0.14194, 0.94663),
     ('H', 0.574, 0.031, 0.9354),
     ('H', 0.863, 0.7548, 0.9991)]
    """

    with open(output_fname, "w") as of:
        of.write("{}".format(len(atom_list)))
        of.write("\n\n")

        for atom in atom_list:
            of.write("{sym:<5} {x:>15.10f} {y:>15.10f} {z:>15.10f}\n".format(sym=atom[0], x=atom[1], y=atom[2], z=atom[3]))


def write_xyz_atmList(filename, atmList):
    """Write atmList as xyz file."""

    lines = []
    for atm in atmList:
        lines.append(f"{atm['sym']} {atm['x']} {atm['y']} {atm['z']}")
    write_xyz_zoe(filename, lines)


def write_xyz_zoe(filename, lines):
    """Write lines to xyz file."""

    with open(filename, 'w') as w:
        w.write(f"{len(lines)}\n\n")
        for line in lines:
            w.write(line+"\n")


def write_central_frag(fragList, atmList, center_ip_id, mx, my, mz):
    """WRITE XYZS, MIDPOINT AND CENTRAL IP TO XYZ"""

    lines = []

    for val, atm in enumerate(atmList):
        # WRITE
        if val in fragList[center_ip_id]['ids']:
            lines.append(f"N {atm['x']} {atm['y']} {atm['z']}")
        else:
            lines.append(f"H {atm['x']} {atm['y']} {atm['z']}")
    lines.append(f"Cl {mx} {my} {mz}")

    write_xyz_zoe("central.xyz", lines)


def write_file(filename, lines):
    """Write any filetype given as list of lines."""

    with open(filename, 'w') as w:
        for line in lines:
            w.write(line + '\n')


### MAIN --------------------------------------------


def main(cif_file, r=100, debug=False, dist_cutoff='smallest', pair_mols="all"):

    timer = Timer(to_print=debug)

    # Step 2: Read into the CIF file and extract data from it into a dictionary
    cif_data = read_cif(cif_file)
    timer.stop("Step 2 - read cif")

    # Step 3: Get threshold values
    Nx, Ny, Nz = cell_lengths(r, cif_data)
    timer.stop("Step 3 - threshold values")

    # Step 4: Get the asymmetric unit from CIF dictionary
    atoms = asym_unit(cif_data)
    timer.stop("Step 4 - asymmetric unit")

    # Step 5: Create the unit cell with symmetry operations
    atoms_uc = unit_cell(atoms, cif_data)
    timer.stop("Step 5 - unit cell")

    # Create a 3x3x3 unit
    atoms_333_f = supercell(atoms_uc, 3, 3, 3)
    atoms_333_m = convert_supercell_tocartes(atoms_333_f, cif_data, 3, 3, 3)
    atoms_333_c = finalise_supercell(atoms_333_m)
    if debug:
        # unit cell
        atoms_uc_m = convert_supercell_tocartes(atoms_uc, cif_data, 1, 1, 1)
        write_xyz("unit_cell.xyz", atoms_uc_m)
        # 3x3x3
        write_xyz("unit_333.xyz", atoms_333_m)
        write_xyz("unit_333_clean.xyz", atoms_333_c)

    timer.stop("Create 3x3x3")

    # Find fragments in 3x3x3
    fragments = pairing_atoms(atoms_333_c)
    timer.stop("Fragment system")

    # Get min of coords
    minu, maxu, minv, maxv, minw, maxw = fract_min_max_from_cart(atoms_333_c, cif_data)
    if debug:
        print("Length of x, y, z in fractional coordinates", maxu-minu, maxv-minv, maxw-minw)

    # Find molecules with centroid in central cell and write to file
    atoms_whole_uc = mol_centroid_in_central_unit_cell(fragments, atoms_333_c, cif_data, minu, minv, minw, atoms_uc)
    write_xyz("atoms_whole_uc.xyz", atoms_whole_uc)

    # Read in system data
    fragList_uc, atmList_uc, totChrg, totMult = systemData("", "atoms_whole_uc.xyz", True)

    # Get midpoint of xyz
    mx, my, mz = coords_midpoint(atmList_uc)
    print("Midpoint:", mx, my, mz)

    # Pair ions
    fragList_uc, center_frag_id = pair_mols_by_type(fragList_uc, atmList_uc, mx, my, mz, pair_mols)
    print("Central fragment ID:", center_frag_id)

    # Translate unit cell
    atmList, fragList = make_sphere_from_whole_unit_cell(fragList_uc, atmList_uc, mx, my, mz, Nx, Ny, Nz, cif_data,
                                                         r, dist_cutoff)

    if debug:
        write_central_frag(fragList, atmList, center_frag_id, mx, my, mz)

    # Create overall json
    nfrag_stop = 500
    json_lines, xyz_lines = make_json_from_frag_ids(list(range(len(fragList))), fragList, atmList, ref_mon=center_frag_id, nfrag_stop=nfrag_stop)
    write_xyz_zoe("sphere.xyz", xyz_lines)
    write_file("sphere.json", json_lines)


if __name__ == "__main__":
    # dist_cutoff: 'smallest' or 'centroid'
    # pair_ions: 'all' or 'central' or 'none'
    try:
        r = float(sys.argv[1])
    except:
        r = 100

    try:
        cif_file = sys.argv[2]
    except:
        cif_file = glob.glob("*cif")[0]

    try:
        pair_mols = sys.argv[3]
    except:
        pair_mols = 'all'

    print("Sphere size:", r)
    print("Cif file:", cif_file)
    main(cif_file, r, debug=True, dist_cutoff='centroid', pair_mols=pair_mols)
