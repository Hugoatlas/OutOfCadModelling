import math
import numpy as np


def reverse_array(array):
    res = np.zeros(np.shape(array))
    n = len(array)
    for i in range(0, n):
        res[i] = array[n - i - 1]
    return res


def slice_solid(y_slice, X, Y, Su, Sl):
    dim = np.shape(Y)
    X_SLICE = np.zeros(np.shape(X[0]))
    Y_SLICE = np.zeros(np.shape(Y[0]))
    Su_SLICE = np.zeros(np.shape(Su[0]))
    Sl_SLICE = np.zeros(np.shape(Sl[0]))
    statements = np.zeros(np.shape(Y[0]))
    for i in range(0, dim[0] - 1):
        between = False
        for j in range(0, dim[1]):
            if (Y[i, j] <= y_slice <= Y[i + 1, j]) or (Y[i, j] >= y_slice >= Y[i + 1, j]):
                between = True
            if between:
                y_0 = Y[i, j]
                y_1 = Y[i + 1, j]
                ratio = (y_slice - y_0) / (y_1 - y_0)
                X_SLICE[j] = X[i, j] + ratio * (X[i + 1, j] - X[i, j])
                Y_SLICE[j] = Y[i, j] + ratio * (Y[i + 1, j] - Y[i, j])
                Su_SLICE[j] = Su[i, j] + ratio * (Su[i + 1, j] - Su[i, j])
                Sl_SLICE[j] = Sl[i, j] + ratio * (Sl[i + 1, j] - Sl[i, j])
                statements[j] = 1
    X_SLICE_FINAL = []
    Y_SLICE_FINAL = []
    Su_SLICE_FINAL = []
    Sl_SLICE_FINAL = []
    for i in range(0, len(statements)):
        if statements[i] == 1:
            X_SLICE_FINAL.append(X_SLICE[i])
            Y_SLICE_FINAL.append(Y_SLICE[i])
            Su_SLICE_FINAL.append(Su_SLICE[i])
            Sl_SLICE_FINAL.append(Sl_SLICE[i])
    return X_SLICE_FINAL, Y_SLICE_FINAL, Su_SLICE_FINAL, Sl_SLICE_FINAL


def check_all_intersections(profiles_x, profiles_z):
    nb_profiles = len(profiles_x)
    intersections_list = []
    for i in range(0, nb_profiles):
        for j in range(i + 1, nb_profiles):
            is_intersecting = check_planar_intersection(profiles_x[i], profiles_z[i], profiles_x[j], profiles_z[j])
            if is_intersecting:
                intersections_list.append([i, j])
    return intersections_list


def check_planar_intersection(X1, Z1, X2, Z2):
    n1 = len(X1)
    n2 = len(X2)
    is_intersecting = False
    if n1 == len(Z1) and n2 == len(Z2):
        for i in range(0, n1 - 1):
            if is_intersecting:
                break
            for j in range(0, n2 - 1):
                is_intersecting, intersection = get_intersection(X1[i], Z1[i], X1[i + 1], Z1[i + 1],
                                                                 X2[j], Z2[j], X2[j + 1], Z2[j + 1])
                if is_intersecting:
                    break
    return is_intersecting


def get_intersection(x1, z1, x2, z2, x3, z3, x4, z4):
    found = False
    px = ((x1 * z2 - z1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * z4 - z3 * x4)) / (
            (x1 - x2) * (z3 - z4) - (z1 - z2) * (x3 - x4))
    pz = ((x1 * z2 - z1 * x2) * (z3 - z4) - (z1 - z2) * (x3 * z4 - z3 * x4)) / (
            (x1 - x2) * (z3 - z4) - (z1 - z2) * (x3 - x4))
    intersection = [px, pz]
    if (x1 <= px <= x2 or x1 >= px >= x2) and (z1 <= pz <= z2 or z1 >= pz >= z2):
        found = True
    return found, intersection


def cleanup_airfoil(x_partial, y_partial, su_partial, sl_partial):
    x_upper = x_partial
    y_upper = y_partial
    z_upper = su_partial
    if len(x_partial) != 0:
        if su_partial[-1] == sl_partial[0]:
            x_upper = x_upper[:len(x_upper) - 1]
            y_upper = y_upper[:len(y_upper) - 1]
            z_upper = z_upper[:len(z_upper) - 1]
        x_lower = reverse_array(x_partial)
        y_lower = reverse_array(y_partial)
        z_lower = reverse_array(sl_partial)
        x_lower = x_lower[:len(x_lower) - 1]
        y_lower = y_lower[:len(y_lower) - 1]
        z_lower = z_lower[:len(z_lower) - 1]
        profile_x = np.append(x_lower, x_upper)
        profile_y = np.append(y_lower, y_upper)
        profile_z = np.append(z_lower, z_upper)
    else:
        profile_x = []
        profile_y = []
        profile_z = []
    return profile_x, profile_y, profile_z


def group_airfoils(x_slice, y_slice, su_slice, sl_slice,
                   x_slices, y_slices, su_slices, sl_slices):
    if len(x_slice) != 0:
        x_slices.append(x_slice)
        y_slices.append(y_slice)
        su_slices.append(su_slice)
        sl_slices.append(sl_slice)
    return x_slices, y_slices, su_slices, sl_slices


def gmsh_export(filename, x_partials, y_partials, su_partials, sl_partials, hs):
    nb_profiles = len(x_partials)
    profiles_x = []
    profiles_y = []
    profiles_z = []
    cord = 0
    quarter_cord_x = 0
    quarter_cord_y = 0
    quarter_cord_z = 0
    for i in range(0, nb_profiles):
        x_partial = x_partials[i]
        y_partial = y_partials[i]
        su_partial = su_partials[i]
        sl_partial = sl_partials[i]
        profile_x, profile_y, profile_z = cleanup_airfoil(x_partial, y_partial, su_partial, sl_partial)
        profiles_x.append(profile_x)
        profiles_y.append(profile_y)
        profiles_z.append(profile_z)

        # Get cord and quarter cord information
        quarter_cord_x_new = round((0.75 * x_partial[0]) + (0.25 * x_partial[len(x_partial) - 1]), 6)
        quarter_cord_y_new = round((0.75 * y_partial[0]) + (0.25 * y_partial[len(x_partial) - 1]), 6)
        quarter_cord_z_new = round((0.75 * (su_partial[0] + sl_partial[0]) / 2) + (
                0.25 * (su_partial[len(su_partial) - 1] + sl_partial[len(sl_partial) - 1]) / 2), 6)
        cord_x = (x_partial[0] - x_partial[len(x_partial) - 1])
        cord_z = (su_partial[0] + sl_partial[0] - su_partial[len(su_partial) - 1] - sl_partial[len(sl_partial) - 1]) / 2
        cord_new = round(math.sqrt(math.pow(cord_x, 2) + math.pow(cord_z, 2)), 6)

        # Keep largest cord
        if cord_new > cord:
            cord = cord_new
            quarter_cord_x = quarter_cord_x_new
            quarter_cord_y = quarter_cord_y_new
            quarter_cord_z = quarter_cord_z_new

    # Print airfoils to .geo file for mesh creation
    print_geo_file(filename, profiles_x, profiles_y, profiles_z, hs,
                   quarter_cord_x, quarter_cord_y, quarter_cord_z, cord)


def print_geo_file(filename, profiles_x, profiles_y, profiles_z, hs,
                   quarter_cord_x, quarter_cord_y, quarter_cord_z, cord):
    f = open(filename, "w")
    # Set factory to OpenCASCADE
    f.write('SetFactory("OpenCASCADE");\n')

    nb_profiles = len(profiles_x)
    current_index = 0
    curve_loops_index = 1
    plane_surface_index = 1
    physical_curve_index = 1
    physical_surface_name = '"internal"'

    h1 = hs[0] * cord
    h2 = hs[1] * cord

    # # Check if airfoils intersect, and which ones do
    # intersections_list = check_all_intersections(profiles_x, profiles_z)

    curve_IDs = []
    surface_IDs = []
    physical_curve_IDs = []
    for i in range(0, nb_profiles):
        profile_x = profiles_x[i]
        profile_y = profiles_y[i]
        profile_z = profiles_z[i]
        if len(profile_x) == len(profile_z):
            nb_points = len(profile_z)

            # Write points to .geo file
            for j in range(current_index, current_index + nb_points):
                string = 'Point(' + str(j + 1) + ') = {' + str(round(profile_x[j - current_index], 6)) + ', ' + str(
                    round(profile_z[j - current_index], 6)) + ', ' + str(
                    round(profile_y[j - current_index], 6)) + ', ' + str(h1) + '};\n'
                f.write(string)

            # Write lines connecting points to .geo file
            for j in range(current_index, current_index + nb_points - 1):
                string = 'Line(' + str(j + 1) + ') = {' + str(j + 1) + ', ' + str(j + 2) + '};\n'
                f.write(string)
            # Add last connecting point to close profile
            string = 'Line(' + str(nb_points + current_index) + ') = {' + str(nb_points + current_index) + ', ' + str(
                1 + current_index) + '};\n'
            f.write(string)
            f.write('\n')

            # Write Curve Loop
            string = '//+\nCurve Loop(' + str(curve_loops_index) + ') = {' + str(1 + current_index)
            for j in range(1 + current_index, current_index + nb_points):
                string += ', ' + str(j + 1)
            string += '};\n'
            f.write(string)
            curve_IDs.append(curve_loops_index)
            curve_loops_index += 1

            # Write Plane Surfaces
            string = '//+\nPlane Surface(' + str(plane_surface_index) + ') = {' + str(curve_loops_index - 1) + '};\n'
            f.write(string)
            surface_IDs.append(plane_surface_index)
            plane_surface_index += 1

            # # Write physical curves
            # string = '//+\nPhysical Curve(' + str(loops_index) + ') = {' + str(1 + current_index)
            # for j in range(1 + current_index, current_index + nb_points):
            #     string += ', ' + str(j + 1)
            # string += '};\n\n'
            # f.write(string)
            # physical_curve_IDs.append(loops_index)
            # loops_index += 1

            # Update current index for next airfoil
            current_index += nb_points

    # Write 3 points for circular domain creation
    index = current_index + 1
    string = '//+\nPoint(' + str(index) + ') = {' + str(quarter_cord_x) + ', ' + str(
        quarter_cord_z) + ', ' + str(quarter_cord_y) + ', ' + str(h2) + '};\n'
    f.write(string)
    string = 'Point(' + str(index + 1) + ') = {' + str(quarter_cord_x) + ', ' + str(
        quarter_cord_z + 50 * cord) + ', ' + str(quarter_cord_y) + ', ' + str(h2) + '};\n'
    f.write(string)
    string = 'Point(' + str(index + 2) + ') = {' + str(quarter_cord_x) + ', ' + str(
        quarter_cord_z - 50 * cord) + ', ' + str(quarter_cord_y) + ', ' + str(h2) + '};\n'
    f.write(string)
    f.write('\n')

    # Write circles
    string = '//+\nCircle(' + str(index) + ') = {' + str(index + 2) + ', ' + str(
        index) + ', ' + str(index + 1) + '};\n'
    f.write(string)
    string = '//+\nCircle(' + str(index + 1) + ') = {' + str(index + 1) + ', ' + str(
        index) + ', ' + str(index + 2) + '};\n'
    f.write(string)

    # Write circle curve loop
    string = '//+\nCurve Loop(' + str(curve_loops_index) + ') = {' + str(index + 1) + ', ' + str(
        index) + '};\n'
    f.write(string)
    domain_curve_ID = curve_loops_index
    curve_loops_index += 1

    # Write Plane Surfaces of circle
    string = '//+\nPlane Surface(' + str(plane_surface_index) + ') = {' + str(domain_curve_ID) + '};\n'
    f.write(string)
    domain_surface_ID = plane_surface_index
    circle_domain_ID = plane_surface_index
    plane_surface_index += 1

    # Write physical curve of circle
    string = '//+\nPhysical Curve(' + str(physical_curve_index) + ') = {' + str(index + 1) + ', ' + str(
        index) + '};\n'
    f.write(string)
    physical_curve_index += 1

    # Fuse airfoils with boolean operations
    nb_surfaces = len(surface_IDs)
    if nb_surfaces > 1:
        string = 'BooleanUnion(' + str(
            plane_surface_index) + ') = { Surface{' + str(surface_IDs[0]) + '}; } { Surface{' + str(surface_IDs[1])
        for i in range(1, nb_surfaces - 1):
            ID = surface_IDs[i]
            string += ', ' + str(int(ID))
        string += '}; };\n'
        f.write(string)
        plane_surface_index += 1

    # Write physical curve domain of airfoils
    string = 'Physical Curve(' + str(physical_curve_index) + ') = { Boundary { Surface{' + \
             str(plane_surface_index - 1) + '};} };\n'
    f.write(string)
    physical_curve_index += 1

    # Remove all airfoils from domain
    nb_surfaces = len(surface_IDs)
    if nb_surfaces != 0:
        string = 'BooleanDifference(' + str(plane_surface_index) + ') = { Surface{' + str(
            domain_surface_ID) + '};}{ Surface{' + str(surface_IDs[0])
        for i in range(1, nb_surfaces):
            ID = surface_IDs[i]
            string += ', ' + str(int(ID))
        string += '}; Delete;};\n'
        f.write(string)
        domain_surface_ID = plane_surface_index
        plane_surface_index += 1

    # Write physical surface
    string = 'Physical Surface(' + str(physical_surface_name) + ') = {' + str(domain_surface_ID) + '};\n'
    f.write(string)

    # # Remove redundant circle domain
    # f.write('\n//+\nRecursive Delete {\n')
    # f.write('    Surface {' + str(circle_domain_ID) + '};\n')
    # f.write('}\n')

    f.close()
