import math
import numpy as np
import matplotlib.pyplot as plt
import cst_geometry as geo


class Mesh:

    class Vortex:
        def __init__(self, vertices, surface, station, row):
            self.type = "Vortex"
            self.id = None
            self.station = station
            self.row = row
            self.surface = surface
            # A vortex always has 4 vertices
            if len(vertices) == 4:
                self.vertices = vertices

    class Doublet:
        def __init__(self, vertices, surface):
            self.type = "Doublet"
            self.id = None
            self.surface = surface
            self.vertices = vertices

    class Surface:
        def __init__(self, surface_id, surface_type, element_count):
            # Two surfaces can have the same ID if they do not have the same type
            self.id = surface_id
            self.type = surface_type
            self.element_count = element_count

    def __init__(self):
        self.vertices = []
        self.elements = []
        self.surfaces = []
        self.vortex_count = 0
        self.doublet_count = 0
        self.surface_count = 0

    def append_vertex(self, vertex):
        if len(vertex) == 3:
            self.vertices.append(vertex)
        else:
            print("ERROR - Incorrect Vertex Dimensions")

    def append_element(self, element):
        # Input element has to be either a Vortex object or a Doublet object. Treat both cases individually
        if element.type == "Vortex":
            # Assign a unique vortex_id to element, then update mesh vortex_count
            element.id = self.vortex_count
            self.vortex_count += 1
            # Append Vortex object to mesh elements list
            self.elements.append(element)

        elif element.type == "Doublet":
            # Assign a unique doublet_id to element, then update mesh doublet_count
            element.id = self.doublet_count
            self.doublet_count += 1
            # Append Doublet object to mesh elements list
            self.elements.append(element)

        else:
            print("ERROR - Incorrect Element Format")

    def update_surfaces(self):
        # Re-initialize element count in mesh surfaces
        for surface in self.surfaces:
            surface.element_count = 0

        for element in self.elements:
            is_new_surface = True
            for surface in self.surfaces:
                if element.surface == surface.id and element.type == surface.type:
                    surface.element_count += 1
                    # Then surface is already initialized, break out of for loop
                    is_new_surface = False
                    break

            if is_new_surface:
                # Create new surface in mesh and initialize element count to 1
                self.surfaces.append(Mesh.Surface(element.surface, element.type, 1))
                self.surface_count += 1


def get_body_camber(su, sl):
    camber = su.copy()
    if len(su) == len(sl):
        for i in range(0, len(su)):
            if len(su[i]) == len(sl[i]):
                for j in range(0, len(su[i])):
                    camber[i][j] = (su[i][j] + sl[i][j]) / 2
    return camber


def cross_product(u, v):
    n = np.zeros(np.shape(u))
    n[0] = u[1] * v[2] - u[2] * v[1]
    n[1] = u[2] * v[0] - u[0] * v[2]
    n[2] = u[0] * v[1] - u[1] * v[0]
    return n


def dot_product(u, v):
    res = 0
    if len(u) == len(v):
        for i in range(0, len(u)):
            res += u[i]*v[i]
    return res


def determinant(a, b, c):
    det = dot_product(a, cross_product(b, c))
    return det


def intersect_line_with_cst_body(p0, p1, param, precision):
    vec_x = np.linspace(p0[0], p1[0], precision)
    vec_y = np.linspace(p0[1], p1[1], precision)
    vec_z = np.linspace(p0[2], p1[2], precision)

    su = None
    sl = None
    res_prev = None
    su_prev = su
    sl_prev = sl
    point_prev_is_inside = -1
    intersect_list = []
    line_is_fully_inside_body = True
    p0_is_inside = False
    p1_is_inside = False

    # Run it for all segments between points
    for i in range(0, precision):

        res = geo.get_domain_position_xy(vec_x[i], vec_y[i], param)

        if res is None:
            # Point is then necessarily outside of body
            point_is_inside = -1
            line_is_fully_inside_body = False

        else:
            _, _, su_, sl_ = geo.build_body_surfaces([res[0]], [res[1]], param)
            su = su_[0][0]
            sl = sl_[0][0]
            # Check if point is located between CST surfaces
            if su >= vec_z[i] >= sl:
                point_is_inside = 1
                if i == 0:
                    p0_is_inside = True
                elif i == precision - 1:
                    p1_is_inside = True
            else:
                point_is_inside = -1
                line_is_fully_inside_body = False

        if (i != 0) and (point_is_inside * point_prev_is_inside < 0):
            # Then the line segment is intersecting the CST body, assuming said body is closed
            # The following code verifies on which surface the intersection manifests, and calculates the coordinates
            found_intersection = False
            t = 0

            # If res or res_prev are None, we need to find the closest boundary point of the CST domain
            if res is None:
                # domain, _ = geo.get_full_domain_boundary(25, 25, param)

                new_res = geo.get_domain_boundary(vec_x[i - 1], vec_y[i - 1],
                                                  vec_x[i] - vec_x[i - 1], vec_y[i] - vec_y[i - 1], param, 50)
                x_, y_, su_, sl_ = geo.build_body_surfaces([new_res[0]], [new_res[1]], param)
                su = su_[0][0]
                sl = sl_[0][0]

            elif res_prev is None:
                # domain, _ = geo.get_full_domain_boundary(25, 25, param)

                new_res = geo.get_domain_boundary(vec_x[i], vec_y[i],
                                                  vec_x[i - 1] - vec_x[i], vec_y[i - 1] - vec_y[i], param, 50)
                x_prev_, y_prev_, su_prev_, sl_prev_ = geo.build_body_surfaces([new_res[0]], [new_res[1]], param)
                su_prev = su_prev_[0][0]
                sl_prev = sl_prev_[0][0]

            # Check if intersection is on upper surface of CST body
            denominator_su = su - vec_z[i] - (su_prev - vec_z[i - 1])
            if denominator_su != 0:
                t = (vec_z[i - 1] - su_prev) / denominator_su
                if 0 <= t <= 1:
                    # Then there is an intersection
                    found_intersection = True

            # Check if intersection is on upper surface of CST body, only if no intersection was found on upper surface
            denominator_sl = sl - vec_z[i] - (sl_prev - vec_z[i - 1])
            if denominator_sl != 0 and not found_intersection:
                t = (vec_z[i - 1] - sl_prev) / denominator_sl
                if 0 <= t <= 1:
                    # Then there is an intersection
                    found_intersection = True

            if found_intersection:
                intersect_x = vec_x[i - 1] + t * (vec_x[i] - vec_x[i - 1])
                intersect_y = vec_y[i - 1] + t * (vec_y[i] - vec_y[i - 1])
                intersect_z = vec_z[i - 1] + t * (vec_z[i] - vec_z[i - 1])
                intersect_list.append([intersect_x, intersect_y, intersect_z])

        # Set res_prev for next iteration
        res_prev = res
        su_prev = su
        sl_prev = sl
        point_prev_is_inside = point_is_inside

    # Return which points are inside of body, if both value is 2
    if p0_is_inside and p1_is_inside:
        point_inside = 2
    elif p1_is_inside:
        point_inside = 1
    elif p0_is_inside:
        point_inside = 0
    else:
        point_inside = None

    return line_is_fully_inside_body, intersect_list, point_inside


def check_intersections_with_cst_body(elem_table, vertex_table, param, precision):

    dim = np.shape(elem_table)
    elem_count = dim[0]
    remove_elements_list = []
    elements_intersecting_body = []
    move_vertices = []
    vertices_destination = []
    for i in range(0, elem_count):
        remove_elem = True

        vertices = elem_table[i, :]
        vertices_count = len(vertices)
        for j in range(0, vertices_count):
            q0 = vertex_table[int(vertices[j]), :]
            if j == vertices_count - 1:
                q1 = vertex_table[int(vertices[0]), :]
            else:
                q1 = vertex_table[int(vertices[j + 1]), :]

            line_is_within_body, intersects_list, point_inside = intersect_line_with_cst_body(q0, q1, param, precision)

            if not line_is_within_body:
                remove_elem = False

            if len(intersects_list) == 1:
                if point_inside == 0:
                    move_vertices.append(int(vertices[j]))
                elif point_inside == 1:
                    if j == vertices_count - 1:
                        move_vertices.append(int(vertices[0]))
                    else:
                        move_vertices.append(int(vertices[j + 1]))
                elements_intersecting_body.append(i)
                vertices_destination.append(intersects_list[0])

        # For all vertices of element, check if it is inside body AND NOT tagged to be moved.
        # If so, the element is to be removed
        #
        #

        if remove_elem:
            remove_elements_list.append(i)

    return elements_intersecting_body, move_vertices, vertices_destination, remove_elements_list


def get_outer_box_dimensions(elem_table, vertex_table):
    p_min = vertex_table[int(elem_table[0, 0]), :].copy()
    p_max = vertex_table[int(elem_table[0, 0]), :].copy()

    for vertex_list in elem_table:
        for vertex_id in vertex_list:
            point = vertex_table[int(vertex_id), :].copy()
            p_min[0] = min(point[0], p_min[0]) - 0.01
            p_min[1] = min(point[1], p_min[1]) - 0.01
            p_min[2] = min(point[2], p_min[2]) - 0.01
            p_max[0] = max(point[0], p_max[0]) + 0.01
            p_max[1] = max(point[1], p_max[1]) + 0.01
            p_max[2] = max(point[2], p_max[2]) + 0.01

    return p_min, p_max


def divide_box(p_min, p_max, int_factor):
    d = p_max - p_min
    d_ijk = d / int_factor
    p_min_list = []
    p_max_list = []
    for i in range(0, int_factor):
        for j in range(0, int_factor):
            for k in range(0, int_factor):
                p_min_ijk = np.zeros(np.shape(p_min))
                p_max_ijk = np.zeros(np.shape(p_max))

                p_min_ijk[0] = p_min[0] + (d_ijk[0] * i)
                p_min_ijk[1] = p_min[1] + (d_ijk[1] * j)
                p_min_ijk[2] = p_min[2] + (d_ijk[2] * k)
                p_max_ijk[0] = p_min_ijk[0] + d_ijk[0]
                p_max_ijk[1] = p_min_ijk[1] + d_ijk[1]
                p_max_ijk[2] = p_min_ijk[2] + d_ijk[2]

                p_min_list.append(p_min_ijk)
                p_max_list.append(p_max_ijk)

    return p_min_list, p_max_list


def generate_box_tree(elem_table, vertex_table, divisions):
    # Initialize box
    p_min, p_max = get_outer_box_dimensions(elem_table, vertex_table)

    box_tree = ["Boxes", [p_min, p_max], []]

    for i in range(0, divisions):
        p_min_list, p_max_list = divide_box(p_min, p_max, int_factor=2)


# def divide_tree_boxes(box_tree, p_min, p_max, int_factor):


def set_point_in_boxes(box_tree, tree_structure, vertex, vertex_id):
    box_tree_type = box_tree[0]
    p_tree_min = box_tree[1][0]
    p_tree_max = box_tree[1][1]
    box_tree_list = box_tree[2]

    if is_point_inside_box(vertex, p_tree_min, p_tree_max):
        if box_tree_type == "Points":
            box_tree_list.append([vertex_id, vertex])
        elif box_tree_type == "Boxes":

            for i in range(box_tree_list):
                box = box_tree_list[i]
                box_tree_list[i] = set_point_in_boxes(box, tree_structure, vertex, vertex_id)

    return box_tree


def set_points_in_box(box_tree, tree_structure, vertex_table):
    # tree_structure = [number of boxes, ]
    # box_tree = [TYPE, [pmin, pmax], [[[pmin, pmax], [list_of_points]], [[pmin, pmax], list_of_points]]]

    for vertex_id in range(0, np.shape(vertex_table)[0]):
        vertex = vertex_table[vertex_id, :]
        box_tree = set_point_in_boxes(box_tree, tree_structure, vertex, vertex_id)

    return box_tree


# def get_elements_in_box(elem_list, elem_table, vertex_table, p_min, p_max):
#     elements_in_box = []
#     for i in elem_list:
#         done_with_elem = False
#         for j in range(0, len(elem_table[i, :])):
#             vertex = vertex_table[int(elem_table[i, j]), :]
#             verif = is_point_inside_box(vertex, p_min, p_max)
#             if verif and not done_with_elem:
#                 elements_in_box.append(i)
#                 done_with_elem = True
#     return elements_in_box


def is_point_inside_box(point, p_min, p_max):
    vertex = point
    u = (vertex[0] - p_min[0]) / (p_max[0] - p_min[0])
    v = (vertex[1] - p_min[1]) / (p_max[1] - p_min[1])
    w = (vertex[2] - p_min[2]) / (p_max[2] - p_min[2])
    if (1 >= u >= 0) and (1 >= v >= 0) and (1 >= w >= 0):
        return True
    else:
        return False


# def set_elements_in_boxes(elem_table, vertex_table, int_factor, n_divisions,
#                           elements_list_by_box, p_min_list, p_max_list):
#     new_elements_list_by_box = []
#     new_p_min_list = []
#     new_p_max_list = []
#
#     if len(p_min_list) == 0:
#         p_min, p_max = get_outer_box_dimensions(elem_table, vertex_table)
#         p_min_list, p_max_list = divide_box(p_min, p_max, int_factor)
#         for i in range(0, len(p_min_list)):
#             elem_list = []
#             for j in range(0, len(elem_table[:, 0])):
#                 elem_list.append(j)
#             elements_in_box = get_elements_in_box(elem_list, elem_table, vertex_table, p_min_list[i], p_max_list[i])
#             elements_list_by_box.append(elements_in_box)
#
#             new_elements_list_by_box = elements_list_by_box
#             new_p_min_list = p_min_list
#             new_p_max_list = p_max_list
#
#     else:
#
#         for i in range(0, len(p_min_list)):
#             p_min_list_i, p_max_list_i = divide_box(p_min_list[i], p_max_list[i], int_factor)
#
#             for j in range(0, len(p_min_list_i)):
#                 elements_in_box_i = get_elements_in_box(elements_list_by_box[i], elem_table, vertex_table,
#                                                         p_min_list_i[j], p_max_list_i[j])
#
#                 new_elements_list_by_box.append(elements_in_box_i)
#                 new_p_min_list.append(p_min_list_i[j])
#                 new_p_max_list.append(p_max_list_i[j])
#
#     n_divisions -= 1
#     if n_divisions >= 0:
#         new_elements_list_by_box, new_p_min_list, new_p_max_list = set_elements_in_boxes(elem_table, vertex_table,
#                                                                                          int_factor, n_divisions,
#                                                                                          new_elements_list_by_box,
#                                                                                          new_p_min_list, new_p_max_list)
#
#     return new_elements_list_by_box, new_p_min_list, new_p_max_list


def get_element_normal(elem_id, elem_table, vertex_table):
    point_a = vertex_table[int(elem_table[elem_id, 0]), :]
    point_b = vertex_table[int(elem_table[elem_id, 1]), :]
    point_c = vertex_table[int(elem_table[elem_id, 2]), :]

    # Evaluate normal vector
    n = cross_product(point_b - point_a, point_c - point_a)

    # Normalize normal vector
    norm = math.sqrt(n[0]**2 + n[1]**2 + n[2]**2)
    if norm > 0:
        n[0] = n[0] / norm
        n[1] = n[1] / norm
        n[2] = n[2] / norm
    else:
        n[0] = 0
        n[1] = 0
        n[2] = 0
    return n


# def get_planar_intersection(n1, p1, n2, p2):
#     # Get direction of intersection line
#     n3 = cross_product(n1, n2)
#
#     # Construct matrix M
#     M = np.zeros([5, 5])
#     M[0, 0] = 2
#     M[1, 1] = 2
#     M[2, 2] = 2
#     M[3, 0] = n1[0]
#     M[3, 1] = n1[1]
#     M[3, 2] = n1[2]
#     M[4, 0] = n2[0]
#     M[4, 1] = n2[1]
#     M[4, 2] = n2[2]
#     M[0, 3] = n1[0]
#     M[1, 3] = n1[1]
#     M[2, 3] = n1[2]
#     M[0, 4] = n2[0]
#     M[1, 4] = n2[1]
#     M[2, 4] = n2[2]
#
#     # Construct b vector
#     b = np.zeros([5, 1])
#     b[0, 0] = 2 * p1[0]
#     b[1, 0] = 2 * p1[1]
#     b[2, 0] = 2 * p1[2]
#     b[3, 0] = dot_product(n1, p1)
#     b[4, 0] = dot_product(n2, p2)
#
#     # Solve system of equations to obtain x
#     x = np.matmul(np.linalg.inv(M), b)
#     p3 = np.zeros(np.shape(n1))
#     p3[0] = x[0, 0]
#     p3[1] = x[1, 0]
#     p3[2] = x[2, 0]
#
#     return n3, p3


# def is_line_within_triangle(n, p, a, b, c):
#     line_is_thru = False
#     p0 = p
#     p1 = p0 + n
#     det2 = determinant(p0, p1, a)
#     det3 = determinant(p0, p1, b)
#     det1 = determinant(p0, p1, c)
#     if det1 * det2 < 0 and det2 * det3 > 0:
#         line_is_thru = True
#     elif det2 * det1 < 0 and det1 * det3 > 0:
#         line_is_thru = True
#     elif det3 * det1 < 0 and det1 * det2 > 0:
#         line_is_thru = True
#
#     return line_is_thru


def intersect_line_triangle(q1, q2, p1, p2, p3):
    def signed_tetra_volume(a, b, c, d):
        return np.sign(np.dot(np.cross(b - a, c - a), d - a) / 6.0)

    s1 = signed_tetra_volume(q1, p1, p2, p3)
    s2 = signed_tetra_volume(q2, p1, p2, p3)

    if s1 != s2:
        s3 = signed_tetra_volume(q1, q2, p1, p2)
        s4 = signed_tetra_volume(q1, q2, p2, p3)
        s5 = signed_tetra_volume(q1, q2, p3, p1)
        if s3 == s4 and s4 == s5:
            n = np.cross(p2 - p1, p3 - p1)
            t = np.dot(p1 - q1, n) / np.dot(q2 - q1, n)
            return q1 + t * (q2 - q1)
    return None


# def check_all_intersections_with_body(elem_table_1, vertex_table_1, elem_table_2, vertex_table_2,
#                                       int_factor, n_divisions):
#
#     elements_list_by_box, p_min_list, p_max_list = set_elements_in_boxes(elem_table_2, vertex_table_2,
#                                                                          int_factor, n_divisions, [], [], [])
#
#     dim = np.shape(elem_table_1)
#     elem_count = dim[0]
#     intersects_list = []
#     intersecting_list = []
#     for i in range(0, elem_count):
#         intersects = check_intersections_with_body(i, elem_table_1, vertex_table_1, elem_table_2, vertex_table_2,
#                                                    elements_list_by_box, p_min_list, p_max_list)
#         intersects_list.append(intersects)
#         if len(intersects) > 0:
#             intersecting_list.append(i)
#     return intersects_list, intersecting_list


# def check_intersections_with_body(k, elem_table_1, vertex_table_1, elem_table_2, vertex_table_2,
#                                   elements_list_by_box, p_min_list, p_max_list):
#     elem_count = len(elem_table_2[:, 0])
#     check_elem_for_intersection = np.zeros([elem_count, 1])
#
#     vertices = elem_table_1[k, :]
#     vertices_count = len(vertices)
#     intersecting_elements = []
#
#     for i in range(0, len(p_min_list)):
#         for m in range(0, vertices_count):
#             # Need to check for lines as well as the vertices
#             point = vertex_table_1[int(vertices[m]), :]
#             point_is_inside_box = is_point_inside_box(point, p_min_list[i], p_max_list[i])
#             if point_is_inside_box:
#                 for elem_id in elements_list_by_box[i]:
#                     check_elem_for_intersection[elem_id] = 1
#
#     for i in range(0, vertices_count):
#         q0 = vertex_table_1[int(vertices[i]), :]
#         if i == vertices_count - 1:
#             q1 = vertex_table_1[int(vertices[0]), :]
#         else:
#             q1 = vertex_table_1[int(vertices[i + 1]), :]
#
#         for j in range(0, elem_count):
#             if check_elem_for_intersection[j] == 1:
#                 p0 = vertex_table_2[int(elem_table_2[j, 0]), :]
#                 p1 = vertex_table_2[int(elem_table_2[j, 1]), :]
#                 p2 = vertex_table_2[int(elem_table_2[j, 2]), :]
#                 intersect_point = intersect_line_triangle(q0, q1, p0, p1, p2)
#                 if intersect_point is None:
#                     pass
#                 else:
#                     intersecting_elements.append(j)
#
#     return intersecting_elements


def get_normals_table(elem_table, vertex_table):
    dim = np.shape(elem_table)
    elem_number = dim[0]
    norm_table = np.zeros([elem_number, 3])
    for i in range(0, elem_number):
        n = get_element_normal(i, elem_table, vertex_table)
        norm_table[i, 0] = n[0]
        norm_table[i, 1] = n[1]
        norm_table[i, 2] = n[2]

    return norm_table


def get_shadow_elements(norm_table):
    # Shadow elements are elements which do not possess a normal vector, as a result of being 1D in 3D space
    elem_number = np.shape(norm_table)[0]

    # Count number of valid elements
    valid_norm_table = np.zeros([elem_number, 1])
    valid_elem_count = 0
    for i in range(0, elem_number):
        norm_length = norm_table[i, 0]**2 + norm_table[i, 1]**2 + norm_table[i, 2]**2
        if norm_length > 0:
            valid_elem_count += 1
            valid_norm_table[i, 0] = 1
    return valid_norm_table


def move_vertices(vertices_to_move, new_locations, vertex_table):
    new_vertex_table = vertex_table.copy()
    for i in range(0, len(vertices_to_move)):
        new_vertex_table[vertices_to_move[i], :] = new_locations[i]
    return new_vertex_table


def remove_mesh_elements(elements_to_keep, norm_table, elem_table, lift_table):
    dim = np.shape(elem_table)
    elem_number = dim[0]
    valid_elem_count = int(np.sum(elements_to_keep))

    # Create new arrays to contain only valid elements
    new_norm_table = np.zeros([valid_elem_count, 3])
    new_elem_table = np.zeros([valid_elem_count, dim[1]])
    new_lift_table = np.zeros([valid_elem_count, 1])
    count = 0
    for i in range(0, elem_number):
        if elements_to_keep[i] == 1:
            new_norm_table[count, :] = norm_table[i, :]
            new_elem_table[count, :] = elem_table[i, :]
            new_lift_table[count, :] = lift_table[i, :]
            count += 1
    return new_norm_table, new_elem_table, new_lift_table


def get_elements_to_keep(elements_to_remove, elem_table):
    dim = np.shape(elem_table)
    elem_number = dim[0]
    elements_to_keep = np.ones([elem_number])
    for elem_id in elements_to_remove:
        elements_to_keep[elem_id] = 0
    return elements_to_keep


def keep_table_lines(lines_to_keep, table):
    dim = np.shape(table)
    lines_number = dim[0]
    valid_lines_count = int(np.sum(lines_to_keep))

    # Create new arrays to contain only valid elements
    new_table = np.zeros([valid_lines_count, dim[1]])
    count = 0
    for i in range(0, lines_number):
        if lines_to_keep[i] == 1:
            new_table[count, :] = table[i, :]
            count += 1
    return new_table


def merge_vertex_tables(vert_table_1, vert_table_2):
    vert_table = np.append(vert_table_1, vert_table_2, axis=0)
    delay = np.shape(vert_table_1)[0]
    return vert_table, delay


def merge_meshes(vert_table_1, elem_table_1, col_table_1, row_table_1, lift_table_1,
                 vert_table_2, elem_table_2, col_table_2, row_table_2, lift_table_2):

    # Merge vertex tables
    vert_table, delay = merge_vertex_tables(vert_table_1, vert_table_2)

    # Transmute vertices of elem_table_2
    transmute_vertex_to = list(range(delay, delay + np.shape(vert_table_2)[0]))
    dim = np.shape(elem_table_2)
    new_elem_table_2 = elem_table_2.copy()
    for i in range(0, dim[0]):
        for j in range(0, dim[1]):
            new_elem_table_2[i, j] = transmute_vertex_to[int(elem_table_2[i, j])]

    surface_1 = np.zeros(np.shape(elem_table_1)[0])
    surface_2 = np.ones(np.shape(elem_table_2)[0])
    surface = np.append(surface_1, surface_2, axis=0)

    # Merge element tables
    elem_table = np.append(elem_table_1, new_elem_table_2, axis=0)
    col_table = np.append(col_table_1, col_table_2, axis=0)
    row_table = np.append(row_table_1, row_table_2, axis=0)
    lift_table = np.append(lift_table_1, lift_table_2, axis=0)

    return vert_table, elem_table, col_table, row_table, surface, lift_table


def get_repeated_elements(vertex_table):
    dim = np.shape(vertex_table)
    vertex_number = dim[0]

    vertices_to_check = list(range(0, vertex_number))
    transmute_vertex_to = list(range(0, vertex_number))
    repeated_elements_count = 0
    # Loop through all vertices in table, except the last one
    for i in range(0, vertex_number - 1):
        if transmute_vertex_to[i] == i:
            vertex_i = vertex_table[i]
            for j in vertices_to_check:
                # Check if vertex i is equal to vertex j
                vertex_j = vertex_table[j]
                # Check if all 3 components are equal. Only implemented for 3D space
                if vertex_i[0] == vertex_j[0] and vertex_i[1] == vertex_j[1] and vertex_i[2] == vertex_j[2]:
                    transmute_vertex_to[j] = i
                    if i != j:
                        repeated_elements_count += 1
    return transmute_vertex_to


def transmute_vertices(transmute_vertex_to, elem_table, vertex_table):
    vertex_number = np.shape(vertex_table)[0]

    # re-index vertices to avoid gaps in ordering
    new_vertices_index = list(range(0, vertex_number))
    current_index = -1
    dummy_transmute_list = transmute_vertex_to.copy()
    keep_vertices_list = []
    for i in range(0, vertex_number):
        try:
            # Try statement, possibility to fail
            dummy_transmute_list.remove(i)
            # Operations to complete if try statement is successful
            keep_vertices_list.append(i)
            current_index += 1
        except ValueError:
            pass
        new_vertices_index[i] = current_index

    # Remove transmuted elements from vertex_table
    new_vertex_table = np.zeros([len(keep_vertices_list), np.shape(vertex_table)[1]])
    for i in range(0, len(keep_vertices_list)):
        new_vertex_table[i, :] = vertex_table[keep_vertices_list[i], :]

    # Update the transmutation vertex list with new vertices index
    for i in range(0, len(transmute_vertex_to)):
        transmute_vertex_to[i] = new_vertices_index[transmute_vertex_to[i]]

    # Transmute vertices in elem_table
    dim = np.shape(elem_table)
    new_elem_table = elem_table.copy()
    for i in range(0, dim[0]):
        for j in range(0, dim[1]):
            new_elem_table[i, j] = transmute_vertex_to[int(elem_table[i, j])]

    return new_elem_table, new_vertex_table


def convert_to_embedded_list(array):
    dim = np.shape(array)
    if len(dim) == 1:
        # Case 1D array
        lists = []
        for i in range(0, dim[0]):
            lists.append(array[i])

    elif len(dim) == 2:
        # Case 2D array
        lists = []
        for i in range(0, dim[0]):
            list_i = []
            for j in range(0, dim[1]):
                list_i.append(array[i, j])
            lists.append(list_i)
    else:
        lists = []
    return lists


def convert_to_numpy_array(list_of_lists):
    # Check dimensions are consistent
    dim_1 = len(list_of_lists)
    dim_2 = len(list_of_lists[0])
    valid_dimensions = True
    for i in range(1, len(list_of_lists)):
        if len(list_of_lists[i]) != dim_2:
            valid_dimensions = False

    # Create array if dimensions are valid
    if valid_dimensions:
        array = np.zeros([dim_1, dim_2])
        for i in range(0, dim_1):
            for j in range(0, dim_2):
                array[i, j] = list_of_lists[i][j]
    else:
        array = []

    # Return resulting array
    return array


def set_lift_property(x_list, lift):
    lift_body = x_list.copy()
    for i in range(0, len(x_list)):
        for j in range(0, len(x_list[i])):
            lift_body[i][j] = lift
    return lift_body


def reverse_array(array):
    res = np.zeros(np.shape(array))
    n = len(array)
    for i in range(0, n):
        res[i] = array[n - i - 1]
    return res


def merge_body_surfaces(x, y, su, sl):
    x_merged = []
    y_merged = []
    z_merged = []

    for i in range(0, len(x)):
        x_partial = x[i]
        y_partial = y[i]
        su_partial = su[i]
        sl_partial = sl[i]

        x_upper = x_partial
        y_upper = y_partial
        z_upper = su_partial
        if len(x_partial) != 0:
            x_lower = reverse_array(x_partial)
            y_lower = reverse_array(y_partial)
            z_lower = reverse_array(sl_partial)
            profile_x = np.append(x_lower, x_upper)
            profile_y = np.append(y_lower, y_upper)
            profile_z = np.append(z_lower, z_upper)
        else:
            profile_x = []
            profile_y = []
            profile_z = []

        if i == 0:
            x_merged = [profile_x]
            y_merged = [profile_y]
            z_merged = [profile_z]
        else:
            x_merged = np.append(x_merged, [profile_x], axis=0)
            y_merged = np.append(y_merged, [profile_y], axis=0)
            z_merged = np.append(z_merged, [profile_z], axis=0)

    # Add loops back to the sides, in case surfaces are not connected
    x_wrapped = [reverse_array(x_merged[0, :])]
    y_wrapped = [reverse_array(y_merged[0, :])]
    z_wrapped = [reverse_array(z_merged[0, :])]
    for i in range(0, len(x_merged)):
        x_wrapped = np.append(x_wrapped, [x_merged[i, :]], axis=0)
        y_wrapped = np.append(y_wrapped, [y_merged[i, :]], axis=0)
        z_wrapped = np.append(z_wrapped, [z_merged[i, :]], axis=0)
    x_wrapped = np.append(x_wrapped, [reverse_array(x_merged[-1, :])], axis=0)
    y_wrapped = np.append(y_wrapped, [reverse_array(y_merged[-1, :])], axis=0)
    z_wrapped = np.append(z_wrapped, [reverse_array(z_merged[-1, :])], axis=0)

    return x_wrapped, y_wrapped, z_wrapped


def get_surface_mesh_tables(x, y, su, sl, is_lift_body, flatten, append_fake_body, triangulate):
    if flatten:
        lift_bodies = set_lift_property(x, is_lift_body)
        z = (su + sl) / 2
        if append_fake_body:
            x, y, z, lift_bodies = append_plane_body_to_wing_surface(x, y, z, lift_bodies)
        # z = set_lift_property(x, 0)
    else:
        x, y, z = merge_body_surfaces(x, y, su, sl)
        lift_bodies = set_lift_property(x, is_lift_body)

    if np.shape(x) == np.shape(y):
        shape = np.shape(x)

        # Flatten x and y vertices into a vertex table
        vertex_table_length = shape[0] * shape[1]
        vertex_table = np.zeros([vertex_table_length, 3])
        for i in range(0, shape[0]):
            for j in range(0, shape[1]):
                vertex_index = i * shape[1] + j
                vertex_table[vertex_index, 0] = x[i, j]
                vertex_table[vertex_index, 1] = y[i, j]
                vertex_table[vertex_index, 2] = z[i, j]

        # Generate elements table
        if triangulate:
            # elem_table_length = 2 * (shape[0] - 1) * (shape[1] - 1)
            elem_table_length = 2 * (shape[0] - 1) * (shape[1] - 1)
            elem_table = np.zeros([elem_table_length, 3])
            lift_table = np.zeros([elem_table_length, 1])
            for i in range(0, shape[0] - 1):

                # Set indexes for i
                i_index = i
                i_plus_index = i + 1

                for j in range(0, shape[1] - 1):
                    elem_1_index = 2 * (i * (shape[1] - 1) + j)
                    elem_2_index = 1 + 2 * (i * (shape[1] - 1) + j)

                    # Set indexes for j
                    j_index = j
                    j_plus_index = j + 1

                    # Set elem_1 and elem_2 vertices in elem_table
                    elem_table[elem_1_index, 0] = int(i_index * shape[1] + j_index)
                    elem_table[elem_1_index, 1] = int(i_index * shape[1] + j_plus_index)
                    elem_table[elem_1_index, 2] = int(i_plus_index * shape[1] + j_plus_index)
                    elem_table[elem_2_index, 0] = int(i_index * shape[1] + j_index)
                    elem_table[elem_2_index, 1] = int(i_plus_index * shape[1] + j_plus_index)
                    elem_table[elem_2_index, 2] = int(i_plus_index * shape[1] + j_index)

                    lift_0 = lift_bodies[i_index][j_index]
                    lift_1 = lift_bodies[i_index][j_plus_index]
                    lift_2 = lift_bodies[i_plus_index][j_plus_index]
                    lift_3 = lift_bodies[i_plus_index][j_index]
                    # Check lift attributes for elem_1 vertices and set in lift_table
                    if lift_0 == lift_1 == lift_2 == 1:
                        lift_table[elem_1_index, 0] = 1
                    # Check lift attributes for elem_2 vertices and set in lift_table
                    if lift_0 == lift_2 == lift_3 == 1:
                        lift_table[elem_2_index, 0] = 1

        else:
            elem_table_length = (shape[0] - 1) * (shape[1] - 1)
            elem_table = np.zeros([elem_table_length, 4])
            lift_table = np.zeros([elem_table_length, 1])
            for i in range(0, shape[0] - 1):
                for j in range(0, shape[1] - 1):
                    elem_index = i * (shape[1] - 1) + j

                    # Set element vertices in elem_table
                    elem_table[elem_index, 0] = int(i * shape[1] + j)
                    elem_table[elem_index, 1] = int((i + 1) * shape[1] + j)
                    elem_table[elem_index, 2] = int((i + 1) * shape[1] + (j + 1))
                    elem_table[elem_index, 3] = int(i * shape[1] + (j + 1))

                    # Check lift attributes of vertices and set in lift_table
                    if lift_bodies[i][j] == lift_bodies[i][j+1] == lift_bodies[i+1][j+1] == lift_bodies[i+1][j] == 1:
                        lift_table[elem_index, 0] = 1

        # Cleanup the mesh
        transmute_vertex_to = get_repeated_elements(vertex_table)
        elem_table, vertex_table = transmute_vertices(transmute_vertex_to, elem_table, vertex_table)
        norm_table = get_normals_table(elem_table, vertex_table)
        elements_to_keep = get_shadow_elements(norm_table)
        norm_table, elem_table, lift_table = remove_mesh_elements(elements_to_keep, norm_table, elem_table, lift_table)

    else:
        vertex_table = []
        elem_table = []
        norm_table = []
        lift_table = []

    return vertex_table, elem_table, norm_table, lift_table


def get_vortex_mesh(x, y, su, sl, is_lift_body, append_fake_body):
    lift_bodies = set_lift_property(x, is_lift_body)
    z = (su + sl) / 2
    if append_fake_body:
        x, y, z, lift_bodies = append_plane_body_to_wing_surface(x, y, z, lift_bodies)

    if np.shape(x) == np.shape(y):
        shape = np.shape(x)

        # Flatten x and y vertices into a vertex table
        vertex_table_length = shape[0] * shape[1]
        vertex_table = np.zeros([vertex_table_length, 3])
        for i in range(0, shape[0]):
            for j in range(0, shape[1]):
                vertex_index = i * shape[1] + j
                vertex_table[vertex_index, 0] = x[i, j]
                vertex_table[vertex_index, 1] = y[i, j]
                vertex_table[vertex_index, 2] = z[i, j]

        # Generate elements table
        elem_table_length = (shape[0] - 1) * (shape[1] - 1)
        elem_table = np.zeros([elem_table_length, 4])
        column_table = np.zeros([elem_table_length, 1])
        row_table = np.zeros([elem_table_length, 1])
        lift_table = np.zeros([elem_table_length, 1])
        for i in range(0, shape[0] - 1):
            for j in range(0, shape[1] - 1):
                elem_index = i * (shape[1] - 1) + j

                elem_table[elem_index, 0] = int(i * shape[1] + j)
                elem_table[elem_index, 1] = int((i + 1) * shape[1] + j)
                elem_table[elem_index, 2] = int((i + 1) * shape[1] + (j + 1))
                elem_table[elem_index, 3] = int(i * shape[1] + (j + 1))

                # Check lift attributes of vertices and set in lift_table
                if lift_bodies[i][j] == lift_bodies[i][j+1] == lift_bodies[i+1][j+1] == lift_bodies[i+1][j] == 1:
                    lift_table[elem_index, 0] = 1

                # Generate column and row element information
                column_table[elem_index, 0] = i
                row_table[elem_index, 0] = j

        # Cleanup the mesh
        # Should not be many repeated elements...
        transmute_vertex_to = get_repeated_elements(vertex_table)
        elem_table, vertex_table = transmute_vertices(transmute_vertex_to, elem_table, vertex_table)
        norm_table = get_normals_table(elem_table, vertex_table)
        elements_to_keep = get_shadow_elements(norm_table)

        # Remove elements
        norm_table = keep_table_lines(elements_to_keep, norm_table)
        elem_table = keep_table_lines(elements_to_keep, elem_table)
        column_table = keep_table_lines(elements_to_keep, column_table)
        row_table = keep_table_lines(elements_to_keep, row_table)
        lift_table = keep_table_lines(elements_to_keep, lift_table)

    else:
        vertex_table = []
        elem_table = []
        column_table = []
        row_table = []
        lift_table = []

    return vertex_table, elem_table, column_table, row_table, lift_table


def append_plane_body_to_wing_surface(x_list, y_list, z_list, lift_bodies):
    y_min = 1000000
    y_min_location = 0
    y_min_is_negative = False
    # Find index i of y values closest to x axis
    for i in range(0, len(y_list)):
        for j in range(0, len(y_list[i])):
            if abs(y_list[i][j]) < y_min:
                y_min_location = i
                y_min = abs(y_list[i][j])
                if y_list[i][j] < 0:
                    y_min_is_negative = True
                else:
                    y_min_is_negative = False

    # Check cases of y position
    new_x_list = []
    new_y_list = []
    new_z_list = []
    new_lift_bodies = []
    # If y_min is negative, the extra panels need to be placed at the position after y_min_location
    y_min_check = y_min_location
    if y_min_is_negative:
        y_min_check += 1

    for i in range(0, len(x_list) + 1):
        x = []
        y = []
        z = []
        lift = []
        if i == y_min_check:
            for j in range(0, len(x_list[y_min_location])):
                x.append(x_list[y_min_location][j])
                y.append(0)
                z.append(z_list[y_min_location][j])
                lift.append(0)
        elif i < y_min_check:
            for j in range(0, len(x_list[i])):
                x.append(x_list[i][j])
                y.append(y_list[i][j])
                z.append(z_list[i][j])
                lift.append(lift_bodies[i][j])
        elif i > y_min_check:
            for j in range(0, len(x_list[i - 1])):
                x.append(x_list[i - 1][j])
                y.append(y_list[i - 1][j])
                z.append(z_list[i - 1][j])
                lift.append(lift_bodies[i - 1][j])
        if i == 0:
            new_x_list = [x]
            new_y_list = [y]
            new_z_list = [z]
            new_lift_bodies = [lift]
        else:
            new_x_list = np.append(new_x_list, [x], axis=0)
            new_y_list = np.append(new_y_list, [y], axis=0)
            new_z_list = np.append(new_z_list, [z], axis=0)
            new_lift_bodies = np.append(new_lift_bodies, [lift], axis=0)

    return new_x_list, new_y_list, new_z_list, new_lift_bodies


def visualize_mesh(vertex_table, elem_table):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    dim = np.shape(elem_table)
    for i in range(0, dim[0]):
        x = []
        y = []
        z = []
        for j in range(0, len(elem_table[i])):
            elem_id = int(elem_table[i, j])
            x.append(vertex_table[elem_id, 0])
            y.append(vertex_table[elem_id, 1])
            z.append(vertex_table[elem_id, 2])
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])

        ax.plot3D(x, y, z)


def visualize_mesh_elements(elements, vertex_table, elem_table, figure):
    if figure is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = figure

    if len(elements) > 0:
        for element_id in elements:
            x = []
            y = []
            z = []
            for j in range(0, len(elem_table[element_id])):
                elem_id = int(elem_table[element_id, j])
                x.append(vertex_table[elem_id, 0])
                y.append(vertex_table[elem_id, 1])
                z.append(vertex_table[elem_id, 2])
            x.append(x[0])
            y.append(y[0])
            z.append(z[0])
            ax.plot3D(x, y, z)

    elif len(elements) == 0:
        for element_id in range(0, np.shape(elem_table)[0]):
            x = []
            y = []
            z = []
            for j in range(0, len(elem_table[element_id])):
                elem_id = int(elem_table[element_id, j])
                x.append(vertex_table[elem_id, 0])
                y.append(vertex_table[elem_id, 1])
                z.append(vertex_table[elem_id, 2])
            x.append(x[0])
            y.append(y[0])
            z.append(z[0])
            ax.plot3D(x, y, z)

    return ax


def visualize_lines(ax, list_of_points):
    x = []
    y = []
    z = []
    for point in list_of_points:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    ax.plot3D(x, y, z)

    return ax


