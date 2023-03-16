import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import cst_geometry as geo
import export_to_gmsh as exp


class PointObject:
    def __init__(self, object_id, list_of_points):
        self.id = object_id
        self.points = list_of_points


class Box:
    def __init__(self, box_id, list_of_points):
        self.id = box_id
        self.points = list_of_points
        self.min = [1e15, 1e15, 1e15]
        self.max = [-1e15, -1e15, -1e15]
        self.set_boundaries()

    def set_boundaries(self):
        for point in self.points:
            # Set min
            self.min[0] = min(self.min[0], point[0])
            self.min[1] = min(self.min[1], point[1])
            self.min[2] = min(self.min[2], point[2])
            # Set max
            self.max[0] = max(self.max[0], point[0])
            self.max[1] = max(self.max[1], point[1])
            self.max[2] = max(self.max[2], point[2])

    def box_touches_box(self, box):
        x_touches = False
        y_touches = False
        z_touches = False
        if not (self.min[0] > box.max[0] or self.max[0] < box.min[0]):
            x_touches = True
        if not (self.min[1] > box.max[1] or self.max[1] < box.min[1]):
            y_touches = True
        if not (self.min[2] > box.max[2] or self.max[2] < box.min[2]):
            z_touches = True

        if x_touches and y_touches and z_touches:
            touches = True
        else:
            touches = False
        return touches

    def is_box_inside_box(self, box):
        is_inside = True
        if self.min[0] == self.max[0] or self.min[1] == self.max[1] or self.min[2] == self.max[2]:
            is_inside = False
        else:
            for point in box.points:
                u = (point[0] - self.min[0]) / (self.max[0] - self.min[0])
                v = (point[1] - self.min[1]) / (self.max[1] - self.min[1])
                w = (point[2] - self.min[2]) / (self.max[2] - self.min[2])
                if (1 >= u >= 0) and (1 >= v >= 0) and (1 >= w >= 0):
                    pass
                else:
                    is_inside = False
            # u1 = (box.min[0] - self.min[0]) / (self.max[0] - self.min[0])
            # v1 = (box.min[1] - self.min[1]) / (self.max[1] - self.min[1])
            # w1 = (box.min[2] - self.min[2]) / (self.max[2] - self.min[2])
            # if (1 >= u1 >= 0) and (1 >= v1 >= 0) and (1 >= w1 >= 0):
            #     pass
            # else:
            #     is_inside = False
            # u2 = (box.max[0] - self.min[0]) / (self.max[0] - self.min[0])
            # v2 = (box.max[1] - self.min[1]) / (self.max[1] - self.min[1])
            # w2 = (box.max[2] - self.min[2]) / (self.max[2] - self.min[2])
            # if (1 >= u2 >= 0) and (1 >= v2 >= 0) and (1 >= w2 >= 0):
            #     pass
            # else:
            #     is_inside = False
        return is_inside

    def is_object_inside_box(self, point_object):
        is_inside = True
        if self.min[0] == self.max[0] or self.min[1] == self.max[1] or self.min[2] == self.max[2]:
            is_inside = False
        else:
            for point in point_object.points:
                u = (point[0] - self.min[0]) / (self.max[0] - self.min[0])
                v = (point[1] - self.min[1]) / (self.max[1] - self.min[1])
                w = (point[2] - self.min[2]) / (self.max[2] - self.min[2])
                if (1 >= u >= 0) and (1 >= v >= 0) and (1 >= w >= 0):
                    pass
                else:
                    is_inside = False
                    break
        return is_inside


class Container:
    def __init__(self, list_of_boxes):
        self.contains = []
        self.box = Box(None, [])
        self.subdivisions = []
        self.nb_subdivisions = 0
        self.set_container_boundaries(list_of_boxes)

    def set_container_boundaries(self, list_of_boxes):
        container_min = [1e15, 1e15, 1e15]
        container_max = [-1e15, -1e15, -1e15]
        for box in list_of_boxes:
            for point in box.points:
                # Set min
                container_min[0] = min(container_min[0], point[0])
                container_min[1] = min(container_min[1], point[1])
                container_min[2] = min(container_min[2], point[2])
                # Set max
                container_max[0] = max(container_max[0], point[0])
                container_max[1] = max(container_max[1], point[1])
                container_max[2] = max(container_max[2], point[2])
        self.box.points = [container_min, container_max]
        self.box.set_boundaries()

    def subdivide(self):
        # Clear subdivisions
        self.subdivisions.clear()
        if self.nb_subdivisions == 0:
            self.nb_subdivisions = 1

        # Fill subdivisions back with the 8 boxes
        d = [self.box.max[0] - self.box.min[0],
             self.box.max[1] - self.box.min[1],
             self.box.max[2] - self.box.min[2]]

        ni = 1
        nj = 1
        nk = 1
        if d[0] > d[1] and d[0] > d[2]:
            ni = 2
        elif d[1] > d[0] and d[1] > d[2]:
            nj = 2
        else:
            nk = 2

        # Divide d by values of ns
        d_ijk = [d[0] / ni, d[1] / nj, d[2] / nk]

        for i in range(0, ni):
            for j in range(0, nj):
                for k in range(0, nk):
                    min_ijk = [self.box.min[0] + i * d_ijk[0],
                               self.box.min[1] + j * d_ijk[1],
                               self.box.min[2] + k * d_ijk[2]]
                    max_ijk = [min_ijk[0] + d_ijk[0],
                               min_ijk[1] + d_ijk[1],
                               min_ijk[2] + d_ijk[2]]
                    # Create new container object and add it to subdivisions list
                    container_dimensions = [Box(None, [min_ijk, max_ijk])]
                    self.subdivisions.append(Container(container_dimensions))

    def place_box_in_container(self, box_object):
        if self.box.is_box_inside_box(box_object):
            self.contains.append(box_object)
            return True
        else:
            return False

    def place_boxes_in_subdivisions(self):
        if len(self.subdivisions) > 0:
            objects_to_remove = []
            for i in range(0, len(self.contains)):
                box_object = self.contains[i]
                for container in self.subdivisions:
                    if container.place_box_in_container(box_object):
                        objects_to_remove.insert(0, i)
                        break
            # Clear the necessary items from the contains list
            if len(objects_to_remove) > 0:
                for object_id in objects_to_remove:
                    self.contains.pop(object_id)

    def is_container_empty(self):
        # Return False if container contains at least one item
        if len(self.contains) > 0:
            return False
        else:
            for subdivision in self.subdivisions:
                # As soon as a subdivision is found containing something, return False
                if not subdivision.is_container_empty():
                    return False
            # If all subdivisions are empty, return True
            return True


def append_box_to_container(container, box_object):
    if container.box.is_box_inside_box(box_object):
        box_fits_in_subdivision = False
        if len(container.subdivisions) > 0:
            for subdivision in container.subdivisions:
                box_fits_in_subdivision = append_box_to_container(subdivision, box_object)
                if box_fits_in_subdivision:
                    break
        # If the object does not fit in any of the subdivisions, then it must be added to the current container
        if not box_fits_in_subdivision:
            container.contains.append(box_object)
        # Return the value True to indicate that the object fits inside this container
        return True
    else:
        # Return the value False to indicate that the object does not fit inside the container
        return False


def generate_containers(container, list_of_boxes, number_of_subdivisions):
    list_of_boxes_copy = list_of_boxes.copy()
    # If container is not previously defined, initialize one with the objects provided
    if container is None:
        container = Container(list_of_boxes)
        # Do not forget to clear the objects from the container, as they will be added back later (but filtered)
        container.nb_subdivisions = number_of_subdivisions

    # Determine if container needs to be subdivided
    subdivide = False
    if number_of_subdivisions > 0:
        subdivide = True

    # Subdivide container if necessary
    if subdivide and len(container.subdivisions) == 0:
        container.subdivide()

    if len(list_of_boxes) > 24 and len(container.subdivisions) > 0:
        for subdivision in container.subdivisions:
            # Get the objects fitting inside the subdivision
            boxes_in_subdivision = []
            boxes_to_remove_from_list = []
            for i in range(0, len(list_of_boxes_copy)):
                box_object = list_of_boxes_copy[i]
                if subdivision.box.is_box_inside_box(box_object):
                    boxes_in_subdivision.append(box_object)
                    boxes_to_remove_from_list.insert(0, i)

            # Clear the necessary items from the contains list
            if len(boxes_to_remove_from_list) > 0:
                for box_index in boxes_to_remove_from_list:
                    list_of_boxes_copy.pop(box_index)

            # Add objects to subdivision
            generate_containers(subdivision, boxes_in_subdivision, number_of_subdivisions - 1)

    # Place leftover objects in container
    # (if the box was initialized from the object list, all objects will end up in a box)
    for box_object in list_of_boxes_copy:
        container.place_box_in_container(box_object)

    # print(number_of_subdivisions, len(list_of_boxes_copy))
    return container


def set_mesh_elements_in_containers(mesh, max_subdivisions):
    # Extract list of Box objects from mesh -> elements
    elements = mesh.elements
    list_of_elements = []
    for i in range(0, len(elements)):
        element = elements[i]
        list_of_points = []
        for vertex in element.vertices:
            list_of_points.append(mesh.vertices[int(vertex)])
        list_of_elements.append(Box(i, list_of_points))

    # Create container structure
    container = generate_containers(None, list_of_elements, max_subdivisions)
    return container


def get_collisions(container_1, container_2, collisions):
    if container_1.box.box_touches_box(container_2.box):
        # Then we need to check for all potential collisions in contains 1 and 2
        for box_object_1 in container_1.contains:
            for box_object_2 in container_2.contains:
                if box_object_1.box_touches_box(box_object_2):
                    # Then there is a potential collision
                    collisions.append([box_object_1.id, box_object_2.id])

        if len(container_1.contains) > 0:
            for subdivision_2 in container_2.subdivisions:
                if not subdivision_2.is_container_empty():
                    collisions = get_collisions(container_1, subdivision_2, collisions)

        for subdivision_1 in container_1.subdivisions:
            if not subdivision_1.is_container_empty():
                collisions = get_collisions(subdivision_1, container_2, collisions)

    return collisions


def get_potential_collisions(container, box_to_test, list_of_collisions):
    add_box_to_container = False
    if container.box.is_box_inside_box(box_to_test):
        add_box_to_container = True

    if container.box.box_touches_box(box_to_test):
        for box_object in container.contains:
            if box_to_test.box_touches_box(box_object):
                list_of_collisions.append([box_to_test.id, box_object.id])

        for subdivision in container.subdivisions:
            if subdivision.box.box_touches_box(box_to_test):
                list_of_collisions, add_box_to_subdivision = get_potential_collisions(subdivision, box_to_test,
                                                                                      list_of_collisions)
                if add_box_to_subdivision:
                    add_box_to_container = False

    if add_box_to_container:
        container.contains.append(box_to_test)

    return list_of_collisions, add_box_to_container


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

    class Wingstation:
        def __init__(self, wingstation_id, surface, elements):
            self.id = wingstation_id
            self.elements = elements
            self.surface = surface

    class Surface:
        def __init__(self, surface_id, surface_type):
            # Two surfaces can have the same ID if they do not have the same type
            self.id = surface_id
            self.type = surface_type
            self.wingstations = []

    def __init__(self):
        self.vertices = []
        self.elements = []
        self.wingstations = []
        self.surfaces = []
        self.vortex_count = 0
        self.doublet_count = 0
        self.wingstation_count = 0
        self.surface_count = 0

    def append_vertex(self, vertex):
        if len(vertex) == 3:
            new_vertex_id = len(self.vertices)
            self.vertices.append(vertex)
            return new_vertex_id
        else:
            print("ERROR - Incorrect Vertex Dimensions")
            return None

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

    def create_new_surface(self, surface_type):
        new_surface_id = self.surface_count
        new_surface = self.Surface(new_surface_id, surface_type)
        self.surfaces.append(new_surface)
        self.surface_count += 1

    def create_new_wingstation(self, surface_id, vortexes):
        new_wingstation_id = self.wingstation_count
        new_wingstation = self.Wingstation(new_wingstation_id, surface_id, vortexes)
        self.wingstations.append(new_wingstation)
        self.wingstation_count += 1
        self.surfaces[surface_id].wingstations.append(new_wingstation_id)

    def append_vortex_to_wingstation(self, vortex, wingstation_id):
        self.wingstations[wingstation_id].elements.append(vortex.id)

    def create_vortex(self, vertices, surface, station, row):
        new_vortex = self.Vortex(vertices, surface, station, row)
        self.append_element(new_vortex)
        if self.wingstation_count > station:
            if self.wingstations[station].surface == surface:
                self.wingstations[station].elements.append(new_vortex.id)
            else:
                print("Error - Wingstation is not on specified surface")
        else:
            print("Error - Wingstation does not yet exist")
            print(station)

    def create_doublet(self, vertices, surface):
        new_doublet = self.Doublet(vertices, surface)
        self.append_element(new_doublet)

    def get_vertices(self, vertex_ids):
        vertex_list = []
        for vertex_id in vertex_ids:
            vertex_list.append(self.vertices[int(vertex_id)])
        return vertex_list

    def get_normals(self):
        normals = []
        for element in self.elements:
            a = np.array(self.vertices[int(element.vertices[0])])
            b = np.array(self.vertices[int(element.vertices[1])])
            c = np.array(self.vertices[int(element.vertices[2])])
            norm = np.cross(b - a, c - a)
            norm_length = math.sqrt(np.dot(norm, norm))
            if norm_length > 0:
                norm = norm / norm_length
            else:
                norm = [0, 0, 0]
            normals.append(norm)
        return normals

    def delete_element(self, element_id):
        # Get element
        type_of_element = self.elements[element_id].type
        # Remove element
        self.elements.pop(element_id)
        if type_of_element == "Vortex":
            # Update vortex count and IDs
            self.vortex_count += -1
            for element in self.elements:
                if element.type == "Vortex" and element.id > element_id:
                    element.id += -1
        elif type_of_element == "Doublet":
            # Update doublet count and IDs
            self.doublet_count += -1
            for element in self.elements:
                if element.type == "Doublet" and element.id > element_id:
                    element.id += -1
        else:
            print("ERROR - Incorrect Element Format")

    def get_unused_vertices(self):
        vertices_usage = np.zeros([len(self.vertices), 1])
        for element in self.elements:
            for vertex in element.vertices:
                vertices_usage[int(vertex)] += 1
        unused_vertices = []
        for i in range(0, len(vertices_usage)):
            if vertices_usage[i] == 0:
                unused_vertices.append(i)
        return unused_vertices

    def get_transmute_vertices(self):
        unused_vertices = self.get_unused_vertices()
        vertex_number = len(self.vertices)
        transmute_vertex_to = list(range(0, vertex_number))
        for i in range(0, vertex_number):
            for unused_vertex in unused_vertices:
                if unused_vertex <= i:
                    transmute_vertex_to[i] += -1
        return transmute_vertex_to

    def delete_unused_vertices(self):
        unused_vertices = self.get_unused_vertices()
        transmute_vertex_to = self.get_transmute_vertices()
        vertex_count = len(self.vertices)

        # re-index vertices to avoid gaps in ordering
        new_vertices_index = list(range(0, vertex_count))
        current_index = -1
        dummy_transmute_list = transmute_vertex_to.copy()
        keep_vertices_list = []
        for i in range(0, vertex_count):
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
        for i in range(0, len(unused_vertices)):
            self.vertices.pop(unused_vertices[len(unused_vertices) - 1 - i])

        # Transmute vertices in elem_table
        dim = len(self.elements)
        for i in range(0, dim):
            for j in range(0, len(self.elements[i].vertices)):
                self.elements[i].vertices[j] = transmute_vertex_to[int(self.elements[i].vertices[j])]

    def convert_used_vertices_to_box_list(self, tol):
        list_of_boxes = []
        for i in range(0, len(self.elements)):
            for j in range(0, len(self.elements[i].vertices)):
                vertex_id = int(self.elements[i].vertices[j])
                point = self.vertices[vertex_id]
                if tol > 0:
                    point_min = [point[0] - tol, point[1] - tol, point[2] - tol]
                    point_max = [point[0] + tol, point[1] + tol, point[2] + tol]
                    list_of_points = [point_min, point_max]
                else:
                    list_of_points = [point]
                list_of_boxes.append(Box([vertex_id, [i, j]], list_of_points))
        return list_of_boxes

    def set_used_vertices_in_containers(self, max_subdivisions, tol):
        # Extract list of PointObjects from mesh -> elements
        list_of_boxes = self.convert_used_vertices_to_box_list(tol)
        # Create container structure
        container = generate_containers(None, list_of_boxes, max_subdivisions)
        return container

    def stitch_used_vertices(self, tol, nb_subdivisions):
        # Get container of used vertices
        print(" Step 1")
        container_1 = self.set_used_vertices_in_containers(nb_subdivisions, 2 * tol)
        container_2 = self.set_used_vertices_in_containers(nb_subdivisions, 2 * tol)

        print(" Step 2")
        list_of_collisions = get_collisions(container_1, container_2, [])

        tol_square = tol ** 2
        print(len(list_of_collisions))
        for collision in list_of_collisions:
            ids_1 = collision[1]
            ids_2 = collision[0]
            address_1 = self.vertices[ids_1[0]]
            address_2 = self.vertices[ids_2[0]]
            vertex_1 = self.vertices[int(address_1)]
            vertex_2 = self.vertices[int(address_2)]

            distance_square = (vertex_1[0] - vertex_2[0]) ** 2 + \
                              (vertex_1[1] - vertex_2[1]) ** 2 + (vertex_1[2] - vertex_2[2]) ** 2
            # If vertices are the same, set vertex_address_2 equal to vertex_address_1
            if distance_square <= tol_square:
                # ids_i contains [vertex_id, [element_id, element_vertex]]
                self.elements[ids_1[1][0]].vertices[ids_1[1][1]] = ids_2

    # def update_surfaces(self):
    #     # Re-initialize element count in mesh surfaces
    #     for surface in self.surfaces:
    #         surface.element_count = 0
    #
    #     for element in self.elements:
    #         is_new_surface = True
    #         for surface in self.surfaces:
    #             if element.surface == surface.id and element.type == surface.type:
    #                 surface.element_count += 1
    #                 # Then surface is already initialized, break out of for loop
    #                 is_new_surface = False
    #                 break
    #
    #         if is_new_surface:
    #             # Create new surface in mesh and initialize element count to 1
    #             self.surfaces.append(Mesh.Surface(element.surface, element.type, 1))
    #             self.surface_count += 1


def convert_doublet_to_mesh_object(mesh, vertex_table, element_table, surface):
    if mesh is None:
        mesh = Mesh()
    for i in range(0, np.shape(vertex_table)[0]):
        # mesh.append_vertex([vertex_table[i, 0], vertex_table[i, 1], vertex_table[i, 2]])
        mesh.append_vertex(vertex_table[i, :])
    for i in range(0, np.shape(element_table)[0]):
        mesh.create_doublet(element_table[i, :], surface)
    return mesh


def convert_vortex_to_mesh_object(mesh, vertex_table, element_table, surface, stations, rows):
    if mesh is None:
        mesh = Mesh()
    for i in range(0, np.shape(vertex_table)[0]):
        # mesh.append_vertex([vertex_table[i, 0], vertex_table[i, 1], vertex_table[i, 2]])
        mesh.append_vertex(vertex_table[i, :])
    for i in range(0, np.shape(element_table)[0]):
        mesh.create_vortex(element_table[i, :], surface, stations[i], rows[i])
    return mesh


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
            _, _, su, sl = geo.build_body_surfaces_complete([res[0]], [res[1]], param)

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
                _, _, su, sl = geo.build_body_surfaces_complete([new_res[0]], [new_res[1]], param)

            elif res_prev is None:
                # domain, _ = geo.get_full_domain_boundary(25, 25, param)

                new_res = geo.get_domain_boundary(vec_x[i], vec_y[i],
                                                  vec_x[i - 1] - vec_x[i], vec_y[i - 1] - vec_y[i], param, 50)
                x_prev, y_prev, su_prev, sl_prev = geo.build_body_surfaces_complete([new_res[0]], [new_res[1]], param)

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
    distance_tol = 1e-9
    # distance_tol = 0

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

                vertex_distance_square = (vertex_i[0] - vertex_j[0])**2 + (vertex_i[1] - vertex_j[1])**2 + (vertex_i[2] - vertex_j[2])**2

                # Check if all 3 components are equal. Only implemented for 3D space
                if vertex_distance_square <= distance_tol:
                # if vertex_i[0] == vertex_j[0] and vertex_i[1] == vertex_j[1] and vertex_i[2] == vertex_j[2]:
                    transmute_vertex_to[j] = i
                    if i != j:
                        repeated_elements_count += 1
    return transmute_vertex_to


def transmute_vertices_mesh(transmute_vertex_to, mesh):
    vertex_table = mesh.vertices
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
    dim = mesh.doublet_count + mesh.vortex_count
    for i in range(0, dim):
        for j in range(0, 3):
            mesh.elements[i].vertices[j] = transmute_vertex_to[int(mesh.elements[i].vertices[j])]

    return mesh


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


def move_first_to_last(array):
    res = np.zeros(np.shape(array))
    if len(np.shape(array)) == 1:
        for i in range(0, len(array)):
            if i == 0:
                res[-1] = array[0]
            else:
                res[i-1] = array[i]
    elif len(np.shape(array)) == 2:
        for i in range(0, np.shape(array)[0]):
            if i == 0:
                res[-1, :] = array[0, :]
            else:
                res[i-1, :] = array[i, :]
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

    # # Add loop at the end of the profiles, to connect
    # x_loop = np.reshape(x_merged[:, 0], (np.shape(x_merged[:, 0])[0], 1))
    # y_loop = np.reshape(y_merged[:, 0], (np.shape(y_merged[:, 0])[0], 1))
    # z_loop = np.reshape(z_merged[:, 0], (np.shape(z_merged[:, 0])[0], 1))
    # x_merged = np.append(x_merged, x_loop, axis=1)
    # y_merged = np.append(y_merged, y_loop, axis=1)
    # z_merged = np.append(z_merged, z_loop, axis=1)

    # # Add loops back to the sides, in case surfaces are not connected
    # x_wrapped = [move_first_to_last(reverse_array(x_merged[0, :]))]
    # y_wrapped = [move_first_to_last(reverse_array(y_merged[0, :]))]
    # z_wrapped = [move_first_to_last(reverse_array(z_merged[0, :]))]
    # for i in range(0, len(x_merged)):
    #     x_wrapped = np.append(x_wrapped, [x_merged[i, :]], axis=0)
    #     y_wrapped = np.append(y_wrapped, [y_merged[i, :]], axis=0)
    #     z_wrapped = np.append(z_wrapped, [z_merged[i, :]], axis=0)
    # x_wrapped = np.append(x_wrapped, [reverse_array(x_merged[-1, :])], axis=0)
    # y_wrapped = np.append(y_wrapped, [reverse_array(y_merged[-1, :])], axis=0)
    # z_wrapped = np.append(z_wrapped, [reverse_array(z_merged[-1, :])], axis=0)

    # return x_wrapped, y_wrapped, z_wrapped
    return x_merged, y_merged, z_merged


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


def get_triangle_intersection(p, q):
    def signed_tetra_volume(a, b, c, d):
        return np.sign(np.dot(np.cross(b - a, c - a), d - a) / 6.0)

    def get_intersection(a, b, x, y, z):
        n = np.cross(y - x, z - x)
        t = np.dot(x - a, n) / np.dot(b - a, n)
        return a + t * (b - a)

    def check_intersection(a, b, x, y, z):
        s3 = signed_tetra_volume(a, b, x, y)
        s4 = signed_tetra_volume(a, b, y, z)
        s5 = signed_tetra_volume(a, b, z, x)
        if s3 == s4 and s4 == s5:
            return True
        else:
            return False

    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    q0 = q[0]
    q1 = q[1]
    q2 = q[2]

    max_intersections = 7
    intersections = []

    # Check vertices of triangle P for intersections
    v0 = signed_tetra_volume(p0, q0, q1, q2)
    v1 = signed_tetra_volume(p1, q0, q1, q2)
    v2 = signed_tetra_volume(p2, q0, q1, q2)
    # Check edge p0-p1
    if len(intersections) < max_intersections and v0 != v1:
        if check_intersection(p0, p1, q0, q1, q2):
            intersections.append([get_intersection(p0, p1, q0, q1, q2), [0, [0, 1]]])
    # Check edge p0-p2
    if len(intersections) < max_intersections and v0 != v2:
        if check_intersection(p0, p2, q0, q1, q2):
            intersections.append([get_intersection(p0, p2, q0, q1, q2), [0, [0, 2]]])
    # Check edge p1-p2
    if len(intersections) < max_intersections and v1 != v2:
        if check_intersection(p1, p2, q0, q1, q2):
            intersections.append([get_intersection(p1, p2, q0, q1, q2), [0, [1, 2]]])

    # Check vertices of triangle Q for intersections
    u0 = signed_tetra_volume(q0, p0, p1, p2)
    u1 = signed_tetra_volume(q1, p0, p1, p2)
    u2 = signed_tetra_volume(q2, p0, p1, p2)
    # Check edge q0-q1
    if len(intersections) < max_intersections and u0 != u1:
        if check_intersection(q0, q1, p0, p1, p2):
            intersections.append([get_intersection(q0, q1, p0, p1, p2), [1, [0, 1]]])
    # Check edge q0-q2
    if len(intersections) < max_intersections and u0 != u2:
        if check_intersection(q0, q2, p0, p1, p2):
            intersections.append([get_intersection(q0, q2, p0, p1, p2), [1, [0, 2]]])
    # Check edge q1-q2
    if len(intersections) < max_intersections and u1 != u2:
        if check_intersection(q1, q2, p0, p1, p2):
            intersections.append([get_intersection(q1, q2, p0, p1, p2), [1, [1, 2]]])

    return intersections


def do_boxes_intersect(box_1, box_2):
    box_1_min = box_1[0]
    box_1_max = box_1[1]
    box_2_min = box_2[0]
    box_2_max = box_2[1]
    a = box_1_max[0] < box_2_min[0] or box_1_max[1] < box_2_min[1] or box_1_max[2] < box_2_min[2]
    b = box_2_max[0] < box_1_min[0] or box_2_max[1] < box_1_min[1] or box_2_max[2] < box_1_min[2]
    if not (a or b):
        return True
    else:
        return False


def get_distance_squared(vertex_1, vertex_2):
    return (vertex_1[0] - vertex_2[0]) ** 2 + (vertex_1[1] - vertex_2[1]) ** 2 + (vertex_1[2] - vertex_2[2]) ** 2


def combine_meshes(mesh1, mesh2, mode):
    def get_element_boundaries(vertices):
        vert_min = vertices[0]
        vert_max = vert_min.copy()
        for vertex in vertices:
            vert_min = [min(vert_min[0], vertex[0]), min(vert_min[1], vertex[1]), min(vert_min[2], vertex[2])]
            vert_max = [max(vert_max[0], vertex[0]), max(vert_max[1], vertex[1]), max(vert_max[2], vertex[2])]
        return [vert_min, vert_max]

    def decode_intersections(intersections_info, vertices_info):
        res = intersections_info.copy()
        for i in range(0, len(intersections_info)):
            res[i][1][1][0] = vertices_info[int(res[i][1][0])][int(res[i][1][1][0])]
            res[i][1][1][1] = vertices_info[int(res[i][1][0])][int(res[i][1][1][1])]
        return res

    def initiate_permutations(vertices_list):
        res = []
        for vertex in vertices_list:
            res.append([vertex.copy()])
        return res

    def append_permutation(permutations, new_permutation):
        location = new_permutation[0]
        for vertex_id in new_permutation[1][1]:
            vertex_id = int(vertex_id)
            permutations[vertex_id].append(location)
        return permutations

    def get_number_of_vertices_with_permutations(permutations):
        permute_number = 0
        for i in range(0, len(permutations)):
            if len(permutations[i]) > 1:
                permute_number += 1
        return permute_number

    def get_closest_location_to_vertex(permutations):
        permute_vertex_to = []
        permuted_vertices = []
        permute_count = 0
        for i in range(0, len(permutations)):
            if len(permutations[i]) == 1:
                permute_vertex_to.append(permutations[i][0])
            else:
                permute_count += 1
                permuted_vertices.append(i)
                min_dist_square = 1e15
                best_permute_yet = 1
                for j in range(1, len(permutations[i])):
                    dist_square = get_distance_squared(permutations[i][0], permutations[i][j])
                    if dist_square <= min_dist_square:
                        min_dist_square = dist_square
                        best_permute_yet = j
                permute_vertex_to.append(permutations[i][best_permute_yet])
        return permute_vertex_to, permuted_vertices, permute_count

    def apply_permutations(mesh, permutations):
        if len(mesh.vertices) == len(permutations):
            mesh.vertices = permutations
        else:
            print("Vertices list length does not match permutations list length")
        return mesh

    # def add_vertices_to_mesh(mesh, vertices):
    #     vertices_index = []
    #     for vertex in vertices:
    #         vertex_id = mesh.append_vertex(vertex)
    #         vertices_index.append(vertex_id)
    #     return mesh, vertices_index
    #
    # def split_element_2in(mesh, element_id, point_1_id, point_2_id):
    #     # Remove 1 element, add 5 elements
    #     replaced_elem = mesh.elements[element_id]
    #     vert_1 = mesh.vertices[int(replaced_elem.vertices[0])]
    #     vert_2 = mesh.vertices[int(replaced_elem.vertices[1])]
    #     vert_3 = mesh.vertices[int(replaced_elem.vertices[2])]
    #     surface = replaced_elem.surface
    #     # Get closest element vertex to point_1
    #     point_1 = mesh.vertices[int(point_1_id)]
    #     point_2 = mesh.vertices[int(point_2_id)]
    #
    #     dist_1_squared = get_distance_squared(point_1, vert_1)
    #     dist_2_squared = get_distance_squared(point_1, vert_2)
    #     dist_3_squared = get_distance_squared(point_1, vert_3)
    #
    #     new_elem_1 = [replaced_elem.vertices[0], replaced_elem.vertices[1], point_1_id]
    #     new_elem_2 = [replaced_elem.vertices[0], point_2_id, point_1_id]
    #     new_elem_3 = [replaced_elem.vertices[1], point_1_id, replaced_elem.vertices[2]]
    #     new_elem_4 = [replaced_elem.vertices[0], replaced_elem.vertices[2], point_2_id]
    #     new_elem_5 = [point_2_id, replaced_elem.vertices[2], point_1_id]
    #
    #     mesh.create_doublet(new_elem_1, surface)
    #     mesh.create_doublet(new_elem_2, surface)
    #     mesh.create_doublet(new_elem_3, surface)
    #     mesh.create_doublet(new_elem_4, surface)
    #     mesh.create_doublet(new_elem_5, surface)
    #     mesh.delete_element(element_id)
    #     return mesh
    #
    # def split_element_2edge(mesh, element_id, point1, point2):
    #     pass
    #
    # def split_element_1in1edge(mesh, element_id, point_in, point_edge):
    #     pass

    intersection_list_1 = []
    intersection_list_2 = []
    permutations_1 = initiate_permutations(mesh1.vertices)
    permutations_2 = initiate_permutations(mesh2.vertices)
    for element1 in mesh1.elements:
        intersects = False
        if len(element1.vertices) == 3:
            p_vertices = mesh1.get_vertices(element1.vertices)
            p_bounds = get_element_boundaries(p_vertices)
            for element2 in mesh2.elements:
                if len(element2.vertices) == 3:
                    q_vertices = mesh2.get_vertices(element2.vertices)
                    q_bounds = get_element_boundaries(q_vertices)

                    # Only run the intersection check if boxes intersect
                    if do_boxes_intersect(p_bounds, q_bounds):
                        intersections = get_triangle_intersection(p_vertices, q_vertices)
                        decoded_intersections = decode_intersections(intersections, [element1.vertices, element2.vertices])
                        for decoded_intersection in decoded_intersections:
                            if decoded_intersection[1][0] == 0:
                                permutations_1 = append_permutation(permutations_1, decoded_intersection)
                            elif decoded_intersection[1][0] == 1:
                                permutations_2 = append_permutation(permutations_2, decoded_intersection)

                        if len(intersections) == 2:
                            intersects = True
                            intersection_list_2.append(element2.id)
                        elif len(intersections) == 1 or len(intersections) > 2:
                            # This should not happen
                            print(len(intersections), "h")

            if intersects:
                intersection_list_1.append(element1.id)

    transmute_1, permuted_1, n1_permutes = get_closest_location_to_vertex(permutations_1)
    transmute_2, permuted_2, n2_permutes = get_closest_location_to_vertex(permutations_2)
    mesh1 = apply_permutations(mesh1, transmute_1)
    mesh2 = apply_permutations(mesh2, transmute_2)

    # n1 = get_number_of_vertices_with_permutations(permutations_1)
    # n2 = get_number_of_vertices_with_permutations(permutations_2)
    # if n1 <= n2:
    #     transmute_1 = get_closest_location_to_vertex(permutations_1)
    #     apply_permutations(mesh1, transmute_1)
    # else:
    #     transmute_2 = get_closest_location_to_vertex(permutations_2)
    #     apply_permutations(mesh2, transmute_2)

    return intersection_list_1, intersection_list_2, mesh1, mesh2, permuted_1, permuted_2


def get_surface_mesh_tables_v2(body, mesh, mode, is_vortex):
    def get_triangle_area_square(vertices):
        point_a = np.array(mesh.vertices[int(vertices[0])])
        point_b = np.array(mesh.vertices[int(vertices[1])])
        point_c = np.array(mesh.vertices[int(vertices[2])])
        n = cross_product(point_b - point_a, point_c - point_a)
        area_square = (n[0] ** 2 + n[1] ** 2 + n[2] ** 2) / 4
        return area_square

    def mesh_rectangle_to_triangles(vertices):
        if np.shape(vertices) == (2, 2):
            # Two ways to cut the shape into triangles
            # Arrangement 1
            triangle_1a = [vertices[0, 0], vertices[0, 1], vertices[1, 1]]
            triangle_1b = [vertices[0, 0], vertices[1, 1], vertices[1, 0]]
            # Measure the "greatness" of these two triangles
            area_min_1 = get_triangle_area_square(triangle_1a)
            area_min_1 = min(area_min_1, get_triangle_area_square(triangle_1b))

            # Arrangement 2
            triangle_2a = [vertices[1, 0], vertices[0, 0], vertices[0, 1]]
            triangle_2b = [vertices[1, 0], vertices[0, 1], vertices[1, 1]]
            # Measure the "greatness" of these two triangles
            area_min_2 = get_triangle_area_square(triangle_2a)
            area_min_2 = min(area_min_2, get_triangle_area_square(triangle_2b))

            # Choose the best arrangement, currently based on triangles area
            if area_min_1 < area_min_2:
                # Arrangement 2 is superior
                if len(np.unique(triangle_2a)) == 3:
                    mesh.create_doublet(triangle_2a, surface_id)
                if len(np.unique(triangle_2b)) == 3:
                    mesh.create_doublet(triangle_2b, surface_id)
            else:
                # Arrangement 1 is superior
                if len(np.unique(triangle_1a)) == 3:
                    mesh.create_doublet(triangle_1a, surface_id)
                if len(np.unique(triangle_1b)) == 3:
                    mesh.create_doublet(triangle_1b, surface_id)

    def mesh_grid(vertex_addressing, grid_offset, clockwise, verify_vertices):
        shape = np.shape(vertex_addressing)
        station_count = 0
        for i in range(0, shape[0] - 1):
            row_count = 0
            wingstation_is_initialized = False
            for j in range(0, shape[1] - 1):
                do_mesh = True
                vertices = vertex_addressing[i:i + 2, j:j + 2]
                # Transpose vertices to reverse the order of vertices if not clockwise
                if clockwise:
                    vertices = vertices.transpose()

                # If specified, make sure the vertices are valid (more than 2 unique vertices)
                if verify_vertices:
                    if len(np.unique(vertices)) < 3:
                        do_mesh = False

                if do_mesh:
                    if mode == "Structured":
                        rectangle = [vertices[0, 0], vertices[0, 1], vertices[1, 1], vertices[1, 0]]
                        # Make sure vortexes are not added if they do not have more than 3 distinct nodes

                        if is_vortex:
                            if len(np.unique(rectangle)) == 3:
                                pass
                            else:
                                # Add to row count
                                row_count += 1

                                # On the last pass, initiate new wingstation
                                if not wingstation_is_initialized:
                                    mesh.create_new_wingstation(surface_id, [])
                                    wingstation_is_initialized = True
                                    # Add to station count
                                    station_count += 1

                                station = station_count - 1 + grid_offset[0]
                                row = row_count - 1 + grid_offset[1]
                                mesh.create_vortex(rectangle, surface_id, station, row)

                        else:
                            mesh.create_doublet(rectangle, surface_id)
                    else:
                        mesh_rectangle_to_triangles(vertices)

    def get_vertices_address(x, y, z):
        shape = np.shape(x)
        vertices_nb = len(mesh.vertices)
        vertex_addressing = np.zeros(np.shape(x))
        for i in range(0, shape[0]):
            for j in range(0, shape[1]):
                vertex_index = i * shape[1] + j + vertices_nb
                vertex_addressing[i, j] = vertex_index
                mesh.append_vertex([x[i, j], y[i, j], z[i, j]])
        return vertex_addressing

    def stitch_addressing(vertex_addressing_1, vertex_addressing_2, tol):
        shape_1 = np.shape(vertex_addressing_1)
        shape_2 = np.shape(vertex_addressing_2)
        tol_square = tol ** 2
        # Run through all addresses in first grid and compare to all elements in second grid
        for i1 in range(0, shape_1[0]):
            for j1 in range(0, shape_1[1]):
                address_1 = vertex_addressing_1[i1, j1]
                vertex_1 = mesh.vertices[int(address_1)]
                for i2 in range(0, shape_2[0]):
                    for j2 in range(0, shape_2[1]):
                        address_2 = vertex_addressing_2[i2, j2]
                        vertex_2 = mesh.vertices[int(address_2)]
                        # Get distance squared
                        distance_square = (vertex_1[0] - vertex_2[0]) ** 2 + \
                                          (vertex_1[1] - vertex_2[1]) ** 2 + (vertex_1[2] - vertex_2[2]) ** 2
                        # If vertices are the same, set vertex_address_2 equal to vertex_address_1
                        if distance_square <= tol_square:
                            vertex_addressing_1[i1, j1] = address_2
        return vertex_addressing_1

    def convert_addressing_to_box_list(vertex_addressing, addressing_id, tol):
        shape = np.shape(vertex_addressing)
        list_of_boxes = []
        for i in range(0, shape[0]):
            for j in range(0, shape[1]):
                vertex = int(vertex_addressing[i, j])
                point = mesh.vertices[vertex]
                if tol > 0:
                    point_min = [point[0] - tol, point[1] - tol, point[2] - tol]
                    point_max = [point[0] + tol, point[1] + tol, point[2] + tol]
                    list_of_points = [point_min, point_max]
                else:
                    list_of_points = [point]
                list_of_boxes.append(Box([addressing_id, [i, j]], list_of_points))
        return list_of_boxes

    def set_addressing_in_containers(vertex_addressing, addressing_id, max_subdivisions, tol):
        # Extract list of PointObjects from mesh -> elements
        list_of_boxes = convert_addressing_to_box_list(vertex_addressing, addressing_id, tol=tol)
        # Create container structure
        container = generate_containers(None, list_of_boxes, max_subdivisions)
        return container

    def stitch_addressing_v1(vertex_addressing_1, vertex_addressing_2, tol, nb_subdivisions):
        # Get container of vertex_addressing_1
        print(" Step 1")
        container_1 = set_addressing_in_containers(vertex_addressing_1, 1, nb_subdivisions, 2 * tol)
        container_2 = set_addressing_in_containers(vertex_addressing_2, 2, nb_subdivisions, 2 * tol)

        print(" Step 2")
        list_of_collisions = get_collisions(container_1, container_2, [])

        tol_square = tol ** 2
        for collision in list_of_collisions:
            ids_1 = collision[1][1]
            ids_2 = collision[0][1]
            address_1 = vertex_addressing_1[ids_1[0], ids_1[1]]
            address_2 = vertex_addressing_2[ids_2[0], ids_2[1]]
            vertex_1 = mesh.vertices[int(address_1)]
            vertex_2 = mesh.vertices[int(address_2)]

            distance_square = (vertex_1[0] - vertex_2[0]) ** 2 + \
                              (vertex_1[1] - vertex_2[1]) ** 2 + (vertex_1[2] - vertex_2[2]) ** 2
            # If vertices are the same, set vertex_address_2 equal to vertex_address_1
            if distance_square <= tol_square:
                vertex_addressing_1[ids_1[0], ids_1[1]] = address_2
        return vertex_addressing_1

    def stitch_addressing_with_mesh(vertex_addressing, tol, nb_subdivisions):
        # Get container of vertex_addressing_1
        print(" Step 1")
        addr_container = set_addressing_in_containers(vertex_addressing, 5, nb_subdivisions, 2 * tol)
        mesh_container = mesh.set_used_vertices_in_containers(nb_subdivisions, 2 * tol)

        print(" Step 2")
        list_of_collisions = get_collisions(addr_container, mesh_container, [])

        tol_square = tol ** 2
        for collision in list_of_collisions:
            addr_ids = collision[0][1]
            mesh_ids = collision[1][0]
            addr_vertex_id = vertex_addressing[addr_ids[0], addr_ids[1]]
            mesh_vertex_id = mesh_ids
            addr_vertex = mesh.vertices[int(addr_vertex_id)]
            mesh_vertex = mesh.vertices[int(mesh_vertex_id)]

            distance_square = (addr_vertex[0] - mesh_vertex[0]) ** 2 + \
                              (addr_vertex[1] - mesh_vertex[1]) ** 2 + (addr_vertex[2] - mesh_vertex[2]) ** 2
            # If vertices are the same, set vertex_address_2 equal to vertex_address_1
            if distance_square <= tol_square:
                vertex_addressing[addr_ids[0], addr_ids[1]] = mesh_vertex_id
        return vertex_addressing

    distance_tol = 1e-14
    # Initiate Mesh object if not already provided
    if mesh is None:
        mesh = Mesh()

    if is_vortex:
        surface_type = "Vortex"
    else:
        surface_type = "Doublet"
    mesh.create_new_surface(surface_type)
    surface_id = mesh.surfaces[len(mesh.surfaces) - 1].id

    # Get body surfaces
    surfaces = body.get_body_surfaces()
    associativity = body.associativity
    vertex_addressing_list = []
    nb = 20

    if is_vortex:
        for stitch in associativity:
            id_1 = stitch[1][0]
            id_2 = stitch[1][1]
            if stitch[0] == "su_sl":
                surface_x = (surfaces[id_1][0] + surfaces[id_2][0]) / 2
                surface_y = (surfaces[id_1][1] + surfaces[id_2][1]) / 2
                surface_z = (surfaces[id_1][2] + surfaces[id_2][2]) / 2

                surface_addressing = get_vertices_address(surface_x, surface_y, surface_z)
                vertex_addressing_list.append(surface_addressing)

        # Stitch every surface together to correct for repeated vertices
        for ik in range(0, len(vertex_addressing_list)):
            vertex_addressing_list[ik] = stitch_addressing_with_mesh(vertex_addressing_list[ik], distance_tol, nb)
            for jk in range(ik + 1, len(vertex_addressing_list)):
                # vertex_addressing_list[ik] = stitch_addressing(vertex_addressing_list[ik],
                #                                                vertex_addressing_list[jk], distance_tol)
                vertex_addressing_list[ik] = stitch_addressing_v1(vertex_addressing_list[ik],
                                                                  vertex_addressing_list[jk], distance_tol, nb)

        # for vertex_addressing_array in vertex_addressing_list:
        current_offset = [mesh.wingstation_count, 0]
        for k in range(0, len(vertex_addressing_list)):
            mesh_grid(vertex_addressing_list[k], current_offset, clockwise=True, verify_vertices=True)
            # Update current grid offset. For now, it is as simple as the sum of previous grids shape[1] values in list
            current_offset[0] += np.shape(vertex_addressing_list[k])[0] - 1

        # Flip wingstations order if body built is marked as flipped
        if body.is_flipped:
            mesh.surfaces[surface_id].wingstations.reverse()

    else:
        for surface in surfaces:
            # Get basic addressing of surface
            surface_addressing = get_vertices_address(surface[0], surface[1], surface[2])
            vertex_addressing_list.append(surface_addressing)

        # Stitch every surface together to correct for repeated vertices
        for ik in range(0, len(vertex_addressing_list)):
            vertex_addressing_list[ik] = stitch_addressing_with_mesh(vertex_addressing_list[ik], distance_tol, nb)
            for jk in range(ik + 1, len(vertex_addressing_list)):
                print(ik, jk)
                # vertex_addressing_list[ik] = stitch_addressing(vertex_addressing_list[ik],
                #                                                vertex_addressing_list[jk], distance_tol)
                vertex_addressing_list[ik] = stitch_addressing_v1(vertex_addressing_list[ik],
                                                                  vertex_addressing_list[jk], distance_tol, nb)

        for stitch in associativity:
            if stitch[0] == "su_sl":
                su_addressing = vertex_addressing_list[stitch[1][0]]
                sl_addressing = vertex_addressing_list[stitch[1][1]]

                su_shape = np.shape(su_addressing)
                sl_shape = np.shape(sl_addressing)
                if su_shape[0] == sl_shape[0] and su_shape[1] == sl_shape[1]:
                    stitch_vertices = np.vstack([su_addressing[0, :], sl_addressing[0, :]])
                    mesh_grid(stitch_vertices, [], clockwise=True, verify_vertices=True)

                    stitch_vertices = np.vstack([su_addressing[su_shape[0] - 1, :], sl_addressing[sl_shape[0] - 1, :]])
                    mesh_grid(stitch_vertices, [], clockwise=True, verify_vertices=True)

                    stitch_vertices = np.vstack([su_addressing[:, 0], sl_addressing[:, 0]])
                    mesh_grid(stitch_vertices, [], clockwise=True, verify_vertices=True)

                    stitch_vertices = np.vstack([su_addressing[:, su_shape[1] - 1], sl_addressing[:, sl_shape[1] - 1]])
                    mesh_grid(stitch_vertices, [], clockwise=True, verify_vertices=True)
                else:
                    print("Error - Su and Sl edges do not match in length")

        # for vertex_addressing_array in vertex_addressing_list:
        for k in range(0, len(vertex_addressing_list)):
            vertex_addressing_array = vertex_addressing_list[k]
            order_vertices_clockwise = True
            if body.surfaces[k].ID == "Lower":
                order_vertices_clockwise = False
            mesh_grid(vertex_addressing_array, [], clockwise=order_vertices_clockwise, verify_vertices=True)

    return mesh


def get_surface_mesh_tables_v1(x, y, su, sl, is_lift_body):
    def get_triangle_area_square(vertices_table, vertices):
        point_a = vertices_table[int(vertices[0])]
        point_b = vertices_table[int(vertices[1])]
        point_c = vertices_table[int(vertices[2])]
        n = cross_product(point_b - point_a, point_c - point_a)
        area_square = (n[0] ** 2 + n[1] ** 2 + n[2] ** 2) / 4
        return area_square

    # def point_inside_triangle(px, py, p0x, p0y, p1x, p1y, p2x, p2y):
    #     Area = 0.5 * (-p1y * p2x + p0y * (-p1x + p2x) + p0x * (p1y - p2y) + p1x * p2y)
    #     if Area != 0:
    #         s = 1 / (2 * Area) * (p0y * p2x - p0x * p2y + (p2y - p0y) * px + (p0x - p2x) * py)
    #         t = 1 / (2 * Area) * (p0x * p1y - p0y * p1x + (p0y - p1y) * px + (p1x - p0x) * py)
    #         if s >= 0 and t >= 0 and 1-s-t >= 0:
    #             return True, s, t
    #         else:
    #             return False, s, t
    #     else:
    #         return False, 0, 0
    #
    # def insert_point(element_table, vertices, point):
    #     # Find element inside which the point resides
    #     for k in range(0, np.shape(element_table)[0]):
    #         p0 = vertices[int(element_table[k, 0]), :]
    #         p1 = vertices[int(element_table[k, 1]), :]
    #         p2 = vertices[int(element_table[k, 2]), :]
    #         is_point_inside, s, t = point_inside_triangle(point[0], point[1], p0[0], p0[1], p1[0], p1[1], p2[0], p2[1])
    #         if is_point_inside:
    #             interpolate = s * (p1[2] - p0[2]) + t * (p2[2] - p0[2]) + (1 - s - t) * (p2[2] - p1[2])
    #             print(p0[0], interpolate)

    def mesh_triangle(element_table, vertices):
        if len(vertices) == 3:
            element_table = np.append(element_table, [vertices], axis=0)
        return element_table

    def mesh_rectangle_to_triangles(vertices_table, element_table, vertices):
        if np.shape(vertices) == (2, 2):
            # Two ways to cut the shape into triangles
            # Arrangement 1
            triangle_1a = [vertices[0, 0], vertices[0, 1], vertices[1, 1]]
            triangle_1b = [vertices[0, 0], vertices[1, 1], vertices[1, 0]]
            # Measure the "greatness" of these two triangles
            area_min_1 = get_triangle_area_square(vertices_table, triangle_1a)
            area_min_1 = min(area_min_1, get_triangle_area_square(vertices_table, triangle_1b))

            # Arrangement 2
            triangle_2a = [vertices[1, 0], vertices[0, 1], vertices[0, 0]]
            triangle_2b = [vertices[1, 0], vertices[1, 1], vertices[0, 1]]
            # Measure the "greatness" of these two triangles
            area_min_2 = get_triangle_area_square(vertices_table, triangle_2a)
            area_min_2 = min(area_min_2, get_triangle_area_square(vertices_table, triangle_2b))

            # Choose the best arrangement, currently based on triangles area
            if area_min_1 < area_min_2:
                # Arrangement 2 is superior
                element_table = mesh_triangle(element_table, triangle_2a)
                element_table = mesh_triangle(element_table, triangle_2b)
            else:
                # Arrangement 1 is superior
                element_table = mesh_triangle(element_table, triangle_1a)
                element_table = mesh_triangle(element_table, triangle_1b)
        return element_table

    distance_tol = 1e-16
    # Initialize elem_table with element pointing to vertices [0, 0, 0], such that it is deleted by the cleanup
    elem_table = np.zeros([1, 3])

    if np.shape(x) == np.shape(y):
        shape = np.shape(x)

        # Flatten x and y vertices into a vertex table
        vertex_table_length = shape[0] * shape[1] * 2
        vertex_table = np.zeros([vertex_table_length, 3])
        vertex_address_su = np.zeros(np.shape(x))
        vertex_address_sl = np.zeros(np.shape(x))
        for i in range(0, shape[0]):
            for j in range(0, shape[1]):
                vertex_index_su = i * shape[1] + j
                vertex_address_su[i, j] = vertex_index_su
                vertex_table[vertex_index_su, 0] = x[i, j]
                vertex_table[vertex_index_su, 1] = y[i, j]
                vertex_table[vertex_index_su, 2] = su[i, j]

                # Check if su and sl vertices are close enough to be considered equal
                # In case not equal, identify to add element on edge
                if abs(su[i, j] - sl[i, j]) <= distance_tol:
                    vertex_index_sl = vertex_index_su
                else:
                    vertex_index_sl = (shape[0] + i) * shape[1] + j
                vertex_address_sl[i, j] = vertex_index_sl
                vertex_table[vertex_index_sl, 0] = x[i, j]
                vertex_table[vertex_index_sl, 1] = y[i, j]
                vertex_table[vertex_index_sl, 2] = sl[i, j]

        for i in range(0, shape[0] - 1):
            for j in range(0, shape[1] - 1):
                vertices_su = vertex_address_su[i:i + 2, j:j + 2]
                elem_table = mesh_rectangle_to_triangles(vertex_table, elem_table, vertices_su)

                vertices_sl = vertex_address_sl[i:i + 2, j:j + 2]
                elem_table = mesh_rectangle_to_triangles(vertex_table, elem_table, vertices_sl)

        for i in range(0, shape[0] - 1):
            for j in [0, shape[1] - 1]:
                # Check if element is needed in between surfaces
                vert_1_check = vertex_address_su[i, j] == vertex_address_sl[i, j]
                vert_2_check = vertex_address_su[i + 1, j] == vertex_address_sl[i + 1, j]
                if vert_1_check and vert_2_check:
                    # Then no element needs to be added
                    pass
                elif vert_1_check:
                    # Add 1 element
                    elem_table = mesh_triangle(elem_table, [vertex_address_su[i, j], vertex_address_su[i + 1, j],
                                                            vertex_address_sl[i + 1, j]])
                elif vert_2_check:
                    elem_table = mesh_triangle(elem_table, [vertex_address_su[i, j], vertex_address_su[i + 1, j],
                                                            vertex_address_sl[i, j]])
                else:
                    # At least 1 element needs to be added
                    side_vertices = np.reshape([vertex_address_su[i, j], vertex_address_su[i + 1, j],
                                                vertex_address_sl[i, j], vertex_address_sl[i + 1, j]], (2, 2))
                    elem_table = mesh_rectangle_to_triangles(vertex_table, elem_table, side_vertices)

        for i in [0, shape[0] - 1]:
            for j in range(0, shape[1] - 1):
                # Check if element is needed in between surfaces
                vert_1_check = vertex_address_su[i, j] == vertex_address_sl[i, j]
                vert_2_check = vertex_address_su[i, j + 1] == vertex_address_sl[i, j + 1]
                if vert_1_check and vert_2_check:
                    # Then no element needs to be added
                    pass
                elif vert_1_check:
                    # Add 1 element
                    elem_table = mesh_triangle(elem_table, [vertex_address_su[i, j], vertex_address_su[i, j + 1],
                                                            vertex_address_sl[i, j + 1]])
                elif vert_2_check:
                    elem_table = mesh_triangle(elem_table, [vertex_address_su[i, j], vertex_address_su[i, j + 1],
                                                            vertex_address_sl[i, j]])
                else:
                    # At least 1 element needs to be added
                    side_vertices = np.reshape([vertex_address_su[i, j], vertex_address_su[i, j + 1],
                                                vertex_address_sl[i, j], vertex_address_sl[i, j + 1]], (2, 2))
                    elem_table = mesh_rectangle_to_triangles(vertex_table, elem_table, side_vertices)

        # Generate lift table
        if is_lift_body:
            lift_table = np.ones([np.shape(elem_table)[0]+1, 1])
        else:
            lift_table = np.zeros([np.shape(elem_table)[0]+1, 1])

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


def get_surface_mesh_tables(x_ini, y_ini, su_ini, sl_ini, is_lift_body, flatten, append_fake_body, triangulate):
    if flatten:
        lift_bodies = set_lift_property(x_ini, is_lift_body)
        x = x_ini
        y = y_ini
        z = (su_ini + sl_ini) / 2
        if append_fake_body:
            x, y, z, lift_bodies = append_plane_body_to_wing_surface(x_ini, y_ini, z, lift_bodies)
        # z = set_lift_property(x, 0)
    else:
        x, y, z = merge_body_surfaces(x_ini, y_ini, su_ini, sl_ini)
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
        norm_table = []

    return vertex_table, elem_table, column_table, row_table, lift_table, norm_table


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


def measure_area_of_triangles(vertex_table, elem_table):
    area_table = np.zeros((np.shape(elem_table)[0], 1))
    for i in range(0, np.shape(elem_table)[0]):
        point_a = vertex_table[int(elem_table[i, 0])]
        point_b = vertex_table[int(elem_table[i, 1])]
        point_c = vertex_table[int(elem_table[i, 2])]
        n = cross_product(point_b - point_a, point_c - point_a)
        area_table[i, 0] = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2) / 2
        if area_table[i, 0] < 1e-8:
            print(i, area_table[i, 0])
    return area_table


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


def visualize_profile_lines(ax, x, y, z, color, linewidth):
    if len(np.shape(x)) == 1:
        ax.plot3D(x, y, z, linewidth=linewidth, color=color)
    elif len(np.shape(x)) == 2:
        for i in range(0, np.shape(x)[0]):
            ax.plot3D(x[i, :], y[i, :], z[i, :], linewidth=linewidth, color=color)
    return ax


def visualize_lines(ax, list_of_points, linewidth):
    color = (0, 0, 0)
    x = []
    y = []
    z = []
    for point in list_of_points:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    ax.plot3D(x, y, z, linewidth=linewidth, color=color)

    return ax


def visualize_points(ax, list_of_points, size):
    color = (0.1, 0.1, 0.1)
    x = []
    y = []
    z = []
    s = []
    for point in list_of_points:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
        s.append(size)
    ax.scatter(x, y, z, s=s, color=color, marker="x")

    return ax


def visualize_profile(ax, x, y, su, sl, color_base, alpha, linewidth):
    N1 = 0
    N2 = np.shape(x)[0]
    for i in range(N1, N2):
        color = (color_base[0]*(N2 - i)/(N2 - N1), color_base[1]*(N2 - i)/(N2 - N1), color_base[2]*(N2 - i)/(N2 - N1))
        x_full, y_full, z_full = exp.cleanup_airfoil(x[i, :], y[i, :], su[i, :], sl[i, :])
        polygon = [list(zip(x_full, y_full, z_full))]
        poly = Poly3DCollection(polygon)
        poly.set_linewidth(0)
        poly.set_alpha(alpha)
        poly.set_color(color)

        ax.plot3D(x_full, y_full, z_full, color=color, linewidth=linewidth)
        ax.add_collection3d(poly)

    return ax


def print_polygons(ax, vertices, elements, normals, color_1, color_2, alpha, line):
    if ax is None:
        fig = plt.figure()
        ax = Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(ax)

    for i in range(0, np.shape(elements)[0]):
        elem = elements[i, :]
        norm = normals[i, :]
        direction = [0, 0, 1]
        luminosity = abs(dot_product(norm, direction))
        color_mix = [0, 0, 0]
        for j in range(0, len(color_mix)):
            color_mix[j] = color_1[j] * luminosity + color_2[j] * (1 - luminosity)
        color = (color_mix[0], color_mix[1], color_mix[2])
        x = []
        y = []
        z = []
        for j in range(0, len(elem)):
            vert_id = int(elem[j])
            x.append(vertices[vert_id, 0])
            y.append(vertices[vert_id, 1])
            z.append(vertices[vert_id, 2])
        triangle = [list(zip(x, y, z))]
        poly = Poly3DCollection(triangle)
        poly.set_linewidth(line)
        poly.set_alpha(alpha)
        poly.set_facecolor(color)
        poly.set_edgecolor((0, 0, 0))
        ax.add_collection3d(poly)

    return ax


def print_mesh_as_polygons(ax, mesh, color_1, color_2, alpha, line):
    if ax is None:
        fig = plt.figure()
        ax = Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(ax)
        x_lim = [1e15, -1e15]
        y_lim = [1e15, -1e15]
        z_lim = [1e15, -1e15]
    else:
        x_lim = [ax.get_xlim()[0], ax.get_xlim()[1]]
        y_lim = [ax.get_ylim()[0], ax.get_ylim()[1]]
        z_lim = [ax.get_zlim()[0], ax.get_zlim()[1]]

    elements = mesh.elements
    vertices = mesh.vertices
    normals = mesh.get_normals()
    for i in range(0, len(elements)):
        elem = mesh.elements[i].vertices
        norm = normals[i]
        direction = [0, 0, 1]
        luminosity = abs(dot_product(norm, direction))
        color_mix = [0, 0, 0]
        for j in range(0, len(color_mix)):
            color_mix[j] = color_1[j] * luminosity + color_2[j] * (1 - luminosity)
        color = (color_mix[0], color_mix[1], color_mix[2])
        x = []
        y = []
        z = []
        for j in range(0, len(elem)):
            vert_id = int(elem[j])
            x.append(vertices[vert_id][0])
            y.append(vertices[vert_id][1])
            z.append(vertices[vert_id][2])

            # Look for larger limits
            x_lim[0] = min(x_lim[0], x[j])
            y_lim[0] = min(y_lim[0], y[j])
            z_lim[0] = min(z_lim[0], z[j])
            x_lim[1] = max(x_lim[1], x[j])
            y_lim[1] = max(y_lim[1], y[j])
            z_lim[1] = max(z_lim[1], z[j])

        # Add polygon to plot
        triangle = [list(zip(x, y, z))]
        poly = Poly3DCollection(triangle)
        poly.set_linewidth(line)
        poly.set_alpha(alpha)
        poly.set_facecolor(color)
        poly.set_edgecolor((0, 0, 0))
        ax.add_collection3d(poly)

    lim_delta_max = max(x_lim[1] - x_lim[0], y_lim[1] - y_lim[0], z_lim[1] - z_lim[0])
    ax.set_xlim((x_lim[0] + x_lim[1] - lim_delta_max) / 2, (x_lim[0] + x_lim[1] + lim_delta_max) / 2)
    ax.set_ylim((y_lim[0] + y_lim[1] - lim_delta_max) / 2, (y_lim[0] + y_lim[1] + lim_delta_max) / 2)
    ax.set_zlim((z_lim[0] + z_lim[1] - lim_delta_max) / 2, (z_lim[0] + z_lim[1] + lim_delta_max) / 2)

    return ax


def visualize_wireframe(ax, x, y, z, color, alpha, linewidth, r, c):
    ax.plot_wireframe(x, y, z, rstride=r, cstride=c, color=color, linewidth=linewidth, alpha=alpha)
    return ax


def visualize_points_as_text(ax, list_of_points, displacement, direction, size, color):
    for i in range(0, np.shape(list_of_points)[0]):
        x = list_of_points[i, 0]
        y = list_of_points[i, 1]
        z = list_of_points[i, 2]
        label = "(%.1f, %.1f, %.1f)" % (x, y, z)
        ax.text(x + displacement[0], y + displacement[1], z + displacement[2], label, zdir=direction, color=color, fontsize=size)
    return ax


def visualize_text(ax, list_of_points, labels, displacement, direction, size, color):
    for i in range(0, np.shape(list_of_points)[0]):
        x = list_of_points[i, 0]
        y = list_of_points[i, 1]
        z = list_of_points[i, 2]
        label = labels[i]
        ax.text(x + displacement[0], y + displacement[1], z + displacement[2], label, zdir=direction, color=color, fontsize=size)
    return ax


def set_plot_properties(ax):
    # ax.set_axis_off()
    # ax.set_xlim(-1.5, 3.5)
    # ax.set_ylim(0, 5)
    # ax.set_zlim(-2.5, 2.5)
    # # ax.view_init(elev=30, azim=-50)
    # # ax.dist = 7
    # ax.view_init(elev=15, azim=-90)
    # ax.dist = 3

    ax.set_axis_off()
    ax.set_xlim(-7.5, 7.5)
    ax.set_ylim(-7.5, 7.5)
    ax.set_zlim(-7.5, 7.5)
    # ax.view_init(elev=30, azim=-50)
    # ax.dist = 7
    ax.view_init(elev=30, azim=-30)
    ax.dist = 5
