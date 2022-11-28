import math
import numpy as np


def print_vertices(f, vertex_table):
    vertex_dim = np.shape(vertex_table)
    for i in range(0, vertex_dim[0]):
        string = "NODE " + str(i)
        for j in range(0, vertex_dim[1]):
            string += " " + str(vertex_table[i, j])
        string += "\n"
        f.write(string)


def print_vortex_elements(f, elem_table, column_table, row_table, surface_id, lift_table):
    # Get a list of the different surfaces in the mesh
    ids_of_surfaces = []
    surface_count = 0

    for i in range(0, np.shape(surface_id)[0]):
        surface_already_found = False
        # Make sure the element is part of a lifting feature
        if lift_table[i] == 1:

            for j in range(0, surface_count):
                if surface_id[i] == ids_of_surfaces[j]:
                    surface_already_found = True
                    break

            if not surface_already_found:
                ids_of_surfaces.append(int(surface_id[i]))
                surface_count += 1

    stations = []
    # For every surface in mesh, get all wing stations and print them to file
    for id_surface in ids_of_surfaces:

        ids_of_stations = []
        station_count = 0
        for i in range(0, np.shape(column_table)[0]):

            # Check if element is in correct surface and is a lifting feature
            if surface_id[i] == id_surface and lift_table[i] == 1:
                station_already_found = False
                for j in range(0, station_count):
                    if column_table[i] == ids_of_stations[j][0]:
                        station_already_found = True
                        break

                # Append station id to list if it is new
                if not station_already_found:
                    ids_of_stations.append([int(column_table[i]), int(station_count), surface_id[i]])
                    station_count += 1

        stations.append(ids_of_stations)

    elem_dim = np.shape(elem_table)
    # Loop through all elements of elem_table to print vortexes
    vortex_count = 0
    for i in range(0, elem_dim[0]):

        # Check if surface element is considered a vortex
        if lift_table[i, 0] == 1:
            station_id = column_table[i, 0]
            for j in range(0, len(stations)):
                station_surface = stations[j]
                for k in range(0, len(station_surface)):
                    # print(station_surface)
                    if station_surface[k][0] == station_id and station_surface[k][2] == surface_id[i]:
                        station_id = station_surface[k][1]
                        break

            string = "VORTEX " + str(vortex_count)
            string += " " + str(int(station_id))
            string += " " + str(int(row_table[i, 0]))
            string += " " + str(int(surface_id[i]))
            for j in range(0, elem_dim[1]):
                string += " " + str(int(elem_table[i, j]))
            string += "\n"
            f.write(string)

            # Update vortex count
            vortex_count += 1


def print_doublet_elements(f, elem_table, column_table, row_table, surface_id, lift_table):
    elem_dim = np.shape(elem_table)

    # Loop through all elements of elem_table to print doublets
    doublet_count = 0
    for i in range(0, elem_dim[0]):

        # Check if surface element is considered a doublet
        if lift_table[i, 0] == 0:
            string = "DOUBLET " + str(doublet_count)
            string += " " + str(int(surface_id[i]))
            for j in range(0, elem_dim[1]):
                string += " " + str(int(elem_table[i, j]))
            string += "\n"
            f.write(string)

            # Update vortex count
            doublet_count += 1


def print_wing_stations(f, column_table, surface_id, lift_table):

    # Get a list of the different surfaces in the mesh
    ids_of_surfaces = []
    surface_count = 0

    for i in range(0, np.shape(surface_id)[0]):
        surface_already_found = False
        # Make sure the element is part of a lifting feature
        if lift_table[i] == 1:

            for j in range(0, surface_count):
                if surface_id[i] == ids_of_surfaces[j]:
                    surface_already_found = True
                    break

            if not surface_already_found:
                ids_of_surfaces.append(int(surface_id[i]))
                surface_count += 1

    # For every surface in mesh, get all wing stations and print them to file
    total_station_count = 0
    for id_surface in ids_of_surfaces:

        ids_of_stations = []
        station_count = 0
        for i in range(0, np.shape(column_table)[0]):

            # Check if element is in correct surface and is a lifting feature
            if surface_id[i] == id_surface and lift_table[i] == 1:
                station_already_found = False
                for j in range(0, station_count):
                    if column_table[i] == ids_of_stations[j][0]:
                        station_already_found = True
                        break

                # Append station id to list if it is new
                if not station_already_found:
                    ids_of_stations.append([int(column_table[i]), station_count])
                    station_count += 1

        # For every station in
        for station in ids_of_stations:
            element_count_in_station = 0

            for i in range(0, np.shape(column_table)[0]):
                # Check if element is in correct surface and correct station
                if surface_id[i] == id_surface and column_table[i] == station[0] and lift_table[i] == 1:
                    element_count_in_station += 1

            string = "WINGSTATION "
            string += str(total_station_count) + " "
            string += str(int(id_surface)) + " "
            string += str(station[1]) + " "
            string += str(element_count_in_station) + "\n"
            f.write(string)

            total_station_count += 1


def print_lifting_surfaces(f, column_table, surface_id, lift_table):
    # Get a list of the different surfaces in the mesh
    ids_of_surfaces = []
    surface_count = 0

    for i in range(0, np.shape(surface_id)[0]):
        surface_already_found = False
        # Make sure the element is part of a lifting feature
        if lift_table[i] == 1:

            for j in range(0, surface_count):
                if surface_id[i] == ids_of_surfaces[j]:
                    surface_already_found = True
                    break

            if not surface_already_found:
                ids_of_surfaces.append(int(surface_id[i]))
                surface_count += 1

    # For every surface in mesh, get all wing stations and print them to file
    total_station_count = 0
    for id_surface in ids_of_surfaces:

        ids_of_stations = []
        station_count = 0
        for i in range(0, np.shape(column_table)[0]):

            # Check if element is in correct surface and is a lifting feature
            if surface_id[i] == id_surface and lift_table[i] == 1:
                station_already_found = False
                for j in range(0, station_count):
                    if column_table[i] == ids_of_stations[j]:
                        station_already_found = True
                        break

                # Append station id to list if it is new
                if not station_already_found:
                    ids_of_stations.append(int(column_table[i]))
                    station_count += 1

        string = "LIFTINGSURFACE "
        string += str(int(id_surface)) + " "
        string += str(int(station_count)) + "\n"
        f.write(string)


def export_vlm_mesh(filename, vertex_table, elem_table, column_table, row_table, surface_id, lift_table):

    f = open(filename, "w")

    print_vertices(f, vertex_table)
    print("Done writing vertices to file")
    print_vortex_elements(f, elem_table, column_table, row_table, surface_id, lift_table)
    print("Done writing vortex elements to file")
    print_doublet_elements(f, elem_table, column_table, row_table, surface_id, lift_table)
    print("Done writing doublet elements to file")
    print_wing_stations(f, column_table, surface_id, lift_table)
    print("Done writing wing stations to file")
    print_lifting_surfaces(f, column_table, surface_id, lift_table)
    print("Done writing lifting surfaces to file")

    f.close()
