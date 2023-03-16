import matplotlib.pyplot as plt
import numpy as np
import math
import time

import vlm_mesh as vlm
import cst_geometry as geo
import export_to_dat as exp
from cst_geometry import Body

WING_RIGHT = Body(body_type="General")
# WING_RIGHT.add_wing_surface_cst(nx=100, ny=50, b=10, d=0,
#                                 c=[1.2, 0.7, 0.35, 0.3, 0.25],
#                                 z_te_half=[0.001],
#                                 r_le=[0.05, 0.05, 0.01],
#                                 beta=[0.4, 0.1],
#                                 x_le=[0, 2], z_n=[0.05, 0.1, 0.4], delta_alpha_t=[0.05, 0.2])
# WING_RIGHT.add_wing_surface_cst(nx=100, ny=2, b=0.1, d=10,
#                                 c=[0.25],
#                                 z_te_half=[0.001],
#                                 r_le=[0.01],
#                                 beta=[0.1],
#                                 x_le=[2, 2.02], z_n=[0.4, 0.5], delta_alpha_t=[0.2])
WING_RIGHT.add_wing_surface_naca(m=[0], p=[0], t=[12],
                                 nx=13, ny=51, b=32, d=0, c=[1],
                                 x_le=[0], z_n=[0], delta_alpha_t=[0])

WING_LEFT = Body(body_type="General")
# WING_LEFT.add_wing_surface_cst(nx=100, ny=50, b=10, d=0,
#                                c=[1.2, 0.7, 0.35, 0.3, 0.25],
#                                z_te_half=[0.001],
#                                r_le=[0.05, 0.05, 0.01],
#                                beta=[0.4, 0.1],
#                                x_le=[0, 2], z_n=[0.05, 0.1, 0.4], delta_alpha_t=[0.05, 0.2])
# WING_LEFT.add_wing_surface_cst(nx=100, ny=2, b=0.1, d=10,
#                                c=[0.25],
#                                z_te_half=[0.001],
#                                r_le=[0.01],
#                                beta=[0.1],
#                                x_le=[2, 2.02], z_n=[0.4, 0.5], delta_alpha_t=[0.2])
WING_LEFT.add_wing_surface_naca(m=[0], p=[0], t=[12],
                                nx=13, ny=51, b=32, d=0, c=[1],
                                x_le=[0], z_n=[0], delta_alpha_t=[0])
WING_LEFT.mirror_body()

BODY = Body(body_type="General")
BODY.add_general_surface_cst(nx=20, ny=20, b=1.5, d=-.75, c=[6],
                             z_te_u=[0.0], z_te_l=[0.0],
                             coefficients_u=[0.30, 0.80, 0.30, 0.30, 0.40, 0.40, 0.40, 0.40, 0.60],
                             coefficients_l=[0.20, 0.30, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10, 0.10],
                             x_le=[-1], z_n=[0.1], delta_alpha_t=[0.0],
                             n1=[0.5], n2=[0.75], nc1=0.5, nc2=0.5, ny1=0.25, ny2=0.25)
BODY.change_all_distributions("sphere", "rounded")

BALL = Body(body_type="General")
BALL.add_general_surface_cst(nx=4, ny=5, b=4, d=0, c=[0.5],
                             z_te_u=[0.0], z_te_l=[0.0],
                             coefficients_u=[2],
                             coefficients_l=[2],
                             x_le=[0], z_n=[0.5], delta_alpha_t=[0.0],
                             n1=[0.5], n2=[0.5], nc1=0.5, nc2=0.5, ny1=0.5, ny2=0.5)
BALL.change_all_distributions("distribute", "cartesian")

REV = Body(body_type="Revolution")
REV.add_general_surface_cst(nx=4, ny=5, b=1, d=0, c=[1],
                            z_te_u=[0.0], z_te_l=[0.0],
                            coefficients_u=[1.0, 0.7],
                            coefficients_l=[1.0, -0.7],
                            x_le=[0], z_n=[1.0], delta_alpha_t=[0.0],
                            n1=[0.5], n2=[1.0], nc1=0, nc2=0, ny1=0, ny2=0)

BALL_REV = Body(body_type="Revolution")
BALL_REV.add_surface_cst(nx=21, ny=11, b=0.25, d=0, c=[4],
                         z_te=[0.0], r_le=[], beta=[],
                         coefficients=[2],
                         x_le=[0.0], z_n=[0.0], delta_alpha_t=[0.0],
                         n1=[0.5], n2=[0.5], nc1=0, nc2=0, ny1=0, ny2=0, identification="Upper")

print("Start mesh")
start_time = time.time()
mesh_W = vlm.get_surface_mesh_tables_v2(WING_RIGHT, mesh=None, mode="Structured", is_vortex=True)
mesh_W = vlm.get_surface_mesh_tables_v2(WING_LEFT, mesh=mesh_W, mode="Structured", is_vortex=True)
# mesh_W = vlm.get_surface_mesh_tables_v2(BODY, mesh=mesh_W, mode="Structured", is_vortex=False)
print("--- %s seconds ---" % (time.time() - start_time))
print("Done mesh")
print(mesh_W.get_unused_vertices())

mesh_W.delete_unused_vertices()
print(mesh_W.get_unused_vertices())
exp.export_mesh("pazzy_wing_vlm.dat", mesh_W)
W_ax = vlm.print_mesh_as_polygons(None, mesh_W, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)

# W_ax = vlm.print_mesh_as_polygons(None, mesh_W, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)
# plt.gca().set_aspect('equal', adjustable='box')

# print("Start to arrange in container")
# container_W = vlm.set_mesh_elements_in_containers(mesh_W, max_subdivisions=4)
# print("Done arranging in container")

# mesh_REV = vlm.get_surface_mesh_tables_v2(BALL, mesh=None, mode="Structured", is_vortex=False)
# mesh_REV.delete_unused_vertices()
# print(mesh_REV.get_unused_vertices())
# REV_ax = vlm.print_mesh_as_polygons(None, mesh_REV, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)
# exp.export_mesh("mesh_revolution.dat", mesh_REV)


# exp.export_mesh("mesh_maybe.dat", mesh_W)
# W_ax = vlm.print_mesh_as_polygons(None, mesh_W, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)
# plt.gca().set_aspect('equal', adjustable='box')

# mesh_BR = vlm.get_surface_mesh_tables_v2(BALL, mesh=None, mode="Structured", is_vortex=True)
# B_ax = vlm.print_mesh_as_polygons(None, mesh_BR, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)
# print(len(mesh_BR.get_unused_vertices()))


# RECT_RIGHT = Body(body_type="General")
# RECT_RIGHT.add_general_surface_cst(nx=15, ny=251, b=32/math.pi, d=-16/math.pi, c=[1],
#                                    z_te_u=[0.0], z_te_l=[0.0],
#                                    coefficients_u=[1],
#                                    coefficients_l=[1],
#                                    x_le=[0.0], z_n=[0.0], delta_alpha_t=[0.0],
#                                    n1=[0.5], n2=[0.5], nc1=0.5, nc2=0.5, ny1=0.5, ny2=0.5)
#
# mesh_RECT = vlm.get_surface_mesh_tables_v2(RECT_RIGHT, mesh=None, mode="Structured", is_vortex=True)
# # mesh_RECT = vlm.get_surface_mesh_tables_v2(RECT_LEFT, mesh=mesh_RECT, mode="Structured", is_vortex=True)
# exp.export_mesh("ellipse_wing.dat", mesh_RECT)
#
# RECT_ax = vlm.print_mesh_as_polygons(None, mesh_RECT, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)
# plt.gca().set_aspect('equal', adjustable='box')




# # Generate VLM mesh
# vertex_table_1, elem_table_1, col_table_1, row_table_1, lift_table_1, _ = vlm.get_vortex_mesh(X1, Y1, Su1, Sl1, 1, True)

# mesh_wing_right = vlm.Mesh()
# for ni in range(0, np.shape(vertex_table_1)[0]):
#     mesh_wing_right.append_vertex(vertex_table_1[ni, :])
#
# for ni in range(0, np.shape(elem_table_1)[0]):
#     if lift_table_1[ni] == 1:
#         element = vlm.Mesh.Vortex(elem_table_1[ni, :], 0, col_table_1[ni], row_table_1[ni])
#     else:
#         element = vlm.Mesh.Doublet(elem_table_1[ni, :], 0)
#     mesh_wing_right.append_element(element)
#
# for ni in range(0, len(mesh_wing_right.vertices)):
#     print(mesh_wing_right.vertices[ni])
#
# for ni in range(0, len(mesh_wing_right.elements)):
#     print(mesh_wing_right.elements[ni].vertices)

# # Get intersection of wings with body
# intersect_1, vert_1, dest_1, remove_1 = vlm.check_intersections_with_cst_body(elem_table_1, vertex_table_1, BODY,
#                                                                               precision=50)
#
# # Move vertices which need moving
# vertex_table_1 = vlm.move_vertices(vert_1, dest_1, vertex_table_1)
#
# # Remove elements of wings which are fully inside body
# keep_1 = vlm.get_elements_to_keep(remove_1, elem_table_1)
# elem_table_1 = vlm.keep_table_lines(keep_1, elem_table_1)
# col_table_1 = vlm.keep_table_lines(keep_1, col_table_1)
# row_table_1 = vlm.keep_table_lines(keep_1, row_table_1)
# lift_table_1 = vlm.keep_table_lines(keep_1, lift_table_1)

# Remove unused vertices


# # Generate VLM mesh
# vertex_table_2, elem_table_2, col_table_2, row_table_2, lift_table_2, _ = vlm.get_vortex_mesh(X2, Y2, Su2, Sl2, 1, True)

# # Get intersection of wings with body
# intersect_2, vert_2, dest_2, remove_2 = vlm.check_intersections_with_cst_body(elem_table_2, vertex_table_2, BODY,
#                                                                               precision=50)
#
# # Move vertices which need moving
# vertex_table_2 = vlm.move_vertices(vert_2, dest_2, vertex_table_2)
#
# # Remove elements of wings which are fully inside body
# keep_2 = vlm.get_elements_to_keep(remove_2, elem_table_2)
# elem_table_2 = vlm.keep_table_lines(keep_2, elem_table_2)
# col_table_2 = vlm.keep_table_lines(keep_2, col_table_2)
# row_table_2 = vlm.keep_table_lines(keep_2, row_table_2)
# lift_table_2 = vlm.keep_table_lines(keep_2, lift_table_2)

# # Merge wings meshes
# vert, elem, col, row, surface, lift = vlm.merge_meshes(vertex_table_1, elem_table_1, col_table_1, row_table_1, lift_table_1,
#                                                        vertex_table_2, elem_table_2, col_table_2, row_table_2, lift_table_2)

# exp.export_vlm_mesh("mesh_maybe.dat", vert, elem, col, row, surface, lift)

# # Treat body mesh
# vertex_table_3, elem_table_3, norm_table_3, lift_bodies_3 = vlm.get_surface_mesh_tables(X3, Y3, Su3, Sl3, 0, False,
#                                                                                         False, True)

# ax2 = vlm.visualize_mesh_elements([], vert, elem, None)
# plt.gca().set_aspect('equal', adjustable='box')

# ax1 = vlm.visualize_mesh_elements([], vertex_table_3, elem_table_3, None)
# ax1 = vlm.visualize_lines(ax1, dest_1, 1)
# plt.gca().set_aspect('equal', adjustable='box')

# vertices_1, elements_1, normals_1, _ = vlm.get_surface_mesh_tables(X1, Y1, Su1, Sl1, 0, False, False, True)
# vertices_1_vlm, elements_1_vlm, columns_1_vlm, rows_1_vlm, _, normals_1_vlm = vlm.get_vortex_mesh(X1, Y1, Su1, Sl1, 1, False)
# centers_1 = geo.get_center_point(X1, Y1, Su1, Sl1)
# X1_full, Y1_full, Z1_full = vlm.merge_body_surfaces(X1, Y1, Su1, Sl1)

# vertices_2, elements_2, normals_2, _ = vlm.get_surface_mesh_tables(X2, Y2, Su2, Sl2, 0, False, False, True)
# X2_full, Y2_full, Z2_full = vlm.merge_body_surfaces(X2, Y2, Su2, Sl2)
# vertices_3, elements_3, normals_3, _ = vlm.get_surface_mesh_tables(X3, Y3, Su3, Sl3, 0, False, False, True)
# X3_full, Y3_full, Z3_full = vlm.merge_body_surfaces(X3, Y3, Su3, Sl3)

# # GRAPHIQUE STRUCTURE
# # Print wing shape
# structure_ax = vlm.print_polygons(None, vertices_1, elements_1, normals_1, color_1=[0.8, 0.8, 1.0], color_2=[0, 0, 0], alpha=0.3, line=0)
# # Print airfoils profiles
# structure_ax = vlm.visualize_profile(structure_ax, X1, Y1, Su1, Sl1, color_base=[0, 0, 1], alpha=0.3, linewidth=0.6)
# # Print centers of profiles
# structure_ax = vlm.visualize_lines(structure_ax, centers_1, 0.8)
# structure_ax = vlm.visualize_points(structure_ax, centers_1, 10)
# structure_ax = vlm.visualize_points_as_text(structure_ax, centers_1, displacement=[-0.6, 0.0, 0.1], direction=None, size=8, color=(0, 0, 0))
# # Print wireframes, body and domain
# structure_ax = vlm.visualize_wireframe(structure_ax, X1_full, Y1_full, Z1_full, color=(0.1, 0.1, 0.1), alpha=0.8, linewidth=0.4, r=1, c=10)
# structure_ax = vlm.visualize_wireframe(structure_ax, X1, Y1, -0.2*np.ones(np.shape(X1)), color=(0.1, 0.1, 0.1), alpha=0.3, linewidth=0.6, r=1, c=10)
# vlm.set_plot_properties(structure_ax)
#
#
# # GRAPHIQUE VLM
# # Print wing shape
# vlm_ax = vlm.print_polygons(None, vertices_1_vlm, elements_1_vlm, normals_1_vlm, color_1=[0.8, 0.8, 1.0], color_2=[0, 0, 0], alpha=1, line=0.4)
# # Print faces and wireframe of mesh
# # vlm_ax = vlm.print_polygons(vlm_ax, vertices_1_vlm, elements_1_vlm, normals_1_vlm, color_base=[0.0, 0.0, 0.1], alpha=1)
# vlm_ax = vlm.visualize_wireframe(vlm_ax, X1, Y1, (Su1 + Sl1) / 2, color=(0.1, 0.1, 0.1), alpha=0.8, linewidth=1.0, r=1, c=10)
# vlm_ax = vlm.visualize_wireframe(vlm_ax, X1, Y1, -0.2*np.ones(np.shape(X1)), color=(0.1, 0.1, 0.1), alpha=0.3, linewidth=0.6, r=1, c=10)
# vlm.set_plot_properties(vlm_ax)
#
#
# #### GEOMETRIE ####
# # DOMAIN
# Z1 = -0.2*np.ones(np.shape(X1))
# id_num = 0
# eps_list = np.zeros([np.shape(X1)[1], 3])
# labels_eps = []
# for ki in range(0, np.shape(X1)[1]):
#     eps_list[ki, 0] = X1[id_num, ki]
#     eps_list[ki, 1] = Y1[id_num, ki]
#     eps_list[ki, 1] = Z1[id_num, ki]
#     labels_eps.append("Eps = %.1f" % (ki / (np.shape(X1)[1] - 1)))
#
# domain_ax = vlm.print_polygons(None, [], [], [], color_1=[0.8, 0.8, 1.0], color_2=[0, 0, 0], alpha=1, line=0)
# domain_ax = vlm.visualize_wireframe(domain_ax, X1, Y1, -0.2*np.ones(np.shape(X1)), color=(0, 0, 0), alpha=0.6, linewidth=1.0, r=1, c=10)
# domain_ax = vlm.visualize_wireframe(domain_ax, X1, Y1, -0.2*np.ones(np.shape(X1)), color=(0, 0, 0), alpha=0.6, linewidth=0.4, r=1, c=1)
# # domain_ax = vlm.visualize_text(domain_ax, eps_list, labels_eps, displacement=[0.0, 0.0, -0.2], direction=None, size=8, color=(0, 0, 0))
# vlm.set_plot_properties(domain_ax)
#
# airfoils_ax = vlm.print_polygons(None, [], [], [], color_1=[0.8, 0.8, 1.0], color_2=[0, 0, 0], alpha=1, line=0)
# airfoils_ax = vlm.visualize_wireframe(airfoils_ax, X1, Y1, Z1, color=(0.1, 0.1, 0.1), alpha=0.3, linewidth=0.6, r=1, c=10)
# airfoils_ax = vlm.visualize_profile(airfoils_ax, X1, Y1, Su1, Sl1, color_base=[0, 0, 1], alpha=0.6, linewidth=0.6)
# airfoils_ax = vlm.visualize_profile_lines(airfoils_ax, X1, Y1, Su1, color=[0.4, 0.4, 0.5], linewidth=1)
# airfoils_ax = vlm.visualize_profile_lines(airfoils_ax, X1, Y1, Sl1, color=[0.2, 0.2, 0.25], linewidth=1)
# vlm.set_plot_properties(airfoils_ax)
#
# surface_ax = vlm.print_polygons(None, vertices_1, elements_1, normals_1, color_1=[0.8, 0.8, 1.0], color_2=[0, 0, 0], alpha=0.8, line=0)
# surface_ax = vlm.print_polygons(surface_ax, vertices_2, elements_2, normals_2, color_1=[0.8, 0.8, 1.0], color_2=[0, 0, 0], alpha=0.8, line=0)
# surface_ax = vlm.print_polygons(surface_ax, vertices_3, elements_3, normals_3, color_1=[0.8, 0.8, 1.0], color_2=[0, 0, 0], alpha=0.8, line=0)
# # surface_ax = vlm.visualize_wireframe(surface_ax, X1, Y1, Z1, color=(0.1, 0.1, 0.1), alpha=0.3, linewidth=0.6, r=1, c=10)
# # surface_ax = vlm.visualize_wireframe(surface_ax, X2, Y2, -0.2*np.ones(np.shape(X2)), color=(0.1, 0.1, 0.1), alpha=0.3, linewidth=0.6, r=1, c=10)
# # surface_ax = vlm.visualize_wireframe(surface_ax, X3, Y3, -0.2*np.ones(np.shape(X3)), color=(0.1, 0.1, 0.1), alpha=0.3, linewidth=0.6, r=1, c=10)
# surface_ax = vlm.visualize_profile(surface_ax, X1, Y1, Su1, Sl1, color_base=[0, 0, 1], alpha=0.3, linewidth=0.6)
# surface_ax = vlm.visualize_profile(surface_ax, X2, Y2, Su2, Sl2, color_base=[0, 0, 1], alpha=0.3, linewidth=0.6)
# # surface_ax = vlm.visualize_profile(surface_ax, X3, Y3, Su3, Sl3, color_base=[0, 0, 1], alpha=0.3, linewidth=0.6)
# surface_ax = vlm.visualize_wireframe(surface_ax, X1_full, Y1_full, Z1_full, color=(0.1, 0.1, 0.1), alpha=0.8, linewidth=0.4, r=1, c=10)
# surface_ax = vlm.visualize_wireframe(surface_ax, X2_full, Y2_full, Z2_full, color=(0.1, 0.1, 0.1), alpha=0.8, linewidth=0.4, r=1, c=10)
# # surface_ax = vlm.visualize_wireframe(surface_ax, X3_full, Y3_full, Z3_full, color=(0.1, 0.1, 0.1), alpha=0.8, linewidth=0.4, r=1, c=10)
# vlm.set_plot_properties(surface_ax)

# GRAPHIQUE BALL
# Print body shape
# vertices_B, elements_B, normals_B, _ = vlm.get_surface_mesh_tables_v1(XB, YB, SuB, SlB, False)
# mesh_B = vlm.convert_doublet_to_mesh_object(None, vertices_B, elements_B, surface=0)

# intersect_B, vert_B, dest_B, remove_B = vlm.check_intersections_with_cst_body(elements_B, vertices_B,
#                                                                               BODY, precision=10)
# # Move vertices which need moving
# vertices_B = vlm.move_vertices(vert_B, dest_B, vertices_B)
# # Remove elements of wings which are fully inside body
# keep_B = vlm.get_elements_to_keep(remove_B, elements_B)
# elements_B = vlm.keep_table_lines(keep_B, elements_B)
# normals_B = vlm.keep_table_lines(keep_B, normals_B)

# vertices_3, elements_3, normals_3, _ = vlm.get_surface_mesh_tables_v1(X3, Y3, Su3, Sl3, False)
# mesh_3 = vlm.convert_doublet_to_mesh_object(None, vertices_3, elements_3, surface=1)

# intersect_3, vert_3, dest_3, remove_3 = vlm.check_intersections_with_cst_body(elements_3, vertices_3,
#                                                                               BALL, precision=10)
# # Move vertices which need moving
# vertices_3 = vlm.move_vertices(vert_3, dest_3, vertices_3)
# # Remove elements of wings which are fully inside body
# keep_3 = vlm.get_elements_to_keep(remove_3, elements_3)
# elements_3 = vlm.keep_table_lines(keep_3, elements_3)
# normals_3 = vlm.keep_table_lines(keep_3, normals_3)

# mesh_B_copy = vlm.convert_doublet_to_mesh_object(None, vertices_B, elements_B, surface=0)
# mesh_3_copy = vlm.convert_doublet_to_mesh_object(None, vertices_3, elements_3, surface=1)

# intersect_B, intersect_3, mesh_B, mesh_3, permuted_B, permuted_3 = vlm.combine_meshes(mesh_B, mesh_3, "dd")
# ax1 = vlm.visualize_mesh_elements(intersect_B, vertices_B, elements_B, None)
# ax1 = vlm.visualize_points(ax1, mesh_B_copy.get_vertices(permuted_B), size=10)
# ax1 = vlm.visualize_mesh_elements(intersect_3, vertices_3, elements_3, ax1)
# ax1 = vlm.visualize_points(ax1, mesh_3_copy.get_vertices(permuted_3), size=10)

# body_ax = vlm.print_mesh_as_polygons(None, mesh_3, color_1=[0.8, 0.8, 1.0], color_2=[0.8, 0.8, 1.0], alpha=0.4, line=0.0)
# body_ax = vlm.print_mesh_as_polygons(None, mesh_3, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)
# body_ax = vlm.print_mesh_as_polygons(None, mesh_B_copy, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)
# body_ax = vlm.print_polygons(None, vertices_B, elements_B, normals_B, color_1=[0.8, 0.8, 1.0], color_2=[0.8, 0.8, 1.0], alpha=0.5, line=0.0)
# body_ax = vlm.print_polygons(body_ax, vertices_3, elements_3, normals_3, color_1=[0.8, 0.8, 1.0], color_2=[0.2, 0.2, 0.25], alpha=1.0, line=0.1)
# plt.gca().set_aspect('equal', adjustable='box')
# print(np.shape(elements_B))

# ax1 = vlm.visualize_mesh_elements([0], vertices_B, elements_B, None)
# ax1 = vlm.visualize_mesh_elements([0], vertices_3, elements_3, ax1)
# ax1 = vlm.visualize_lines(ax1, dest_B, 1)


# print(area_3)
# ax2 = vlm.visualize_mesh_elements([20,22,23,24,25,26,27,29], vertices_3, elements_3, None)

# domain_ax_3 = vlm.print_polygons(None, [], [], [], color_1=[0.8, 0.8, 1.0], color_2=[0, 0, 0], alpha=1, line=0)
# domain_ax_3 = vlm.visualize_wireframe(domain_ax_3, XB, YB, -0.2*np.ones(np.shape(XB)), color=(0, 0, 0), alpha=0.6, linewidth=1.0, r=1, c=10)
# domain_ax_3 = vlm.visualize_wireframe(domain_ax_3, XB, YB, -0.2*np.ones(np.shape(XB)), color=(0, 0, 0), alpha=0.6, linewidth=0.4, r=1, c=1)
# plt.gca().set_aspect('equal', adjustable='box')

# fig = plt.figure()
# ax = plt.axes(projection='3d')

# ax.plot_surface(SB[0][0], SB[0][2], SB[0][4], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(SB[0][1], SB[0][3], SB[0][5], rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# ax.plot_surface(SR[0][0], SR[0][2], SR[0][4], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(SR[0][1], SR[0][3], SR[0][5], rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# ax.plot_surface(S1[0][0], S1[0][2], S1[0][4], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(S1[0][1], S1[0][3], S1[0][5], rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# ax.plot_surface(S2[0][0], S2[0][2], S2[0][4], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(S2[0][1], S2[0][3], S2[0][5], rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# ax.plot_surface(S3[0][0], S3[0][2], S3[0][4], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(S3[0][1], S3[0][3], S3[0][5], rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# # Plot Right Wings
# ax.plot_surface(X1, Y1, Su1, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(X1, Y1, Sl1, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
#
# # Plot Left Wings
# ax.plot_surface(X2, Y2, Su2, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(X2, Y2, Sl2, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
#
# # Plot Body
# ax.plot_surface(X3, Y3, Su3, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(X3, Y3, Sl3, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
#
# # Plot Engine Right
# ax.plot_surface(X4, Y4, Su4, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(X4, Y4, Sl4, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
#
# # Plot Engine Left
# ax.plot_surface(X5, Y5, Su5, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
# ax.plot_surface(X5, Y5, Sl5, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

plt.gca().set_aspect('equal', adjustable='box')
plt.show()
