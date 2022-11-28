import matplotlib.pyplot as plt
import numpy as np

import export_to_gmsh
import vlm_mesh as vlm
import cst_geometry as geo
import export_to_dat as exp
from cst_geometry import Param

WING_RIGHT = Param(nx=20, ny=50, b=10, d=.5,
                   c=[1.2, 0.7, 0.35, 0.3, 0.25],
                   scale_u=[0.4, 0.1],
                   scale_l=[0.0, 0.0],
                   z_te=[0.002, 0.002, 0.0],
                   r_le_u=[0.05, 0.05, 0.01],
                   r_le_l=[0.025, 0.025, 0.005],
                   beta_u=[0.4, 0.1],
                   beta_l=[-0.2, -0.05],
                   coeffs_u=[], coeffs_l=[],
                   x_le=[0, 2],
                   z_n=[0.05, 0.1, 0.4],
                   delta_alpha_t=[0.05, 0.2],
                   n1=[0.5],
                   n2=[1.0],
                   nc1=0.0, nc2=0.5, ny1=0.0, ny2=0.0, order=2)

WING_LEFT = geo.mirror_body_xz(WING_RIGHT)

BODY = Param(nx=15, ny=15, b=1.5, d=-.75,
             c=[6],
             scale_u=[0.25],
             scale_l=[.15],
             z_te=[0],
             r_le_u=[0.5],
             r_le_l=[0.25],
             beta_u=[.5],
             beta_l=[.2],
             coeffs_u=[0.30, 0.80, 0.30, 0.30, 0.40, 0.40, 0.40, 0.40, 0.60],
             coeffs_l=[0.20, 0.30, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10, 0.10],
             x_le=[-1],
             z_n=[0.1],
             delta_alpha_t=[0],
             n1=[0.50],
             n2=[0.75],
             nc1=0.5, nc2=0.5, ny1=0.25, ny2=0.25, order=4)

ENGINE_RIGHT = Param(nx=20, ny=20, b=1, d=3.5,
                     c=[0.5],
                     scale_u=[1],
                     scale_l=[1],
                     z_te=[0],
                     r_le_u=[0.5],
                     r_le_l=[0.5],
                     beta_u=[.5],
                     beta_l=[.2],
                     coeffs_u=[], coeffs_l=[],
                     x_le=[.75],
                     z_n=[-0.2],
                     delta_alpha_t=[0],
                     n1=[0.01],
                     n2=[0.01],
                     nc1=0.5, nc2=0.5, ny1=0.0, ny2=0.0, order=0)

ENGINE_LEFT = geo.mirror_body_xz(ENGINE_RIGHT)

TAIL = Param(nx=50, ny=50, b=1.5, d=-.75,
             c=[6],
             scale_u=[.25],
             scale_l=[.15],
             z_te=[0],
             r_le_u=[0.5],
             r_le_l=[0.25],
             beta_u=[.5],
             beta_l=[.2],
             coeffs_u=[], coeffs_l=[],
             x_le=[-1],
             z_n=[0.1],
             delta_alpha_t=[0],
             n1=[0.35],
             n2=[0.75],
             nc1=0.5, nc2=0.5, ny1=0.25, ny2=0.25, order=4)

X1, Y1, Su1, Sl1 = geo.build_body_surfaces(geo.gen_base(WING_RIGHT, "eps"), geo.gen_base(WING_RIGHT, "eta"), WING_RIGHT)
X2, Y2, Su2, Sl2 = geo.build_body_surfaces(geo.gen_base(WING_LEFT, "eps"), geo.gen_base(WING_LEFT, "eta"), WING_LEFT)
X3, Y3, Su3, Sl3 = geo.build_body_surfaces(geo.gen_base(BODY, "eps"), geo.gen_base(BODY, "eta"), BODY)
X4, Y4, Su4, Sl4 = geo.build_body_surfaces(geo.gen_base(ENGINE_RIGHT, "eps"), geo.gen_base(ENGINE_RIGHT, "eta"),
                                           ENGINE_RIGHT)
X5, Y5, Su5, Sl5 = geo.build_body_surfaces(geo.gen_base(ENGINE_LEFT, "eps"), geo.gen_base(ENGINE_LEFT, "eta"),
                                           ENGINE_LEFT)
X6, Y6, Su6, Sl6 = geo.build_body_surfaces(geo.gen_base(TAIL, "eps"), geo.gen_base(TAIL, "eta"), TAIL)

X_slices = []
Y_slices = []
Su_slices = []
Sl_slices = []

filename = 'airfoil.geo'
slice_id = 0
X_slice = vlm.convert_to_embedded_list(X1[slice_id, :])
Y_slice = vlm.convert_to_embedded_list(Y1[slice_id, :])
Su_slice = vlm.convert_to_embedded_list(Su1[slice_id, :])
Sl_slice = vlm.convert_to_embedded_list(Sl1[slice_id, :])
X_slices, Y_slices, Su_slices, Sl_slices = export_to_gmsh.group_airfoils(X_slice, Y_slice,
                                                                         Su_slice, Sl_slice, X_slices, Y_slices,
                                                                         Su_slices, Sl_slices)
export_to_gmsh.gmsh_export(filename, [X_slice], [Y_slice], [Su_slice], [Sl_slice], [0.003, 3])

# Generate VLM mesh
vertex_table_1, elem_table_1, col_table_1, row_table_1, lift_table_1 = vlm.get_vortex_mesh(X1, Y1, Su1, Sl1, 1, True)

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


# Generate VLM mesh
vertex_table_2, elem_table_2, col_table_2, row_table_2, lift_table_2 = vlm.get_vortex_mesh(X2, Y2, Su2, Sl2, 1, True)

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

# Merge wings meshes
vert, elem, col, row, surface, lift = vlm.merge_meshes(vertex_table_1, elem_table_1, col_table_1, row_table_1, lift_table_1,
                                                       vertex_table_2, elem_table_2, col_table_2, row_table_2, lift_table_2)

exp.export_vlm_mesh("mesh_maybe.dat", vert, elem, col, row, surface, lift)

# Treat body mesh
vertex_table_3, elem_table_3, norm_table_3, lift_bodies_3 = vlm.get_surface_mesh_tables(X3, Y3, Su3, Sl3, 0, False,
                                                                                        False, True)

# ax2 = vlm.visualize_mesh_elements([], vert, elem, None)
# plt.gca().set_aspect('equal', adjustable='box')

# ax1 = vlm.visualize_mesh_elements([], vertex_table_3, elem_table_3, None)
# ax1 = vlm.visualize_lines(ax1, dest_1)
# plt.gca().set_aspect('equal', adjustable='box')

fig = plt.figure()
ax = plt.axes(projection='3d')

# Plot Right Wings
ax.plot_surface(X1, Y1, Su1, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.plot_surface(X1, Y1, Sl1, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# Plot Left Wings
ax.plot_surface(X2, Y2, Su2, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.plot_surface(X2, Y2, Sl2, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# Plot Body
ax.plot_surface(X3, Y3, Su3, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.plot_surface(X3, Y3, Sl3, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# Plot Engine Right
ax.plot_surface(X4, Y4, Su4, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.plot_surface(X4, Y4, Sl4, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# Plot Engine Left
ax.plot_surface(X5, Y5, Su5, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.plot_surface(X5, Y5, Sl5, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

plt.gca().set_aspect('equal', adjustable='box')
plt.show()
