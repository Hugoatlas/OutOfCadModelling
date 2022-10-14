import math
import numpy as np
import matplotlib.pyplot as plt
import export_to_gmsh


class Param:
    def __init__(self, nx, ny, b, d, c, scale_u, scale_l,
                 z_te, r_le_u, r_le_l, beta_u, beta_l,
                 x_le, z_n, delta_alpha_t,
                 n1, n2, nc1, nc2, ny1, ny2, order):
        self.Nx = nx
        self.Ny = ny

        self.B = b  # Largeur totale des deux ailes combinées
        self.D = d  # Diamètre du fuselage
        self.C_ = c  # Largeur locale de l'aile
        self.SCALE_Upper_ = scale_u  # Scale général du profile supérieur
        self.SCALE_Lower_ = scale_l  # Scale général du profile inférieur

        self.Z_TE_ = z_te  # Espacement au bord de fuite
        self.R_LE_Upper_ = r_le_u  # Rayon de courbure de la surface supérieure
        self.R_LE_Lower_ = r_le_l  # Rayon de courbure de la surface inférieure
        self.BETA_Upper_ = beta_u  # Angle de la surface supérieure au bord de fuite
        self.BETA_Lower_ = beta_l  # Angle de la surface inférieure au bord de fuite

        self.X_LE_ = x_le  # Translation en x de l'aile
        self.Z_N_ = z_n  # Translation en z de l'aile
        self.DELTA_ALPHA_T_ = delta_alpha_t  # Inclinaison de l'aile

        # Class function coefficients
        self.N1_ = n1  # Coefficient de classe du bord d'attaque
        self.N2_ = n2  # Coefficient de classe du bord de fuite
        self.NC1 = nc1  # Coefficient de classe de l'extrémité -Y de l'aile
        self.NC2 = nc2  # Coefficient de classe de l'extrémité +Y de l'aile
        self.NY1 = ny1
        self.NY2 = ny2
        self.ORDER = order  # Ordre des polynômes de Bernstein


def vector_multiply(scalar, vector):
    y = []
    for elem in vector:
        y.append(scalar * elem)
    return y


def vector_divide(scalar, vector):
    y = []
    for elem in vector:
        y.append(elem / scalar)
    return y


def normalize(array):
    res = np.zeros(np.shape(array))
    minimum = min(array)
    for i in range(0, len(array)):
        res[i] = array[i] + minimum
    maximum = max(res)
    if maximum != 0:
        for i in range(0, len(array)):
            res[i] = res[i] / maximum
    return res


def distribute_points(start, end, nb):
    points_init = np.linspace(-1, 1, nb)
    points_final = []
    for point in points_init:
        point_new = point / (1 + math.pow(point, 2))
        point_new = start + (end - start) * (0.5 + point_new)
        points_final.append(point_new)
    return points_final


def c_function(eps, ns):
    y = []
    for elem in eps:
        if ns[0] == 0 and elem == 0:
            y.append(1)
        elif ns[1] == 0 and elem == 1:
            y.append(1)
        else:
            y.append(math.pow(elem, ns[0]) * math.pow(1 - elem, ns[1]))
    return y


def bernstein_function(eps, coefficients):
    res = np.zeros(np.shape(eps))
    order = len(coefficients)
    for i in range(0, order):
        bern_i = math.factorial(order - 1) / (math.factorial(i) * math.factorial(order - 1 - i))
        for j in range(0, len(eps)):
            res[j] += bern_i * coefficients[i] * math.pow(eps[j], i) * math.pow(1 - eps[j], order - 1 - i)
    return res


def bernstein_coefficients(z_te, r_le, beta, c, order, scale):
    if order == 0:
        a = [scale]
    else:
        a = np.ones([order, 1]) * scale
        a[0] = (math.sqrt(2 * r_le / c))
        a[order - 1] = math.tan(beta) + (z_te / c)
    return a


# def derivative(x, func):
#     res = np.zeros(np.shape(x))
#     for i in range(0, len(x)):
#         if i == 0:
#             res[i] = (func[i + 1] - func[i]) / (x[i + 1] - x[i])
#         elif i == len(x) - 1:
#             res[i] = (func[i] - func[i - 1]) / (x[i] - x[i - 1])
#         else:
#             res[i] = (func[i + 1] - func[i - 1]) / (x[i + 1] - x[i - 1])
#     return res
#
#
# def find_roots(x, func):
#     f_current = func[0]
#     x_current = x[0]
#     roots = []
#     for i in range(0, len(func)):
#         if f_current * func[i] <= 0:
#             if f_current != func[i]:
#                 root = x_current + ((x[i] - x_current) / (func[i] - f_current)) * (0 - f_current)
#                 roots.append(root)
#         f_current = func[i]
#         x_current = x[i]
#     return roots
#
#
# def bernstein_coefficients_optimization(canter, eps, zeta_t, shape_coefficients, class_coefficients):
#     n = 10
#     final_shape_coefficients = shape_coefficients.copy
#     roots = []
#     for i in range(0, n):
#         # z = bernstein_function(eps, shape_coefficients)
#         z = airfoil_profile(eps, zeta_t, shape_coefficients, class_coefficients)
#         z_dot = derivative(eps, z)
#         roots = find_roots(eps, z_dot)
#         new_shape_coefficients = shape_coefficients.copy
#         if len(roots) > 0:
#             canter_position = roots[0]
#             err = canter - canter_position
#
#     return roots


def airfoil_profile(eps, zeta_n, zeta_t, delta_alpha_t, shape_coefficients, class_coefficients, class_eta):
    y = []
    shape_function = bernstein_function(eps, shape_coefficients)
    class_function = c_function(eps, class_coefficients)
    for i in range(0, len(eps)):
        p = zeta_n + (shape_function[i] * class_function[i] * class_eta) + eps[i] * (zeta_t - math.tan(delta_alpha_t))
        y.append(p)
    return y


def build_airfoil(param):
    eps = distribute_points(0, 1, param.Nx)
    eta = distribute_points(0, 1, param.Ny)

    b = param.B
    d = param.D
    c = bernstein_function(eta, param.C_)
    scale_upper = bernstein_function(eta, param.SCALE_Upper_)
    scale_lower = bernstein_function(eta, param.SCALE_Lower_)

    z_te = bernstein_function(eta, param.Z_TE_)
    r_le_upper = bernstein_function(eta, param.R_LE_Upper_)
    r_le_lower = bernstein_function(eta, param.R_LE_Lower_)
    beta_upper = bernstein_function(eta, param.BETA_Upper_)
    beta_lower = bernstein_function(eta, param.BETA_Lower_)

    x_le = bernstein_function(eta, param.X_LE_)
    z_n = bernstein_function(eta, param.Z_N_)
    delta_alpha_t = bernstein_function(eta, param.DELTA_ALPHA_T_)

    # Class function coefficients
    n1 = bernstein_function(eta, param.N1_)
    n2 = bernstein_function(eta, param.N2_)
    nc1 = param.NC1
    nc2 = param.NC2
    ny1 = param.NY1
    ny2 = param.NY2
    order = param.ORDER

    # Adjust BETA with DELTA_ALPHA
    beta_upper = beta_upper + delta_alpha_t
    beta_lower = beta_lower - delta_alpha_t

    # Adjust C with DELTA_ALPHA
    c = np.multiply(np.cos(delta_alpha_t), c)

    x = []
    y = []
    su = []
    sl = []
    class_eta = c_function(eta, [nc1, nc2])
    class_eta_y = normalize(c_function(eta, [ny1, ny2]))

    for k in range(0, len(eta)):

        x_partial = vector_multiply(c[k] * class_eta_y[k], eps) + np.ones(np.shape(eps)) * (
                x_le[k] + (c[k] / 2) * (1 - class_eta_y[k]))
        y_partial = ((eta[k] + d / b) * np.ones(np.shape(eps))) * (b / 2)

        if k == 0:
            x = [x_partial]
            y = [y_partial]
        else:
            x = np.append(x, [x_partial], axis=0)
            y = np.append(y, [y_partial], axis=0)

        zeta_t = z_te[k] / c[k]
        zeta_n = z_n[k] / c[k]

        a_upper = bernstein_coefficients(z_te[k] / 2, r_le_upper[k], beta_upper[k], c[k], order, scale_upper[k])
        a_lower = bernstein_coefficients(z_te[k] / 2, r_le_lower[k], beta_lower[k], c[k], order, scale_lower[k])

        su_partial = vector_multiply(c[k],
                                     airfoil_profile(eps, zeta_n, zeta_t / 2, delta_alpha_t[k], a_upper, [n1[k], n2[k]],
                                                     class_eta[k]))
        sl_partial = vector_multiply(-c[k],
                                     airfoil_profile(eps, -zeta_n, zeta_t / 2, -delta_alpha_t[k], a_lower,
                                                     [n1[k], n2[k]],
                                                     class_eta[k]))

        if k == 0:
            su = [su_partial]
            sl = [sl_partial]
        else:
            su = np.append(su, [su_partial], axis=0)
            sl = np.append(sl, [sl_partial], axis=0)

    return x, y, su, sl


WING_RIGHT = Param(nx=25, ny=10, b=10, d=.5,
                   c=[1.2, 0.7, 0.35, 0.3, 0.25],
                   scale_u=[0.4, 0.1],
                   scale_l=[0.0, 0.0],
                   z_te=[0.002, 0.002, 0.0],
                   r_le_u=[0.05, 0.05, 0.01],
                   r_le_l=[0.025, 0.025, 0.005],
                   beta_u=[0.4, 0.1],
                   beta_l=[-0.2, -0.05],
                   x_le=[0, 2],
                   z_n=[0.05, 0.1, 0.4],
                   delta_alpha_t=[0.05, 0.2],
                   n1=[0.5],
                   n2=[1.0],
                   nc1=0.0, nc2=0.5, ny1=0.0, ny2=0.0, order=2)

WING_LEFT = Param(nx=25, ny=10, b=-10, d=-.5,
                  c=[1.2, 0.7, 0.35, 0.3, 0.25],
                  scale_u=[0.4, 0.1],
                  scale_l=[0.0, 0.0],
                  z_te=[0.002, 0.002, 0.0],
                  r_le_u=[0.05, 0.05, 0.01],
                  r_le_l=[0.025, 0.025, 0.005],
                  beta_u=[0.4, 0.1],
                  beta_l=[-0.2, -0.05],
                  x_le=[0, 2],
                  z_n=[0.05, 0.1, 0.4],
                  delta_alpha_t=[0.05, 0.2],
                  n1=[0.5],
                  n2=[1.0],
                  nc1=0.0, nc2=0.5, ny1=0.0, ny2=0.0, order=2)

BODY = Param(nx=50, ny=50, b=1.5, d=-.75,
             c=[6],
             scale_u=[.25],
             scale_l=[.15],
             z_te=[0],
             r_le_u=[0.5],
             r_le_l=[0.25],
             beta_u=[.5],
             beta_l=[.2],
             x_le=[-1],
             z_n=[0.1],
             delta_alpha_t=[0],
             n1=[0.35],
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
                     x_le=[.75],
                     z_n=[-0.2],
                     delta_alpha_t=[0],
                     n1=[0.01],
                     n2=[0.01],
                     nc1=0.5, nc2=0.5, ny1=0.0, ny2=0.0, order=0)

ENGINE_LEFT = Param(nx=20, ny=20, b=1, d=-4.5,
                    c=[0.5],
                    scale_u=[1],
                    scale_l=[1],
                    z_te=[0],
                    r_le_u=[0.5],
                    r_le_l=[0.5],
                    beta_u=[.5],
                    beta_l=[.2],
                    x_le=[.75],
                    z_n=[-0.2],
                    delta_alpha_t=[0],
                    n1=[0.01],
                    n2=[0.01],
                    nc1=0.5, nc2=0.5, ny1=0.0, ny2=0.0, order=0)

X1, Y1, Su1, Sl1 = build_airfoil(WING_RIGHT)
X2, Y2, Su2, Sl2 = build_airfoil(WING_LEFT)
X3, Y3, Su3, Sl3 = build_airfoil(BODY)
X4, Y4, Su4, Sl4 = build_airfoil(ENGINE_RIGHT)
X5, Y5, Su5, Sl5 = build_airfoil(ENGINE_LEFT)

X_slices = []
Y_slices = []
Su_slices = []
Sl_slices = []

y_slice = 0.345
# y_slice = 2.00
# y_slice = -1.0
X1_slice, Y1_slice, Su1_slice, Sl1_slice = export_to_gmsh.slice_solid(y_slice, X1, Y1, Su1, Sl1)
X_slices, Y_slices, Su_slices, Sl_slices = export_to_gmsh.group_airfoils(X1_slice, Y1_slice,
                                                                         Su1_slice, Sl1_slice, X_slices, Y_slices,
                                                                         Su_slices, Sl_slices)

X2_slice, Y2_slice, Su2_slice, Sl2_slice = export_to_gmsh.slice_solid(y_slice, X2, Y2, Su2, Sl2)
X_slices, Y_slices, Su_slices, Sl_slices = export_to_gmsh.group_airfoils(X2_slice, Y2_slice,
                                                                         Su2_slice, Sl2_slice, X_slices, Y_slices,
                                                                         Su_slices, Sl_slices)

X3_slice, Y3_slice, Su3_slice, Sl3_slice = export_to_gmsh.slice_solid(y_slice, X3, Y3, Su3, Sl3)
X_slices, Y_slices, Su_slices, Sl_slices = export_to_gmsh.group_airfoils(X3_slice, Y3_slice,
                                                                         Su3_slice, Sl3_slice, X_slices, Y_slices,
                                                                         Su_slices, Sl_slices)

X4_slice, Y4_slice, Su4_slice, Sl4_slice = export_to_gmsh.slice_solid(y_slice, X4, Y4, Su4, Sl4)
X_slices, Y_slices, Su_slices, Sl_slices = export_to_gmsh.group_airfoils(X4_slice, Y4_slice,
                                                                         Su4_slice, Sl4_slice, X_slices, Y_slices,
                                                                         Su_slices, Sl_slices)

X5_slice, Y5_slice, Su5_slice, Sl5_slice = export_to_gmsh.slice_solid(y_slice, X5, Y5, Su5, Sl5)
X_slices, Y_slices, Su_slices, Sl_slices = export_to_gmsh.group_airfoils(X5_slice, Y5_slice,
                                                                         Su5_slice, Sl5_slice, X_slices, Y_slices,
                                                                         Su_slices, Sl_slices)
plt.plot(X1_slice, Su1_slice)
plt.plot(X1_slice, Sl1_slice)
plt.plot(X2_slice, Su2_slice)
plt.plot(X2_slice, Sl2_slice)
plt.plot(X3_slice, Su3_slice)
plt.plot(X3_slice, Sl3_slice)
plt.plot(X4_slice, Su4_slice)
plt.plot(X4_slice, Sl4_slice)
plt.plot(X5_slice, Su5_slice)
plt.plot(X5_slice, Sl5_slice)

filename = 'airfoil.geo'
export_to_gmsh.gmsh_export(filename, X_slices, Y_slices, Su_slices, Sl_slices, [0.01, 1])

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
