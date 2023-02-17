import math
import numpy as np
import export_to_gmsh
import vlm_mesh as vlm
import matplotlib.pyplot as plt


class Body:
    class CSTAirfoil:
        def __init__(self, c, z_te, r_le, beta, coefficients, n1, n2, nc1, nc2):
            self.type = "CST"
            self.C = list(c)
            self.Z_TE = list(z_te)
            self.R_LE = list(r_le)
            self.BETA = list(beta)
            self.COEFFICIENTS = coefficients

            # Class function coefficients
            self.N1 = list(n1)
            self.N2 = list(n2)
            self.NC1 = nc1
            self.NC2 = nc2

        def build_normalized_surface(self, eps, eta):
            def get_airfoil_profile(eps_list, zeta_t, shape_coefficients, class_coefficients, eta_class):
                y = []
                shape_function = bernstein_function(eps_list, shape_coefficients)
                class_function = c_function(eps_list, class_coefficients)
                for k in range(0, len(eps_list)):
                    p = (shape_function[k] * class_function[k] * eta_class) + eps_list[k] * zeta_t
                    y.append(p)
                if len(eps_list) == 1:
                    res = y[0]
                else:
                    res = y
                return res

            # eps and eta are 2D numpy arrays
            dim_eps = np.shape(eps)
            dim_eta = np.shape(eta)

            if dim_eps[0] == dim_eta[0]:
                x = []
                s = []

                if len(dim_eps) == 1:
                    num = 1
                else:
                    num = dim_eps[0]

                for i in range(0, num):
                    if len(dim_eps) == 1:
                        eps_i = eps
                        eta_i = eta
                    else:
                        eps_i = eps[i, :]
                        eta_i = eta[i, :]

                    c = bernstein_function(eta_i, self.C)
                    z_te = bernstein_function(eta_i, self.Z_TE)
                    r_le = bernstein_function(eta_i, self.R_LE)
                    beta = bernstein_function(eta_i, self.BETA)
                    n1 = bernstein_function(eta_i, self.N1)
                    n2 = bernstein_function(eta_i, self.N2)
                    class_eta = c_function(eta_i, [self.NC1, self.NC2])

                    if len(eps_i) == len(eta_i):
                        x_i = []
                        s_i = []
                        for j in range(0, len(eta_i)):
                            eps_ij = eps_i[j]

                            # If COEFFICIENTS attribute is empty, use r_le and beta
                            if len(self.COEFFICIENTS) == 0:
                                a = bernstein_coefficients(z_te[j], r_le[j], beta[j], c[j], 2, 0)
                            else:
                                a = self.COEFFICIENTS

                            x_ij = c[j] * eps_ij
                            s_ij = c[j] * get_airfoil_profile([eps_ij], z_te[j] / c[j], a, [n1[j], n2[j]], class_eta[j])

                            x_i = append_to_np_array(x_i, x_ij)
                            s_i = append_to_np_array(s_i, s_ij)
                        x = append_to_np_array(x, x_i)
                        s = append_to_np_array(s, s_i)
                return x, s

        def flip(self):
            self.C.reverse()
            self.Z_TE.reverse()
            self.R_LE.reverse()
            self.BETA.reverse()
            self.N1.reverse()
            self.N2.reverse()
            nc1 = self.NC2
            nc2 = self.NC1
            self.NC1 = nc1
            self.NC2 = nc2

    class NACAAirfoil:
        def __init__(self, m, p, t):
            self.type = "NACA"
            self.M = list(m)
            self.P = list(p)
            self.T = list(t)

        def flip(self):
            self.M.reverse()
            self.P.reverse()
            self.T.reverse()

    class Surface:
        def __init__(self, airfoil_data, nx, ny, b, d, c, x_le, z_n, delta_alpha_t, nc1, nc2, ny1, ny2, identification):
            self.Nx = nx
            self.Ny = ny
            if airfoil_data.type == "CST" or airfoil_data.type == "NACA":
                self.airfoil_data = airfoil_data
            else:
                print("ERROR - airfoil_data is not of the correct type")

            self.B = b
            self.D = d
            self.C = list(c)
            self.X_LE = list(x_le)
            self.Z_N = list(z_n)
            self.DELTA_ALPHA_T = list(delta_alpha_t)
            self.NC1 = nc1
            self.NC2 = nc2
            self.NY1 = ny1
            self.NY2 = ny2
            self.ID = identification

            # Get adequate distributions
            self.x_distribution = "cartesian"
            self.y_distribution = "cartesian"
            self.eps = np.zeros(0)
            self.eta = np.zeros(0)
            # Set eps and eta distributions
            self.set_distributions()

        def flip_on_xz_plane(self):
            self.C.reverse()
            self.X_LE.reverse()
            self.Z_N.reverse()
            self.DELTA_ALPHA_T.reverse()
            self.D = - (self.B + self.D)
            self.airfoil_data.flip()
            nc1 = self.NC2
            nc2 = self.NC1
            self.NC1 = nc1
            self.NC2 = nc2

        def change_distributions(self, x_dist, y_dist):
            self.x_distribution = x_dist
            self.y_distribution = y_dist
            self.set_distributions()

        def build_complete_surface(self, surface_type):
            eps = self.eps
            eta = self.eta
            u, r = self.airfoil_data.build_normalized_surface(eps, eta)
            if self.ID == "Lower":
                r = -r

            dim_eps = np.shape(eps)
            x = []
            y = []
            s = []
            for i in range(0, dim_eps[0]):
                x_i = []
                y_i = []
                s_i = []
                c = bernstein_function(eta[i, :], self.C)
                x_le = bernstein_function(eta[i, :], self.X_LE)
                z_n = bernstein_function(eta[i, :], self.Z_N)
                delta_alpha_t = bernstein_function(eta[i, :], self.DELTA_ALPHA_T)
                class_eta_y = vector_divide(c_function_max([self.NY1, self.NY2]),
                                            c_function(eta[i, :], [self.NY1, self.NY2]))

                for j in range(0, dim_eps[1]):
                    # Rotate with alpha_delta_t and z_n
                    x_ij = + math.cos(math.pi * delta_alpha_t[j]) * u[i, j] + math.sin(math.pi * delta_alpha_t[j]) * r[i, j]
                    r_ij = - math.sin(math.pi * delta_alpha_t[j]) * u[i, j] + math.cos(math.pi * delta_alpha_t[j]) * r[i, j] + z_n[j]

                    if surface_type == "Revolution":
                        x_ij = x_ij + x_le[j]
                        y_ij = (self.D + self.B * math.cos(2 * math.pi * eta[i, j]) * r_ij) / 2
                        s_ij = (self.D + self.B * math.sin(2 * math.pi * eta[i, j]) * r_ij) / 2
                    else:
                        x_ij = (class_eta_y[j] * x_ij) + (x_le[j] + (c[j] / 2) * (1 - class_eta_y[j]))
                        y_ij = (self.D + self.B * eta[i, j]) / 2
                        s_ij = r_ij

                    x_i = append_to_np_array(x_i, x_ij)
                    y_i = append_to_np_array(y_i, y_ij)
                    s_i = append_to_np_array(s_i, s_ij)
                x = append_to_np_array(x, x_i)
                y = append_to_np_array(y, y_i)
                s = append_to_np_array(s, s_i)
            return x, y, s

        def set_distributions(self):
            mode = self.x_distribution
            shape = self.y_distribution

            if mode == "distribute":
                eps_base = distribute_points(0, 1, self.Nx, mode)
                eta_base = distribute_points(0, 1, self.Ny, mode)
            elif mode == "partial":
                eps_base = distribute_points(0, 1, self.Nx, mode)
                eta_base = np.linspace(0, 1, self.Ny)
            elif mode == "crushed":
                delta = 0.005
                eps_base = distribute_points(0, 1, self.Nx, mode)
                eta_base = np.linspace(0 + delta, 1 - delta, self.Ny)
            elif mode == "sphere":
                eps_base = distribute_points(0, 1, self.Nx, mode)
                eta_base = distribute_points(0, 1, self.Ny, mode)
            else:
                eps_base = np.linspace(0, 1, self.Nx)
                eta_base = np.linspace(0, 1, self.Ny)
            eps = []
            eta = []

            for i in range(0, len(eta_base)):
                eps_i = eps_base
                eta_i = eta_base[i] * np.ones(np.shape(eps_base))

                if i == 0:
                    eps = [eps_i]
                    eta = [eta_i]
                else:
                    eps = np.append(eps, [eps_i], axis=0)
                    eta = np.append(eta, [eta_i], axis=0)

            if shape == "rounded":
                for i in range(0, np.shape(eps)[0]):
                    for j in range(0, np.shape(eps)[1]):
                        y_ij = eta[i, j] - 0.5
                        x_ij = eps[i, j] - 0.5

                        # Correct positions for equal distance on circumference
                        if abs(x_ij) < abs(y_ij):
                            if y_ij != 0:
                                correction = abs(math.tan(y_ij))
                                if y_ij - x_ij * correction != 0:
                                    delta_x = correction * (y_ij ** 2 - x_ij ** 2) / (y_ij - x_ij * correction)
                                    delta_x = (x_ij / y_ij) * delta_x
                                    x_ij = x_ij - delta_x
                        else:
                            if x_ij != 0:
                                correction = abs(math.tan(x_ij))
                                if x_ij - y_ij * correction != 0:
                                    delta_y = correction * (x_ij ** 2 - y_ij ** 2) / (x_ij - y_ij * correction)
                                    delta_y = (y_ij / x_ij) * delta_y
                                    y_ij = y_ij - delta_y

                        if x_ij ** 2 + y_ij ** 2 == 0:
                            x_ij_new = 0
                            y_ij_new = 0
                        else:
                            y_ij_new = y_ij * max(abs(x_ij), abs(y_ij)) / math.sqrt(x_ij ** 2 + y_ij ** 2)
                            if abs(y_ij_new) == 0.5:
                                x_ij_new = 0
                            else:
                                x_ij_new = x_ij * max(abs(x_ij), abs(y_ij)) / math.sqrt(
                                    (x_ij ** 2 + y_ij ** 2) * (1 - 4 * (y_ij_new ** 2)))

                        eta_ij = y_ij_new + 0.5
                        eps_ij = x_ij_new + 0.5

                        eta[i, j] = max(0, min(1, eta_ij))
                        eps[i, j] = max(0, min(1, eps_ij))

            self.eps = eps
            self.eta = eta

    def __init__(self, body_type):
        self.type = body_type
        self.surfaces = []
        self.associativity = []

    def set_to_general_solid(self):
        self.type = "General"

    def set_to_revolution_solid(self):
        self.type = "Revolution"

    def add_surface_cst(self, nx, ny, b, d, c, z_te, r_le, beta, coefficients, x_le, z_n, delta_alpha_t,
                        n1, n2, nc1, nc2, ny1, ny2, identification):
        airfoil_data = self.CSTAirfoil(c, z_te, r_le, beta, coefficients, n1, n2, nc1, nc2)
        self.surfaces.append(self.Surface(airfoil_data, nx, ny, b, d, c, x_le, z_n, delta_alpha_t, nc1, nc2, ny1, ny2, identification))

    def add_surface_naca(self, m, p, t, nx, ny, b, d, c, x_le, z_n, delta_alpha_t, nc1, nc2, ny1, ny2, identification):
        airfoil_data = self.NACAAirfoil(m, p, t)
        self.surfaces.append(self.Surface(airfoil_data, nx, ny, b, d, c, x_le, z_n, delta_alpha_t, nc1, nc2, ny1, ny2, identification))

    def add_wing_surface_cst(self, nx, ny, b, d, c, z_te_half, r_le, beta, x_le, z_n, delta_alpha_t):
        # Get current number of surfaces in Body
        # New surfaces will be stored at indices nb (upper) and nb+1 (lower)
        nb = len(self.surfaces)
        # Set associativity of surfaces
        self.associativity.append(["su_sl", [nb, nb+1]])

        # Append surfaces to list
        self.add_surface_cst(nx, ny, b, d, c, z_te_half, r_le, beta, [],
                             x_le, z_n, delta_alpha_t, [0.5], [1.0], 0, 0, 0, 0, "Upper")
        self.add_surface_cst(nx, ny, b, d, c, -1*z_te_half, r_le, -1*beta, [],
                             x_le, z_n, delta_alpha_t, [0.5], [1.0], 0, 0, 0, 0, "Lower")

    def add_general_surface_cst(self, nx, ny, b, d, c, z_te_u, z_te_l, coefficients_u, coefficients_l,
                                x_le, z_n, delta_alpha_t, n1, n2, nc1, nc2, ny1, ny2):
        # Get current number of surfaces in Body
        # New surfaces will be stored at indices nb (upper) and nb+1 (lower)
        nb = len(self.surfaces)
        # Set associativity of surfaces
        self.associativity.append(["su_sl", [nb, nb+1]])

        # Append surfaces to list
        self.add_surface_cst(nx, ny, b, d, c, z_te_u, [], [], coefficients_u,
                             x_le, z_n, delta_alpha_t, n1, n2, nc1, nc2, ny1, ny2, "Upper")
        self.add_surface_cst(nx, ny, b, d, c, z_te_l, [], [], coefficients_l,
                             x_le, z_n, delta_alpha_t, n1, n2, nc1, nc2, ny1, ny2, "Lower")

    def add_wing_surface_naca(self, m, p, t, nx, ny, b, d, c, x_le, z_n, delta_alpha_t):
        # Get current number of surfaces in Body
        # New surfaces will be stored at indices nb (upper) and nb+1 (lower)
        nb = len(self.surfaces)
        # Set associativity of surfaces
        self.associativity.append(["su_sl", [nb, nb + 1]])

        # Append surfaces to list
        self.add_surface_naca(m, p, t, nx, ny, b, d, c, x_le, z_n, delta_alpha_t, 0, 0, 0, 0, "Upper")
        self.add_surface_naca(m, p, t, nx, ny, b, d, c, x_le, z_n, delta_alpha_t, 0, 0, 0, 0, "Lower")

    def mirror_body(self):
        for surface in self.surfaces:
            surface.flip_on_xz_plane()

    def change_all_distributions(self, x_dist, y_dist):
        for surface in self.surfaces:
            surface.change_distributions(x_dist, y_dist)

    def get_paired_body_surfaces(self):
        surfaces = []
        for pair in self.associativity:
            su_id = int(pair[1][0])
            sl_id = int(pair[1][1])
            surface_type = self.type
            Xu, Yu, Zu = self.surfaces[su_id].build_complete_surface(surface_type)
            Xl, Yl, Zl = self.surfaces[sl_id].build_complete_surface(surface_type)
            surfaces.append([Xu, Xl, Yu, Yl, Zu, Zl])
        return surfaces

    def get_body_surfaces(self):
        surfaces = []
        for surface in self.surfaces:
            surface_type = self.type
            X, Y, Z = surface.build_complete_surface(surface_type)
            surfaces.append([X, Y, Z])
        return surfaces


# class Param:
#
#     def __init__(self, nx, ny, b, d, c, scale_u, scale_l,
#                  z_te, r_le_u, r_le_l, beta_u, beta_l, coeffs_u, coeffs_l,
#                  x_le, z_n, delta_alpha_t,
#                  n1, n2, nc1, nc2, ny1, ny2, order):
#         self.Nx = nx
#         self.Ny = ny
#
#         self.B = b  # Largeur totale des deux ailes combinées
#         self.D = d  # Diamètre du fuselage
#         self.C_ = list(c)  # Largeur locale de l'aile
#         self.SCALE_Upper_ = list(scale_u)  # Scale général du profile supérieur
#         self.SCALE_Lower_ = list(scale_l)  # Scale général du profile inférieur
#
#         self.Z_TE_ = list(z_te)  # Espacement au bord de fuite
#         self.R_LE_Upper_ = list(r_le_u)  # Rayon de courbure de la surface supérieure
#         self.R_LE_Lower_ = list(r_le_l)  # Rayon de courbure de la surface inférieure
#         self.BETA_Upper_ = list(beta_u)  # Angle de la surface supérieure au bord de fuite
#         self.BETA_Lower_ = list(beta_l)  # Angle de la surface inférieure au bord de fuite
#         self.COEFFICIENTS_Upper = list(coeffs_u)    # Coefficients du polynôme de Bernstein de la surface supérieure
#         self.COEFFICIENTS_Lower = list(coeffs_l)  # Coefficients du polynôme de Bernstein de la surface supérieure
#
#         self.X_LE_ = list(x_le)  # Translation en x de l'aile
#         self.Z_N_ = list(z_n)  # Translation en z de l'aile
#         self.DELTA_ALPHA_T_ = list(delta_alpha_t)  # Inclinaison de l'aile
#
#         # Class function coefficients
#         self.N1_ = list(n1)  # Coefficient de classe du bord d'attaque
#         self.N2_ = list(n2)  # Coefficient de classe du bord de fuite
#         self.NC1 = nc1  # Coefficient de classe de l'extrémité -Y de l'aile
#         self.NC2 = nc2  # Coefficient de classe de l'extrémité +Y de l'aile
#         self.NY1 = ny1
#         self.NY2 = ny2
#         self.ORDER = order  # Ordre des polynômes de Bernstein
#
#
# def mirror_body_xz(param):
#     D = - (param.B + param.D)
#     C_ = vlm.reverse_array(param.C_)
#     SCALE_Upper_ = vlm.reverse_array(param.SCALE_Upper_)
#     SCALE_Lower_ = vlm.reverse_array(param.SCALE_Lower_)
#     Z_TE_ = vlm.reverse_array(param.Z_TE_)
#     R_LE_Upper_ = vlm.reverse_array(param.R_LE_Upper_)
#     R_LE_Lower_ = vlm.reverse_array(param.R_LE_Lower_)
#     BETA_Upper_ = vlm.reverse_array(param.BETA_Upper_)
#     BETA_Lower_ = vlm.reverse_array(param.BETA_Lower_)
#     COEFFICIENTS_Upper = param.COEFFICIENTS_Upper
#     COEFFICIENTS_Lower = param.COEFFICIENTS_Lower
#     X_LE_ = vlm.reverse_array(param.X_LE_)
#     Z_N_ = vlm.reverse_array(param.Z_N_)
#     DELTA_ALPHA_T_ = vlm.reverse_array(param.DELTA_ALPHA_T_)
#     N1_ = vlm.reverse_array(param.N1_)
#     N2_ = vlm.reverse_array(param.N2_)
#     NC1 = param.NC2
#     NC2 = param.NC1
#
#     param_mirror = Param(nx=param.Nx, ny=param.Ny, b=param.B, d=D, c=C_, scale_u=SCALE_Upper_, scale_l=SCALE_Lower_,
#                          z_te=Z_TE_, r_le_u=R_LE_Upper_, r_le_l=R_LE_Lower_, beta_u=BETA_Upper_, beta_l=BETA_Lower_,
#                          coeffs_u=COEFFICIENTS_Upper, coeffs_l=COEFFICIENTS_Lower, x_le=X_LE_, z_n=Z_N_,
#                          delta_alpha_t=DELTA_ALPHA_T_, n1=N1_, n2=N2_, nc1=NC1, nc2=NC2,
#                          ny1=param.NY1, ny2=param.NY2, order=param.ORDER)
#     return param_mirror
#
#
def append_to_np_array(array, element):
    if len(array) == 0:
        res = [element]
    else:
        res = np.append(array, [element], axis=0)
    return res
#
#
# def vector_multiply(scalar, vector):
#     y = []
#     for elem in vector:
#         y.append(scalar * elem)
#     return y
#
#
def vector_divide(scalar, vector):
    y = []
    for elem in vector:
        y.append(elem / scalar)
    return y
#
#
# def normalize(array):
#     res = np.zeros(np.shape(array))
#     minimum = min(array)
#     for i in range(0, len(array)):
#         res[i] = array[i] + minimum
#     maximum = max(res)
#     if maximum != 0:
#         for i in range(0, len(array)):
#             res[i] = res[i] / maximum
#     return res
#
#
def c_function(eps, ns):
    y = []
    for elem in eps:
        if ns[0] == 0 and elem == 0:
            y.append(1)
        elif ns[1] == 0 and elem == 1:
            y.append(1)
        else:
            y.append(math.pow(elem, ns[0]) * math.pow((1 - elem), ns[1]))
    return y
#
#
def c_function_max(ns):
    eps_max = 0.5
    y_max = math.pow(eps_max, ns[0]) * math.pow((1 - eps_max), ns[1])
    return y_max
#
#
def bernstein_function(eps, coefficients):
    res = np.zeros(np.shape(eps))
    order = len(coefficients)
    for i in range(0, order):
        bern_i = math.factorial(order - 1) / (math.factorial(i) * math.factorial(order - 1 - i))
        for j in range(0, len(eps)):
            res[j] += bern_i * coefficients[i] * math.pow(eps[j], i) * math.pow(1 - eps[j], order - 1 - i)
    return res
#
#
# def is_point_out_of_bounds(value, minimum, maximum):
#     if value >= maximum or value <= minimum:
#         return True
#     else:
#         return False
#
#
def bernstein_coefficients(z_te, r_le, beta, c, order, scale):
    if order == 0:
        a = [scale]
    else:
        a = np.ones([order, 1]) * scale
        a[0] = (math.sqrt(2 * r_le / c))
        a[order - 1] = math.tan(beta) + (z_te / c)
    return a
#
#
# def derivative(x, func):
#     res = np.zeros(np.shape(x))
#     for i in range(0, len(x)):
#         if i == 0:
#             if x[i + 1] == x[i]:
#                 res[i] = None
#             else:
#                 res[i] = (func[i + 1] - func[i]) / (x[i + 1] - x[i])
#         elif i == len(x) - 1:
#             if x[i] == x[i - 1]:
#                 res[i] = None
#             else:
#                 res[i] = (func[i] - func[i - 1]) / (x[i] - x[i - 1])
#         else:
#             if x[i + 1] == x[i - 1]:
#                 res[i] = None
#             else:
#                 res[i] = (func[i + 1] - func[i - 1]) / (x[i + 1] - x[i - 1])
#     return res
#
#
# # def find_roots(x, func):
# #     f_current = func[0]
# #     x_current = x[0]
# #     roots = []
# #     for i in range(0, len(func)):
# #         if f_current * func[i] <= 0:
# #             if f_current != func[i]:
# #                 root = x_current + ((x[i] - x_current) / (func[i] - f_current)) * (0 - f_current)
# #                 roots.append(root)
# #         f_current = func[i]
# #         x_current = x[i]
# #     return roots
# #
# #
# # def bernstein_coefficients_optimization(canter, eps, zeta_t, shape_coefficients, class_coefficients):
# #     n = 10
# #     final_shape_coefficients = shape_coefficients.copy
# #     roots = []
# #     for i in range(0, n):
# #         # z = bernstein_function(eps, shape_coefficients)
# #         z = airfoil_profile(eps, zeta_t, shape_coefficients, class_coefficients)
# #         z_dot = derivative(eps, z)
# #         roots = find_roots(eps, z_dot)
# #         new_shape_coefficients = shape_coefficients.copy
# #         if len(roots) > 0:
# #             canter_position = roots[0]
# #             err = canter - canter_position
# #
# #     return roots
#
#
# def airfoil_profile(eps, zeta_n, zeta_t, delta_alpha_t, shape_coefficients, class_coefficients, class_eta):
#     y = []
#     shape_function = bernstein_function(eps, shape_coefficients)
#     class_function = c_function(eps, class_coefficients)
#     for i in range(0, len(eps)):
#         p = zeta_n + (shape_function[i] * class_function[i] * class_eta) + eps[i] * (zeta_t - math.tan(delta_alpha_t))
#         y.append(p)
#
#     if len(eps) == 1:
#         res = y[0]
#     else:
#         res = y
#     return res
#
#
# def get_domain_derivative(eps, eta, body):
#     dim_eps = np.shape(eps)
#     dim_eta = np.shape(eta)
#
#     Xd_eps = np.zeros(np.shape(eps))
#     Xd_eta = np.zeros(np.shape(eps))
#     Yd_eps = np.zeros(np.shape(eps))
#     Yd_eta = np.zeros(np.shape(eps))
#
#     if dim_eps == dim_eta:
#         X, Y, Su, Sl = build_body_surfaces_complete(eps, eta, body)
#
#         for i in range(0, dim_eps[0]):
#             Xd_eps[i, :] = derivative(eps[i, :], X[i, :])
#             Yd_eps[i, :] = derivative(eps[i, :], Y[i, :])
#         for j in range(0, dim_eps[1]):
#             Xd_eta[:, j] = derivative(eta[:, j], X[:, j])
#             Yd_eta[:, j] = derivative(eta[:, j], Y[:, j])
#
#     return Xd_eps, Xd_eta, Yd_eps, Yd_eta
#
#
# def get_gradient_norm_square(grad_x, grad_y):
#     dim_x = np.shape(grad_x)
#     dim_y = np.shape(grad_y)
#     norm = np.zeros(np.shape(grad_x))
#     if dim_x == dim_y:
#         for i in range(0, dim_x[0]):
#             for j in range(0, dim_x[1]):
#                 norm[i, j] = math.sqrt(grad_x[i, j]**2 + grad_y[i, j]**2)
#     return norm
#
#
# def get_full_domain_boundary(nx, ny, body):
#     list_of_points = []
#     list_of_coords = []
#     eta = 0
#     eps = 0
#     for i in range(1, nx):
#         eps = i / (nx - 1)
#         x_, y_, _, _ = build_body_surfaces_complete([eps], [eta], body)
#         list_of_points.append([x_, y_])
#         list_of_coords.append([eps, eta])
#
#     for i in range(1, ny):
#         eta = i / (ny - 1)
#         x_, y_, _, _ = build_body_surfaces_complete([eps], [eta], body)
#         list_of_points.append([x_, y_])
#         list_of_coords.append([eps, eta])
#
#     for i in range(1, nx):
#         eps = 1 - (i / (nx - 1))
#         x_, y_, _, _ = build_body_surfaces_complete([eps], [eta], body)
#         list_of_points.append([x_, y_])
#         list_of_coords.append([eps, eta])
#
#     for i in range(1, ny):
#         eta = 1 - (i / (ny - 1))
#         x_, y_, _, _ = build_body_surfaces_complete([eps], [eta], body)
#         list_of_points.append([x_, y_])
#         list_of_coords.append([eps, eta])
#
#     list_of_points.append(list_of_points[0])
#     list_of_coords.append(list_of_coords[0])
#     return list_of_points, list_of_coords
#
#
# def get_domain_boundary(x, y, dir_x, dir_y, body, max_attempts):
#     res = get_domain_position_xy(x, y, body)
#     eps = res[0]
#     eta = res[1]
#     next_res = None
#     for i in range(0, max_attempts):
#         new_x = x + (dir_x / (i + 1))
#         new_y = y + (dir_y / (i + 1))
#         next_res = get_domain_position_xy(new_x, new_y, body)
#         if next_res is not None:
#             break
#
#     if next_res is not None:
#         next_eps = next_res[0]
#         next_eta = next_res[1]
#
#         if eps != 0 and eps != 1 and eta != 0 and eta != 1:
#             # Get t for eps
#             if next_eps > eps:
#                 t_eps = (1 - eps) / (next_eps - eps)
#             elif next_eps < eps:
#                 t_eps = eps / (eps - next_eps)
#             else:
#                 t_eps = 0
#             # Get t for eta
#             if next_eta > eta:
#                 t_eta = (1 - eta) / (next_eta - eta)
#             elif next_eta < eta:
#                 t_eta = eta / (eta - next_eta)
#             else:
#                 t_eta = 0
#
#             # Get boundary point, for the smaller t
#             t = min(t_eps, t_eta)
#             boundary_eps = min(max(eps + t * (next_eps - eps), 0), 1)
#             boundary_eta = min(max(eta + t * (next_eta - eta), 0), 1)
#             return [boundary_eps, boundary_eta]
#         else:
#             return [eps, eta]
#     return None
#
#
# def get_domain_position_xy(x, y, body):
#     surface = body.surface_u
#
#     # Get trivial coefficients
#     b = surface.B
#     d = surface.D
#     ny1 = surface.NY1
#     ny2 = surface.NY2
#     if b != 0:
#         # Get eta
#         eta = (2 * y - d) / b
#
#         # Check eta is not outside of domain
#         if min(1, max(0, eta)) == eta:
#             # Calculate missing coefficients
#             c = bernstein_function([eta], surface.C)
#             class_eta_y = c_function([eta], [ny1, ny2])[0] / c_function_max([ny1, ny2])
#             x_le = bernstein_function([eta], surface.X_LE)
#
#             # Get eps
#             if (c * class_eta_y) != 0:
#                 eps = (x - (x_le[0] + (c[0] / 2) * (1 - class_eta_y))) / (c[0] * class_eta_y)
#             else:
#                 # Return the average value for eps if division by zero
#                 eps = 0.5
#             if min(1, max(0, eps)) == eps:
#                 return [eps, eta]
#     else:
#         return None
#
#
# def moment(x, f, power):
#     mmnt = 0
#     for i in range(0, np.shape(x)[0]):
#         mmnt += math.pow(x[i], power) * f[i]
#     return mmnt
#
#
# def get_center_point(x, y, su, sl):
#     if np.shape(x) == np.shape(y) == np.shape(su) == np.shape(sl):
#         dim = np.shape(x)
#         center = np.zeros([dim[0], 3])
#         for i in range(0, dim[0]):
#             xi = x[i, :]
#             yi = y[i, :]
#             sui = su[i, :]
#             sli = sl[i, :]
#
#             ax = moment(xi, sui - sli, 0)
#             if ax != 0:
#                 qx = moment(xi, sui - sli, 1)
#             else:
#                 ax = moment(xi, np.ones(np.shape(xi)), 0)
#                 qx = moment(xi, np.ones(np.shape(xi)), 1)
#
#             ay = moment(yi, sui - sli, 0)
#             if ay != 0:
#                 qy = moment(yi, sui - sli, 1)
#             else:
#                 ay = moment(yi, np.ones(np.shape(yi)), 0)
#                 qy = moment(yi, np.ones(np.shape(yi)), 1)
#
#             z_moy_i = (sui + sli) / 2
#             az = moment(z_moy_i, sui - sli, 0)
#             if az != 0:
#                 qz = moment(z_moy_i, sui - sli, 1)
#             else:
#                 az = moment(z_moy_i, np.ones(np.shape(z_moy_i)), 0)
#                 qz = moment(z_moy_i, np.ones(np.shape(z_moy_i)), 1)
#
#             center[i, 0] = qx / ax
#             center[i, 1] = qy / ay
#             center[i, 2] = qz / az
#     else:
#         center = []
#     return center
#
#
def distribute_points(start, end, nb, mode):
    points_init = np.linspace(-1, 1, nb)
    points_final = []
    for point in points_init:
        if mode == "sphere":
            point_new = start + (end - start) * (0.5 * (1 + math.sin(math.pi * point / 2)))
        else:
            point_new = point / (1 + math.pow(point, 2))
            point_new = start + (end - start) * (0.5 + point_new)
        # Make sure new points are bounded
        point_new = max(start, min(end, point_new))
        points_final.append(point_new)
    return points_final
#
#
# def gen_base(body, output, mode, shape):
#
#     if mode == "distribute":
#         eps_base = distribute_points(0, 1, body.Nx, mode)
#         eta_base = distribute_points(0, 1, body.Ny, mode)
#     elif mode == "partial":
#         eps_base = distribute_points(0, 1, body.Nx, mode)
#         eta_base = np.linspace(0, 1, body.Ny)
#     elif mode == "crushed":
#         delta = 0.005
#         eps_base = distribute_points(0, 1, body.Nx, mode)
#         eta_base = np.linspace(0+delta, 1-delta, body.Ny)
#     elif mode == "sphere":
#         eps_base = distribute_points(0, 1, body.Nx, mode)
#         eta_base = distribute_points(0, 1, body.Ny, mode)
#     else:
#         eps_base = np.linspace(0, 1, body.Nx)
#         eta_base = np.linspace(0, 1, body.Ny)
#     eps = []
#     eta = []
#
#     for i in range(0, len(eta_base)):
#         eps_i = eps_base
#         eta_i = eta_base[i] * np.ones(np.shape(eps_base))
#
#         if i == 0:
#             eps = [eps_i]
#             eta = [eta_i]
#         else:
#             eps = np.append(eps, [eps_i], axis=0)
#             eta = np.append(eta, [eta_i], axis=0)
#
#     if shape == "rounded":
#         for i in range(0, np.shape(eps)[0]):
#             for j in range(0, np.shape(eps)[1]):
#                 y_ij = eta[i, j] - 0.5
#                 x_ij = eps[i, j] - 0.5
#
#                 # Correct positions for equal distance on circumference
#                 if abs(x_ij) < abs(y_ij):
#                     if y_ij != 0:
#                         correction = abs(math.tan(y_ij))
#                         if y_ij - x_ij * correction != 0:
#                             delta_x = correction * (y_ij ** 2 - x_ij ** 2) / (y_ij - x_ij * correction)
#                             delta_x = (x_ij / y_ij) * delta_x
#                             x_ij = x_ij - delta_x
#                 else:
#                     if x_ij != 0:
#                         correction = abs(math.tan(x_ij))
#                         if x_ij - y_ij * correction != 0:
#                             delta_y = correction * (x_ij**2 - y_ij**2) / (x_ij - y_ij * correction)
#                             delta_y = (y_ij / x_ij) * delta_y
#                             y_ij = y_ij - delta_y
#
#                 if x_ij ** 2 + y_ij ** 2 == 0:
#                     x_ij_new = 0
#                     y_ij_new = 0
#                 else:
#                     y_ij_new = y_ij * max(abs(x_ij), abs(y_ij)) / math.sqrt(x_ij ** 2 + y_ij ** 2)
#                     if abs(y_ij_new) == 0.5:
#                         x_ij_new = 0
#                     else:
#                         x_ij_new = x_ij * max(abs(x_ij), abs(y_ij)) / math.sqrt(
#                             (x_ij ** 2 + y_ij ** 2) * (1 - 4 * (y_ij_new ** 2)))
#
#                 eta_ij = y_ij_new + 0.5
#                 eps_ij = x_ij_new + 0.5
#
#                 eta[i, j] = max(0, min(1, eta_ij))
#                 eps[i, j] = max(0, min(1, eps_ij))
#
#     if output == "eps":
#         return eps
#     elif output == "eta":
#         return eta
#     else:
#         return eps, eta
#
#
# # def transform_eps_eta(eps, eta, N):
# #     eps = eps - 0.5
# #     eta = eta - 0.5
# #
# #     R = max(abs(eps), abs(eta))
# #
# #     if eps == 0 or N == 0:
# #         eps_new = eps
# #     else:
# #         eps_new = R / math.pow(1 + math.pow(abs(eta / eps), 1/N), N)
# #         if eps < 0:
# #             eps_new = -eps_new
# #
# #     if eta == 0 or N == 0:
# #         eta_new = eta
# #     else:
# #         eta_new = R / math.pow(1 + math.pow(abs(eps / eta), 1/N), N)
# #         if eta < 0:
# #             eta_new = -eta_new
# #
# #     eps_new = eps_new + 0.5
# #     eta_new = eta_new + 0.5
# #
# #     return eps_new, eta_new
#
#
# def is_point_inside_body(point, body):
#     tolerance = 1e-10
#     x = point[0]
#     y = point[1]
#     z = point[2]
#     coord = get_domain_position_xy(x, y, body)
#     if coord is not None:
#         eps = coord[0]
#         eta = coord[1]
#         # x_body_list, y_body_list, su_body_list, sl_body_list = build_body_surfaces([eps], [eta], param)
#         x_body, y_body, su_body, sl_body = build_body_surfaces_complete([eps], [eta], body)
#
#         if (x_body - x) ** 2 + (y_body - y) ** 2 < tolerance:
#             # Good to go
#             if su_body > z > sl_body:
#                 return True
#             else:
#                 return False
#         else:
#             print("Returned point is incorrect")
#             return False
#     else:
#         return False
#
#
# # def map_domain_to_bodies(param_main, param_list, intersections_list):
# #     # Set eps_domain and eta_domain to
# #
# #     eps_list = []
# #     eta_list = []
# #     Xu_list = []
# #     Xl_list = []
# #     Yu_list = []
# #     Yl_list = []
# #     Su_list = []
# #     Sl_list = []
# #     map_u_list = []
# #     map_l_list = []
# #     for i in range(0, len(param_list)):
# #         eps_i, eta_i = gen_base(param_main, "both", "cartesian", "square")
# #         X_i, Y_i, Su_i, Sl_i = build_body_surfaces(eps_i, eta_i, param_list[i])
# #         inside_body_i = np.ones(np.shape(X_i)) * i
# #         eps_list.append(eps_i)
# #         eta_list.append(eta_i)
# #         Xu_list.append(X_i)
# #         Xl_list.append(X_i)
# #         Yu_list.append(Y_i)
# #         Yl_list.append(Y_i)
# #         Su_list.append(Su_i)
# #         Sl_list.append(Sl_i)
# #         map_u_list.append(inside_body_i)
# #         map_l_list.append(inside_body_i)
# #
# #     for n in range(0, len(param_list)):
# #         # Run through Su vertices first, then run through Sl vertices
# #         for i in range(0, np.shape(Xu_list[n])[0]):
# #             for j in range(0, np.shape(Xu_list[n])[1]):
# #                 point_u = [Xu_list[n][i, j], Yu_list[n][i, j], Su_list[n][i, j]]
# #                 point_l = [Xl_list[n][i, j], Yl_list[n][i, j], Sl_list[n][i, j]]
# #
# #                 for k in range(0, len(param_list)):
# #                     # Only first param in param_list is checked here
# #                     if is_point_inside_body(point_u, param_list[k]):
# #                         map_u_list[n][i, j] = k
# #                     if is_point_inside_body(point_l, param_list[k]):
# #                         map_l_list[n][i, j] = k
#
#
# def build_body_surfaces(eps, eta, param):
#     # eps and eta are 2D numpy arrays
#     dim_eps = np.shape(eps)
#     dim_eta = np.shape(eta)
#
#     if dim_eps[0] == dim_eta[0]:
#         x = []
#         y = []
#         su = []
#         sl = []
#
#         # Set constant parameters
#         b = param.B
#         d = param.D
#         nc1 = param.NC1
#         nc2 = param.NC2
#         ny1 = param.NY1
#         ny2 = param.NY2
#         order = param.ORDER
#
#         a_upper = param.COEFFICIENTS_Upper
#         a_lower = param.COEFFICIENTS_Lower
#
#         if len(dim_eps) == 1:
#             num = 1
#         else:
#             num = dim_eps[0]
#
#         for i in range(0, num):
#             if len(dim_eps) == 1:
#                 eps_i = eps
#                 eta_i = eta
#             else:
#                 eps_i = eps[i, :]
#                 eta_i = eta[i, :]
#
#             c = bernstein_function(eta_i, param.C_)
#             scale_upper = bernstein_function(eta_i, param.SCALE_Upper_)
#             scale_lower = bernstein_function(eta_i, param.SCALE_Lower_)
#
#             z_te = bernstein_function(eta_i, param.Z_TE_)
#             r_le_upper = bernstein_function(eta_i, param.R_LE_Upper_)
#             r_le_lower = bernstein_function(eta_i, param.R_LE_Lower_)
#             beta_upper = bernstein_function(eta_i, param.BETA_Upper_)
#             beta_lower = bernstein_function(eta_i, param.BETA_Lower_)
#
#             x_le = bernstein_function(eta_i, param.X_LE_)
#             z_n = bernstein_function(eta_i, param.Z_N_)
#             delta_alpha_t = bernstein_function(eta_i, param.DELTA_ALPHA_T_)
#
#             # Class function coefficients
#             n1 = bernstein_function(eta_i, param.N1_)
#             n2 = bernstein_function(eta_i, param.N2_)
#
#             # Adjust BETA with DELTA_ALPHA
#             beta_upper = beta_upper + delta_alpha_t
#             beta_lower = beta_lower - delta_alpha_t
#
#             # Adjust C with DELTA_ALPHA
#             c = np.multiply(np.cos(delta_alpha_t), c)
#
#             class_eta = c_function(eta_i, [nc1, nc2])
#             class_eta_y = vector_divide(c_function_max([ny1, ny2]), c_function(eta_i, [ny1, ny2]))
#
#             if len(eps_i) == len(eta_i):
#                 x_i = []
#                 y_i = []
#                 su_i = []
#                 sl_i = []
#
#                 for j in range(0, len(eta_i)):
#                     eps_ij = eps_i[j]
#                     eta_ij = eta_i[j]
#
#                     zeta_t = z_te[j] / c[j]
#                     zeta_n = z_n[j] / c[j]
#
#                     if len(param.COEFFICIENTS_Upper) == 0:
#                         a_upper = bernstein_coefficients(z_te[j] / 2, r_le_upper[j], beta_upper[j],
#                                                          c[j], order, scale_upper[j])
#                     if len(param.COEFFICIENTS_Lower) == 0:
#                         a_lower = bernstein_coefficients(z_te[j] / 2, r_le_lower[j], beta_lower[j],
#                                                          c[j], order, scale_lower[j])
#
#                         # x_ij = (c[j] * class_eta_y[j] * eps_ij) + (x_le[j] + (c[j] / 2) * (1 - class_eta_y[j]))
#                         # ru_ij = c[j] * airfoil_profile([eps_ij], zeta_n, zeta_t / 2, delta_alpha_t[j],
#                         #                                a_upper, [n1[j], n2[j]], class_eta[j])
#                         # rl_ij = c[j] * airfoil_profile([eps_ij], zeta_n, zeta_t / 2, delta_alpha_t[j],
#                         #                                a_lower, [n1[j], n2[j]], class_eta[j])
#                         #
#                         # yu_ij = np.cos(eta_ij) * ru_ij + (d / 2) + (b / 4)
#                         # su_ij = np.sin(eta_ij) * ru_ij
#                         # yl_ij = np.cos(eta_ij) * rl_ij + (d / 2) + (b / 4)
#                         # sl_ij = np.sin(eta_ij) * rl_ij
#
#                     x_ij = (c[j] * class_eta_y[j] * eps_ij) + (x_le[j] + (c[j] / 2) * (1 - class_eta_y[j]))
#                     y_ij = (d + b * eta_ij) / 2
#                     su_ij = c[j] * airfoil_profile([eps_ij], zeta_n, zeta_t / 2, delta_alpha_t[j],
#                                                    a_upper, [n1[j], n2[j]], class_eta[j])
#                     sl_ij = -c[j] * airfoil_profile([eps_ij], -zeta_n, zeta_t / 2, -delta_alpha_t[j],
#                                                     a_lower, [n1[j], n2[j]], class_eta[j])
#
#                     if j == 0:
#                         x_i = [x_ij]
#                         y_i = [y_ij]
#                         su_i = [su_ij]
#                         sl_i = [sl_ij]
#                     else:
#                         x_i = np.append(x_i, [x_ij], axis=0)
#                         y_i = np.append(y_i, [y_ij], axis=0)
#                         su_i = np.append(su_i, [su_ij], axis=0)
#                         sl_i = np.append(sl_i, [sl_ij], axis=0)
#
#                 if i == 0:
#                     x = [x_i]
#                     y = [y_i]
#                     su = [su_i]
#                     sl = [sl_i]
#                 else:
#                     x = np.append(x, [x_i], axis=0)
#                     y = np.append(y, [y_i], axis=0)
#                     su = np.append(su, [su_i], axis=0)
#                     sl = np.append(sl, [sl_i], axis=0)
#         return x, y, su, sl
#     else:
#         return None
#
#
# def build_body_surfaces_complete(eps_p, eta_p, body):
#     def build_surface(eps, eta, surface, system):
#         # eps and eta are 2D numpy arrays
#         dim_eps = np.shape(eps)
#         dim_eta = np.shape(eta)
#
#         if dim_eps[0] == dim_eta[0]:
#             x = []
#             y = []
#             s = []
#
#             # Set constant parameters
#             b = surface.B
#             d = surface.D
#             nc1 = surface.NC1
#             nc2 = surface.NC2
#             ny1 = surface.NY1
#             ny2 = surface.NY2
#             order = surface.ORDER
#
#             a = surface.COEFFICIENTS
#
#             if len(dim_eps) == 1:
#                 num = 1
#             else:
#                 num = dim_eps[0]
#
#             for i in range(0, num):
#                 if len(dim_eps) == 1:
#                     eps_i = eps
#                     eta_i = eta
#                 else:
#                     eps_i = eps[i, :]
#                     eta_i = eta[i, :]
#
#                 c = bernstein_function(eta_i, surface.C)
#                 scale = bernstein_function(eta_i, surface.SCALE)
#
#                 z_te = bernstein_function(eta_i, surface.Z_TE)
#                 r_le = bernstein_function(eta_i, surface.R_LE)
#                 beta = bernstein_function(eta_i, surface.BETA)
#
#                 x_le = bernstein_function(eta_i, surface.X_LE)
#                 z_n = bernstein_function(eta_i, surface.Z_N)
#                 delta_alpha_t = bernstein_function(eta_i, surface.DELTA_ALPHA_T)
#
#                 if surface.ID == "Lower":
#                     z_n = -z_n
#                     delta_alpha_t = -delta_alpha_t
#
#                 # Class function coefficients
#                 n1 = bernstein_function(eta_i, surface.N1)
#                 n2 = bernstein_function(eta_i, surface.N2)
#
#                 # Adjust BETA with DELTA_ALPHA
#                 beta = beta + delta_alpha_t
#
#                 # Adjust C with DELTA_ALPHA
#                 c = np.multiply(np.cos(delta_alpha_t), c)
#
#                 class_eta = c_function(eta_i, [nc1, nc2])
#                 class_eta_y = vector_divide(c_function_max([ny1, ny2]), c_function(eta_i, [ny1, ny2]))
#
#                 if len(eps_i) == len(eta_i):
#                     x_i = []
#                     y_i = []
#                     s_i = []
#
#                     for j in range(0, len(eta_i)):
#                         eps_ij = eps_i[j]
#                         eta_ij = eta_i[j]
#
#                         zeta_t = z_te[j] / c[j]
#                         zeta_n = z_n[j] / c[j]
#
#                         if len(surface.COEFFICIENTS) == 0:
#                             a = bernstein_coefficients(z_te[j], r_le[j], beta[j], c[j], order, scale[j])
#
#                         r_ij = c[j] * airfoil_profile([eps_ij], zeta_n, zeta_t, delta_alpha_t[j],
#                                                       a, [n1[j], n2[j]], class_eta[j])
#                         if surface.ID == "Lower":
#                             r_ij = -r_ij
#
#                         if system == "Revolution":
#                             x_ij = (c[j] * eps_ij) + x_le[j]
#                             y_ij = (d + b * math.cos(2 * math.pi * eta_ij) * r_ij) / 2
#                             s_ij = (d + b * math.sin(2 * math.pi * eta_ij) * r_ij) / 2
#
#                         else:
#                             x_ij = (c[j] * class_eta_y[j] * eps_ij) + (x_le[j] + (c[j] / 2) * (1 - class_eta_y[j]))
#                             y_ij = (d + b * eta_ij) / 2
#                             s_ij = r_ij
#
#                         if j == 0:
#                             x_i = [x_ij]
#                             y_i = [y_ij]
#                             s_i = [s_ij]
#                         else:
#                             x_i = np.append(x_i, [x_ij], axis=0)
#                             y_i = np.append(y_i, [y_ij], axis=0)
#                             s_i = np.append(s_i, [s_ij], axis=0)
#
#                     if i == 0:
#                         x = [x_i]
#                         y = [y_i]
#                         s = [s_i]
#                     else:
#                         x = np.append(x, [x_i], axis=0)
#                         y = np.append(y, [y_i], axis=0)
#                         s = np.append(s, [s_i], axis=0)
#             return x, y, s
#         else:
#             return None
#
#     if len(eps_p) == 0 and len(eta_p) == 0:
#         eps_u, eta_u = gen_base(body.surface_u, "both", body.x_distribution, body.y_distribution)
#         eps_l, eta_l = gen_base(body.surface_l, "both", body.x_distribution, body.y_distribution)
#
#         Xu, Yu, Zu = build_surface(eps_u, eta_u, body.surface_u, body.type)
#         Xl, Yl, Zl = build_surface(eps_l, eta_l, body.surface_l, body.type)
#         return [Xu, Xl, Yu, Yl, Zu, Zl]
#
#     else:
#         Xu, Yu, Zu = build_surface(eps_p, eta_p, body.surface_u, body.type)
#         Xl, Yl, Zl = build_surface(eps_p, eta_p, body.surface_l, body.type)
#
#         if len(eps_p) == 1 and len(eta_p) == 1:
#             return Xu[0][0], Yu[0][0], Zu[0][0], Zl[0][0]
#         else:
#             return Xu, Yu, Zu, Zl

