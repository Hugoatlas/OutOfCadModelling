import math
import numpy as np
import export_to_gmsh
import vlm_mesh as vlm
import matplotlib.pyplot as plt


class Param:
    def __init__(self, nx, ny, b, d, c, scale_u, scale_l,
                 z_te, r_le_u, r_le_l, beta_u, beta_l, coeffs_u, coeffs_l,
                 x_le, z_n, delta_alpha_t,
                 n1, n2, nc1, nc2, ny1, ny2, order):
        self.Nx = nx
        self.Ny = ny

        self.B = b  # Largeur totale des deux ailes combinées
        self.D = d  # Diamètre du fuselage
        self.C_ = list(c)  # Largeur locale de l'aile
        self.SCALE_Upper_ = list(scale_u)  # Scale général du profile supérieur
        self.SCALE_Lower_ = list(scale_l)  # Scale général du profile inférieur

        self.Z_TE_ = list(z_te)  # Espacement au bord de fuite
        self.R_LE_Upper_ = list(r_le_u)  # Rayon de courbure de la surface supérieure
        self.R_LE_Lower_ = list(r_le_l)  # Rayon de courbure de la surface inférieure
        self.BETA_Upper_ = list(beta_u)  # Angle de la surface supérieure au bord de fuite
        self.BETA_Lower_ = list(beta_l)  # Angle de la surface inférieure au bord de fuite
        self.COEFFICIENTS_Upper = list(coeffs_u)    # Coefficients du polynôme de Bernstein de la surface supérieure
        self.COEFFICIENTS_Lower = list(coeffs_l)  # Coefficients du polynôme de Bernstein de la surface supérieure

        self.X_LE_ = list(x_le)  # Translation en x de l'aile
        self.Z_N_ = list(z_n)  # Translation en z de l'aile
        self.DELTA_ALPHA_T_ = list(delta_alpha_t)  # Inclinaison de l'aile

        # Class function coefficients
        self.N1_ = list(n1)  # Coefficient de classe du bord d'attaque
        self.N2_ = list(n2)  # Coefficient de classe du bord de fuite
        self.NC1 = nc1  # Coefficient de classe de l'extrémité -Y de l'aile
        self.NC2 = nc2  # Coefficient de classe de l'extrémité +Y de l'aile
        self.NY1 = ny1
        self.NY2 = ny2
        self.ORDER = order  # Ordre des polynômes de Bernstein


def mirror_body_xz(param):
    D = - (param.B + param.D)
    C_ = vlm.reverse_array(param.C_)
    SCALE_Upper_ = vlm.reverse_array(param.SCALE_Upper_)
    SCALE_Lower_ = vlm.reverse_array(param.SCALE_Lower_)
    Z_TE_ = vlm.reverse_array(param.Z_TE_)
    R_LE_Upper_ = vlm.reverse_array(param.R_LE_Upper_)
    R_LE_Lower_ = vlm.reverse_array(param.R_LE_Lower_)
    BETA_Upper_ = vlm.reverse_array(param.BETA_Upper_)
    BETA_Lower_ = vlm.reverse_array(param.BETA_Lower_)
    COEFFICIENTS_Upper = param.COEFFICIENTS_Upper
    COEFFICIENTS_Lower = param.COEFFICIENTS_Lower
    X_LE_ = vlm.reverse_array(param.X_LE_)
    Z_N_ = vlm.reverse_array(param.Z_N_)
    DELTA_ALPHA_T_ = vlm.reverse_array(param.DELTA_ALPHA_T_)
    N1_ = vlm.reverse_array(param.N1_)
    N2_ = vlm.reverse_array(param.N2_)
    NC1 = param.NC2
    NC2 = param.NC1

    param_mirror = Param(nx=param.Nx, ny=param.Ny, b=param.B, d=D, c=C_, scale_u=SCALE_Upper_, scale_l=SCALE_Lower_,
                         z_te=Z_TE_, r_le_u=R_LE_Upper_, r_le_l=R_LE_Lower_, beta_u=BETA_Upper_, beta_l=BETA_Lower_,
                         coeffs_u=COEFFICIENTS_Upper, coeffs_l=COEFFICIENTS_Lower, x_le=X_LE_, z_n=Z_N_,
                         delta_alpha_t=DELTA_ALPHA_T_, n1=N1_, n2=N2_, nc1=NC1, nc2=NC2,
                         ny1=param.NY1, ny2=param.NY2, order=param.ORDER)
    return param_mirror


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


def c_function_max(ns):
    eps_max = 0.5
    y_max = math.pow(eps_max, ns[0]) * math.pow(1 - eps_max, ns[1])
    return y_max


def bernstein_function(eps, coefficients):
    res = np.zeros(np.shape(eps))
    order = len(coefficients)
    for i in range(0, order):
        bern_i = math.factorial(order - 1) / (math.factorial(i) * math.factorial(order - 1 - i))
        for j in range(0, len(eps)):
            res[j] += bern_i * coefficients[i] * math.pow(eps[j], i) * math.pow(1 - eps[j], order - 1 - i)
    return res


def is_point_out_of_bounds(value, minimum, maximum):
    if value >= maximum or value <= minimum:
        return True
    else:
        return False


def bernstein_coefficients(z_te, r_le, beta, c, order, scale):
    if order == 0:
        a = [scale]
    else:
        a = np.ones([order, 1]) * scale
        a[0] = (math.sqrt(2 * r_le / c))
        a[order - 1] = math.tan(beta) + (z_te / c)
    return a


def derivative(x, func):
    res = np.zeros(np.shape(x))
    for i in range(0, len(x)):
        if i == 0:
            if x[i + 1] == x[i]:
                res[i] = None
            else:
                res[i] = (func[i + 1] - func[i]) / (x[i + 1] - x[i])
        elif i == len(x) - 1:
            if x[i] == x[i - 1]:
                res[i] = None
            else:
                res[i] = (func[i] - func[i - 1]) / (x[i] - x[i - 1])
        else:
            if x[i + 1] == x[i - 1]:
                res[i] = None
            else:
                res[i] = (func[i + 1] - func[i - 1]) / (x[i + 1] - x[i - 1])
    return res


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

    if len(eps) == 1:
        res = y[0]
    else:
        res = y
    return res


def get_domain_derivative(eps, eta, param):
    dim_eps = np.shape(eps)
    dim_eta = np.shape(eta)

    Xd_eps = np.zeros(np.shape(eps))
    Xd_eta = np.zeros(np.shape(eps))
    Yd_eps = np.zeros(np.shape(eps))
    Yd_eta = np.zeros(np.shape(eps))

    if dim_eps == dim_eta:
        X, Y, Su, Sl = build_body_surfaces(eps, eta, param)

        for i in range(0, dim_eps[0]):
            Xd_eps[i, :] = derivative(eps[i, :], X[i, :])
            Yd_eps[i, :] = derivative(eps[i, :], Y[i, :])

        for j in range(0, dim_eps[1]):
            Xd_eta[:, j] = derivative(eta[:, j], X[:, j])
            Yd_eta[:, j] = derivative(eta[:, j], Y[:, j])

    return Xd_eps, Xd_eta, Yd_eps, Yd_eta


def get_gradient_norm_square(grad_x, grad_y):
    dim_x = np.shape(grad_x)
    dim_y = np.shape(grad_y)
    norm = np.zeros(np.shape(grad_x))
    if dim_x == dim_y:
        for i in range(0, dim_x[0]):
            for j in range(0, dim_x[1]):
                norm[i, j] = math.sqrt(grad_x[i, j]**2 + grad_y[i, j]**2)
    return norm


def get_full_domain_boundary(nx, ny, param):
    list_of_points = []
    list_of_coords = []
    eta = 0
    eps = 0
    for i in range(1, nx):
        eps = i / (nx - 1)
        x_, y_, _, _ = build_body_surfaces([eps], [eta], param)
        list_of_points.append([x_[0][0], y_[0][0]])
        list_of_coords.append([eps, eta])

    for i in range(1, ny):
        eta = i / (ny - 1)
        x_, y_, _, _ = build_body_surfaces([eps], [eta], param)
        list_of_points.append([x_[0][0], y_[0][0]])
        list_of_coords.append([eps, eta])

    for i in range(1, nx):
        eps = 1 - (i / (nx - 1))
        x_, y_, _, _ = build_body_surfaces([eps], [eta], param)
        list_of_points.append([x_[0][0], y_[0][0]])
        list_of_coords.append([eps, eta])

    for i in range(1, ny):
        eta = 1 - (i / (ny - 1))
        x_, y_, _, _ = build_body_surfaces([eps], [eta], param)
        list_of_points.append([x_[0][0], y_[0][0]])
        list_of_coords.append([eps, eta])

    list_of_points.append(list_of_points[0])
    list_of_coords.append(list_of_coords[0])
    return list_of_points, list_of_coords


def get_domain_boundary(x, y, dir_x, dir_y, param, max_attempts):
    res = get_domain_position_xy(x, y, param)
    eps = res[0]
    eta = res[1]
    next_res = None
    for i in range(0, max_attempts):
        new_x = x + (dir_x / (i + 1))
        new_y = y + (dir_y / (i + 1))
        next_res = get_domain_position_xy(new_x, new_y, param)
        if next_res is not None:
            break

    if next_res is not None:
        next_eps = next_res[0]
        next_eta = next_res[1]

        if eps != 0 and eps != 1 and eta != 0 and eta != 1:
            # Get t for eps
            if next_eps > eps:
                t_eps = (1 - eps) / (next_eps - eps)
            elif next_eps < eps:
                t_eps = eps / (eps - next_eps)
            else:
                t_eps = 0
            # Get t for eta
            if next_eta > eta:
                t_eta = (1 - eta) / (next_eta - eta)
            elif next_eta < eta:
                t_eta = eta / (eta - next_eta)
            else:
                t_eta = 0

            # Get boundary point, for the smaller t
            t = min(t_eps, t_eta)
            boundary_eps = min(max(eps + t * (next_eps - eps), 0), 1)
            boundary_eta = min(max(eta + t * (next_eta - eta), 0), 1)
            return [boundary_eps, boundary_eta]
        else:
            return [eps, eta]
    return None


def get_domain_position_xy(x, y, param):
    # Get trivial coefficients
    b = param.B
    d = param.D
    ny1 = param.NY1
    ny2 = param.NY2
    if b != 0:
        # Get eta
        eta = (2 * y - d) / b

        # Check eta is not outside of domain
        if min(1, max(0, eta)) == eta:
            # Calculate missing coefficients
            c = bernstein_function([eta], param.C_)
            class_eta_y = c_function([eta], [ny1, ny2])[0] / c_function_max([ny1, ny2])
            x_le = bernstein_function([eta], param.X_LE_)

            # Get eps
            if (c * class_eta_y) != 0:
                eps = (x - (x_le[0] + (c[0] / 2) * (1 - class_eta_y))) / (c[0] * class_eta_y)
            else:
                # Return the average value for eps if division by zero
                eps = 0.5
            if min(1, max(0, eps)) == eps:
                return [eps, eta]
    else:
        return None


def gen_base(param, output):

    mode = "cartesian"
    # mode = "rounded"

    if mode == "distribute":
        eps_base = distribute_points(0, 1, param.Nx)
        eta_base = distribute_points(0, 1, param.Ny)
    else:
        eps_base = np.linspace(0, 1, param.Nx)
        eta_base = np.linspace(0, 1, param.Ny)
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

    if mode == "rounded":
        eps = 1/2 + 4 * ((eps - 1/2)**3)
        eta = 1 / 4 + (eta / 2)

        print(eps)

    if output == "eps":
        return eps
    elif output == "eta":
        return eta
    else:
        return eps, eta


def build_body_surfaces(eps, eta, param):
    # eps and eta are 2D numpy arrays
    dim_eps = np.shape(eps)
    dim_eta = np.shape(eta)

    if dim_eps[0] == dim_eta[0]:
        x = []
        y = []
        su = []
        sl = []

        # Set constant parameters
        b = param.B
        d = param.D
        nc1 = param.NC1
        nc2 = param.NC2
        ny1 = param.NY1
        ny2 = param.NY2
        order = param.ORDER

        a_upper = param.COEFFICIENTS_Upper
        a_lower = param.COEFFICIENTS_Lower

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

            c = bernstein_function(eta_i, param.C_)
            scale_upper = bernstein_function(eta_i, param.SCALE_Upper_)
            scale_lower = bernstein_function(eta_i, param.SCALE_Lower_)

            z_te = bernstein_function(eta_i, param.Z_TE_)
            r_le_upper = bernstein_function(eta_i, param.R_LE_Upper_)
            r_le_lower = bernstein_function(eta_i, param.R_LE_Lower_)
            beta_upper = bernstein_function(eta_i, param.BETA_Upper_)
            beta_lower = bernstein_function(eta_i, param.BETA_Lower_)

            x_le = bernstein_function(eta_i, param.X_LE_)
            z_n = bernstein_function(eta_i, param.Z_N_)
            delta_alpha_t = bernstein_function(eta_i, param.DELTA_ALPHA_T_)

            # Class function coefficients
            n1 = bernstein_function(eta_i, param.N1_)
            n2 = bernstein_function(eta_i, param.N2_)

            # Adjust BETA with DELTA_ALPHA
            beta_upper = beta_upper + delta_alpha_t
            beta_lower = beta_lower - delta_alpha_t

            # Adjust C with DELTA_ALPHA
            c = np.multiply(np.cos(delta_alpha_t), c)

            class_eta = c_function(eta_i, [nc1, nc2])
            class_eta_y = vector_divide(c_function_max([ny1, ny2]), c_function(eta_i, [ny1, ny2]))

            if len(eps_i) == len(eta_i):
                x_i = []
                y_i = []
                su_i = []
                sl_i = []

                for j in range(0, len(eta_i)):
                    eps_ij = eps_i[j]
                    eta_ij = eta_i[j]
                    x_ij = (c[j] * class_eta_y[j] * eps_ij) + (x_le[j] + (c[j] / 2) * (1 - class_eta_y[j]))
                    y_ij = (d + b * eta_ij) / 2

                    zeta_t = z_te[j] / c[j]
                    zeta_n = z_n[j] / c[j]

                    if len(param.COEFFICIENTS_Upper) == 0:
                        a_upper = bernstein_coefficients(z_te[j] / 2, r_le_upper[j], beta_upper[j],
                                                         c[j], order, scale_upper[j])
                    if len(param.COEFFICIENTS_Lower) == 0:
                        a_lower = bernstein_coefficients(z_te[j] / 2, r_le_lower[j], beta_lower[j],
                                                         c[j], order, scale_lower[j])

                    su_ij = c[j] * airfoil_profile([eps_ij], zeta_n, zeta_t / 2, delta_alpha_t[j],
                                                   a_upper, [n1[j], n2[j]], class_eta[j])
                    sl_ij = -c[j] * airfoil_profile([eps_ij], -zeta_n, zeta_t / 2, -delta_alpha_t[j],
                                                    a_lower, [n1[j], n2[j]], class_eta[j])

                    if j == 0:
                        x_i = [x_ij]
                        y_i = [y_ij]
                        su_i = [su_ij]
                        sl_i = [sl_ij]
                    else:
                        x_i = np.append(x_i, [x_ij], axis=0)
                        y_i = np.append(y_i, [y_ij], axis=0)
                        su_i = np.append(su_i, [su_ij], axis=0)
                        sl_i = np.append(sl_i, [sl_ij], axis=0)

                if i == 0:
                    x = [x_i]
                    y = [y_i]
                    su = [su_i]
                    sl = [sl_i]
                else:
                    x = np.append(x, [x_i], axis=0)
                    y = np.append(y, [y_i], axis=0)
                    su = np.append(su, [su_i], axis=0)
                    sl = np.append(sl, [sl_i], axis=0)
        return x, y, su, sl
    else:
        return None

