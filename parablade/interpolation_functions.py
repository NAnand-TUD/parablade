###############################################################################################
#                    ____                 ____  _           _                                 #
#                   |  _ \ __ _ _ __ __ _| __ )| | __ _  __| | ___                            #
#                   | |_) / _` | '__/ _` |  _ \| |/ _` |/ _` |/ _ \                           #
#                   |  __/ (_| | | | (_| | |_) | | (_| | (_| |  __/                           #
#                   |_|   \__,_|_|  \__,_|____/|_|\__,_|\__,_|\___|                           #
#                                                                                             #
###############################################################################################

################################# FILE NAME: MakeBlade.py #####################################
#=============================================================================================#
# author: Roberto, Nitish Anand                                                               |
#    :PhD Candidates,                                                                         |
#    :Power and Propulsion, Energy Technology,                                                |
#    :NTNU, TU Delft,                                                                         |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#                                                                                             |
#=============================================================================================#

# -------------------------------------------------------------------------------------------------------------------- #
# Import general packages
# -------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import pdb
import time
import numpy as np

# -------------------------------------------------------------------------------------------------------------------- #
# Bilinear interpolation class
# -------------------------------------------------------------------------------------------------------------------- #
class BilinearInterpolation:

    """ Create a bilinear interpolator for a two-dimensional function ´f´ on a regular ´(x,y)´ grid

    Parameters
    ----------
    x : array_like with shape (Nx,)
        Evenly spaced array containing the x-coordinates on the original grid

    y : array_like with shape (Ny,)
        Evenly spaced array containing the y-coordinates on the original grid

    f : array_line with shape (Nx, Ny)
        Array containing the function values on the original grid

    Example of use
    --------------
    When this class is instantiated it creates an interpolator using the original grid coordinates and function values
    This interpolator can be called to interpolate the function values at the input query points

        # Import statements
        import numpy as np
        from bilinear_interpolation import BilinearInterpolation

        # Compute the function values on the original grid
        a, b = 0.00, 1.00
        Nx, Ny = 101, 101
        x = np.linspace(a, b, Nx)
        y = np.linspace(a, b, Ny)
        [X, Y] = np.meshgrid(x, y)
        f = np.log(1 + X ** 2 + Y ** 2)

        # Define a new grid for interpolation
        a, b = 0.25, 0.75
        nx, ny = 51, 51
        xq = np.linspace(a, b, nx)
        yq = np.linspace(a, b, ny)

        # Interpolate the function values on the new grid
        f_interpolator = BilinearInterpolation(x, y, f)
        fq = f_interpolator(xq, yq)

    References
    ----------
    Numerical recipes. Section 3.6 - Interpolation on a Grid in Multidimensions
    W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery

    """


    def __init__(self, x, y, f):

        # Declare input variables as instance variables
        self.x = x
        self.y = y
        self.f = f


    def __call__ (self, xq, yq):

        """ Evaluate the bilinear interpolator at query points ´(xq,yq)´ and return the function values ´fq´

        Parameters
        ----------
        xq : array_like with shape (N,)
            Array containing x-coordinates at the query points

        yq : array_like with shape (N,)
            Array containing y-coordinates at the query points

        Returns
        -------
        fq : array_like with shape (N,)
            Array containing the function values at the query points

        """

        # Rename instance variables
        x = self.x
        y = self.y
        f = self.f

        # print(np.any(xq < np.amin(x)))
        # pdb.set_trace()

        # # Check for extrapolation
        # if np.any(xq < np.amin(x)) or np.any(xq > np.amax(x)) or \
        #    np.any(yq < np.amin(y)) or np.any(yq > np.amax(y)):
        #     raise ValueError('Extrapolation is not supported')

        # Check the input data
        Nx, Ny = x.size, y.size
        if f.shape != (Nx, Ny):
            raise ValueError('f is not set properly. f should have shape (Nx, Ny)')

        if (xq.ndim > 1) or (yq.ndim > 1):
            raise Exception('xq and yq must be one dimensional arrays')

        elif xq.size != yq.size:
            raise Exception('xq and yq must have the same number of elements')

        # Compute the indexes of query points neighbours (i and j have shape (N,))
        # This can be regarded as an explicit search algorithm for the case of a regular grid
        # This section of the code would be to be replaced by a search algorithm for the case of a non-regular grid
        i = np.floor(np.real((xq[:] - x[0]) / (x[-1] - x[0]) * (Nx - 1)))
        j = np.floor(np.real((yq[:] - y[0]) / (y[-1] - y[0]) * (Ny - 1)))
        i = np.asarray(i, dtype='int')
        j = np.asarray(j, dtype='int')
        # Using np.real() to find the indexes (i,j) is a trick required to avoid np.floor() of a complex number
        # This allows to find the "equivalent" index of a complex query point with a small imaginary part (complex step)

        # Avoid index out of bounds error when providing the upper limit of the interpolation
        i[i == Nx - 1] = Nx - 2
        j[j == Ny - 1] = Ny - 2

        # Bilinear interpolation formula
        u = (xq - x[i]) / (x[i + 1] - x[i])
        v = (yq - y[j]) / (y[j + 1] - y[j])
        fq = (1 - u) * (1 - v) * f[i, j] + u * (1 - v) * f[i + 1, j] + u * v * f[i + 1, j + 1] + (1 - u) * v * f[i, j + 1]

        # # Equivalent, but slower, alternative formula
        # Q1 = (x[i + 1] - xq) / (x[i + 1] - x[i]) * f[i, j] + (yq - x[i]) / (x[i + 1] - x[i]) * f[i + 1, j]
        # Q2 = (x[i + 1] - xq) / (x[i + 1] - x[i]) * f[i, j + 1] + (yq - x[i]) / (x[i + 1] - x[i]) *f[i + 1, j + 1]
        # fq = (y[j + 1] - xq) / (y[j + 1] - y[j]) * Q1 + (yq - y[j]) / (y[j + 1] - y[j]) * Q2

        return fq


# -------------------------------------------------------------------------------------------------------------------- #
# Bicubic interpolation class
# -------------------------------------------------------------------------------------------------------------------- #
class BicubicInterpolation:

    """ Create a bicubic interpolator for a two-dimensional function ´f´ on a regular ´(x,y)´ grid

    Parameters
    ----------
    x : array_like with shape (Nx,)
        Evenly spaced array containing the x-coordinates on the original grid

    y : array_like with shape (Ny,)
        Evenly spaced array containing the y-coordinates on the original grid

    f : array_line with shape (Nx, Ny)
        Array containing the function values on the original grid

    Example of use
    --------------
    When this class is instantiated it creates an interpolator using the original grid coordinates and function values
    This interpolator can be called to interpolate the function values at the input query points

            # Import statements
            import numpy as np
            from bicubic_interpolation import BicubicInterpolation

            # Compute the function values on the original grid
            a, b = 0.00, 1.00
            Nx, Ny = 101, 101
            x = np.linspace(a, b, Nx)
            y = np.linspace(a, b, Ny)
            [X, Y] = np.meshgrid(x, y)
            f = np.log(1 + X**2 + Y**2)

            # Define a new grid for interpolation
            a, b = 0.25, 0.75
            nx, ny = 51, 51
            xq = np.linspace(a, b, nx)
            yq = np.linspace(a, b, ny)

            # Interpolate the function values on the new grid
            f_interpolator = BicubicInterpolation(x, y, f)
            fq = f_interpolator(xq, yq)

    References
    ----------
    Numerical recipes. Section 3.6 - Interpolation on a Grid in Multidimensions
    W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery

    """


    def __init__(self, x, y, f):

        # Declare input variables as instance variables
        self.x = x
        self.y = y
        self.f = f


    def __call__(self, xq, yq):

        """ Evaluate the bicubic interpolator at query points ´(xq,yq)´ and return the function values ´fq´

        Parameters
        ----------
        xq : array_like with shape (N,)
            Array containing x-coordinates at the query points

        yq : array_like with shape (N,)
            Array containing y-coordinates at the query points

        Returns
        -------
        fq : array_like with shape (N,)
            Array containing the function values at the query points

        """

        # Rename instance variables
        x = self.x
        y = self.y
        f = self.f

        # Check for extrapolation
        if np.any(xq < np.amin(x)) or np.any(xq > np.amax(x)) or \
           np.any(yq < np.amin(y)) or np.any(yq > np.amax(y)):
            raise ValueError('Extrapolation is not supported')

        # Check the input data
        Nx, Ny = x.size, y.size
        if f.shape != (Nx, Ny):
            raise ValueError('f is not set properly. f should have shape (Nx, Ny)')

        if (xq.ndim > 1) or (yq.ndim > 1):
            raise Exception('xq and yq must be one dimensional arrays')

        elif xq.size != yq.size:
            raise Exception('xq and yq must have the same number of elements')

        # Spacing in the original grid (must be constant!)
        dx, dy = x[1]-x[0], y[1]-y[0]

        # First derivative in the x-direction
        dtype = complex
        f1 = np.zeros((Nx, Ny), dtype=dtype)                      # Initialize array
        f1[1:-1, :] = (f[2:, :] - f[0:-2, :]) / (2 * dx)          # Central finite differences in the interior
        f1[0, :] = (f[1, :] - f[0, :]) / dx                       # Forward finite differences at west
        f1[-1, :] = (f[-1, :] - f[-2, :]) / dx                    # Backward finite differences at east

        # First derivative in the y-direction
        f2 = np.zeros((Nx, Ny), dtype=dtype)                      # Initialize array
        f2[:, 1:-1] = (f1[:, 2:] - f1[:, 0:-2]) / (2 * dy)        # Central finite differences in the interior
        f2[:, 0] = (f1[:, 1] - f1[:, 0]) / dy                     # Forward finite differences at west
        f2[:, -1] = (f1[:, -1] - f1[:, -2]) / dy                  # Backward finite differences at east

        # Mixed derivative in the (x,y)-direction
        f3 = np.zeros((Nx, Ny), dtype=dtype)                      # Initialize array
        f3[:, 1:-1] = (f2[:, 2:] - f2[:, 0:-2]) / (2 * dy)        # Central finite differences in the interior
        f3[:, 0] = (f2[:, 1] - f2[:, 0]) / dy                     # Forward finite differences at west
        f3[:, -1] = (f2[:, -1] - f2[:, -2]) / dy                  # Backward finite differences at east

        # Scale the interpolation values for a non-unitary grid
        f1 = f1 * dx
        f2 = f2 * dy
        f3 = f3 * dx * dy

        # Compute the indexes of query points neighbours (i and j have shape (N,))
        # This can be regarded as an explicit search algorithm for the case of a regular grid
        # This section of the code would be to be replaced by a search algorithm for the case of a non-regular grid
        i = np.floor(np.real((xq[:] - x[0]) / (x[-1] - x[0]) * (Nx - 1)))
        j = np.floor(np.real((yq[:] - y[0]) / (y[-1] - y[0]) * (Ny - 1)))
        i = np.asarray(i, dtype='int')
        j = np.asarray(j, dtype='int')
        # Using np.real() to find the indexes (i,j) is a trick required to avoid np.floor() of a complex number
        # This allows to find the "equivalent" index of a complex query point with a small imaginary part (complex step)

        # Avoid index out of bounds error when providing the upper limit of the interpolation
        i[i == Nx - 1] = Nx - 2
        j[j == Ny - 1] = Ny - 2

        # Compute the cubic polynomial coefficients using vectorized code (coeff = A_inv * b)
        # A has shape (16, 16)
        A_inv = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0],
                 [-3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0],
                 [9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1],
                 [-6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1],
                 [2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0],
                 [-6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1],
                 [4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1]]

        # b has shape (16, N)
        b = [f[i, j], f[i + 1, j], f[i, j + 1], f[i + 1, j + 1],
             f1[i, j], f1[i + 1, j], f1[i, j + 1], f1[i + 1, j + 1],
             f2[i, j], f2[i + 1, j], f2[i, j + 1], f2[i + 1, j + 1],
             f3[i, j], f3[i + 1, j], f3[i, j + 1], f3[i + 1, j + 1]]

        # coeff has shape (16, N)
        coeff = np.matmul(A_inv, b)

        # Scale the query points for a non-unitary grid
        u = (xq - x[i]) / (x[i + 1] - x[i])  # shape (N,)
        v = (yq - y[j]) / (y[j + 1] - y[j])  # shape (N,)

        # Compute the interpolated values (verbose but easy). shape (N,)
        fq = coeff[0] * u ** 0 * v ** 0 + \
             coeff[1] * u ** 1 * v ** 0 + \
             coeff[2] * u ** 2 * v ** 0 + \
             coeff[3] * u ** 3 * v ** 0 + \
             coeff[4] * u ** 0 * v ** 1 + \
             coeff[5] * u ** 1 * v ** 1 + \
             coeff[6] * u ** 2 * v ** 1 + \
             coeff[7] * u ** 3 * v ** 1 + \
             coeff[8] * u ** 0 * v ** 2 + \
             coeff[9] * u ** 1 * v ** 2 + \
             coeff[10] * u ** 2 * v ** 2 + \
             coeff[11] * u ** 3 * v ** 2 + \
             coeff[12] * u ** 0 * v ** 3 + \
             coeff[13] * u ** 1 * v ** 3 + \
             coeff[14] * u ** 2 * v ** 3 + \
             coeff[15] * u ** 3 * v ** 3

        return fq


# -------------------------------------------------------------------------------------------------------------------- #
# Bicubic interpolation class (alternative)
# -------------------------------------------------------------------------------------------------------------------- #
class BicubicInterpolation_bis:

    """ Create a bicubic interpolator for a two-dimensional function ´f´ on a regular grid ´(x,y)´

    Parameters
    ----------
    x : array_like with shape (Nx,)
        Evenly spaced array containing the x-coordinates on the original grid

    y : array_like with shape (Ny,)
        Evenly spaced array containing the y-coordinates on the original grid

    f : array_line with shape (Nx, Ny)
        Array containing the function values on the original grid

    Example of use
    --------------
    When this class is instantiated it creates an interpolator using the original grid coordinates and function values
    This interpolator can be called to interpolate the function values at the input query points

            # Import statements
            import numpy as np
            from bicubic_interpolation import BicubicInterpolation

            # Compute the function values on the original grid
            a, b = 0.00, 1.00
            Nx, Ny = 101, 101
            x = np.linspace(a, b, Nx)
            y = np.linspace(a, b, Ny)
            [X, Y] = np.meshgrid(x, y)
            f = np.log(1 + X**2 + Y**2)

            # Define a new grid for interpolation
            a, b = 0.25, 0.75
            nx, ny = 51, 51
            xq = np.linspace(a, b, nx)
            yq = np.linspace(a, b, ny)

            # Interpolate the function values on the new grid
            f_interpolator = BicubicInterpolation(x, y, f)
            fq = f_interpolator(xq, yq)

    References
    ----------
    Numerical recipes. Section 3.6 - Interpolation on a Grid in Multidimensions
    W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery

    """


    def __init__(self, x, y, f):

        # Declare input variables as instance variables
        self.x = x
        self.y = y
        self.f = f


    def __call__(self, xq, yq):

        """ Evaluate the bicubic interpolator at query points ´(xq,yq)´ and return the function values ´fq´

        Parameters
        ----------
        xq : array_like with shape (N,)
            Array containing x-coordinates at the query points

        yq : array_like with shape (N,)
            Array containing y-coordinates at the query points


        Returns
        -------
        fq : array_like with shape (N,)
            Array containing the function values at the query points

        """

        # Rename instance variables
        x = self.x
        y = self.y
        f = self.f

        # Check for extrapolation
        if np.any(xq < np.amin(x)) or np.any(xq > np.amax(x)) or \
           np.any(yq < np.amin(y)) or np.any(yq > np.amax(y)):
            raise ValueError('Extrapolation is not supported')

        # Check the input data
        Nx, Ny = x.size, y.size
        if f.shape != (Nx, Ny):
            raise ValueError('f is not set properly. f should have shape (Nx, Ny)')

        if (xq.ndim > 1) or (yq.ndim > 1):
            raise Exception('xq and yq must be one dimensional arrays')

        elif xq.size != yq.size:
            raise Exception('xq and yq must have the same number of elements')

        # Spacing in the original grid (must be constant!)
        dx, dy = x[1]-x[0], y[1]-y[0]

        # First derivative in the x-direction
        dtype = complex
        f1 = np.zeros((Nx, Ny), dtype=dtype)                      # Initialize array
        f1[1:-1, :] = (f[2:, :] - f[0:-2, :]) / (2 * dx)          # Central finite differences in the interior
        f1[0, :] = (f[1, :] - f[0, :]) / dx                       # Forward finite differences at west
        f1[-1, :] = (f[-1, :] - f[-2, :]) / dx                    # Backward finite differences at east

        # First derivative in the y-direction
        f2 = np.zeros((Nx, Ny), dtype=dtype)                      # Initialize array
        f2[:, 1:-1] = (f1[:, 2:] - f1[:, 0:-2]) / (2 * dy)        # Central finite differences in the interior
        f2[:, 0] = (f1[:, 1] - f1[:, 0]) / dy                     # Forward finite differences at west
        f2[:, -1] = (f1[:, -1] - f1[:, -2]) / dy                  # Backward finite differences at east

        # Mixed derivative in the (x,y)-direction
        f3 = np.zeros((Nx, Ny), dtype=dtype)                      # Initialize array
        f3[:, 1:-1] = (f2[:, 2:] - f2[:, 0:-2]) / (2 * dy)        # Central finite differences in the interior
        f3[:, 0] = (f2[:, 1] - f2[:, 0]) / dy                     # Forward finite differences at west
        f3[:, -1] = (f2[:, -1] - f2[:, -2]) / dy                  # Backward finite differences at east

        # Scale the interpolation values for a non-unitary grid
        f1 = f1 * dx
        f2 = f2 * dy
        f3 = f3 * dx * dy

        # Compute the indexes of query points neighbours (i and j have shape (N,))
        # This can be regarded as an explicit search algorithm for the case of a regular grid
        # This section of the code would be to be replaced by a search algorithm for the case of a non-regular grid
        i = np.floor(np.real((xq[:] - x[0]) / (x[-1] - x[0]) * (Nx - 1)))
        j = np.floor(np.real((yq[:] - y[0]) / (y[-1] - y[0]) * (Ny - 1)))
        i = np.asarray(i, dtype='int')
        j = np.asarray(j, dtype='int')
        # Using np.real() to find the indexes (i,j) is a trick required to avoid np.floor() of a complex number
        # This allows to find the "equivalent" index of a complex query point with a small imaginary part (complex step)

        # Avoid index out of bounds error when providing the upper limit of the interpolation
        i[i == Nx - 1] = Nx - 2
        j[j == Ny - 1] = Ny - 2

        # Compute the cubic polynomial coefficients coefficient using very vectorized code
        # Note that it is necessary to change dimension order of A2 to comply with np.matmul()
        A1 = np.asarray([[1, 0, 0, 0],
                         [0, 0, 1, 0],
                         [-3, 3, -2, -1],
                         [2, -2, 1, 1]])  # shape (4, 4)

        A3 = np.asarray([[1, 0, -3, 2],
                         [0, 0, 3, -2],
                         [0, 1, -2, 1],
                         [0, 0, -1, 1]])  # shape 4 x 4

        A2 = np.asarray([[f[i, j], f[i, j + 1], f2[i, j], f2[i, j + 1]],
                         [f[i + 1, j], f[i + 1, j + 1], f2[i + 1, j], f2[i + 1, j + 1]],
                         [f1[i, j], f1[i, j + 1], f3[i, j], f3[i, j + 1]],
                         [f1[i + 1, j], f1[i + 1, j + 1], f3[i + 1, j], f3[i + 1, j + 1]]])  # shape (4, 4, N)

        A2 = np.moveaxis(A2, [0, 1, 2], [1, 2, 0])  # shape (N, 4, 4)
        coeff = np.matmul(A1, np.matmul(A2, A3))  # shape (N, 4, 4)

        # Scale the query points for a non-unitary grid
        u = (xq - x[i]) / (x[i + 1] - x[i])  # shape (N,)
        v = (yq - y[j]) / (y[j + 1] - y[j])  # shape (N,)

        # Prepare the arrays with the polynomial bases
        x_array = np.asarray([u ** 0, u ** 1, u ** 2, u ** 3])[np.newaxis]  # shape (1, 4, N)
        y_array = np.asarray([v ** 0, v ** 1, v ** 2, v ** 3])[np.newaxis]  # shape (1, 4, N)

        # Reorder the arrays for matrix multiplication
        x_array = np.moveaxis(x_array, [0, 1, 2], [1, 2, 0])  # shape (N, 1, 4)
        y_array = np.moveaxis(y_array, [0, 1, 2], [2, 1, 0])  # shape (N, 4, 1)

        # Compute the interpolated values
        fq = np.matmul(x_array, np.matmul(coeff, y_array)).squeeze()  # shape (N,)

        return fq


# -------------------------------------------------------------------------------------------------------------------- #
# Transfinite interpolation class
# -------------------------------------------------------------------------------------------------------------------- #
class TransfiniteInterpolation:

    """ Create a transfinite interpolator for a surface with a boundary described by four parametric curves

    Parameters
    ----------
    C1_func : function returning ndarray with shape (ndim, Nv)
        Function of ´v´ that returns an array containing the coordinates of the west boundary

    C2_func : function returning ndarray with shape (ndim, Nu)
        Function of ´u´ that returns an array containing the coordinates of the south boundary

    C3_func : function returning ndarray with shape (ndim, Nv)
        Function of ´v´ that returns an array containing the coordinates of the east boundary

    C4_func : function returning ndarray with shape (ndim, Nu)
        Function of ´u´ that returns an array containing of the north boundary

    P12 : ndarray with shape (ndim,)
        1D array containing the coordinates of the point connecting the west and south boundaries

    P23 : ndarray with shape (ndim,)
        1D array containing the coordinates of the point connecting the south and east boundaries

    P34 : ndarray with shape (ndim,)
        1D array containing the coordinates of the point connecting the east and north boundaries

    P41 : ndarray with shape (ndim,)
        1D array containing the coordinates of the point connecting the north and west boundaries

    Notes
    -----
    A surface generated by transfinite interpolation is also known as a bilinearly blended Coons patch

    References
    ----------
    Construction of Curvilinear Co-ordinate Systems and Applications to Mesh Generation
    W. J. Gordon and C. A. Hall
    International Journal for Numerical Methods in Engineering (1973)

    The NURBS book. Section 10.6 - Coons Surfaces
    L. Piegl and W. Tiller
    Springer, second edition

    """


    def __init__(self, C1_func, C2_func, C3_func, C4_func, P12, P23, P34, P41):

        # Declare input variables as instance variables
        self.C1_func = C1_func
        self.C2_func = C2_func
        self.C3_func = C3_func
        self.C4_func = C4_func
        self.P12 = P12
        self.P23 = P23
        self.P34 = P34
        self.P41 = P41


    def __call__(self, u, v):

        """ Evaluate the transfinite interpolator for input parameters ´(u,v)´ and return the surface coordinates ´S´

        Parameters
        ----------

        u: scalar or ndarray with shape (N,)
        Parameter in the west-east direction

        v: scalar or ndarray with shape (N,)
        Parameter in the south-north direction

        Returns
        -------
        S : ndarray with shape (ndim, N)
            Array containing the coordinates computed by transfinite interpolation | S(u,v)

        Notes
        -----
        The following formats are allowed for the input ´u´ and ´v´:

            - Both ´u´ and ´v´ are scalars
            - Both ´u´ and ´v´ are one-dimensional arrays with shape (N,)
            - ´u´ is a scalar and ´v´ is a one-dimensional array with shape (N,)
            - ´v´ is a scalar and ´u´ is a one-dimensional array with shape (N,)

        This function does not support ´u´ or ´v´ as two-dimensional arrays
        Use np.flatten() before calling this function if the (u,v) parametrization is in meshgrid format

        """

        # Rename instance variables
        C1_func = self.C1_func
        C2_func = self.C2_func
        C3_func = self.C3_func
        C4_func = self.C4_func
        P12 = self.P12
        P23 = self.P23
        P34 = self.P34
        P41 = self.P41

        # Adjust the data format depending on the shape of the (u,v) input
        if np.isscalar(u) and np.isscalar(v):
            N = 1
            u = np.asarray([u])
            v = np.asarray([v])

        elif np.isscalar(u) and (v.ndim == 1):
            N = np.shape(v)[0]
            u = np.asarray([u])
            u = np.repeat(u, repeats=N, axis=0)

        elif np.isscalar(v) and (u.ndim == 1):
            N = np.shape(u)[0]
            v = np.asarray([v])
            v = np.repeat(v, repeats=N, axis=0)

        elif (u.ndim == 1) and (v.ndim == 1):
            N = np.shape(u)[0]
            Nv = np.shape(v)[0]
            assert (N == Nv), 'u and v must have the same size when they are one-dimensional arrays'

        elif (u.ndim > 1) or (v.ndim > 1):
            raise Exception('u or v arrays with more than one dimension are not supported')

        else:
            raise Exception('Te format of the u or v input is not supported')

        # Evaluate the functions that define the boundaries (before reshaping u and v!)
        C1 = C1_func(v)
        C2 = C2_func(u)
        C3 = C3_func(v)
        C4 = C4_func(u)

        # Number of coordinate dimensions of the problem
        if np.isscalar(P12):
            P12 = np.asarray([P12])
            P23 = np.asarray([P23])
            P34 = np.asarray([P34])
            P41 = np.asarray([P41])

        n_dim = np.shape(P12)[0]

        # Reshape variables so that they are conformable for matrix multiplication
        u = np.repeat(u[np.newaxis, :], repeats=n_dim, axis=0)
        v = np.repeat(v[np.newaxis, :], repeats=n_dim, axis=0)
        P12 = np.repeat(P12[:, np.newaxis], repeats=N, axis=1)
        P23 = np.repeat(P23[:, np.newaxis], repeats=N, axis=1)
        P34 = np.repeat(P34[:, np.newaxis], repeats=N, axis=1)
        P41 = np.repeat(P41[:, np.newaxis], repeats=N, axis=1)

        # Linear interpolation in v between the C1(v) and C3(v) curves
        term_1a = (1 - u) * C1 + u * C3

        # Linear interpolation in u between the C2(v) and C4(v) curves
        term_1b = (1 - v) * C2 + v * C4

        # Bilinear interpolation between the four corners of the domain (P12, P23, P34, P41)
        term_2 = (1 - u) * (1 - v) * P12 + u * v * P34 + (1 - u) * v * P41 + u * (1 - v) * P23

        # Transfinite interpolation formula (also known as bilinearly blended Coons patch)
        # Note that the order of (u,v) is inverted with respect to the original formula
        S = term_1a + term_1b - term_2

        return S
