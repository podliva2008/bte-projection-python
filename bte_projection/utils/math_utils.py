import math

ROOT3 = math.sqrt(3)
"""Square root of 3"""

TAU = 2 * math.pi
"""Two times pi"""

def geo2Spherical(geo: list[float]) -> list[float]:
    """
    Converts geographic latitude and longitude coordinates
    to spherical coordinates on a sphere of radius 1.

    :param list geo: - geographic coordinates as a double array of length 2,
    {longitude, latitude}, in degrees
    :return: the corresponding spherical coordinates in radians:
    {longitude, colatitude}
    """
    lambda_ = math.radians(geo[0])
    phi = math.radians(90 - geo[1])

    return [lambda_, phi]

def spherical2Geo(spherical: list[float]) -> list[float]:
    """
    Converts spherical coordinates to geographic coordinates
    on a sphere of radius 1.

    :param list spherical: - spherical coordinates in radians as a double
    array of length 2: {longitude, colatitude}
    :return: the corresponding geographic coordinates in degrees:
    {longitude, latitude}
    """
    lon = math.degrees(spherical[0])
    lat = 90 - math.degrees(spherical[1])

    return [lon, lat]

def spherical2Cartesian(spherical: list[float]) -> list[float]:
    """
    Converts spherical coordinates to Cartesian
    coordinates on a sphere of radius 1.

    :param list spherical: - spherical coordinates in radians
    as a double array of length 2: {longitude, colatitude}
    :return: the corresponding Cartesian coordinates: {x, y, z}
    """
    sinphi = math.sin(spherical[1])
    x = sinphi * math.cos(spherical[0])
    y = sinphi * math.sin(spherical[0])
    z = math.cos(spherical[1])

    return [x, y, z]

def cartesian2Spherical(cartesian: list[float]) -> list[float]:
    """
    Converts Cartesian coordinates to spherical coordinates on a sphere of radius 1.

    :param cartesian: - Cartesian coordinates as double array of length 3: {x, y, z}
    :return: the spherical coordinates of the corresponding normalized vector
    """
    lambda_ = math.atan2(cartesian[1], cartesian[0])
    phi = math.atan2(math.sqrt(cartesian[0]**2 + cartesian[1]**2), cartesian[2])

    return [lambda_, phi]

def produceZYZRotationMatrix(a: float, b: float, c: float) -> list[list[float]]:
    # TODO: make a docstring
    sina = math.sin(a)
    cosa = math.cos(a)
    sinb = math.sin(b)
    cosb = math.cos(b)
    sinc = math.sin(c)
    cosc = math.cos(c)
    
    mat = [[None] * 3] * 3
    mat[0][0] = cosa * cosb * cosc - sinc * sina
    mat[0][1] = -sina * cosb * cosc - sinc * cosa
    mat[0][2] = cosc * sinb

    mat[1][0] = sinc * cosb * cosa + cosc * sina
    mat[1][1] = cosc * cosa - sinc * cosb * sina
    mat[1][2] = sinc * sinb

    mat[2][0] = -sinb * cosa
    mat[2][1] = sinb * sina
    mat[2][2] = cosb
    
    return mat

def matVecProdD(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """
    Multiples the given matrix with the given vector.
    The matrix is assumed to be square and the vector is assumed
    to be of the same dimension as the matrix.

    :param list matrix: - the matrix as a n*n double array
    :param list vector: - the vector as double array of length n
    :return: the result of the multiplication as an array of double on length n
    """
    result = [None] * len(vector)
    for i in range(len(result)):
        for j in range(len(matrix[i])):
            result[i] = matrix[i][j] * vector[j]

    return result

def toRadians(arr: list[float]) -> None:
    """
    Converts all values in a double array from degrees to radians

    :param list arr: - array to work on
    """
    for i in len(arr):
        arr[i] = math.radians(i)

def toDegrees(arr: list[float]) -> None:
    """
    Converts all values in a double array from radians to degrees

    :param list arr: - array to work on
    """
    for i in len(arr):
        arr[i] = math.degrees(i)

def safeDirectionalShift(val: int, shift: int) -> int:
    """
    Right-shifts the given value by the given number of bits, safely handling negative shifts and checking for overflow.

    :param int val: - the value
    :param int shift: - the number of bits to shift by
    :return: the shifted value
    """
    if shift == 0:
        res = val
    elif shift > 0:
        res = val << shift
        assert (res >> shift) == val, f"numeric overflow: val: {val}, shift: {shift}"
    else:
        res = val >> -shift
    return res
