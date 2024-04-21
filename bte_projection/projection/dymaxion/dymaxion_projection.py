import math
import sys

from ..geographic_projection import GeographicProjection
from ..oob import OutOfProjectionBoundsException
from ...utils import math_utils

class DymaxionProjection(GeographicProjection):
    """
    Implementation of the Dynmaxion projection.
    Also known as Airocean or Fuller projection.

    :see: Wikipedia's article on the Dynmaxion projection: https://clck.ru/39fCEa
    """
    ARC = 2 * math.asin(math.sqrt(5 - math.sqrt(5)) / math.sqrt(10))
    Z = math.sqrt(5 + 2 * math.sqrt(5)) / math.sqrt(15)
    EL = math.sqrt(8) / math.sqrt(5 + math.sqrt(5))
    EL6 = EL / 6
    DVE = math.sqrt(3 + math.sqrt(5)) / math.sqrt(5 + math.sqrt(5))
    R = -3 * EL6 / DVE

    NEWTON = 5
    """Number of iterations for Newton's method"""

    VERTICES = [
        [10.536199, 64.700000],
        [-5.245390, 2.300882],
        [58.157706, 10.447378],
        [122.300000, 39.100000],
        [-143.478490, 50.103201],
        [-67.132330, 23.717925],
        [36.521510, -50.103200],
        [112.867673, -23.717930],
        [174.754610, -2.300882],
        [-121.842290, -10.447350],
        [-57.700000, -39.100000],
        [-169.463800, -64.700000],
    ]
    """
    This contains the vertices of the icosahedron,
    identified by their geographic longitude and latitude in degrees.
    When the class is loaded, a static block below converts all these coordinates
    to the equivalent spherical coordinates (longitude and colatitude), in radians.
    
    :see: Wikipedia: https://clck.ru/39fCMk
    """

    ISO = [
        [2, 1, 6],
        [1, 0, 2],
        [0, 1, 5],
        [1, 5, 10],
        [1, 6, 10],
        [7, 2, 6],
        [2, 3, 7],
        [3, 0, 2],
        [0, 3, 4],
        [4, 0, 5],  # 9, qubec
        [5, 4, 9],
        [9, 5, 10],
        [10, 9, 11],
        [11, 6, 10],
        [6, 7, 11],
        [8, 3, 7],
        [8, 3, 4],
        [8, 4, 9],
        [9, 8, 11],
        [7, 8, 11],
        [11, 6, 7],  # child of 14
        [3, 7, 8]  # child of 15
    ]
    """
    Indicates the vertices forming each face of the icosahedron.
    Each entry refers to the index of a vertex in VERTICES
    """

    CENTER_MAP = [
        [-3, 7],
        [-2, 5],
        [-1, 7],
        [2, 5],
        [4, 5],
        [-4, 1],
        [-3, -1],
        [-2, 1],
        [-1, -1],
        [0, 1],
        [1, -1],
        [2, 1],
        [3, -1],
        [4, 1],
        [5, -1],  # 14, left side, right to be cut
        [-3, -5],
        [-1, -5],
        [1, -5],
        [2, -7],
        [-4, -7],
        [-5, -5],  # 20, pseudo triangle, child of 14
        [-2, -7]  # 21, pseudo triangle, child of 15
    ]

    FLIP_TRIANGLE = [
        True, False, True, False, False,
        True, False, True, False, True, False, True, False, True, False,
        True, True, True, False, False,
        True, False
    ]
    """Indicates for each face if it needs to be flipped after projecting"""

    CENTROIDS = [[None] * 22] * 3
    """
    This contains the Cartesian coordinates the centroid
    of each face of the icosahedron.
    """

    ROTATION_MATRICES = [[[None] * 22] * 3] * 3
    """
    Rotation matrices to move the triangles to the reference
    coordinates from the original positions.
    Indexed by the face's indices.
    """

    INVERSE_ROTATION_MATRICES = [[[None] * 22] * 3] * 3
    """
    Rotation matrices to move the triangles to the reference
    coordinates from the original positions.
    Indexed by the face's indices.
    """

    FACE_ON_GRID = [
        -1, -1, 0, 1, 2, -1, -1, 3, -1, 4, -1,
        -1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
        20, 19, 15, 21, 16, -1, 17, 18, -1, -1, -1,
    ]

    @staticmethod
    def __init__(self):
        for i in range(22):
            self.CENTER_MAP[i][0] *= 0.5 * self.ARC
            self.CENTER_MAP[i][1] *= self.ARC * math_utils.ROOT3 / 12

        # Will contain the list of vertices in Cartesian coordinates
        verticesCartesian = [[None] * len(self.VERTICES)] * 3
        print(verticesCartesian)

        # Convert the geographic vertices to spherical in radians
        for i in range(len(self.VERTICES)):
            print(i)
            vertexSpherical = math_utils.geo2Spherical(self.VERTICES[i])
            vertex = math_utils.spherical2Cartesian(vertexSpherical)
            verticesCartesian[i] = vertex
            print(verticesCartesian)
            self.VERTICES[i] = vertexSpherical

        for i in range(22):
            # Vertices of the current face
            vec1 = verticesCartesian[self.ISO[i][0]]
            vec2 = verticesCartesian[self.ISO[i][1]]
            vec3 = verticesCartesian[self.ISO[i][2]]

            # Find the centroid's projection onto the sphere
            xsum = vec1[0] + vec2[0] + vec3[0]
            ysum = vec1[1] + vec2[1] + vec3[1]
            zsum = vec1[2] + vec2[2] + vec3[2]
            mag = math.sqrt(xsum * xsum + ysum * ysum + zsum * zsum)
            self.CENTROIDS[i] = [xsum / mag, ysum / mag, zsum / mag]

            centroidSpherical = math_utils.cartesian2Spherical(self.CENTROIDS[i])
            centroidLambda = centroidSpherical[0]
            centroidPhi = centroidSpherical[1]

            vertex = self.VERTICES[self.ISO[i][0]]
            v = [vertex[0] - centroidLambda, vertex[1]]
            v = self.yRot(v, -centroidPhi)

            self.ROTATION_MATRICES[i] = math_utils.produceZYZRotationMatrix(-centroidLambda, -centroidPhi, (math.pi / 2) - v[0])
            self.INVERSE_ROTATION_MATRICES[i] = math_utils.produceZYZRotationMatrix(v[0] - (math.pi / 2), centroidPhi, centroidLambda)

    @staticmethod
    def findTriangleGrid(self, x: float, y: float) -> int:
        xp = x / self.ARC
        yp = y / (self.ARC * math_utils.ROOT3)

        row = 0
        if yp > -0.25:
            if yp < 0.25:  # middle
                row = 1
            elif yp <= 0.75:  # top
                row = 0
                yp = 0.5 - yp  #translate to middle and flip
            else:
                return -1
        elif yp >= -0.75:  # bottom
            row = 2
            yp = -yp - 0.5  # translate to middle and flip
        else:
            return -1

        yp += 0.25  # change origin to vertex 4, to allow grids to align

        # rotate coords 45 degrees so left and right sides of the triangle become the x/y axies (also side lengths are now 1)
        xr = xp - yp
        yr = xp + yp

        # assign a order to what grid along the y=x line it is
        gx = int(math.floor(xr))
        gy = int(math.floor(yr))
        col = 2 * gx + (1 if gy != gx else 0) + 6

        # out of bounds
        if col < 0 or col >= 11:
            return -1
        
        return self.FACE_ON_GRID[row * 11 + col]  # get face at this position

    @staticmethod
    def yRot(self, spherical: list[float], rot: float) -> list[float]:
        c = math_utils.spherical2Cartesian(spherical)

        x = c[0]
        c[0] = c[2] * math.sin(rot) + x * math.cos(rot)
        c[2] = c[2] * math.cos(rot) - x * math.sin(rot)

        mag = math.sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
        c[0] /= mag
        c[1] /= mag
        c[2] /= mag

        return [
            math.atan2(c[1], c[0]),
            math.atan2(math.sqrt(c[0] * c[0] + c[1] * c[1]), c[2])
        ]

    def findTriangle(self, vector):
        """
        Finds the face of the icosahedron on which to project a point.
        In practice, it works by finding the face
        with the closest centroid to the point.

        :param list vector: - position vector as double array of
        length 3, using Cartesian coordinates
        :return: an integer identifying the face on which to project the point
        """
        min = sys.float_info.max
        face = 0

        for i in range(20):
            xd = self.CENTROIDS[i][0] - vector[0]
            yd = self.CENTROIDS[i][1] - vector[1]
            zd = self.CENTROIDS[i][2] - vector[2]

            dissq = xd * xd + yd * yd + zd * zd
            if dissq < min:
                if dissq < 0.1:  # TODO: enlarge radius
                    return i
                face = i
                min = dissq
        return face

    def triangleTransform(self, vec: list[float]) -> list[float]:
        S = self.Z / vec[2]

        xp = S * vec[0]
        yp = S * vec[1]

        a = math.atan((2 * yp / math_utils.ROOT3 - self.EL6) / self.DVE)  # ARC/2 terms cancel
        b = math.atan((xp - yp / math_utils.ROOT3 - self.EL6) / self.DVE)
        c = math.atan((-xp - yp / math_utils.ROOT3 - self.EL6) / self.DVE)

        return [0.5 * (b - c), (2 * a - b - c) / (2 * math_utils.ROOT3)]

    def inverseTriangleTransformNewton(self, xpp: float, ypp: float) -> list[float]:
        # a & b are linearly related to c, so using the tan of sum formula we know:
        # tan(c+off) = (tanc + tanoff)/(1-tanc*tanoff)
        tanaoff = math.tan(math_utils.ROOT3 * ypp + xpp)  # a = c + root3*y'' + x''
        tanboff = math.tan(2 * xpp)  # b = c + 2x''

        anumer = tanaoff * tanaoff + 1
        bnumer = tanboff * tanboff + 1

        # we will be solving for tanc, starting at t=0, tan(0) = 0
        tana = tanaoff
        tanb = tanboff
        tanc = 0

        adenom = 1
        bdenom = 1

        # double fp = anumer + bnumer + 1; //derivative relative to tanc

        for _ in range(self.NEWTON):
            f = tana + tanb + tanc - self.R  # R = tana + tanb + tanc
            fp = anumer * adenom * adenom + bnumer * bdenom * bdenom + 1  # derivative relative to tanc

            # TODO: fp could be simplified on first loop: 1 + anumer + bnumer

            tanc -= f / fp

            adenom = 1 / (1 - tanc * tanaoff)
            bdenom = 1 / (1 - tanc * tanboff)

            tana = (tanc + tanaoff) * adenom
            tanb = (tanc + tanboff) * bdenom

        # simple reversal algebra based on tan values
        yp = math_utils.ROOT3 * (self.DVE * tana + self.EL6) / 2
        xp = self.DVE * tanb + yp / math_utils.ROOT3 + self.EL6

        # x = z*xp/Z, y = z*yp/Z, x^2 + y^2 + z^2 = 1
        xpoZ = xp / self.Z
        ypoZ = yp / self.Z

        z = 1 / math.sqrt(1 + xpoZ * xpoZ + ypoZ * ypoZ)

        return [z * xpoZ, z * ypoZ, z]

    def inverseTriangleTransform(self, x: float, y: float) -> list[float]:
        return self.inverseTriangleTransformNewton(x, y)

    def fromGeo(self, longitude: float, latitude: float) -> list[float]:
        OutOfProjectionBoundsException.checkLongitudeLatitudeInRange(longitude, latitude)

        vector = math_utils.spherical2Cartesian(math_utils.geo2Spherical([longitude, latitude]))

        face = self.findTriangle(vector)

        # apply rotation matrix (move triangle onto template triangle)
        pvec = math_utils.matVecProdD(self.ROTATION_MATRICES[face], vector)
        projectedVec = self.triangleTransform(pvec)

        # flip triangle to correct orientation
        if self.FLIP_TRIANGLE[face]:
            projectedVec[0] = -projectedVec[0]
            projectedVec[1] = -projectedVec[1]
        
        vector[0] = projectedVec[0]
        # deal with special snowflakes (child faces 20, 21)
        if ((face == 15 and vector[0] > projectedVec[1] * math_utils.ROOT3) or face == 14) and vector[0] > 0:
            projectedVec[0] = 0.5 * vector[0] - 0.5 * math_utils.ROOT3 * projectedVec[1]
            projectedVec[1] = 0.5 * math_utils.ROOT3 * vector[0] + 0.5 * projectedVec[1]
            face += 6  # shift 14->20 & 15->21
        
        projectedVec[0] += self.CENTER_MAP[face][0]
        projectedVec[1] += self.CENTER_MAP[face][1]

        return projectedVec

    def toGeo(self, x: float, y: float) -> list[float]:
        face = self.findTriangleGrid(x, y)
        
        if face == -1:
            raise OutOfProjectionBoundsException

        x -= self.CENTER_MAP[face][0]
        y -= self.CENTER_MAP[face][1]

        # deal with bounds of special snowflakes
        if face == 14:
            if x > 0:
                raise OutOfProjectionBoundsException
        elif face == 20:
            if -y * math_utils.ROOT3 > x:
                raise OutOfProjectionBoundsException
        elif face == 15:
            if x > 0 and x > y * math_utils.ROOT3:
                raise OutOfProjectionBoundsException
        elif face == 21:
            if x < 0 or -y * math_utils.ROOT3 > x:
                raise OutOfProjectionBoundsException

        # flip triangle to upright orientation (if not already)
        if self.FLIP_TRIANGLE[face]:
            x = -x
            y = -y

        # invert triangle transform
        c = self.inverseTriangleTransform(x, y)
        x = c[0]
        y = c[1]
        z = c[2]

        vec = [x, y, z]
        # apply inverse rotation matrix (move triangle from template triangle to correct position on globe)
        vecp = math_utils.matVecProdD(self.INVERSE_ROTATION_MATRICES[face], vec)

        # convert back to geo coordinates
        return math_utils.spherical2Geo(math_utils.cartesian2Spherical(vecp))

    def bounds(self) -> list[float]:
        return [-3 * self.ARC, -0.75 * self.ARC * math_utils.ROOT3, 2.5 * self.ARC, 0.75 * self.ARC * math_utils.ROOT3]

    def upright(self) -> bool:
        return False

    def metersPerUnit(self) -> float:
        return math.sqrt(510100000000000.0 / (20 * math_utils.ROOT3 * self.ARC * self.ARC / 4))

    def __str__(self):
        return "Dymaxion"
