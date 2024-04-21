import lzma
import json
import math
import io

from .dymaxion_projection import DymaxionProjection
from ...utils import math_utils
from ...data.conformal_data import CONFORMAL_DATA

class ConformalDynmaxionProjection(DymaxionProjection):
    """
    Implementation of the Dynmaxion like conformal projection.
    Slightly modifies the Dynmaxion projection to make it (almost) conformal.

    :see: DymaxionProjection
    """
    VECTOR_SCALE_FACTOR = 1.0 / 1.1473979730192934
    SIDE_LENGTH = 256
    INVERSE_CACHE = None
    
    def __init__(self):
        super().__init__(self)

        xs = [None] * (self.SIDE_LENGTH + 1)
        ys = [None] * len(xs)

        for u in range(len(xs)):
            px = [None] * (len(xs) - u)
            py = [None] * (len(xs) - u)
            xs[u] = px
            ys[u] = py
            
        for v in range(len(xs)):
            for u in range(len(xs)):
                split = CONFORMAL_DATA[u]
                xs[u][v] = split[0] * self.VECTOR_SCALE_FACTOR
                ys[u][v] = split[1] * self.VECTOR_SCALE_FACTOR

        self.inverse = InvertableVectorField(xs, ys)
    
    def triangleTransform(self, vec: list[float]) -> list[float]:
        c = super().triangleTransform(vec)

        x = c[0]
        y = c[1]

        c[0] /= self.ARC
        c[1] /= self.ARC

        c[0] += 0.5
        c[1] += math_utils.ROOT3 / 6
        
        c = self.inverse.applyNewtonsMethod(x, y, c[0], c[1], 5)
        c[0] -= 0.5
        c[1] -= math.sqrt(3) / 6
        c[0] *= self.ARC
        c[1] *= self.ARC
        return c
    
    def inverseTriangleTransform(self, x: float, y: float) -> list[float]:
        x /= self.ARC
        y /= self.ARC

        x += 0.5
        y += math_utils.ROOT3 / 6

        c = self.inverse.getInterpolatedVector(x, y)
        return super().inverseTriangleTransform(c[0], c[1])
    
    def metersPerUnit(self) -> float:
        return (40075017.0 / (2.0 * math.pi)) / self.VECTOR_SCALE_FACTOR
    
    def __str__(self) -> str:
        return "Conformal Dymaxion"

class InvertableVectorField:
    def __init__(self, vx: list[list[float]], vy: list[list[float]]):
        self.vx = vx
        self.vy = vy

    def getInterpolatedVector(self, x: float, y: float) -> list[float]:
        # scale up triangle to be triangleSize across
        x *= ConformalDynmaxionProjection.SIDE_LENGTH
        y *= ConformalDynmaxionProjection.SIDE_LENGTH
        
        # convert to triangle units
        v = 2 * y / math_utils.ROOT3
        u = x - v * 0.5

        u1 = int(u)
        v1 = int(v)

        if u1 < 0:
            u1 = 0
        elif u1 >= ConformalDynmaxionProjection.SIDE_LENGTH:
            u1 = ConformalDynmaxionProjection.SIDE_LENGTH - 1
        
        if v1 < 0:
            v1 = 0
        elif v1 >= ConformalDynmaxionProjection.SIDE_LENGTH - u1:
            v1 = ConformalDynmaxionProjection.SIDE_LENGTH - u1 - 1

        flip = 1
        
        if (y < -math_utils.ROOT3 * (x - u1 - v1 - 1)
            or v1 == ConformalDynmaxionProjection.SIDE_LENGTH - u1 - 1):
            valx1 = self.vx[u1][v1]
            valy1 = self.vy[u1][v1]
            valx2 = self.vx[u1][v1 + 1]
            valy2 = self.vy[u1][v1 + 1]
            valx3 = self.vx[u1 + 1][v1]
            valy3 = self.vy[u1 + 1][v1]

            y3 = 0.5 * math_utils.ROOT3 * v1
            x3 = (u1 + 1) + 0.5 * v1
        else:
            valx1 = self.vx[u1][v1 + 1]
            valy1 = self.vy[u1][v1 + 1]
            valx2 = self.vx[u1 + 1][v1]
            valy2 = self.vy[u1 + 1][v1]
            valx3 = self.vx[u1 + 1][v1 + 1]
            valy3 = self.vy[u1 + 1][v1 + 1]

            flip = -1
            y = -y

            y3 = -(0.5 * math_utils.ROOT3 * (v1 + 1))
            x3 = (u1 + 1) + 0.5 * (v1 + 1)
        
        # TODO: not sure if weights are right
        # (but weirdly mirrors stuff so there may be simplifcation yet)
        w1 = -(y - y3) / math_utils.ROOT3 - (x - x3)
        w2 = 2 * (y - y3) / math_utils.ROOT3
        w3 = 1 - w1 - w2

        return [
            valx1 * w1 + valx2 * w2 + valx3 * w3,
            valy1 * w1 + valy2 * w2 + valy3 * w3,
            (valx3 - valx1) * ConformalDynmaxionProjection.SIDE_LENGTH,
            ConformalDynmaxionProjection.SIDE_LENGTH * flip * (2 * valx2 - valx1 - valx3) / math_utils.ROOT3,
            (valy3 - valy1) * ConformalDynmaxionProjection.SIDE_LENGTH,
            ConformalDynmaxionProjection.SIDE_LENGTH * flip * (2 * valy2 - valy1 - valy3) / math_utils.ROOT3
        ]
    
    def applyNewtonsMethod(self, expectedf: float, expectedg: float, xest: float, yest: float, iter: int) -> list[float]:
        for i in range(iter):
            c = self.getInterpolatedVector(xest, yest)

            f = c[0] - expectedf
            g = c[1] - expectedg
            dfdx = c[2]
            dfdy = c[3]
            dgdx = c[4]
            dgdy = c[5]

            determinant = 1 / (dfdx * dgdy - dfdy * dgdx)

            xest -= determinant * (dgdy * f - dfdy * g)
            yest -= determinant * (-dgdx * f + dfdx * g)

        return [xest, yest]
