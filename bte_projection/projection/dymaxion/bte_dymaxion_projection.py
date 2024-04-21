import math

from ...utils import math_utils
from ..oob import OutOfProjectionBoundsException
from .conformal_dynmaxion_projection import ConformalDynmaxionProjection

class BTEDymaxionProjection(ConformalDynmaxionProjection):
    """
    Implementation of the BTE modified Dynmaxion projection.

    :see: DymaxionProjection
    :see: ConformalDynmaxionProjection
    """

    THETA = math.radians(-150)
    SIN_THETA = math.sin(THETA)
    COS_THETA = math.cos(THETA)
    BERING_X = -0.3420420960118339  # -0.3282152608138795
    BERING_Y = -0.322211064085279  # -0.3281491467713469
    ARCTIC_Y = -0.2  # -0.3281491467713469
    ARCTIC_M = (ARCTIC_Y - math_utils.ROOT3 * ConformalDynmaxionProjection.ARC / 4) / (BERING_X - -0.5 * ConformalDynmaxionProjection.ARC)
    ARCTIC_B = ARCTIC_Y - ARCTIC_M * BERING_X
    ALEUTIAN_Y = -0.5000446805492526  # -0.5127463765943157
    ALEUTIAN_XL = -0.5149231279757507  # -0.4957832938238718
    ALEUTIAN_XR = -0.45
    ALEUTIAN_M = (BERING_Y - ALEUTIAN_Y) / (BERING_X - ALEUTIAN_XR)
    ALEUTIAN_B = BERING_Y - ALEUTIAN_M * BERING_X

    def fromGeo(self, longitude: float, latitude: float) -> list[float]:
        c = super().fromGeo(longitude, latitude)
        x = c[0]
        y = c[1]

        easia = self.isEurasianPart(x, y)

        y -= 0.75 * self.ARC * math_utils.ROOT3

        if easia:
            x += self.ARC

            t = x
            x = self.COS_THETA * x - self.SIN_THETA * y
            y = self.SIN_THETA * t + self.COS_THETA * y
        else:
            x -= self.ARC

        c[0] = y
        c[1] = -x
        return c
    
    def toGeo(self, x: float, y: float) -> list[float]:
        if y < 0:
            easia = x > 0
        elif y > self.ARC / 2:
            easia = x > -math_utils.ROOT3 * self.ARC / 2
        else:
            easia = y * -math_utils.ROOT3 < x

        t = x
        x = -y
        y = t

        if easia:
            t = x;
            x = self.COS_THETA * x + self.SIN_THETA * y
            y = self.COS_THETA * y - self.SIN_THETA * t
            x -= self.ARC
        else:
            x += self.ARC

        y += 0.75 * self.ARC * math_utils.ROOT3

        # check to make sure still in right part
        if easia != self.isEurasianPart(x, y):
            raise OutOfProjectionBoundsException
        
        return super().toGeo(x, y)
    
    def isEurasianPart(self, x: float, y: float) -> bool:
        # catch vast majority of cases in not near boundary
        if x > 0:
            return False
        if x < 0.5 * self.ARC:
            return True
        
        if y > math_utils.ROOT3 * self.ARC / 4:  # above arctic ocean
            return x < 0
        
        if y < self.ALEUTIAN_Y:
            return y < (self.ALEUTIAN_Y + self.ALEUTIAN_XL) - x
        
        if y > self.BERING_Y:  # boundary across arctic ocean
            if y < self.ARCTIC_Y:
                return x < self.BERING_X  # in strait
            return y < self.ARCTIC_M * x + self.ARCTIC_B  # above strait
        return y > self.ALEUTIAN_M * x + self.ALEUTIAN_B
    
    def bounds(self) -> list[float]:
        return [
            -1.5 * self.ARC * math_utils.ROOT3,
            -1.5 * self.ARC, 3 * self.ARC,
            math_utils.ROOT3 * self.ARC
        ]
    
    def __str__() -> str:
        return "BuildTheEarth Conformal Dymaxion"
