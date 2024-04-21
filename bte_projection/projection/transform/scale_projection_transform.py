import math

from .projection_transform import ProjectionTransform
from ..geographic_projection import GeographicProjection
from ..oob import OutOfProjectionBoundsException

class ScaleProjectionTransform(ProjectionTransform):
    """
    Scales the warps projection's projected space up or down.
    More specifically, it multiplies x and y by there respective scale factors.
    """

    def __init__(self, delegate: GeographicProjection, x: float, y: float):
        """
        Creates a new ScaleProjection with different scale factors for the x and y axis.

        :param GeographicProjection delegate: - projection to transform
        :param float x: - scaling to apply along the x axis
        :param float y: - scaling to apply along the y axis
        """
        super().__init__(delegate)

        assert math.isfinite(x) and math.isfinite(y), "Projection scales should be finite"
        assert x != 0 and y != 0, "Projection scale cannot be 0!"

        self.x = x
        self.y = y

    def toGeo(self, x: float, y: float) -> list[float]:
        return super().delegate.toGeo(x / self.x, y / self.y)
    
    def fromGeo(self, longitude: float, latitude: float) -> list[float]:
        p = super().delegate.fromGeo(longitude, latitude)
        p[0] *= self.x
        p[1] *= self.y
        return p
    
    def upright(self) -> bool:
        return (self.y < 0) ^ self.delegate.upright()
    
    def bounds(self) -> list[float]:
        b = super().delegate.bounds()
        b[0] *= self.x
        b[1] *= self.y
        b[2] *= self.x
        b[3] *= self.y
        return b
    
    def metersPerUnit(self) -> float:
        return super().delegate.metersPerUnit() / math.sqrt((self.x * self.x + self.y * self.y) / 2)  # TODO: better transform
    
    def __str__(self) -> str:
        return "Scale (" + super().delegate + ") by " + self.x + ", " + self.y
