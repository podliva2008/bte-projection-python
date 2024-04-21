from .projection_transform import ProjectionTransform
from ..geographic_projection import GeographicProjection
from ..oob import OutOfProjectionBoundsException

class FlipVerticalProjectionTransform(ProjectionTransform):
    """
    Mirrors the warped projection vertically.
    I.E. x' = x and y' = -y
    """

    def __init__(self, delegate: GeographicProjection):
        """:param GeographicProjection delegate: - projection to transform"""
        super().__init__(delegate)

    def toGeo(self, x: float, y: float) -> list[float]:
        return self.delegate.toGeo(x, -y)
    
    def fromGeo(self, longitude: float, latitude: float) -> list[float]:
        p = self.delegate.fromGeo(longitude, latitude)
        p[1] = -p[1]
        return p
    
    def upright(self) -> bool:
        return not self.delegate.upright()
    
    def bounds(self) -> list[float]:
        b = self.delegate.bounds()
        return [b[0], -b[3], b[2], -b[1]]
    
    def __str__(self) -> str:
        return "Vertical Flip (" + super().delegate + ')'
