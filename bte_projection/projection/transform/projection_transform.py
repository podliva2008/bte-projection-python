from ..geographic_projection import GeographicProjection

class ProjectionTransform(GeographicProjection):
    """Warps a Geographic projection and applies a transformation to it."""

    def __init__(self, delegate: GeographicProjection):
        """:param GeographicProjection delegate: - projection to transform"""
        self.delegate = delegate()

    def upright(self) -> bool:
        return self.delegate.upright()
    
    def bounds(self) -> list[float]:
        return self.delegate.bounds()
    
    def metersPerUnit(self) -> float:
        return self.delegate.metersPerUnit()
