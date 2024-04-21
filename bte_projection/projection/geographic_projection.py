import math
from abc import ABC, abstractmethod

from ..constants import EARTH_CIRCUMFERENCE, EARTH_POLAR_CIRCUMFERENCE
from ..utils import math_utils
from .oob import OutOfProjectionBoundsException

class GeographicProjection(ABC):
    """
    Support for various projection types.

    The geographic space is the surface of the earth, parameterized by the
    usual spherical coordinates system of latitude and longitude. The
    projected space is a plane on to which the geographic space is being
    projected, and is parameterized by a 2D Cartesian coordinate system
    (x and y).

    A projection as defined here is something that projects a point in the
    geographic space to a point of the projected space (and vice versa).

    All geographic coordinates are in degrees.
    """
    
    @abstractmethod
    def toGeo(self, x: float, y: float) -> list[float]:
        """
        Converts map coordinates to geographic coordinates

        :param float x: - x map coordinate
        :param float y: - y map coordinate
        :return: {longitude, latitude} in degrees
        :throws OutOfProjectionBoundsException: if the specified point on
        the projected space cannot be mapped to a point of the geographic space
        """
        
        raise NotImplementedError
    
    @abstractmethod
    def fromGeo(self, longitude: float, latitude: float) -> list[float]:
        """
        Converts geographic coordinates to map coordinates

        :param float longitude: - longitude, in degrees
        :param float latitude: - latitude, in degrees
        :return: {x, y} map coordinates
        :throws OutOfProjectionBoundsException: if the specified point on
        the geographic space cannot be mapped to a point of the projected space
        """
        raise NotImplementedError
    
    @abstractmethod
    def metersPerUnit(self) -> float:
        """
        Gives an estimation of the scale of this projection.
        This is just an estimation, as distortion is inevitable when
        projecting a sphere onto a flat surface, so this value varies from
        places to places in reality.
        
        :return: an estimation of the scale of this projection
        """
        raise NotImplementedError
    
    def bounds(self) -> list[float]:
        """
        Indicates the minimum and maximum X and Y coordinates on the projected space.

        :return: {minimum X, minimum Y, maximum X, maximum Y}
        """
        try:
            # get max in by using extreme coordinates
            bounds = [
                self.fromGeo(-180, 0)[0],
                self.fromGeo(0, -90)[1],
                self.fromGeo(180, 0)[0],
                self.fromGeo(0, 90)[1]
            ]

            if bounds[0] > bounds[2]:
                bounds[0], bounds[2] = bounds[2], bounds[0]
            
            if bounds[1] > bounds[3]:
                bounds[1], bounds[3] = bounds[3], bounds[1]
            
            return bounds
        except OutOfProjectionBoundsException:
            return [0, 0, 1, 1]
    
    def upright(self) -> bool:
        """
        Indicates whether or not the north pole is projected to the north of
        the south pole on the projected space, assuming Minecraft's coordinate
        system cardinal directions for the projected space (north is negative Z).
     
        :return: north pole Z <= south pole Z
        """
        try:
            return self.fromGeo(0, 90)[1] <= self.fromGeo(0, -90)[1]
        except OutOfProjectionBoundsException:
            return False
    
    def vector(self, x: float, y: float, north: float, east: float) -> list[float]:
        """
        Calculates the vector that goes a given distance north and a given
        distance east from the given point in the projected space.
        This is useful to get a direction in the projected space, e.g. it
        is used to calculate the north vector used when sending eyes of ender.

        :param float x: - x coordinate in the projected space
        :param float y: - y coordinate in the projected space
        :param float north: - how far north to go, in meters on the geographic space
        :param float south: - how far south to go, in meters on the geographic space
        :return: {distance x, distance y} on the projected space
        """
        geo = self.toGeo(x, y)

        # TODO: east may be slightly off because earth not a sphere
        off = self.fromGeo(
            geo[0] + east * 360.0 / (math.cos(math.radians(geo[1])) * EARTH_CIRCUMFERENCE),
            geo[1] + north * 360.0 / EARTH_POLAR_CIRCUMFERENCE
        )
        return [off[0] - x, off[1] - y]
    
    def tissot(self, longitude: float, latitude: float) -> list[float]:
        """
        Computes the Tissot's indicatrix of this projection
        at the given point (i.e. the distortion).

        :param float longitude: - longitude, in degrees
        :param float latitude: - latitude, in degrees
        :return: {area inflation, maximum angular distortion, maximum scale
        factor, minimum scale factor}
        :see: Wikipedia's article on Tissot's indicatrix: https://clck.ru/39ex3x
        """
        R = EARTH_CIRCUMFERENCE / (2 * math.pi)
        D = 1E-7


        ddeg = math.degrees(D)

        base = self.fromGeo(longitude, latitude)
        lonoff = self.fromGeo(longitude + ddeg, latitude)
        latoff = self.fromGeo(longitude, latitude + ddeg)

        dxdl = (lonoff[0] - base[0]) / D
        dxdp = (latoff[0] - base[0]) / D
        dydl = (lonoff[1] - base[1]) / D
        dydp = (latoff[1] - base[1]) / D

        cosp = math.cos(math.radians(latitude))

        h = math.sqrt(dxdp * dxdp + dydp * dydp) / R
        k = math.sqrt(dxdl * dxdl + dydl * dydl) / (cosp * R)

        sint = abs(dydp * dxdl - dxdp * dydl) / (R * R * cosp * h * k)
        ap = math.sqrt(h * h + k * k + 2 * h * k * sint)
        bp = math.sqrt(h * h + k * k - 2 * h * k * sint)

        a = (ap + bp) / 2
        b = (ap - bp) / 2

        return [h * k * sint, 2 * math.asin(bp / ap), a, b]
    
    def azimuth(self, x: float, y: float, angle: float) -> float:
        """
        Converts an angle in the projected space to an azimuth in the
        geographic space, at a specific point. This is useful to get the
        direction an entity is looking at. With conformal projections,
        this should be equivalent to using
        `GeographicProjection.vector(double, double, double, double)`
        and computing the facing azimuth in the projected space, but on
        non-conformal projections angles are not preserved when projecting and
        this will be right when using
        `GeographicProjection.vector(double, double, double, double)`
        is likely to be wrong.

        :param float x: - x coordinate of the point in the projected space
        :param float y: - y coordinate of the point in the projected space
        :param float angle: - the angle to convert, in degrees, in
        minecraft's coordinate system (angular origin at the positive
        side of the Z axis, positive clockwise)
        :return: the corresponding azimuth, in degrees, counted positively
        clockwise, between 0° and 360°.
        :throws OutOfProjectionBoundsException: if the given point is
        outside the projection domain
        """
        D = 1E-5

        x2 = x - D * math.sin(math.radians(angle))
        y2 = y + D * math.cos(math.radians(angle))

        geo1 = self.toGeo(x, y)
        geo2 = self.toGeo(x2, y2)

        math_utils.toRadians(geo1)
        math_utils.toRadians(geo2)

        dlon = geo2[0] - geo1[0]
        dlat = geo2[1] - geo1[1]

        a = math.degrees(math.atan2(dlat, dlon * math.cos(geo1[1])))
        a = 90 - a        
        if a < 0:
            a += 360

        return a
    
    def properties(self) -> dict[str, object]:
        return {}
