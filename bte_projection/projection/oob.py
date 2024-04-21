class OutOfProjectionBoundsException(Exception):
    @staticmethod
    def checkInRange(x: float, y: float, maxX: float, maxY: float) -> None:
        """
        :param float x:
        :param float y:
        :param float maxX:
        :param float maxY:
        :raises OutOfProjectionBoundsException: if `abs(x) > maxX or abs(y) > maxY`
        """
        if abs(x) > maxX or abs(y) > maxY:
            raise OutOfProjectionBoundsException

    @staticmethod
    def checkLongitudeLatitudeInRange(longitude: float, latitude: float) -> None:
        """
        :param float longitude:
        :param float latitude:
        :raises OutOfProjectionBoundsException: if `abs(longitude) > 180 or abs(latitude) > 90`
        """
        OutOfProjectionBoundsException.checkInRange(longitude, latitude, 180, 90)
