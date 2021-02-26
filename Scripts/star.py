"""
Class that stores data about a star.
"""

class StarData:
    """Class that stores data about a star.

    Parameters
    ----------
    star : str
        ID of the star usually in a standardised format i.e. HD
    main_id : str
        the main ID of the star usually with a more human-readable name i.e. 70 Vir
    tic : str
        the TIC ID of the star
    sectors : list
        the TESS sectors the star was observed in
    spectral_type : str
        the spectral type of the star
    diameter : float
        the angular diameter of the star as measured by interferometry
    diameter_error : float
        the error in the angular diameter
    instrument : str
        the interferometry instrument used to measure the angular diameter

    Attributes
    ----------
    star : str
        ID of the star usually in a standardised format i.e. HD
    main_id : str
        the main ID of the star usually with a more human-readable name i.e. 70 Vir
    tic : str
        the TIC ID of the star
    sectors : list
        the TESS sectors the star was observed in
    spectral_type : str
        the spectral type of the star
    diameter : float
        the angular diameter of the star as measured by interferometry
    diameter_error : float
        the error in the angular diameter
    instrument : str
        the interferometry instrument used to measure the angular diameter
    """

    def __init__(
            self,
            star: str,
            main_id: str,
            tic: str,
            sectors: list,
            spectral_type: str,
            diameter: float,
            diameter_error: float,
            instrument: str
    ) -> None:
        self.star = star
        self.main_id = main_id
        self.tic = tic
        self.sectors = sectors
        self.spectral_type = spectral_type
        self.diameter = diameter
        self.diameter_error = diameter_error
        self.instrument = instrument

    def __repr__(self):
        """Returns the string representation of StarData

        Returns
        -------
        result : str
            the string representation of StarData
        """

        result = f"Star: {self.star}\n"
        result += f"Main ID: {self.main_id}\n"
        result += f"TIC: {self.tic}\n"

        sectors = self.sectors.split(", ")
        if len(sectors) > 4:
            sectors = sectors[:3] + ["..."]
        sectors = ", ".join(sectors)
        result += f"TESS Sectors: {sectors}\n"

        result += f"Spectral Type: {self.spectral_type}\n"
        result += f"Diameter: {self.diameter:.4f}\n"
        result += f"Diameter error: {self.diameter_error:.3f}\n"
        result += f"Instrument: {self.instrument}"

        return result

    def __str__(self):
        """Returns the string representation of StarData

        Returns
        -------
        result : str
            the string representation of StarData
        """

        result = f"Star: {self.star}\n"
        result += f"Main ID: {self.main_id}\n"
        result += f"TIC: {self.tic}\n"

        sectors = self.sectors.split(", ")
        if len(sectors) > 4:
            sectors = sectors[:3] + ["..."]
        sectors = ", ".join(sectors)
        result += f"TESS Sectors: {sectors}\n"

        result += f"Spectral Type: {self.spectral_type}\n"
        result += f"Diameter: {self.diameter:.4f}\n"
        result += f"Diameter error: {self.diameter_error:.3f}\n"
        result += f"Instrument: {self.instrument}"

        return result
