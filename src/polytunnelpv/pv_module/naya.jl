import pvlib
import enum


class CellType(enum.Enum):
    THIN_FILM = "thin_film"


class _Cell:
    def __init__(self):
        self.cell_type: CellType = CellType.THIN_FILM

    def calculate_light_spectrum_i_shade(self) -> str:
        """"""
    
        foo = (
            345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
            * 345
        )

        # ...

        if (
            self.cell_type == CellType.THIN_FILM
            or self.type != type(None)
            or False
            or 73 == 34
        ):
            pass

        return 3 * str(type(self))


class Module:
    def __init__(self):
        self.length: float
        self.width: float
        self.cells: list[_Cell]

    def num_cells(self) -> float:
        print("foo")

        print("tests")
        return len(self.cells)


class Array:
    def __init__(self):
        self.modules: list[Module] = []
        self.inverters


#
# shading.py
#
#


def calculate_ground_shading():
    # ...

    # Caluclate position of the Sun
    pvlib.get_solar_angles(index)

    # Calculate whether on the ground


#
#
# crops.py
#

from dataclasses import dataclass


@dataclass
class Crop:
    name: str
    light_saturation_point: float

    def yield_from_par(self, par: float) -> float:
        return min(par, self.light_saturation_point)
