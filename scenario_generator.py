#!/usr/bin/python3
"""
scenario_generator.py

File which generates a scenarios.json file containing all the different bypass-diode
configurations along with a pv_modules.json file.

"""

import dataclasses
import json
import os
import sys

from typing import Any

from tqdm import tqdm

# ALLOWED_DIODE_LENGTHS:
#    Allowed lengths of diode.
ALLOWED_DIODE_LENGTHS: list[int] = [33, 30, 25, 20, 15, 10, 5, 2]

# MIN_DIODE_LENGTH:
#   The minimum length for bypass diodes.
MIN_DIODE_LENGTH: int = 2


@dataclasses.dataclass
class BypassDiode:
    start: int
    end: int

    @property
    def length(self) -> int:
        """Return the length of the diode."""

        return self.end - self.start

    def as_output(self) -> dict[str, float | int]:
        """
        Return an output which is savable.

        :returns: a nice-looking output which can be saved.

        """

        return {
            "start_index": self.start,
            "end_index": self.end,
            "bypass_voltage": -0.5,
        }


class ConfigurationError(Exception):
    """Raised when the configuration is invalid."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__("Too many configurations were attempted.")


BASE_PANEL_CONFIG: dict[str, Any] = {
    "name": "mia_sole_flex_03_280nl_bd_n_{config_index}",
    "cell_type": "Miasole_FLEX_03_280NL",
    "type": "thin_film",
    "cell_breakdown_voltage": -4,
    "cell_length": 0.04473,  # [m] - The dimension parallel to the module
    "cell_spacing": 0.0,  # [m] - The spacing between the cells
    "cell_width": 0.348,  # [m] - The dimension perpendicular to the module
    "module_centre_offset": 0,
    "n_cells": 132,
    "offset_angle": 90,
    "polytunnel_curve": "tilted_kent_circle",
}

BYPASS_DIODE_KEY: str = "bypass_diodes"

BASE_SCENARIO: dict[str, Any] = {
    "name": "circular_polysolar_kent_bd_n_{config_index}",
    "location": "polysolar_kent",
}

PV_MODULE_KEY: str = "pv_module"


def _recursive_bypassing(
    remaining_cells: list[int],
    remaining_diodes: int,
    bypass_diode_configuration: list[BypassDiode],
    bypass_diode_configurations: list[list[BypassDiode]],
) -> None:
    """
    Recursively bypass.

    From the remaining cells, iterate through all combinations that a bypass diode could
    bypass, and add the current configuration to the list of possible configurations.

    :param: remaining_cells
        The remaining cells which can be bypassed.

    :param: remaining_dioes
        The number of remaining bypass diodes which can be assigned.

    :param: bypass_diode_configuration
        The list of possible bypass-diode configurations.

    """

    # Remove one diode to work on.
    remaining_diodes -= 1

    # Iterate through all possible start positions.
    min_diode_start = min(remaining_cells)
    max_diode_end = max(remaining_cells) + 1 - MIN_DIODE_LENGTH * remaining_diodes
    # for diode_start in tqdm(range(
    #     min_diode_start := min(remaining_cells),
    #     max_diode_end := ,
    # ), desc="diode start", leave=False):

    # Iterate through all possible diode lengths.
    for diode_end in tqdm(
        range(min_diode_start + MIN_DIODE_LENGTH, max_diode_end),
        desc="diode end",
        leave=False,
    ):
        # Skip configurations which are too short.
        if (diode_length := diode_end - min_diode_start) not in ALLOWED_DIODE_LENGTHS:
            continue

        # Skip configurations where the diode length is increasing towards the centre.
        if len(bypass_diode_configuration) > 0:
            if diode_length > bypass_diode_configuration[-1].length:
                continue

        # If this is the final diode being considered, then add the current config.
        if remaining_diodes == 0:
            bypass_diode_configurations.append(
                bypass_diode_configuration + [BypassDiode(min_diode_start, diode_end)]
            )

            if len(bypass_diode_configurations) > 10000:
                raise ConfigurationError()

            continue

        # Otherwise, call the function again with the remaining diodes.
        # Add the current diode configuration to the config.
        bypass_diode_configuration.append(BypassDiode(min_diode_start, diode_end))

        # Allow for diodes only to be bypassed which are more central.
        _recursive_bypassing(
            list(range(diode_end, max(remaining_cells) + 1)),
            remaining_diodes,
            bypass_diode_configuration,
            bypass_diode_configurations,
        )

        bypass_diode_configuration.pop()


def main(args: list[Any]) -> None:
    """
    Main function.

    :param: args
        Un-parsed command-line arguments.

    """

    bypass_diode_configurations: list[list[BypassDiode]] = []

    try:
        for num_diodes in tqdm(
            range(1, 67), desc="number of bypass dioes", unit="dioes"
        ):
            _recursive_bypassing(
                list(range(67)), num_diodes, [], bypass_diode_configurations
            )
    except ConfigurationError:
        pass

    print(len(bypass_diode_configurations))

    # Mirror and combine the configurations.
    mirrored_configurations: list[list[BypassDiode]] = [
        [
            BypassDiode(start=132 - diode.end, end=132 - diode.start)
            for diode in configuration
        ]
        for configuration in bypass_diode_configurations
    ]
    combined_configurations: list[list[BypassDiode]] = [
        configuration + mirrored_configurations[index]
        for index, configuration in enumerate(bypass_diode_configurations)
    ]

    # Generate a set of output panels and matching scenarios.
    output_panels: list[dict[str, Any]] = []
    scenarios: list[dict[str, str]] = []
    for index, configuration in tqdm(
        enumerate(combined_configurations), desc="generating output panels"
    ):

        # Generate the panel data
        this_panel = BASE_PANEL_CONFIG.copy()
        this_panel["name"] = this_panel["name"].format(config_index=index)
        this_panel[BYPASS_DIODE_KEY] = [entry.as_output() for entry in configuration]
        output_panels.append(this_panel)

        # Do the same for the scenario.
        this_scenario = BASE_SCENARIO.copy()
        this_scenario[PV_MODULE_KEY] = this_panel["name"]
        this_scenario["name"] = this_scenario["name"].format(config_index=index)
        scenarios.append(this_scenario)

    # Save to the files.
    with open(
        os.path.join((input_directory := "input_data"), "pv_modules.json"), "w"
    ) as panel_file:
        json.dump(output_panels, panel_file)

    with open(os.path.join(input_directory, "scenarios.json"), "w") as scenarios_file:
        json.dump(scenarios, scenarios_file)

    with open("scenarios.txt", "w") as scenarios_text_file:
        scenarios_text_file.write(
            "\n".join([scenario["name"] for scenario in scenarios])
        )


if __name__ == "__main__":
    main(sys.argv[1:])

#!/usr/bin/python3
"""
scenario_generator.py

File which generates a scenarios.json file containing all the different bypass-diode
configurations along with a pv_modules.json file.

"""

import dataclasses
import json
import os
import sys

from typing import Any

from tqdm import tqdm

# ALLOWED_DIODE_LENGTHS:
#    Allowed lengths of diode.
ALLOWED_DIODE_LENGTHS: list[int] = [33, 30, 25, 20, 15, 10, 5, 2]

# MIN_DIODE_LENGTH:
#   The minimum length for bypass diodes.
MIN_DIODE_LENGTH: int = 2


@dataclasses.dataclass
class BypassDiode:
    start: int
    end: int

    @property
    def length(self) -> int:
        """Return the length of the diode."""

        return self.end - self.start

    def as_output(self) -> dict[str, float | int]:
        """
        Return an output which is savable.

        :returns: a nice-looking output which can be saved.

        """

        return {
            "start_index": self.start,
            "end_index": self.end,
            "bypass_voltage": -0.5,
        }


class ConfigurationError(Exception):
    """Raised when the configuration is invalid."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__("Too many configurations were attempted.")


BASE_PANEL_CONFIG: dict[str, Any] = {
    "name": "mia_sole_flex_03_280nl_bd_n_{config_index}",
    "cell_type": "Miasole_FLEX_03_280NL",
    "type": "thin_film",
    "cell_breakdown_voltage": -4,
    "cell_length": 0.04473,  # [m] - The dimension parallel to the module
    "cell_spacing": 0.0,  # [m] - The spacing between the cells
    "cell_width": 0.348,  # [m] - The dimension perpendicular to the module
    "module_centre_offset": 0,
    "n_cells": 132,
    "offset_angle": 90,
    "polytunnel_curve": "tilted_kent_circle",
}

BYPASS_DIODE_KEY: str = "bypass_diodes"

BASE_SCENARIO: dict[str, Any] = {
    "name": "circular_polysolar_kent_bd_n_{config_index}",
    "location": "polysolar_kent",
}

PV_MODULE_KEY: str = "pv_module"


def _recursive_bypassing(
    remaining_cells: list[int],
    remaining_diodes: int,
    bypass_diode_configuration: list[BypassDiode],
    bypass_diode_configurations: list[list[BypassDiode]],
) -> None:
    """
    Recursively bypass.

    From the remaining cells, iterate through all combinations that a bypass diode could
    bypass, and add the current configuration to the list of possible configurations.

    :param: remaining_cells
        The remaining cells which can be bypassed.

    :param: remaining_dioes
        The number of remaining bypass diodes which can be assigned.

    :param: bypass_diode_configuration
        The list of possible bypass-diode configurations.

    """

    # Remove one diode to work on.
    remaining_diodes -= 1

    # Iterate through all possible start positions.
    min_diode_start = min(remaining_cells)
    max_diode_end = max(remaining_cells) + 1 - MIN_DIODE_LENGTH * remaining_diodes
    # for diode_start in tqdm(range(
    #     min_diode_start := min(remaining_cells),
    #     max_diode_end := ,
    # ), desc="diode start", leave=False):

    # Iterate through all possible diode lengths.
    for diode_end in tqdm(
        range(min_diode_start + MIN_DIODE_LENGTH, max_diode_end),
        desc="diode end",
        leave=False,
    ):
        # Skip configurations which are too short.
        if (diode_length := diode_end - min_diode_start) not in ALLOWED_DIODE_LENGTHS:
            continue

        # Skip configurations where the diode length is increasing towards the centre.
        if len(bypass_diode_configuration) > 0:
            if diode_length > bypass_diode_configuration[-1].length:
                continue

        # If this is the final diode being considered, then add the current config.
        if remaining_diodes == 0:
            bypass_diode_configurations.append(
                bypass_diode_configuration + [BypassDiode(min_diode_start, diode_end)]
            )

            if len(bypass_diode_configurations) > 10000:
                raise ConfigurationError()

            continue

        # Otherwise, call the function again with the remaining diodes.
        # Add the current diode configuration to the config.
        bypass_diode_configuration.append(BypassDiode(min_diode_start, diode_end))

        # Allow for diodes only to be bypassed which are more central.
        _recursive_bypassing(
            list(range(diode_end, max(remaining_cells) + 1)),
            remaining_diodes,
            bypass_diode_configuration,
            bypass_diode_configurations,
        )

        bypass_diode_configuration.pop()


def main(args: list[Any]) -> None:
    """
    Main function.

    :param: args
        Un-parsed command-line arguments.

    """

    bypass_diode_configurations: list[list[BypassDiode]] = []

    try:
        for num_diodes in tqdm(
            range(1, 67), desc="number of bypass dioes", unit="dioes"
        ):
            _recursive_bypassing(
                list(range(67)), num_diodes, [], bypass_diode_configurations
            )
    except ConfigurationError:
        pass

    print(len(bypass_diode_configurations))

    # Mirror and combine the configurations.
    mirrored_configurations: list[list[BypassDiode]] = [
        [
            BypassDiode(start=132 - diode.end, end=132 - diode.start)
            for diode in configuration
        ]
        for configuration in bypass_diode_configurations
    ]
    combined_configurations: list[list[BypassDiode]] = [
        configuration + mirrored_configurations[index]
        for index, configuration in enumerate(bypass_diode_configurations)
    ]

    # Generate a set of output panels and matching scenarios.
    output_panels: list[dict[str, Any]] = []
    scenarios: list[dict[str, str]] = []
    for index, configuration in tqdm(
        enumerate(combined_configurations), desc="generating output panels"
    ):

        # Generate the panel data
        this_panel = BASE_PANEL_CONFIG.copy()
        this_panel["name"] = this_panel["name"].format(config_index=index)
        this_panel[BYPASS_DIODE_KEY] = [entry.as_output() for entry in configuration]
        output_panels.append(this_panel)

        # Do the same for the scenario.
        this_scenario = BASE_SCENARIO.copy()
        this_scenario[PV_MODULE_KEY] = this_panel["name"]
        this_scenario["name"] = this_scenario["name"].format(config_index=index)
        scenarios.append(this_scenario)

    # Save to the files.
    with open(
        os.path.join((input_directory := "input_data"), "pv_modules.json"), "w"
    ) as panel_file:
        json.dump(output_panels, panel_file)

    with open(os.path.join(input_directory, "scenarios.json"), "w") as scenarios_file:
        json.dump(scenarios, scenarios_file)

    with open("scenarios.txt", "w") as scenarios_text_file:
        scenarios_text_file.write(
            "\n".join([scenario["name"] for scenario in scenarios])
        )


if __name__ == "__main__":
    main(sys.argv[1:])
