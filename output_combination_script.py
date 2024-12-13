#!/usr/bin/python3
"""
output_combination_script.py - Responsible for combining the outputs from HPC runs.

"""

import os

import json
import pandas as pd

from tqdm import tqdm


def main() -> None:
    """Main script for the HPC combination work for the polytunnel-PV outputs."""

    # Check the script is being executed correctly and that all exepcted files are
    # present.
    if os.path.basename(os.getcwd()) != "Polytunnel-PV":
        raise Exception("Script must be launched from the polytunnel root directory.")

    if not os.path.isdir(output_dirname := "outputs_1"):
        raise FileNotFoundError(f"Could not find outputs directory: {output_dirname}.")

    if not os.path.isfile(scenarios_filename := "scenarios.txt"):
        raise FileNotFoundError(
            f"Could not find scenarios.txt file: {scenarios_filename}"
        )

    # Go through each file in the outputs directory and extract the total output data
    # at each time.
    filename_to_july_1st_data_map: dict[str, list[float]] = {}
    filename_to_july_2nd_data_map: dict[str, list[float]] = {}
    for filename in tqdm(
        os.listdir(output_dirname), desc="Parsing files", leave=True, unit="file"
    ):
        with open(
            os.path.join(output_dirname, filename), "r", encoding="UTF-8"
        ) as output_file:
            filedata = json.load(output_file)

        if "4344" in filename:
            filename_to_july_1st_data_map[
                filename.split("hourly_mpp_")[1].split("_4344_to")[0]
            ] = [entry[2] for entry in filedata]
            continue

        if "4368" in filename:
            filename_to_july_2nd_data_map[
                filename.split("hourly_mpp_")[1].split("_4368_to")[0]
            ] = [entry[2] for entry in filedata]
            continue

    time_series = [entry[1] for entry in filedata]
    filename_to_july_1st_data_map["hour"] = time_series
    filename_to_july_2nd_data_map["hour"] = time_series

    combined_july_1st_frame = pd.DataFrame(filename_to_july_1st_data_map)
    with open(
        "combined_bypass_diode_investigation_july_1st.csv", "w", encoding="UTF-8"
    ) as combined_file:
        combined_july_1st_frame.to_csv(combined_file)

    combined_july_2nd_frame = pd.DataFrame(filename_to_july_2nd_data_map)
    with open(
        "combined_bypass_diode_investigation_july_2nd.csv", "w", encoding="UTF-8"
    ) as combined_file:
        combined_july_2nd_frame.to_csv(combined_file)


if __name__ == "__main__":
    main()
