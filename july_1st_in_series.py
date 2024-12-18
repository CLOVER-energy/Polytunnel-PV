"""Run for July 1st"""

import subprocess

# import warnings

from tqdm import tqdm

BASE_COMMAND: str = (
    "python -m src.polytunnelpv --scenario {scenario} -st 4344 -i 24 --operating-mode hourly_mpp --timestamps-file hourly_december_data_martyn.csv"
)

# Filter warnings
# warnings.filterwarnings("ignore")


def main() -> None:
    """Main function."""

    # Open the scenarios file and parse the data.
    with open("scenarios.txt", "r", encoding="UTF-8") as scenarios_file:
        scenarios = scenarios_file.readlines()

    try:
        for scenario in tqdm(scenarios, desc="Scenarios", leave=True):
            subprocess.run(BASE_COMMAND.format(scenario=scenario.strip()).split(" "))
    except BaseException:
        print("FAILED")
        print(f"Parsing of scenario {scenario} failed.")
        raise


if __name__ == "__main__":
    main()
