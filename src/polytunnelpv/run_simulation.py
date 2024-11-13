import os
import yaml
from typing import Any
from .__main__ import main as polytunnelpv_main

# HPC Job Number:
#   Name of the environment variable for the HPC job number.
HPC_JOB_NUMBER: str = "PBS_ARRAY_INDEX"


def load_scenarios(file_path: str) -> list:
    with open(file_path, "r") as file:
        scenarios = yaml.safe_load(file)
    return scenarios


def main(args: list[Any]) -> None:
    """
    Wrapper around Polytunnel-PV when run on the HPC.
    """
    # Determine the run that is to be carried out.
    try:
        hpc_job_number = int(os.getenv(HPC_JOB_NUMBER))  # type: ignore
    except (ValueError, TypeError) as e:
        print(
            f"HPC environmental variable {HPC_JOB_NUMBER} was not of type int or not set.",
            e,
        )
        raise

    # Determine the run.
    run_number: int = hpc_job_number - 1

    # Load scenarios from the file.
    scenarios = load_scenarios("input_data/scenarios.yaml")

    # Get the scenario name based on the job number.
    try:
        scenario = scenarios[run_number]["name"]
    except IndexError:
        print(
            f"Invalid job number {run_number}. Check the number of scenarios available."
        )
        raise

    # Setup the arguments to pass to Polytunnel-PV.
    polytunnel_pv_arguments = [
        "--scenario",
        scenario,
        "--operating-mode",
        "hourly_mpp",
        "-st",
        "0",
        "-i",
        "8760",
    ]

    # Call the main function from Polytunnel-PV with the scenario arguments.
    polytunnelpv_main(polytunnel_pv_arguments)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
