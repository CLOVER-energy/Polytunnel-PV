"""Run for July 1st"""

import subprocess

# import warnings

from tqdm import tqdm

# April 20th 
BASE_COMMAND: str = (
    "python -m src.polytunnelpv --scenario {scenario} --operating-mode {operating_mode} -st {start_time} -i {iteration_length} "
)

# Filter warnings
# warnings.filterwarnings("ignore")


def main() -> None:
    """Main function."""

    # # Open the scenarios file and parse the data.
    # with open("scenarios.txt", "r", encoding="UTF-8") as scenarios_file:
    #     scenarios = scenarios_file.readlines()
    start_time = 0
    iteration_length = [24,48,72,96,120]
    
    
    scenarios = ['circular_polysolar_kent']
                # 'circular_polysolar_kent_bypassed_2',
                # 'circular_polysolar_kent_bypassed_10']

    modes = ['hourly_mpp']#['hourly_mpp','hourly_mpp_interpolation']#
    try:
        for scenario in tqdm(scenarios, desc="Scenarios", leave=True):
            for mode in modes:
                for it in iteration_length:
                    command = BASE_COMMAND.format(scenario=scenario,operating_mode=mode,start_time=start_time, iteration_length=it)
                    print(command)
                    subprocess.run(command)
    except BaseException:
        print("FAILED")
        print(f"Parsing of scenario {scenario} failed.")
        raise


if __name__ == "__main__":
    main()
