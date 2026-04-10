#!/usr/bin/python3
# Script will update the qsub validation scrpit and submit.

import os
import subprocess

# HOURS_PER_DAY:
#   The number of hours per day.
HOURS_PER_DAY: int = 24

# SCRIPTS_DIR:
#   The directory in which scripts are saved.
SCRIPTS_DIR: str = "bin"

# START_TIME:
#   The start time to use.
START_TIME: int = 2184

# VALIDATION_FILENAME:
#   The name of the validation script to use.
VALIDATION_FILENAME: str = "validation_april.sh"

def main() -> None:
    """
    Main script that will loop through the hours and qsub.

    """

    with open(os.path.join(SCRIPTS_DIR, VALIDATION_FILENAME), "r"), encoding="UTF-8") as validation_file:
        base_validation_data: str = validation_file.read()

    for start_date in range(31):
        start_time = START_TIME + start_date * hours_per_day

        # Update the script information and write to a file.
        with open(_temp_launchscript:=os.path.join(SCRIPTS_DIR, f"temp_april_{start_date+1}_{VALIDATION_FILENAME}"), "w", encoding="UTF-8") as temp_validation_file:
            temp_validation_file.write(base_validation_data.format(START_TIME=start_time))

        # Run the file, then remove.
        try:
            subprocess.run(f"qsub _temp_launchscript".split(" "))
        except:
            print("Failed to launch script.")
            raise
        else:
            print("Launch script successfully launched.")

        try:
            os.remove(_temp_launchscript)
        except:
            pass
        else:
            print("Temp launchscript removed")


if __name__ == "__main__":
    main()

