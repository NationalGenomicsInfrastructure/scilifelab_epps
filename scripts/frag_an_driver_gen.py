#!/usr/bin/env python

from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

DESC = """EPP used to create csv files to drive the fragment analyzer"""


def main(lims, args):
    currentStep = Process(lims, id=args.pid)
    driver_file_out = None
    driver = []
    ar_driver = {}
    valid_cols = set()
    for output in currentStep.all_outputs():
        if output.name == "Driver File":
            driver_file_out = output
        else:
            if output.location[1]:
                location_ar = output.location[1].split(":")
                valid_cols.add(location_ar[0])
                # idx = (ord(location_ar[0])-65)*12 + int(location_ar[1])-1
                # ar_driver[output.location[1].replace(":", "")] = output.samples[0].name
                ar_driver[output.location[1].replace(":", "")] = output.name.split(
                    "Fragment Analyzer "
                )[1]

    col_idx = -1
    for column in sorted(list(valid_cols)):
        col_idx += 1
        for i in range(1, 13):
            location = f"{column}{i}"
            driver.append(
                (
                    col_idx * 12 + i,
                    location,
                    ar_driver.get(location, "ladder" if i == 12 else ""),
                )
            )

    with open("frag_an_driver.csv", "w") as f:
        for line in driver:
            f.write(f"{line[0]},{line[1]},{line[2]}\n")

    lims.upload_new_file(driver_file_out, "frag_an_driver.csv")


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
