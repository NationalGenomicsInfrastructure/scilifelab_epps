#!/usr/bin/env python

DESC = """EPP for parsing VC100 output CSV file
Volumes will be filled in by matching well positions
A warning message will be given for the wells that give a more than
threshold volume but should not have been used
Author: Chuan Wang, Science for Life Laboratory, Stockholm, Sweden
"""

import csv
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from scilifelab_epps.epp import EppLogger, set_field

VOL_WARNING_THRESHOLD = 5


def get_vc100_file(process, log):
    for outart in process.all_outputs():
        # get the right output artifact
        if outart.type == "ResultFile" and outart.name == "VC100 CSV File":
            try:
                fid = outart.files[0].id
                content = lims.get_file_contents(id=fid)
            except:
                log.append("Cannot parse VC100 output file")
            break
    return content, log


def get_data(content):
    data = dict()
    headers = dict()
    dialect = csv.Sniffer().sniff(content)
    pf = csv.reader(content.splitlines(), dialect=dialect)
    for line in pf:
        # this is the header row
        if "TUBE" in line:
            for item in line:
                headers[item] = line.index(item)
        else:
            well = line[headers["TUBE"]]
            row = well[0]
            col = str(int(well[1:]))
            new_well = row + ":" + col
            volume = line[headers["VOLAVG"]]
            data[new_well] = volume
    return data


def parse_vc100_results(process):
    # strings returned to the EPP user
    log = []
    # get file contents by parsing lims artifacts
    (content, log) = get_vc100_file(process, log)
    # parse the file and get the interesting data out
    data = get_data(content)
    used_wells = []
    # Fill in LIMS field Volume (ul)
    for target_file in process.result_files():
        well = target_file.samples[0].artifact.location[1]
        used_wells.append(well)
        if well in data:
            # Set to 0 for negative values
            if float(data[well]) > 0:
                target_file.udf["Volume (ul)"] = float(data[well])
            else:
                target_file.udf["Volume (ul)"] = 0
        else:
            log.append(f"Cannot find volume for well {well} in the VC100 CSV file.")
        target_file.put()
        set_field(target_file)

    # Give warning messages if a well supposed to be empty give volume
    wells_with_warning = []
    for k, v in data.items():
        if float(v) > VOL_WARNING_THRESHOLD and k not in used_wells:
            wells_with_warning.append(k)
    log.append(
        f"The following wells are supposed to be empty but give a volume higher than {VOL_WARNING_THRESHOLD}: {','.join(wells_with_warning)}"
    )

    # Throw warnings and errors
    if log:
        sys.stderr.write("; ".join(log))


def main(lims, pid, epp_logger):
    process = Process(lims, id=pid)
    parse_vc100_results(process)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument(
        "--pid", default="24-38458", dest="pid", help="Lims id for current Process"
    )
    parser.add_argument(
        "--log",
        dest="log",
        help=("File name for standard log file, for runtime information and problems."),
    )

    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()

    with EppLogger(log_file=args.log, lims=lims, prepend=True) as epp_logger:
        main(lims, args.pid, epp_logger)
