#!/usr/bin/env python

import re
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

DESC = """EPP diluting the concentration of samples in the fragment analyzer step"""


def parse_range_bp(range_value):
    bp_values = re.findall(r"\d+(?:\.\d+)?", range_value)
    if len(bp_values) < 2:
        raise ValueError("Could not parse Range value '{}'".format(range_value))
    return int(round(float(bp_values[0]))), int(round(float(bp_values[1])))


def main(lims, args):
    conc_is_local = True
    process = Process(lims, id=args.pid)
    log_art = None
    log = []
    fid = None
    # first, read the fragment analyzer results
    for o in process.all_outputs():
        if o.name == "CSV Result File":
            try:
                fid = o.files[0].id
            except:
                sys.exit("Please upload a CSV result file.")
        if o.name == "Calculation Log":
            log_art = o

    file_contents = lims.get_file_contents(id=fid)
    if isinstance(file_contents, bytes):
        file_contents = file_contents.decode("utf-8")
    frag_data = {}
    keys = []
    for line in file_contents.splitlines():
        if not keys:
            keys = line.split(",")
        else:
            values = line.split(",")
            frag_data[values[0]] = {}
            for i in range(1, len(values)):
                frag_data[values[0]][keys[i]] = values[i]
    # Then, read the concentration from the step defined in the process udf
    try:
        conc_process_name = process.udf["Concentration Source"]
        conc_is_local = False
    except KeyError:
        conc_is_local = True

    for io in process.input_output_maps:
        if (
            "Fragment Analyzer" in io[1]["uri"].name
            and io[1]["output-generation-type"] == "PerInput"
        ):
            base_concentration = None
            base_conc_unit = None
            well = io[1]["uri"].location[1].replace(":", "")
            if conc_is_local:
                base_concentration = float(frag_data[well]["ng/uL"])
                base_conc_unit = "ng/uL"
            else:
                try:
                    concentration_step = lims.get_processes(
                        type=conc_process_name, inputartifactlimsid=io[0]["limsid"]
                    )[0]
                except IndexError:
                    log.append(
                        "Cannot find a {} step starting with {}".format(
                            conc_process_name, io[0]["limsid"]
                        )
                    )
                else:
                    for io2 in concentration_step.input_output_maps:
                        if (
                            io2[0]["limsid"] == io[0]["limsid"]
                            and "Concentration" in io2[1]["uri"].udf
                        ):
                            base_concentration = io2[1]["uri"].udf["Concentration"]
                            base_conc_unit = io2[1]["uri"].udf["Conc. Units"]
            try:
                min_size_bp, max_size_bp = parse_range_bp(frag_data[well]["Range"])
                io[1]["uri"].udf["Min Size (bp)"] = min_size_bp
                io[1]["uri"].udf["Max Size (bp)"] = max_size_bp
                if "Ratio (%)" not in io[1]["uri"].udf:
                    io[1]["uri"].udf["Ratio (%)"] = float(frag_data[well]["% Total"])
                io[1]["uri"].udf["Size (bp)"] = int(frag_data[well]["Avg. Size"])
                io[1]["uri"].put()

                if base_concentration and base_conc_unit:
                    if conc_is_local:
                        io[1]["uri"].udf["Concentration"] = base_concentration
                    else:
                        io[1]["uri"].udf["Concentration"] = base_concentration * (
                            float(io[1]["uri"].udf["Ratio (%)"]) / 100.0
                        )
                    io[1]["uri"].udf["Conc. Units"] = base_conc_unit
                    io[1]["uri"].put()
                    log.append("Updated values for output {}".format(io[1]["uri"].name))
                else:
                    log.append(
                        "Failed to update the concentration of output {}".format(
                            io[1]["uri"].name
                        )
                    )

            except Exception as e:
                log.append(
                    "Error updating {} with fragment analyzer data : {}".format(
                        io[1]["uri"].name, e
                    )
                )

        if log:
            with open(f"{log_art.id}_frag_analyzer.log", "w") as logContext:
                logContext.write("\n".join(log))


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    parser.add_argument(
        "--read",
        dest="read",
        action="store_true",
        help="reads the output csv and populates the fields",
    )
    parser.add_argument(
        "--calc",
        dest="calc",
        action="store_true",
        help="recalculates the concentration",
    )
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
