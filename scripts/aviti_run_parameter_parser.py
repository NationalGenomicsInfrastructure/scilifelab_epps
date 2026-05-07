#!/usr/bin/env python

import glob
import json
import os
import sys
from argparse import ArgumentParser
from datetime import datetime

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

DESC = """EPP for parsing run paramters for AVITI runs
Author: Chuan Wang, Science for Life Laboratory, Stockholm, Sweden
"""


def fetch_fc(process):
    fc_id = ""
    if "Load to Flowcell (AVITI)" in process.parent_processes()[0].type.name:
        fc_id = process.parent_processes()[0].output_containers()[0].name
    else:
        sys.stderr.write("No associated parent step can be found.")
        sys.exit(2)
    return fc_id


def fetch_rundir(fc_id):
    run_dir = ""
    metadata_dir = "ngi-nas-ns"
    data_dir = "Aviti_data"
    run_dir_path = os.path.join(os.sep, "srv", metadata_dir, data_dir, f"*{fc_id}")
    if len(glob.glob(run_dir_path)) == 1:
        run_dir = glob.glob(run_dir_path)[0]
    elif len(glob.glob(run_dir_path)) == 0:
        sys.stderr.write(f"No run dir can be found for FC {fc_id}")
        sys.exit(2)
    else:
        sys.stderr.write(f"Multiple run dirs found for FC {fc_id}")
        sys.exit(2)
    return run_dir


def attach_json_files(process, run_dir):
    for outart in process.all_outputs():
        if outart.type == "ResultFile" and outart.name == "Run Parameters":
            try:
                lims.upload_new_file(outart, f"{run_dir}/RunParameters.json")
            except OSError:
                sys.stderr.write("No RunParameters.json found")
                sys.exit(2)
        elif outart.type == "ResultFile" and outart.name == "Run Stats":
            try:
                lims.upload_new_file(outart, f"{run_dir}/AvitiRunStats.json")
            except OSError:
                sys.stderr.write("No AvitiRunStats.json found")


def parse_run_parameters(run_dir):
    if os.path.exists(f"{run_dir}/RunParameters.json"):
        with open(f"{run_dir}/RunParameters.json") as run_parameters_json:
            run_parameters = json.load(run_parameters_json)
        return run_parameters
    else:
        sys.stderr.write(f"No RunParameters.json found in path {run_dir}")
        sys.exit(2)


def set_step_udfs(process, run_dir):
    run_parameters = parse_run_parameters(run_dir)

    process.udf["Run ID"] = run_parameters.get("RunFolderName")
    process.udf["Flow Cell ID"] = run_parameters.get("FlowcellID")
    process.udf["Side"] = run_parameters.get("Side")
    process.udf["Run Series"] = run_parameters.get("RunID")

    if run_parameters.get("Cycles"):
        process.udf["Read 1 Cycles"] = run_parameters["Cycles"].get("R1", 0)
        process.udf["Read 2 Cycles"] = run_parameters["Cycles"].get("R2", 0)
        process.udf["Index Read 1"] = run_parameters["Cycles"].get("I1", 0)
        process.udf["Index Read 2"] = run_parameters["Cycles"].get("I2", 0)

    process.udf["Read Order"] = run_parameters.get("ReadOrder")
    process.udf["Throughput Selection"] = run_parameters.get("ThroughputSelection")
    process.udf["Kit Configuration"] = run_parameters.get("KitConfiguration")
    process.udf["Preparation Workflow"] = run_parameters.get("PreparationWorkflow")
    process.udf["Chemistry Version"] = run_parameters.get("ChemistryVersion")
    process.udf["Low Diversity"] = str(run_parameters.get("LowDiversity"))
    process.udf["Platform Version"] = run_parameters.get("PlatformVersion")
    process.udf["Analysis Lanes"] = run_parameters.get("AnalysisLanes")
    process.udf["Library Type"] = run_parameters.get("LibraryType")

    if run_parameters.get("Consumables"):
        if run_parameters["Consumables"].get("Flowcell"):
            process.udf["Flowcell Serial Number"] = run_parameters["Consumables"][
                "Flowcell"
            ].get("SerialNumber")
            process.udf["Flowcell Part Number"] = run_parameters["Consumables"][
                "Flowcell"
            ].get("PartNumber")
            process.udf["Flowcell Lot Number"] = run_parameters["Consumables"][
                "Flowcell"
            ].get("LotNumber")
            process.udf["Flowcell Expiration Date"] = datetime.strptime(
                run_parameters["Consumables"]["Flowcell"].get("Expiration")[0:10],
                "%Y-%m-%d",
            ).date()
        if run_parameters["Consumables"].get("SequencingCartridge"):
            process.udf["Sequencing Cartridge Serial Number"] = run_parameters[
                "Consumables"
            ]["SequencingCartridge"].get("SerialNumber")
            process.udf["Sequencing Cartridge Part Number"] = run_parameters[
                "Consumables"
            ]["SequencingCartridge"].get("PartNumber")
            process.udf["Sequencing Cartridge Lot Number"] = run_parameters[
                "Consumables"
            ]["SequencingCartridge"].get("LotNumber")
            process.udf["Sequencing Cartridge Expiration Date"] = datetime.strptime(
                run_parameters["Consumables"]["SequencingCartridge"].get("Expiration")[
                    0:10
                ],
                "%Y-%m-%d",
            ).date()
        if run_parameters["Consumables"].get("Buffer"):
            process.udf["Buffer Serial Number"] = run_parameters["Consumables"][
                "Buffer"
            ].get("SerialNumber")
            process.udf["Buffer Part Number"] = run_parameters["Consumables"][
                "Buffer"
            ].get("PartNumber")
            process.udf["Buffer Lot Number"] = run_parameters["Consumables"][
                "Buffer"
            ].get("LotNumber")
            process.udf["Buffer Expiration Date"] = datetime.strptime(
                run_parameters["Consumables"]["Buffer"].get("Expiration")[0:10],
                "%Y-%m-%d",
            ).date()

    process.put()


def parse_run_stats(run_dir):
    if os.path.exists(f"{run_dir}/AvitiRunStats.json"):
        with open(f"{run_dir}/AvitiRunStats.json") as run_stats_json:
            run_stats = json.load(run_stats_json)
        return run_stats
    else:
        sys.stderr.write(f"No AvitiRunStats.json found in path {run_dir}")
        sys.exit(2)


def calculate_mean(input_list, key):
    values = [d[key] for d in input_list if key in d and d[key] > 0]
    mean_value = sum(values) / len(values) if values else ""
    return mean_value


def set_run_stats(process, run_dir):
    run_stats = parse_run_stats(run_dir)
    for art in process.all_outputs():
        if "Lane" in art.name:
            lane_nbr = int(art.name.split(" ")[1])
            lanes = [d["Lane"] for d in run_stats["LaneStats"]]
            # When there is no runmanifest provided, the Lane number will be displayed as 0
            # In this case we have to parse the lanes in order
            if lane_nbr not in lanes:
                if lane_nbr <= len(run_stats["LaneStats"]):
                    lane_stats = run_stats["LaneStats"][lane_nbr - 1]
                else:
                    sys.stderr.write("Inconsistent lane number detected!")
                    sys.exit(2)
            else:
                lane_stats = next(
                    d for d in run_stats["LaneStats"] if d["Lane"] == lane_nbr
                )
            for read in lane_stats["Reads"]:
                read_key = read["Read"]
                art.udf[f"Reads PF (M) {read_key}"] = lane_stats["PFCount"] / 1000000
                art.udf[f"%PF {read_key}"] = lane_stats["PercentPF"]
                art.udf[f"Yield PF (Gb) {read_key}"] = (
                    lane_stats["TotalYield"] / 1000000000
                )
                art.udf[f"% Aligned {read_key}"] = read["PhiXAlignmentRate"]
                art.udf[f"% Bases >=Q30 {read_key}"] = calculate_mean(
                    read["Cycles"], "PercentQ30"
                )
                art.udf[f"% Bases >=Q40 {read_key}"] = calculate_mean(
                    read["Cycles"], "PercentQ40"
                )
                art.udf[f"Avg Q Score {read_key}"] = calculate_mean(
                    read["Cycles"], "AverageQScore"
                )
                art.udf[f"% Error Rate {read_key}"] = calculate_mean(
                    read["Cycles"], "PercentPhixErrorRate"
                )
            art.put()
    process.put()


def main(lims, args):
    process = Process(lims, id=args.pid)
    # Fetch FC ID
    fc_id = fetch_fc(process)
    # Fetch run dir
    run_dir = fetch_rundir(fc_id)
    # Attach json files
    attach_json_files(process, run_dir)
    # Set step UDFs
    set_step_udfs(process, run_dir)
    # Set run stats
    set_run_stats(process, run_dir)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
