#!/usr/bin/env python
DESC = """
This file together with manage_demux_stats_thresholds.py performs the "bclconversion" step of LIMS workflow.
In common tongue, it:

Fetches info from the sequencing process (RunID, FCID; derives instrument and data type)
Assigns (Q30, Clust per Lane) thresholds to the process (workflow step)
Reformats laneBarcode.html to "demuxstats_FCID.csv" for usage of other applications
Assigns a lot of info from laneBarcode.html to individual samples of the process (e.g. %PF)
Flags samples as QC PASSED/FAILED based on thresholds

Written by Isak Sylvin; isak.sylvin@scilifelab.se"""

import csv
import logging
import os
import re
import sys
from argparse import ArgumentParser
from shutil import move

import flowcell_parser.classes as classes
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from manage_demux_stats_thresholds import Thresholds

logger = logging.getLogger("demux_logger")


def my_float(value):
    if value == "":
        return 0.0
    else:
        return float(value)


def problem_handler(type, message):
    if type == "exit":
        logger.error(message)
        sys.exit(message)
    elif type == "warning":
        logger.warning(message)
        sys.stderr.write(message)
    else:
        logger.info(message)


def get_process_stats(demux_process):
    """Fetches overarching process info"""
    seq_processes = {
        "MiSeq Run (MiSeq) 4.0",
        "AUTOMATED - NovaSeq Run (NovaSeq 6000 v2.0)",
        "Illumina Sequencing (NextSeq) v1.0",
        "NovaSeqXPlus Run v1.0",
        "AVITI Run v1.0",
        "Illumina Sequencing (MiSeq i100) v1.0",
    }
    try:
        # Query LIMS for all steps containing the first input artifact of this step and match to the set of sequencing steps
        seq_process = lims.get_processes(
            inputartifactlimsid=demux_process.all_inputs()[0].id, type=seq_processes
        )[0]
    except Exception as e:
        problem_handler("exit", f"Undefined prior workflow step (run type): {str(e)}")
    # Copies LIMS sequencing step content
    proc_stats = dict(list(seq_process.udf.items()))
    # Instrument is denoted the way it is since it is also used to find
    # the folder of the laneBarcode.html file
    if "MiSeq Run (MiSeq) 4.0" == seq_process.type.name:
        if seq_process.udf["Run Type"] == "null":
            proc_stats["Chemistry"] = "MiSeq"
        else:
            proc_stats["Chemistry"] = seq_process.udf["Run Type"]
        proc_stats["Instrument"] = "miseq"
        proc_stats["Read Length"] = (
            max(seq_process.udf["Read 1 Cycles"], seq_process.udf["Read 2 Cycles"])
            if seq_process.udf.get("Read 2 Cycles")
            else seq_process.udf["Read 1 Cycles"]
        )
        proc_stats["Paired"] = True if seq_process.udf.get("Read 2 Cycles") else False

    elif "AUTOMATED - NovaSeq Run (NovaSeq 6000 v2.0)" == seq_process.type.name:
        try:
            proc_stats["Chemistry"] = seq_process.udf["Flow Cell Mode"]
        except Exception as e:
            problem_handler(
                "exit", f"No flowcell version set in sequencing step: {str(e)}"
            )
        proc_stats["Instrument"] = "NovaSeq"
        proc_stats["Read Length"] = (
            max(seq_process.udf["Read 1 Cycles"], seq_process.udf["Read 2 Cycles"])
            if seq_process.udf.get("Read 2 Cycles")
            else seq_process.udf["Read 1 Cycles"]
        )
        proc_stats["Paired"] = True if seq_process.udf.get("Read 2 Cycles") else False

    elif "NovaSeqXPlus Run" in seq_process.type.name:
        try:
            proc_stats["Chemistry"] = seq_process.udf["Flow Cell Mode"]
        except Exception as e:
            problem_handler(
                "exit", f"No flowcell version set in sequencing step: {str(e)}"
            )
        proc_stats["Instrument"] = "NovaSeqXPlus"
        proc_stats["Read Length"] = (
            max(seq_process.udf["Read 1 Cycles"], seq_process.udf["Read 2 Cycles"])
            if seq_process.udf.get("Read 2 Cycles")
            else seq_process.udf["Read 1 Cycles"]
        )
        proc_stats["Paired"] = True if seq_process.udf.get("Read 2 Cycles") else False

    elif "Illumina Sequencing (NextSeq) v1.0" == seq_process.type.name:
        try:
            proc_stats["Chemistry"] = seq_process.udf["Chemistry"]
        except Exception as e:
            problem_handler("exit", f"No run type set in sequencing step: {str(e)}")
        proc_stats["Instrument"] = "NextSeq"
        proc_stats["Read Length"] = (
            max(seq_process.udf["Read 1 Cycles"], seq_process.udf["Read 2 Cycles"])
            if seq_process.udf.get("Read 2 Cycles")
            else seq_process.udf["Read 1 Cycles"]
        )
        proc_stats["Paired"] = True if seq_process.udf.get("Read 2 Cycles") else False

    elif "Illumina Sequencing (MiSeq i100) v1.0" == seq_process.type.name:
        raw_chem = seq_process.udf["Chemistry"]

        # Extract 5M, 25M, 50M, 100M if present
        match = re.search(r"(5M|25M|50M|100M)", raw_chem)
        if match:
            proc_stats["Chemistry"] = match.group(1)
        else:
            proc_stats["Chemistry"] = raw_chem  # fallback
        proc_stats["Instrument"] = "MiSeqi100"
        proc_stats["Read Length"] = (
            max(seq_process.udf["Read 1 Cycles"], seq_process.udf["Read 2 Cycles"])
            if seq_process.udf.get("Read 2 Cycles")
            else seq_process.udf["Read 1 Cycles"]
        )

        proc_stats["Paired"] = True if seq_process.udf.get("Read 2 Cycles") else False

    elif "AVITI Run" in seq_process.type.name:
        try:
            proc_stats["Chemistry"] = (
                "AVITI" + " " + seq_process.udf["Throughput Selection"]
            )
        except Exception as e:
            problem_handler(
                "exit", f"No flowcell version set in sequencing step: {str(e)}"
            )
        proc_stats["Instrument"] = "Aviti"
        proc_stats["Read Length"] = (
            max(seq_process.udf["Read 1 Cycles"], seq_process.udf["Read 2 Cycles"])
            if seq_process.udf.get("Read 2 Cycles")
            else seq_process.udf["Read 1 Cycles"]
        )
        proc_stats["Paired"] = True if seq_process.udf.get("Read 2 Cycles") else False

    else:
        problem_handler("exit", "Unhandled workflow step (run type)")
    logger.info("Run type/chemistry set to {}".format(proc_stats["Chemistry"]))
    logger.info("Instrument set to {}".format(proc_stats["Instrument"]))

    try:
        proc_stats["Paired"] = proc_stats.get("Paired", False)
    except Exception as e:
        problem_handler("exit", f"Unable to fetch workflow information: {str(e)}")
    if "Read 2 Cycles" in proc_stats:
        proc_stats["Paired"] = True
    logger.info("Paired libraries: {}".format(proc_stats["Paired"]))

    # Assignment to make usage more explicit
    try:
        proc_stats["Read Length"] = (
            proc_stats["Read Length"]
            if proc_stats.get("Read Length", None)
            else max(proc_stats["Read 1 Cycles"], proc_stats["Read 2 Cycles"])
        )
    except Exception as e:
        problem_handler(
            "exit",
            f"Read Cycles not found. Unable to read Read Length: {str(e)}",
        )

    logger.info("Read length set to {}".format(proc_stats["Read Length"]))
    return proc_stats


def fill_process_fields(demux_process, process_stats):
    """Sets run thresholds"""
    thresholds = Thresholds(
        process_stats["Instrument"],
        process_stats["Chemistry"],
        process_stats["Paired"],
        process_stats["Read Length"],
    )

    if "Threshold for % bases >= Q30" not in demux_process.udf:
        thresholds.set_Q30()
        try:
            demux_process.udf["Threshold for % bases >= Q30"] = thresholds.Q30
            logger.info(
                "Q30 threshold set to {}".format(
                    demux_process.udf["Threshold for % bases >= Q30"]
                )
            )
        except Exception as e:
            problem_handler(
                "exit",
                f"Udf improperly formatted. Unable to set Q30 threshold: {str(e)}",
            )

    # Would REALLY prefer "Minimum Reads per Lane" over "Threshold for # Reads"
    if "Minimum Reads per Lane" not in demux_process.udf:
        thresholds.set_exp_lane_clust()
        try:
            demux_process.udf["Minimum Reads per Lane"] = thresholds.exp_lane_clust
            logger.info(
                "Minimum clusters per lane set to {}".format(
                    demux_process.udf["Minimum Reads per Lane"]
                )
            )
        except Exception as e:
            problem_handler(
                "exit",
                f"Udf improperly formatted. Unable to set # Reads threshold: {str(e)}",
            )

    # Would REALLY prefer "Maximum % Undetermined Reads per Lane" over "Threshold for Undemultiplexed Index Yield"
    if "Maximum % Undetermined Reads per Lane" not in demux_process.udf:
        thresholds.set_undet_indexes_perc()
        try:
            demux_process.udf["Maximum % Undetermined Reads per Lane"] = (
                thresholds.undet_indexes_perc
            )
            logger.info(
                "Maximum percentage of undetermined per lane set to {} %".format(
                    demux_process.udf["Maximum % Undetermined Reads per Lane"]
                )
            )
        except Exception as e:
            problem_handler(
                "exit",
                f"Udf improperly formatted. Unable to set Undemultiplexed Index Yield threshold: {str(e)}",
            )

    # Sets Run ID if not already exists:
    if "Run ID" not in demux_process.udf:
        try:
            demux_process.udf["Run ID"] = process_stats["Run ID"]
        except Exception as e:
            logger.info(f"Unable to automatically regenerate Run ID: {str(e)}")

    # Checks for document version
    if "Document Version" not in demux_process.udf:
        problem_handler("exit", "No Document Version set. Please set one.")

    try:
        demux_process.put()
    except Exception as e:
        problem_handler("exit", f"Failed to apply process thresholds to LIMS: {str(e)}")


def set_sample_values(demux_process, parser_struct, process_stats):
    """Sets artifact = sample values"""

    thresholds = Thresholds(
        process_stats["Instrument"],
        process_stats["Chemistry"],
        process_stats["Paired"],
        process_stats["Read Length"],
    )
    failed_entries = 0
    undet_included = False
    noIndex = False
    undet_lanes = list()
    proj_pattern = re.compile(r"(P\w+_\d+)")

    # Necessary for noindexruns, should always resolve
    try:
        seq_processes = {
            "MiSeq Run (MiSeq) 4.0",
            "AUTOMATED - NovaSeq Run (NovaSeq 6000 v2.0)",
            "Illumina Sequencing (NextSeq) v1.0",
            "NovaSeqXPlus Run v1.0",
            "AVITI Run v1.0",
            "Illumina Sequencing (MiSeq i100) v1.0",
        }
        seq_process = lims.get_processes(
            inputartifactlimsid=demux_process.all_inputs()[0].id, type=seq_processes
        )[0]
    except Exception as e:
        problem_handler("exit", f"Undefined prior workflow step (run type): {str(e)}")

    if "Lanes to include undetermined" in demux_process.udf:
        try:
            undet_lanes = re.split(
                "[ ,.]", demux_process.udf["Lanes to include undetermined"]
            )
            undet_lanes = [int(i) for i in undet_lanes]
        except:
            problem_handler(
                "exit",
                "Unable to typecast included undetermined lanes. Possibly non-number in list",
            )

    # Prevent multiple different flowcells in the same step
    container_names = [pool.container.name for pool in demux_process.all_inputs()]
    assert len(set(container_names)) == 1, (
        f"All input pools must be in the same flowcell, found: {container_names}"
    )
    run_id = process_stats["Run ID"]
    if process_stats["Instrument"] in ["NextSeq", "miseq"]:
        run_id = process_stats["Reagent Cartridge ID"]
    assert container_names[0] in run_id, (
        f"Flowcell name {container_names[0]} seems unrelated to run {run_id}"
    )

    for pool in demux_process.all_inputs():
        undet_reads = 0
        lane_reads = 0
        undet_lane_reads = 0
        samplesum = dict()

        try:
            outarts_per_lane = []
            for art_tuple in demux_process.input_output_maps:
                if (
                    art_tuple[0]["uri"].id == pool.id
                    and art_tuple[1]["output-generation-type"] == "PerReagentLabel"
                ):
                    outarts_per_lane.append(art_tuple[1]["uri"])
        except Exception as e:
            problem_handler("exit", f"Unable to fetch artifacts of process: {str(e)}")
        if process_stats["Instrument"] in ("miseq", "MiSeqi100"):
            lane_no = "1"

        else:
            try:
                lane_no = pool.location[1][0]
            except Exception as e:
                problem_handler(
                    "exit",
                    f"Unable to determine lane number. Incorrect location variable in process: {str(e)}",
                )
        logger.info(f"Lane number set to {lane_no}")
        try:
            correction_factor = (
                thresholds.correction_factor_for_sample_in_pool
                if len(outarts_per_lane) > 1
                else 1
            )
            exp_smp_per_lne = round(
                demux_process.udf["Minimum Reads per Lane"]
                / my_float(len(outarts_per_lane))
                * correction_factor,
                0,
            )
        except ZeroDivisionError as e:
            problem_handler(
                "exit",
                f"Faulty LIMS setup. Pool in lane {lane_no} has no samples: {e}",
            )
        logger.info(f"Expected sample clusters for this lane: {exp_smp_per_lne}")

        # Artifacts in each lane
        for target_file in outarts_per_lane:
            try:
                # This block addresses a LIMS bug in which multiple samples are tied to the same demux artifact
                # In this case, try to find a single sample name matching the name of the demux artifact
                if len(target_file.samples) > 1:
                    logging.warning(
                        f"Multiple samples {[i.name for i in target_file.samples]} linked to demux artifact '{target_file.name}' ({target_file.id})."
                    )
                    matching_names = [
                        sample.name
                        for sample in target_file.samples
                        if sample.name in target_file.name
                    ]
                    if len(matching_names) > 1:
                        raise AssertionError(
                            f"Multiple linked samples '{matching_names}' matches demux artifact name. Cannot tell which is which."
                        )
                    elif len(matching_names) == 0:
                        raise AssertionError(
                            "No linked sample matches demux artifact name."
                        )
                    else:
                        current_name = matching_names[0]
                else:
                    current_name = target_file.samples[0].name
            except Exception as e:
                problem_handler(
                    "exit",
                    f"Unable to determine sample name. Incorrect sample variable in process: {str(e)}",
                )
            for entry in parser_struct:
                if lane_no == entry["Lane"]:
                    sample = entry["Sample"]
                    # Finds name subset "P Anything Underscore Digits"
                    if sample != "Undetermined":
                        try:
                            sample = proj_pattern.search(sample).group(0)
                        # PhiX cases for AVITI
                        except AttributeError:
                            pass

                    if (
                        entry["Barcode sequence"] == "unknown"
                        and sample != "Undetermined"
                    ):
                        noIndex = True
                        if undet_included:
                            problem_handler(
                                "error",
                                "Logical error, undetermined cannot be included for a noIndex lane!",
                            )

                    # Bracket for adding undetermined to results
                    if not sample == "Undetermined" and int(lane_no) in undet_lanes:
                        undet_included = True
                        # Sanity check for including undetermined
                        # Next entry is undetermined and previous is for a different lane
                        current_index = parser_struct.index(entry)
                        undet = parser_struct[current_index + 1]
                        if (
                            undet["Sample"] == "Undetermined"
                            and parser_struct[current_index - 1]["Lane"] != lane_no
                        ):
                            try:
                                clusterType = None
                                if "PF Clusters" in undet:
                                    clusterType = "PF Clusters"
                                else:
                                    clusterType = "Clusters"
                                # Paired runs are divided by two within flowcell parser
                                if process_stats["Paired"]:
                                    undet_reads = (
                                        int(undet[clusterType].replace(",", "")) * 2
                                    )
                                # Since a single ended run has no pairs, pairs is set to equal reads
                                else:
                                    undet_reads = int(
                                        undet[clusterType].replace(",", "")
                                    )
                                logger.info(
                                    f"Included undetermined for lane number {lane_no}"
                                )
                            except Exception as e:
                                problem_handler(
                                    "exit",
                                    f"Unable to set values for undetermined #Reads and #Read Pairs: {str(e)}",
                                )
                        else:
                            problem_handler(
                                "exit",
                                f"Undetermined for lane {lane_no} requested, which has more than one sample",
                            )

                    # Bracket for adding typical sample info
                    if sample == current_name:
                        # Sample samplesum construction
                        if sample not in samplesum:
                            samplesum[sample] = dict()
                            samplesum[sample]["count"] = 1
                        else:
                            samplesum[sample]["count"] += 1

                        try:
                            def_atr = {
                                "% of thelane": "% of Raw Clusters Per Lane",
                                "% Perfectbarcode": "% Perfect Index Read",
                                "% One mismatchbarcode": "% One Mismatch Reads (Index)",
                                "Yield (Mbases)": "Yield PF (Gb)",
                                "% PFClusters": "%PF",
                                "Mean QualityScore": "Ave Q Score",
                                "% >= Q30bases": "% Bases >=Q30",
                            }
                            for old_attr, attr in def_atr.items():
                                # Sets default value for unwritten fields
                                if old_attr in entry.keys():
                                    if (
                                        entry[old_attr] == ""
                                        or entry[old_attr] == "NaN"
                                    ):
                                        if old_attr == "% of Raw Clusters Per Lane":
                                            default_value = 100.0
                                        else:
                                            default_value = 0.0

                                        samplesum[sample][attr] = (
                                            default_value
                                            if attr not in samplesum[sample]
                                            else samplesum[sample][attr] + default_value
                                        )
                                        logger.info(
                                            f"{attr} field not found. Setting default value: {default_value}"
                                        )

                                    else:
                                        # Yields needs division by 1K, is also non-percentage
                                        if old_attr == "Yield (Mbases)":
                                            samplesum[sample][attr] = (
                                                my_float(
                                                    entry[old_attr].replace(",", "")
                                                )
                                                / 1000
                                                if attr not in samplesum[sample]
                                                else samplesum[sample][attr]
                                                + my_float(
                                                    entry[old_attr].replace(",", "")
                                                )
                                                / 1000
                                            )
                                        else:
                                            samplesum[sample][attr] = (
                                                my_float(entry[old_attr])
                                                if attr not in samplesum[sample]
                                                else samplesum[sample][attr]
                                                + my_float(entry[old_attr])
                                            )

                        except Exception as e:
                            problem_handler(
                                "exit",
                                f"Unable to set artifact values. Check laneBarcode.html for odd values: {str(e)}",
                            )

                        # Fetches clusters from laneBarcode.html file
                        if noIndex:
                            # For the case of NovaSeq run, parse lane yield from the ResultsFile of all_outputs.
                            if seq_process.type.name in [
                                "AUTOMATED - NovaSeq Run (NovaSeq 6000 v2.0)",
                                "Illumina Sequencing (NextSeq) v1.0",
                                "NovaSeqXPlus Run v1.0",
                                "AVITI Run v1.0",
                                "Illumina Sequencing (MiSeq i100) v1.0",
                            ]:
                                try:
                                    for inp in seq_process.all_outputs():
                                        if (
                                            inp.output_type == "ResultFile"
                                            and inp.name.split(" ")[1] == lane_no
                                            and "Reads PF (M) R1" in inp.udf
                                        ):
                                            if process_stats["Paired"]:
                                                target_file.udf["# Reads"] = (
                                                    inp.udf["Reads PF (M) R1"]
                                                    * 1000000
                                                    * 2
                                                )
                                                target_file.udf["# Read Pairs"] = (
                                                    target_file.udf["# Reads"] / 2
                                                )
                                            else:
                                                target_file.udf["# Reads"] = (
                                                    inp.udf["Reads PF (M) R1"] * 1000000
                                                )
                                                target_file.udf["# Read Pairs"] = (
                                                    target_file.udf["# Reads"]
                                                )
                                    logger.info(
                                        "{}# Reads".format(target_file.udf["# Reads"])
                                    )
                                    logger.info(
                                        "{}# Read Pairs".format(
                                            target_file.udf["# Read Pairs"]
                                        )
                                    )
                                except Exception as e:
                                    problem_handler(
                                        "exit",
                                        f"Unable to set values for #Reads and #Read Pairs for perceived noIndex lane: {str(e)}",
                                    )
                            # For all other cases, parse lane yield from all_inputs
                            else:
                                try:
                                    for inp in seq_process.all_inputs():
                                        # If reads in seq step, and the lane is equal to the current lane
                                        # Handle special case for MiSeq with noIndex case:
                                        inp_location = (
                                            "1"
                                            if inp.location[1][0] == "A"
                                            else inp.location[1][0]
                                        )
                                        if (
                                            inp_location == lane_no
                                            and "Clusters PF R1" in inp.udf
                                        ):
                                            if process_stats["Paired"]:
                                                target_file.udf["# Reads"] = (
                                                    inp.udf["Clusters PF R1"] * 2
                                                )
                                                target_file.udf["# Read Pairs"] = (
                                                    target_file.udf["# Reads"] / 2
                                                )
                                            else:
                                                target_file.udf["# Reads"] = inp.udf[
                                                    "Clusters PF R1"
                                                ]
                                                target_file.udf["# Read Pairs"] = (
                                                    target_file.udf["# Reads"]
                                                )
                                    logger.info(
                                        "{}# Reads".format(target_file.udf["# Reads"])
                                    )
                                    logger.info(
                                        "{}# Read Pairs".format(
                                            target_file.udf["# Read Pairs"]
                                        )
                                    )
                                except Exception as e:
                                    problem_handler(
                                        "exit",
                                        f"Unable to set values for #Reads and #Read Pairs for perceived noIndex lane: {str(e)}",
                                    )

                        elif not noIndex:
                            try:
                                clusterType = None
                                if "PF Clusters" in entry:
                                    clusterType = "PF Clusters"
                                else:
                                    clusterType = "Clusters"
                                # Paired runs are divided by two within flowcell parser
                                basenumber = int(entry[clusterType].replace(",", ""))
                                if process_stats["Paired"]:
                                    # Undet always 0 unless manually included
                                    samplesum[sample]["# Reads"] = (
                                        basenumber * 2 + undet_reads
                                        if "# Reads" not in samplesum[sample]
                                        else samplesum[sample]["# Reads"]
                                        + basenumber * 2
                                        + undet_reads
                                    )

                                    samplesum[sample]["# Read Pairs"] = (
                                        basenumber + undet_reads / 2
                                        if "# Read Pairs" not in samplesum[sample]
                                        else samplesum[sample]["# Read Pairs"]
                                        + basenumber
                                        + undet_reads / 2
                                    )
                                # Since a single ended run has no pairs, pairs is set to equal reads
                                else:
                                    # Undet always 0 unless manually included
                                    samplesum[sample]["# Reads"] = (
                                        basenumber + undet_reads
                                        if "# Reads" not in samplesum[sample]
                                        else samplesum[sample]["# Reads"]
                                        + basenumber
                                        + undet_reads
                                    )

                                    samplesum[sample]["# Read Pairs"] = (
                                        samplesum[sample]["# Reads"]
                                        if "# Read Pairs" not in samplesum[sample]
                                        else samplesum[sample]["# Read Pairs"]
                                        + samplesum[sample]["# Reads"]
                                    )
                            except Exception as e:
                                problem_handler(
                                    "exit",
                                    f"Unable to set values for #Reads and #Read Pairs: {str(e)}",
                                )

                        # Spools samplesum into samples
                        try:
                            if samplesum[sample]["count"] > 1:
                                logger.info("Iteratively pooling samples in same lane.")
                            for thing in samplesum:
                                for k, v in samplesum[thing].items():
                                    if thing == sample and thing == current_name:
                                        if k == "count":
                                            logger.info(
                                                f"Setting values for sample {thing} of lane {lane_no}"
                                            )
                                        # Average for percentages
                                        elif k in [
                                            "% One Mismatch Reads (Index)",
                                            "% Perfect Index Read",
                                            "Ave Q Score",
                                            "%PF",
                                            "% of Raw Clusters Per Lane",
                                            "% Bases >=Q30",
                                        ]:
                                            target_file.udf[k] = (
                                                v / samplesum[thing]["count"]
                                            )
                                        elif k != "count":
                                            target_file.udf[k] = samplesum[thing][k]
                                        if samplesum[sample]["count"] > 1:
                                            logger.info(
                                                f"Pooled total for {k} of sample {thing} is set to {v}"
                                            )
                                        else:
                                            logger.info(
                                                f"Attribute {k} of sample {thing} is set to {v}"
                                            )
                        except Exception as e:
                            problem_handler(
                                "exit",
                                f"Unable to set artifact values. Check laneBarcode.html for odd values: {str(e)}",
                            )

                        # Applies thresholds to samples
                        try:
                            if (
                                demux_process.udf["Threshold for % bases >= Q30"]
                                <= my_float(entry["% >= Q30bases"])
                                and int(exp_smp_per_lne)
                                <= target_file.udf["# Read Pairs"]
                            ):
                                target_file.udf["Include reads"] = "YES"
                                target_file.qc_flag = "PASSED"
                            else:
                                target_file.udf["Include reads"] = "NO"
                                target_file.qc_flag = "FAILED"
                                failed_entries = failed_entries + 1
                            logger.info(
                                "Q30 %: {}% found, minimum at {}%".format(
                                    my_float(entry["% >= Q30bases"]),
                                    demux_process.udf["Threshold for % bases >= Q30"],
                                )
                            )
                            logger.info(
                                "Expected reads: {} found, minimum at {}".format(
                                    target_file.udf["# Read Pairs"],
                                    int(exp_smp_per_lne),
                                )
                            )
                            logger.info(
                                f"Sample QC status set to {target_file.qc_flag}"
                            )
                        except Exception as e:
                            problem_handler(
                                "exit",
                                f"Unable to set QC status for sample: {str(e)}",
                            )

                        lane_reads = lane_reads + target_file.udf["# Reads"]

                    # Counts undetermined
                    elif sample == "Undetermined":
                        if "PF Clusters" in entry:
                            clusterType = "PF Clusters"
                        else:
                            clusterType = "Clusters"

                        if process_stats["Paired"]:
                            undet_lane_reads = (
                                int(entry[clusterType].replace(",", "")) * 2
                            )
                        else:
                            undet_lane_reads = int(entry[clusterType].replace(",", ""))

            if list(target_file.udf.items()) == [] and current_name != "Undetermined":
                problem_handler(
                    "exit",
                    f'Lanebarcode mismatch. Expected sample "{current_name}" of lane "{lane_no}", found "{sample}"',
                )

            # Push lane into lims
            try:
                target_file.put()
            except Exception as e:
                problem_handler(
                    "exit",
                    f"Failed to apply artifact data to LIMS. Possibly due to data in laneBarcode.html; {str(e)}",
                )

        # Counts undetermined per lane
        if not undet_included:
            try:
                found_undet = round(
                    my_float(undet_lane_reads) / (lane_reads + undet_lane_reads) * 100,
                    2,
                )

            # Only plausible error situation. Avoids zero division
            except Exception:
                problem_handler(
                    "error",
                    f"BCLConverter parsing error. No reads detected for lane {lane_no}.",
                )

            # If undetermined reads are greater than threshold*reads_in_lane
            if not noIndex:
                if (
                    found_undet
                    > demux_process.udf["Maximum % Undetermined Reads per Lane"]
                ):
                    problem_handler(
                        "warning",
                        f"Undemultiplexed reads for lane {lane_no} was {undet_lane_reads} ({found_undet})% thus exceeding defined limit.",
                    )
                else:
                    logger.info(
                        f"Found {undet_lane_reads} ({found_undet}%) undemultiplexed reads for lane {lane_no}."
                    )

    if undet_included:
        problem_handler("warning", "Undetermined reads included in read count!")

    if failed_entries > 0:
        problem_handler("warning", f"{failed_entries} entries failed automatic QC")


def write_demuxfile(process_stats, demux_id):
    """Creates demux_{FCID}.csv and attaches it to process"""
    # Includes windows drive letter support

    metadata_dir_name = "ngi-nas-ns"
    instrument_dir_name = "{}_data".format(process_stats["Instrument"])

    lanebc_path = os.path.join(
        os.sep,
        "srv",
        metadata_dir_name,
        instrument_dir_name,
        process_stats["Run ID"],
        "laneBarcode.html",
    )
    try:
        laneBC = classes.LaneBarcodeParser(lanebc_path)
    except Exception as e:
        problem_handler(
            "exit",
            f"Unable to fetch laneBarcode.html from {lanebc_path}: {str(e)}",
        )
    fname = "{}_demuxstats_{}.csv".format(demux_id, process_stats["Flow Cell ID"])

    # Writes less undetermined info than undemultiplex_index.py. May cause problems downstreams
    with open(fname, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Project",
                "Sample ID",
                "Lane",
                "# Reads",
                "Index",
                "Index name",
                "% of >= Q30 Bases (PF)",
            ]
        )
        for entry in laneBC.sample_data:
            index_name = ""
            if "PF Clusters" in entry:
                reads = entry["PF Clusters"]
            else:
                reads = entry["Clusters"]

            if process_stats["Paired"]:
                reads = int(reads.replace(",", "")) * 2
            else:
                reads = int(reads.replace(",", ""))

            try:
                writer.writerow(
                    [
                        entry["Project"],
                        entry["Sample"],
                        entry["Lane"],
                        reads,
                        entry["Barcode sequence"],
                        index_name,
                        entry["% >= Q30bases"],
                    ]
                )
            except Exception as e:
                problem_handler(
                    "exit",
                    f"Flowcell parser is unable to fetch all necessary fields for demux file: {str(e)}",
                )
    return laneBC.sample_data


def write_demuxfile_aviti(process_stats, demux_id):
    """Creates demux_{FCID}.csv and attaches it to process"""
    # Includes windows drive letter support

    metadata_dir_name = "ngi-nas-ns"
    instrument_dir_name = "{}_data".format(process_stats["Instrument"])

    lanebc_path = os.path.join(
        os.sep,
        "srv",
        metadata_dir_name,
        instrument_dir_name,
        process_stats["Run ID"],
        "IndexAssignment.csv",
    )

    try:
        laneBC = {}
        laneBC["sample_data"] = []
        with open(lanebc_path) as lanebc_file:
            reader = csv.DictReader(lanebc_file)
            for row in reader:
                if "+" not in row["Lane"]:
                    index = row.get("I1", "")
                    if row.get("I2"):
                        index += "-"
                        index += row["I2"]

                    laneBC["sample_data"].append(
                        {
                            "Lane": row.get("Lane", ""),
                            "Sample": row.get("SampleName", ""),
                            "Project": row.get("Project", ""),
                            "Barcode sequence": index,
                            "PF Clusters": str(row.get("NumPoloniesAssigned", "0")),
                            "% of thelane": float(
                                row.get("PercentPoloniesAssigned", "0")
                            ),
                            "% >= Q30bases": float(row.get("PercentQ30", "0") or 0),
                            "Mean QualityScore": float(
                                row.get("QualityScoreMean", "0") or 0
                            ),
                            "% Perfectbarcode": 100
                            - float(row.get("PercentMismatch", "0") or 0),
                            "% One mismatchbarcode": float(
                                row.get("PercentMismatch", "0") or 0
                            ),
                            "Yield (Mbases)": str(
                                float(row.get("Yield(Gb)", "0")) * 1000
                            ),
                        }
                    )
    except Exception as e:
        problem_handler(
            "exit",
            f"Unable to fetch IndexAssignment.csv from {lanebc_path}: {str(e)}",
        )

    fname = "{}_demuxstats_{}.csv".format(demux_id, process_stats["Flow Cell ID"])

    # Writes less undetermined info than undemultiplex_index.py. May cause problems downstreams
    with open(fname, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Project",
                "Sample ID",
                "Lane",
                "# Reads",
                "Index",
                "% of >= Q30 Bases (PF)",
                "Mean QualityScore",
                "% Perfectbarcode",
                "% One mismatchbarcode",
                "Yield (Mbases)",
            ]
        )
        for entry in laneBC["sample_data"]:
            reads = entry["PF Clusters"]

            if process_stats["Paired"]:
                reads = int(reads.replace(",", "")) * 2
            else:
                reads = int(reads.replace(",", ""))

            try:
                writer.writerow(
                    [
                        entry["Project"],
                        entry["Sample"],
                        entry["Lane"],
                        reads,
                        entry["Barcode sequence"],
                        entry["% >= Q30bases"],
                        entry["Mean QualityScore"],
                        entry["% Perfectbarcode"],
                        entry["% One mismatchbarcode"],
                        entry["Yield (Mbases)"],
                    ]
                )
            except Exception as e:
                problem_handler(
                    "exit",
                    f"Flowcell parser is unable to fetch all necessary fields for demux file: {str(e)}",
                )
    return laneBC["sample_data"]


def main(process_lims_id, demux_id, log_id):
    # Sets up logger
    basic_name = f"{log_id}_logfile.txt"
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(basic_name)
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.addHandler(fh)

    logger.info(
        f"--process_lims_id {process_lims_id} --demux_id {demux_id} --log_id {log_id}"
    )

    demux_process = Process(lims, id=process_lims_id)

    # Fetches info on "workflow" level
    process_stats = get_process_stats(demux_process)

    # Sets up the process values
    fill_process_fields(demux_process, process_stats)

    # Create the demux output file
    if "AVITI" in demux_process.type.name:
        parser_struct = write_demuxfile_aviti(process_stats, demux_id)
    else:
        parser_struct = write_demuxfile(process_stats, demux_id)

    # Alters artifacts
    set_sample_values(demux_process, parser_struct, process_stats)

    # Changing log file name, can't do this step earlier since proc_stats is made during runtime.
    new_name = "{}_logfile_{}.txt".format(log_id, process_stats["Flow Cell ID"])
    move(basic_name, new_name)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument(
        "--process_lims_id",
        required=True,
        dest="process_lims_id",
        help="Lims ID of process. Example:24-92373",
    )
    parser.add_argument(
        "--demux_id",
        required=True,
        dest="demux_id",
        help=("Id prefix for demux output."),
    )
    parser.add_argument(
        "--log_id", required=True, dest="log_id", help=("Id prefix for logfile")
    )
    args = parser.parse_args()
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(args.process_lims_id, args.demux_id, args.log_id)
