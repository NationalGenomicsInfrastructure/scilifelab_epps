#!/usr/bin/env python
DESC = """EPP script to aggregate the number of reads from different demultiplexing runs,
based on the flag 'include reads' located at the same level as '# reads'

Denis Moreno, Science for Life Laboratory, Stockholm, Sweden
"""
import logging
import os
from argparse import ArgumentParser
from datetime import datetime as dt

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from scilifelab_epps.epp import attach_file
from scilifelab_epps.wrapper import epp_decorator

TIMESTAMP: str = dt.now().strftime("%y%m%d_%H%M%S")

# Master step IDs and names
DEMULTIPLEX = {
    "13": "Bcl Conversion & Demultiplexing (Illumina SBS) 4.0",
    "3205": "ONT Finish Sequencing v3",  # TODO Update this to reflect prod
    # TODO Add AVITI
}
SUMMARY = {
    "356": "Project Summary 1.3",
}
SEQUENCING = {
    "38": "Illumina Sequencing (Illumina SBS) 4.0",
    "46": "MiSeq Run (MiSeq) 4.0",
    "714": "Illumina Sequencing (HiSeq X) 1.0",
    "1454": "AUTOMATED - NovaSeq Run (NovaSeq 6000 v2.0)",
    "1908": "Illumina Sequencing (NextSeq) v1.0",
    "2612": "NovaSeqXPlus Run v1.0",
    "2955": "ONT Start Sequencing v3.0",  # TODO Update this to reflect prod
}


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    """This should be run at project summary level"""
    process = Process(lims, id=args.pid)

    sample_counter = 0
    error_counter = 0

    summary = {}  # { sample_name : { flowcell : { lane1, lane2, ... } } }  # dict -> dict -> set
    log_artifact = [
        art
        for art in process.all_outputs()
        if art.type == "ResultFile" and art.name == "AggregationSummary"
    ][0]

    # Iterate across output analytes
    arts_out = [art for art in process.all_outputs() if art.type == "Analyte"]
    for art_out in arts_out:
        assert (
            len(art_out.samples) == 1
        ), f"Found {len(art_out.samples)} samples for the output analyte {art_out.id}, that should not happen"

        sample = art_out.samples[0]
        sample_counter += 1

        # Traceback and calculate the total number of reads
        total_reads = sum_reads(sample, summary)

        # Set total reads for sample and artifact UDFs
        sample.udf["Total Reads (M)"] = total_reads
        art_out.udf["Set Total Reads"] = total_reads
        logging.info(f"Total reads is {total_reads} for sample '{sample.name}'")

        # Set min reads sample UDF from project UDF
        min_reads = sample.project.udf.get("Reads Min", 0) / 1e6
        logging.info(f"Updating '{sample.name}' UDF 'Reads Min' to {min_reads}")
        sample.udf["Reads Min"] = min_reads

        # Set sample UDFs for status and sequencing QC based on min reads and total reads
        if total_reads <= min_reads:
            logging.info(
                "Total reads is below minimum, setting status to 'In Progress' and 'Passed Sequencing QC' to 'False'"
            )
            sample.udf["Status (auto)"] = "In Progress"
            sample.udf["Passed Sequencing QC"] = "False"
        elif total_reads > min_reads:
            logging.info(
                "Total reads is above minimum, setting status to 'Finished' and 'Passed Sequencing QC' to 'True'"
            )
            sample.udf["Passed Sequencing QC"] = "True"
            sample.udf["Status (auto)"] = "Finished"

        # Commit changes to sample and sample artifact
        sample.put()
        art_out.put()

    # Write the csv file, separated by pipes, no cell delimiter
    with open("AggregationSummary.csv", "w") as f:
        f.write("sep=,\n")
        f.write(
            "sample name,number of flowcells,number of lanes,flowcell1:lane1|lane2;flowcell2:lane1|lane2|lane3 ...\n"
        )
        for sample in summary:
            view = []
            n_flowcells = len(summary[sample])
            n_lanes = 0
            for fc in summary[sample]:
                view.append("{}:{}".format(fc, "|".join(summary[sample][fc])))
                n_lanes += len(summary[sample][fc])
            f.write(
                "{},{},{},{}\n".format(sample, n_flowcells, n_lanes, ";".join(view))
            )

    # Upload the log file
    try:
        attach_file(os.path.join(os.getcwd(), "AggregationSummary.csv"), log_artifact)
        logging.info(f"Updated {sample_counter} samples with {error_counter} errors")
    except AttributeError:
        # Happens if the log artifact does not exist, if the step has been started before the configuration changes
        logging.info("Could not upload the log file")


def sum_reads(sample, summary):
    """For a given submitted sample and it's summary,
    calculate the total number of reads and append to the summary.
    """

    logging.info(f"Aggregating reads of sample '{sample.name}'...")

    # Append to summary
    if sample.name not in summary:
        summary[sample.name] = {}

    expected_name = f"{sample.name} (FASTQ reads)"
    # Look for artifacts matching the sample name and expected analyte name in the demultiplexing processeses
    demux_arts = lims.get_artifacts(
        sample_name=sample.name,
        process_type=list(DEMULTIPLEX.values()),
        name=expected_name,
    )

    # Iterate across found demux artifacts to aggregate reads and collect flowcell information
    tot_reads = 0
    flowcell_lane_list = []
    for demux_art in sorted(
        demux_arts, key=lambda art: art.parent_process.date_run, reverse=True
    ):
        logging.info(
            f"Looking at '{demux_art.name}' ({demux_art.id}) of step '{demux_art.parent_process.type.name}' ({demux_art.parent_process.id})..."
        )

        # Evaluate skip conditions
        if "# Reads" not in demux_art.udf:
            logging.warning("Missing or unpopulated UDF '# Reads', skipping.")
            continue

        if "Include reads" not in demux_art.udf:
            logging.warning(
                "Missing or unpopulated UDF 'Include_reads' filled, skipping."
            )
            continue

        if demux_art.udf["Include reads"] == "NO":
            logging.info("UDF 'Include reads' is set to 'NO', skipping.")
            continue

        assert demux_art.udf["Include reads"] == "YES"

        # From the demux artifact, find the parent analyte from the actual sequencing step
        demux_art_parents = [
            parent
            for parent in get_parent_inputs(demux_art)
            if sample in parent.samples
        ]
        assert len(demux_art_parents) == 1
        demux_art_parent = demux_art_parents[0]

        # Check whether we are dealing with dual reads
        dual_reads = False
        try:
            seq_process = lims.get_processes(
                type=list(SEQUENCING.values()),
                inputartifactlimsid=demux_art_parent.id,
            )[0]
        except TypeError:
            logging.error(
                f"Did not manage to get sequencing process for artifact '{demux_art_parent.name}' ({demux_art_parent.id})"
            )
        else:
            if (
                "Read 2 Cycles" in seq_process.udf
                and seq_process.udf["Read 2 Cycles"] is not None
            ):
                dual_reads = True

        # Gather flowcell information
        if "ONT flow cell ID" in demux_art_parent.udf:
            # ONT
            ont_flowcell = demux_art_parent.udf["ONT flow cell ID"]
            if ont_flowcell not in flowcell_lane_list:
                flowcell_lane_list.append(ont_flowcell)
            if ont_flowcell not in summary[sample.name]:
                summary[sample.name][ont_flowcell] = set()
        else:
            # Illumina
            flowcell_and_lane = "{}:{}".format(
                demux_art_parent.location[0].name,
                demux_art_parent.location[1].split(":")[0],
            )
            if flowcell_and_lane not in flowcell_lane_list:
                flowcell_lane_list.append(flowcell_and_lane)
            if demux_art_parent.location[0].name in summary[sample.name]:
                summary[sample.name][demux_art_parent.location[0].name].add(
                    demux_art_parent.location[1].split(":")[0]
                )
            else:
                summary[sample.name][demux_art_parent.location[0].name] = set(
                    demux_art_parent.location[1].split(":")[0]
                )

        # Aggregate reads to total
        if dual_reads:
            tot_reads += float(demux_art.udf["# Reads"]) / 2
        else:
            tot_reads += float(demux_art.udf["# Reads"])

    # Total is displayed as millions
    tot_reads_m = tot_reads / 1e6
    return tot_reads_m


def get_parent_inputs(art):
    input_arts = set()
    for input_output_tuple in art.parent_process.input_output_maps:
        if input_output_tuple[1]["uri"].id == art.id:
            input_arts.add(input_output_tuple[0]["uri"])

    return input_arts


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", type=str, help="Lims ID for current Process")
    parser.add_argument("--log", type=str, help="Which log file slot to use")

    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)

    main(args)
