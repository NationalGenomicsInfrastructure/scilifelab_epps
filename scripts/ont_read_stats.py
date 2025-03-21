#!/usr/bin/env python

import logging
from argparse import ArgumentParser
from datetime import datetime as dt

from couchdb.client import Database, Row, ViewResults
from generate_minknow_samplesheet import get_pool_sample_label_mapping
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Artifact, Process, Sample
from genologics.lims import Lims
from ont_send_reloading_info_to_db import get_ONT_db

from data.ONT_barcodes import ont_label2dict
from scilifelab_epps.wrapper import epp_decorator

DESC = """Assign metrics from ONT run onto samples in a LIMS demultiplexing step.
"""

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    # Connect to database
    db: Database = get_ONT_db()
    barcodes_view: ViewResults = db.view("info/barcodes")
    stats_view: ViewResults = db.view("info/all_stats")

    for library in process.all_inputs():
        # Use ONT run name to navigate database
        run_name: str = library.udf.get("ONT run name")
        logging.info(
            f"Processing library '{library.name}' ({library.id}) of run {run_name}."
        )

        # Get sample-barcode linkage
        is_barcoded: bool = len(library.samples) > 1
        sample2label: dict[str, str] = (
            get_pool_sample_label_mapping(library) if is_barcoded else None
        )

        # For both views, get the row corresponding to the current run
        stats_row: Row = [row for row in stats_view.rows if run_name == row.key][0]
        if is_barcoded:
            barcodes_row: Row = [
                row for row in barcodes_view.rows if run_name == row.key
            ][0]

        # Get the demultiplexing artifacts of the current run
        demux_arts: list[Artifact] = [
            io[1]["uri"]
            for io in process.input_output_maps
            if io[0]["uri"].id == library.id
        ]
        logging.info(f"Handling {len(demux_arts)} demultiplexing artifacts.")

        for demux_art in demux_arts:
            sample: Sample = demux_art.samples[0]
            logging.info(
                f"Processing demultiplexing artifact '{demux_art.name}' ({demux_art.id})."
            )

            # Get dict containing sequencing run metrics
            if is_barcoded:
                label: str = sample2label[sample.name]
                barcode_info: dict = ont_label2dict[label]
                barcode_num: int = barcode_info["num"]
                barcode_generic_name: str = f"barcode{str(barcode_num).zfill(2)}"
                logging.info(f"Using LIMS label '{label}' as '{barcode_generic_name}'.")

                metrics: dict | None = barcodes_row.value.get(barcode_generic_name)
                if not metrics:
                    logging.warning(f"Barcode {label} not found in database. Skipping.")
                    continue
            else:
                metrics: dict = stats_row.value

            # Parse and calculate values
            passed_reads = int(metrics.get("basecalled_pass_read_count", "0"))
            total_reads = int(metrics.get("read_count", "0"))
            passed_gb = float(metrics.get("basecalled_pass_bases", "0")) / 1e9
            passfail_pc: float = (passed_reads / total_reads) * 100
            avg_len: int = round(passed_gb / passed_reads * 1e9)

            # Assign UDFs
            demux_art.udf["# Reads"] = passed_reads
            demux_art.udf["Yield PF (Gb)"] = passed_gb
            demux_art.udf["%PF"] = passfail_pc
            demux_art.udf["Avg. Read Length"] = avg_len

            # Publish metrics
            demux_art.put()


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument(
        "--pid",
        required=True,
        type=str,
        help="Lims ID for current Process",
    )
    parser.add_argument(
        "--log",
        required=True,
        type=str,
        help="Which log file slot to use",
    )
    args = parser.parse_args()

    main(args)
