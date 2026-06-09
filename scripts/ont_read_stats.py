#!/usr/bin/env python

import logging
from argparse import ArgumentParser
from datetime import datetime as dt

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Artifact, Process
from genologics.lims import Lims
from ibmcloudant import cloudant_v1
from ont_send_reloading_info_to_db import get_ONT_db

from data.ONT_barcodes import get_barcode_info
from scilifelab_epps.wrapper import epp_decorator

DESC = """Assign metrics from ONT run onto samples in a LIMS demultiplexing step.
"""

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    # Connect to database
    client: cloudant_v1.CloudantV1
    client, db_name = get_ONT_db()

    for library in process.all_inputs():
        # Use ONT run name to navigate database
        run_name: str = library.udf.get("ONT run name")
        logging.info(
            f"Processing library '{library.name}' ({library.id}) of run {run_name}."
        )

        # Get sample-barcode linkage
        is_barcoded: bool = len(library.samples) > 1

        # For both views, get the row corresponding to the current run
        stats_rows: list[dict] = client.post_view(
            db=db_name,
            ddoc="info",
            view="all_stats",
            include_docs=False,
            key=run_name,
        ).get_result()["rows"]
        if not stats_rows:
            logging.error(
                f"No stats entry found in database for run '{run_name}'. Skipping."
            )
            continue
        stats_row: dict = stats_rows[0]

        if is_barcoded:
            logging.info("Library seems barcoded.")
            barcodes_rows: list[dict] = client.post_view(
                db=db_name,
                ddoc="info",
                view="barcodes",
                include_docs=False,
                key=run_name,
            ).get_result()["rows"]
            if not barcodes_rows:
                logging.error(
                    "Run database entry does not appear to contain barcodes. Skipping."
                )
                continue
            barcodes_row: dict = barcodes_rows[0]
            if not barcodes_row.get("value"):
                logging.error(
                    "Run database entry does not appear to contain barcodes. Skipping."
                )
                continue
        else:
            logging.info("Library seems unlabeled.")

        # Get the demultiplexing artifacts of the current run
        demux_arts: list[Artifact] = [
            io[1]["uri"]
            for io in process.input_output_maps
            if io[0]["uri"].id == library.id
            and io[1]["output-generation-type"] == "PerReagentLabel"
        ]
        logging.info(f"Handling {len(demux_arts)} demultiplexing artifacts.")

        for demux_art in demux_arts:
            logging.info(
                f"Processing demultiplexing artifact '{demux_art.name}' ({demux_art.id})."
            )

            # Get dict containing sequencing run metrics
            if is_barcoded:
                label: str = demux_art.reagent_labels[0]
                barcode_info: dict = get_barcode_info(label)
                barcode_num: int = barcode_info["num"]
                barcode_generic_name: str = f"barcode{str(barcode_num).zfill(2)}"
                logging.info(f"Using LIMS label '{label}' as '{barcode_generic_name}'.")

                if not barcodes_row["value"].get(barcode_generic_name):
                    logging.error(
                        f"Run database entry does not appear to contain data for {barcode_generic_name}. Skipping."
                    )
                    continue
                metrics: dict = barcodes_row["value"][barcode_generic_name]

            else:
                metrics: dict = stats_row["value"]

            # Parse and calculate values
            reads_pass = int(metrics.get("basecalled_pass_read_count", "0"))
            reads_tot = int(metrics.get("read_count", "0"))
            gb_pass = float(metrics.get("basecalled_pass_bases", "0")) / 1e9
            pf_pc: float = (reads_pass / reads_tot) * 100 if reads_tot != 0 else 0
            avg_len: int = round(gb_pass / reads_pass * 1e9) if reads_pass != 0 else 0

            # Assign UDFs
            demux_art.udf["# Reads"] = reads_pass
            demux_art.udf["Yield PF (Gb)"] = gb_pass
            demux_art.udf["%PF"] = pf_pc
            demux_art.udf["Avg. Read Length"] = avg_len
            demux_art.udf["Include reads"] = "YES"

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
