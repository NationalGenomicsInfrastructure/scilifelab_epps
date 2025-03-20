#!/usr/bin/env python

import logging
import os
import re
from argparse import ArgumentParser, Namespace
from datetime import datetime as dt

from couchdb.client import Database, Document, Row, ViewResults
from generate_minknow_samplesheet import (
    generate_MinKNOW_samplesheet,
    get_ont_library_contents,
    sanitize_string,
)
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Artifact, Process
from genologics.lims import Lims
from ont_send_reloading_info_to_db import get_ONT_db

from scilifelab_epps.utils import udf_tools
from scilifelab_epps.wrapper import epp_decorator

DESC = """Script for finishing the step to start ONT sequencing in LIMS.

Makes sure there are no discrepancies in the provided information and
appends LIMS-specific information to the ONT run in StatusDB.
"""

TIMESTAMP: str = dt.now().strftime("%y%m%d_%H%M%S")


def assert_samplesheet(process: Process, args: Namespace, lims: Lims):
    """Check that there isn't a loaded samplesheet that is contradicted by the step UDFs."""

    # Get ONT samplesheet file artifact
    file_art: Artifact = [
        op for op in process.all_outputs() if op.name == args.samplesheet
    ][0]

    # If samplesheet file is loaded
    if file_art.files:
        logging.info("Detected samplesheet.")
        samplesheet_contents: str = lims.get_file_contents(uri=file_art.files[0].uri)

        logging.info("Checking that the loaded samplesheet is up to date...")
    else:
        logging.info("No samplesheet file loaded.")
        return True

    # Generate new samplesheet from step, then read it and remove the file
    new_samplesheet_path = generate_MinKNOW_samplesheet(process)
    new_samplesheet_contents = open(new_samplesheet_path).read()
    os.remove(new_samplesheet_path)

    # Check step samplesheet is up-to-date with step UDFs
    if samplesheet_contents != new_samplesheet_contents:
        logging.error(
            f"The current sample sheet doesn't correspond to the current UDFs.\nCurrent sample sheet:\n{samplesheet_contents}\nNew sample sheet:\n{new_samplesheet_contents}"
        )
        raise AssertionError(
            "The current sample sheet doesn't correspond to the current UDFs."
        )
    else:
        logging.info("Samplesheet is up to date.")
        return True


def udfs_matches_run_name(art: Artifact) -> bool:
    """Check that artifact run name is not contradicted by other UDFs."""

    matches: bool = True
    _yyyymmdd, _hhmm, pos, fc_id, _hash = art.udf["ONT run name"].split("_")

    if art.udf["ONT flow cell ID"] != fc_id:
        matches = False
        msg = f"Mismatch between UDFs 'ONT flow cell ID': '{art.udf['ONT flow cell ID']}' and 'ONT run name': '{art.udf['ONT run name']}'"
        logging.error(msg)
    if (
        art.udf["ONT flow cell position"] != "None"
        and art.udf["ONT flow cell position"] != pos
    ):
        matches = False
        msg = f"Mismatch between UDFs 'ONT flow cell position': '{art.udf['ONT flow cell position']}' and 'ONT run name': '{art.udf['ONT run name']}'"
        logging.error(msg)

    return matches


def get_matching_db_rows(
    art: Artifact,
    process: Process,
    view: ViewResults,
    run_name: str | None,
) -> list[Row]:
    """Find the rows in the database view that match the given artifact."""
    matching_rows = []

    qc = True if "QC" in process.type.name else False
    logging.info(f"QC run: {qc}")

    # If run name is not supplied, try to find it in the database, assuming it follows the samplesheet naming convention
    if run_name is None:
        # Define query pattern
        experiment_name = process.id if not qc else f"QC_{process.id}"
        sample_name = (
            sanitize_string(art.name) if not qc else f"QC_{sanitize_string(art.name)}"
        )
        position = (
            art.udf["ONT flow cell position"]
            if art.udf["ONT flow cell position"] != "None"
            else "MN19414"
        )

        pattern = rf"{experiment_name}/{sample_name}/[^/]*_{position}_{art.udf['ONT flow cell ID']}_[^/]*"

        logging.info(
            f"No run name supplied. Quering the database for run with path pattern '{pattern}'."
        )

        for row in view.rows:
            query = row.value["TACA_run_path"]
            if re.match(pattern, query):
                matching_rows.append(row)

    # If the run name is supplied, query the database directly
    else:
        logging.info(
            f"Full run name supplied. Quering the database for run '{run_name}'."
        )
        for row in view.rows:
            if run_name == row.key:
                matching_rows.append(row)

    return matching_rows


def write_to_doc(
    doc: Document, db: Database, process: Process, art: Artifact, args: Namespace
):
    """Update a given document with the given artifact's loading information."""

    library_df = get_ont_library_contents(
        ont_library=art,
        list_contents=True,
        print_dataframe=True,
    )

    # Info to add to the db doc
    dict_to_add = {
        "step_name": process.type.name,
        "step_id": process.id,
        "timestamp": TIMESTAMP,
        "operator": process.technician.name,
        "load_fmol": art.udf["ONT flow cell loading amount (fmol)"],
        "load_vol": art.udf["Volume to take (uL)"],
        "library_name": art.name,
        "library_id": art.id,
        "sample_data": library_df.to_dict(orient="records"),
    }

    if "lims" not in doc:
        doc["lims"] = {}
    if "loading" not in doc["lims"]:
        doc["lims"]["loading"] = []
    doc["lims"]["loading"].append(dict_to_add)

    db[doc.id] = doc


def sync_runs_to_db(process: Process, args: Namespace, lims: Lims):
    """Executive script, called once."""

    # Assert samplesheet, if any, makes sense in the context of the step UDFs
    assert_samplesheet(process=process, args=args, lims=lims)

    arts: list[Artifact] = [
        art for art in process.all_outputs() if art.type == "Analyte"
    ]

    # Assert that only one input is provided for QC runs
    if "QC" in process.type.name:
        assert len(arts) == 1, (
            "When starting QC sequencing runs, only one input is allowed."
        )

    # Keep track of which artifacts were successfully updated
    arts_successful = []

    nanopore_runs_db: Database = get_ONT_db()
    view: ViewResults = nanopore_runs_db.view("info/all_stats")

    for art in arts:
        logging.info(f"Processing '{art.name}'...")

        run_name: str | None = udf_tools.fetch(art, "ONT run name", on_fail=None)
        if run_name is not None:
            # Assert run name is not contradicted by other UDFs
            if not udfs_matches_run_name(art):
                logging.warning("Run name contradicted by other UDFs. Skipping.")
                continue

        # Get matching run docs
        matching_rows: list[Row] = get_matching_db_rows(art, process, view, run_name)

        if len(matching_rows) == 0:
            logging.warning("Run was not found in the database. Skipping.")
            continue

        elif len(matching_rows) > 1:
            matching_run_names = [row.key for row in matching_rows]
            logging.warning("Query was found in multiple instances in the database: ")
            for matching_run_name in matching_run_names:
                logging.warning(f"Matching run name: '{matching_run_name}'.")
            logging.warning("Contact a database administrator. Skipping.")
            continue

        doc_run_name: str = matching_rows[0].key
        doc_id: str = matching_rows[0].id
        doc: Document = nanopore_runs_db[doc_id]

        logging.info(f"Found matching run '{doc_run_name}' in the database.")

        if run_name is not None and run_name != doc_run_name:
            logging.error(
                f"UDF Run name '{run_name}' contradicted by database run name '{doc_run_name}'. Skipping."
            )
            continue

        logging.info(f"Assigning UDF 'ONT run name': '{doc_run_name}'.")
        udf_tools.put(art, "ONT run name", doc_run_name)

        write_to_doc(doc, nanopore_runs_db, process, art, args)
        logging.info(f"'{doc_run_name}' was found and updated successfully.")
        arts_successful.append(art)

    if len(arts_successful) < len(arts):
        raise AssertionError(
            f"Only {len(arts_successful)} out of {len(arts)} artifacts were successfully updated. Check log."
        )


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    # Set up LIMS
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    process = Process(lims, id=args.pid)

    sync_runs_to_db(process=process, lims=lims, args=args)


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument(
        "--pid",
        type=str,
        help="Lims ID for current Process",
    )
    parser.add_argument(
        "--log",
        required=True,
        type=str,
        help="Which log file slot to use",
    )
    parser.add_argument(
        "--samplesheet",
        required=True,
        type=str,
        help="Which samplesheet file slot to use",
    )
    args: Namespace = parser.parse_args()

    main(args)
