#!/usr/bin/env python

import logging
import os
import re
import shutil
from argparse import ArgumentParser
from datetime import datetime as dt

from generate_minknow_samplesheet import get_ont_library_contents
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from data.Chromium_10X_indexes import Chromium_10X_indexes
from data.ONT_barcodes import ont_label2dict
from scilifelab_epps.epp import upload_file
from scilifelab_epps.wrapper import epp_decorator

DESC = """Script to generate Anglerfish samplesheet for ONT runs.
"""

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")


def generate_anglerfish_samplesheet(process):
    """Generate an Anglerfish samplesheet.

    The samplesheet is a headerless .csv-file in which the columns correspond to:
    'sample_name', 'adaptors', 'index', 'fastq_path'
    """

    ont_libraries = [art for art in process.all_outputs() if art.type == "Analyte"]
    assert len(ont_libraries) == 1, (
        "Samplesheet can only be generated for a single sequencing library."
    )
    ont_library = ont_libraries[0]

    df = get_ont_library_contents(
        ont_library=ont_library,
        print_dataframe=True,
        list_contents=True,
    )

    # Add columns pertaining to barcode properties
    if "ont_barcode" in df.columns:
        for i in ["num", "well", "seq"]:
            df[f"ont_barcode_{i}"] = df["ont_barcode"].apply(
                lambda barcode_label: ont_label2dict[barcode_label][i]
            )

        df["fastq_path"] = df["ont_barcode_num"].apply(
            lambda num: f"./fastq_pass/barcode{str(num).zfill(2)}/*.fastq.gz"
        )
    else:
        df["fastq_path"] = "./fastq_pass/*.fastq.gz"

    # Extract index sequence and adaptor type
    df["index_seq"] = df["illumina_index"].apply(lambda x: extract_sequence(x))
    df["adaptor_type"] = df["illumina_index"].apply(lambda x: get_adaptor_name(x))

    # Subset columns
    df_anglerfish = df[["sample_name", "adaptor_type", "index_seq", "fastq_path"]]

    # Get run name from LIMS field
    run_name = ont_library.udf.get("ONT run name")

    file_name = f"anglerfish_samplesheet_{run_name}_{TIMESTAMP}.csv"
    df_anglerfish.to_csv(
        file_name,
        header=False,
        index=False,
    )

    return file_name


def extract_sequence(reagent_label: str) -> str | None:
    """Extract sequence from string."""

    index_pattern = re.compile("([ACTG]{4,})-?([ACTG]{4,})?")
    index_search = re.search(index_pattern, reagent_label)

    if index_search:
        return index_search.group()
    else:
        return None


def get_adaptor_name(reagent_label: str) -> str | list[str]:
    """Derive adaptor name from reagent label."""

    seq = extract_sequence(reagent_label)

    if seq:
        if "-" in seq:
            return "truseq_dual"
        else:
            return "truseq"

    elif reagent_label in Chromium_10X_indexes.keys():
        matching_10x_indices = Chromium_10X_indexes[reagent_label]

        if len(matching_10x_indices) == 2:
            # Return i7-i5
            return "-".join(matching_10x_indices)

        elif len(matching_10x_indices) == 4:
            # Return list of combination i7 indices
            return matching_10x_indices

        else:
            raise AssertionError(
                f"Could not determine adaptor of reagent label {reagent_label}"
            )

    else:
        raise AssertionError(
            f"Could not determine adaptor of reagent label {reagent_label}"
        )


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    file_name = generate_anglerfish_samplesheet(process)

    logging.info("Uploading samplesheet to LIMS...")
    upload_file(
        file_name,
        args.file,
        process,
        lims,
    )

    logging.info("Moving samplesheet to ngi-nas-ns...")
    try:
        shutil.copyfile(
            file_name,
            f"/srv/ngi-nas-ns/samplesheets/anglerfish/{dt.now().year}/{file_name}",
        )
        os.remove(file_name)
    except:
        logging.error("Failed to move samplesheet to ngi-nas-ns.")
    else:
        logging.info("Samplesheet moved to ngi-nas-ns.")


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument(
        "--pid",
        required=True,
        type=str,
        help="Lims ID for current Process.",
    )
    parser.add_argument(
        "--log",
        required=True,
        type=str,
        help="Which file slot to use for the script log.",
    )
    parser.add_argument(
        "--file",
        required=True,
        type=str,
        help="Which file slot to use for the samplesheet.",
    )
    args = parser.parse_args()

    main()
