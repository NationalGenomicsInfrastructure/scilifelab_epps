#!/usr/bin/env python

import logging
import os
import re
import shutil
import sys
from argparse import ArgumentParser, Namespace
from datetime import datetime as dt

import pandas as pd
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Artifact, Process
from genologics.lims import Lims
from tabulate import tabulate

from data.ONT_barcodes import ONT_BARCODE_LABEL_PATTERN, ONT_BARCODES
from epp_utils.udf_tools import fetch
from scilifelab_epps.epp import traceback_to_step, upload_file

DESC = """ Script to generate MinKNOW samplesheet for starting ONT runs.
"""

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")
SCRIPT_NAME: str = os.path.basename(__file__).split(".")[0]


def get_ont_library_contents(
    ont_library: Artifact,
    ont_pooling_step_name: str,
    list_contents: bool = False,
    print_dataframe: bool = False,
) -> pd.DataFrame:
    """For an ONT sequencing library, compile a dataframe with sample-level information.

    Will backtrack the library to previous pooling step (if any) to elucidate
    sample and index information and decide whether to demultiplex at the level of
    ONT barcodes, Illumina indices, both or neither.

    The logic is a bit messy to accommodate both the QC and non-QC workflows.

    """

    ont_barcode_well2label: dict[str:str] = {
        ont_barcode["well"]: ont_barcode["label"] for ont_barcode in ONT_BARCODES
    }

    logging.info(
        f"Compiling sample-level information for library '{ont_library.name}'..."
    )
    pool_contents_msg = f"ONT sequencing library '{ont_library.name}' consists of:"

    # Instantiate list to collect dataframe rows
    rows = []

    # Remaining possibilities:
    # (1) ONT-barcodes and Illumina indexes
    # (2) ONT-barcodes only
    # (3) Illumina-indexes only
    # (4) No labels

    # See if library can be backtracked to an ONT pooling step
    ont_pooling_traceback = traceback_to_step(
        ont_library, ont_pooling_step_name, allow_multiple_inputs=True
    )
    if ont_pooling_traceback is not None:
        # Remaining possibilities:
        # (1) ONT-barcodes and Illumina indexes
        # (2) ONT-barcodes only

        _ont_pooling_step, ont_pooling_inputs, ont_pooling_output = (
            ont_pooling_traceback
        )

        pool_contents_msg += f"\n - '{ont_pooling_output.name}': ONT-barcoded pool"

        # Iterate across ONT pooling inputs
        for ont_pooling_input in ont_pooling_inputs:
            if len(ont_pooling_input.samples) > 1:
                # Remaining possibilities:
                # (1) ONT-barcodes and Illumina indexes

                # ONT barcodes on Illumina libraries are assigned via UDFs rather than reagent labels, since LIMS can't handle double demultiplexing
                udf_ont_barcode_well = fetch(
                    ont_pooling_input, "ONT Barcode Well", on_fail=None
                )
                assert udf_ont_barcode_well, f"Pooling input '{ont_pooling_input.name}' consists of multiple samples, but has not been assigned an ONT barcode."
                ont_barcode = ont_barcode_well2label[
                    udf_ont_barcode_well.upper().replace(":", "")
                ]

                pool_contents_msg += f"\n\t - '{ont_pooling_input.name}': Illumina indexed pool with ONT-barcode '{ont_barcode}'"

                for sample, illumina_index in zip(
                    ont_pooling_input.samples, ont_pooling_input.reagent_labels
                ):
                    pool_contents_msg += f"\n\t\t - '{sample.name}': Illumina sample with index '{illumina_index}'."
                    rows.append(
                        {
                            "sample_name": sample.name,
                            "sample_id": sample.id,
                            "project_name": sample.project.name,
                            "project_id": sample.project.id,
                            "illumina_index": illumina_index,
                            "illumina_pool_name": ont_pooling_input.name,
                            "illumina_pool_id": ont_pooling_input.id,
                            "ont_barcode": ont_barcode,
                            "ont_pool_name": ont_pooling_output.name,
                            "ont_pool_id": ont_pooling_output.id,
                        }
                    )

            elif len(ont_pooling_input.samples) == 1:
                # Remaining possibilities:
                # (2) ONT-barcodes only
                assert (
                    len(ont_pooling_input.reagent_labels) == 1
                ), f"ONT-pooling input '{ont_pooling_input.name}' lacks any reagent labels."

                # ONT barcode-level demultiplexing
                for ont_sample, ont_barcode in zip(
                    ont_pooling_input.samples, ont_pooling_input.reagent_labels
                ):
                    pool_contents_msg += f"\n\t - '{ont_pooling_input.name}: ONT sample with barcode '{ont_barcode}'"
                    rows.append(
                        {
                            "sample_name": ont_sample.name,
                            "sample_id": ont_sample.id,
                            "project_name": ont_sample.project.name,
                            "project_id": ont_sample.project.id,
                            "ont_barcode": ont_barcode,
                            "ont_pool_name": ont_pooling_output.name,
                            "ont_pool_id": ont_pooling_output.id,
                        }
                    )

            else:
                raise AssertionError(
                    f"ONT-pooling input '{ont_pooling_input.name}' lacks any samples."
                )

    else:
        # Remaining possibilities:
        # (3) Illumina-indexes only
        # (4) No labels

        if len(ont_library.reagent_labels) > 0:
            # Remaining possibilities:
            # (3) Illumina-indexes only
            for sample, illumina_index in zip(
                ont_library.samples, ont_library.reagent_labels
            ):
                pool_contents_msg += f"\n - '{sample.name}': Illumina sample with index '{illumina_index}'."
                rows.append(
                    {
                        "sample_name": sample.name,
                        "sample_id": sample.id,
                        "project_name": sample.project.name,
                        "project_id": sample.project.id,
                        "illumina_index": illumina_index,
                        "illumina_pool_name": ont_library.name,
                        "illumina_pool_id": ont_library.id,
                    }
                )

        else:
            # Remaining possibilities:
            # (4) No labels
            sample = ont_library.samples[0]

            pool_contents_msg += f"\n - {sample.name}: Non-labeled sample"
            rows.append(
                {
                    "sample_name": sample.name,
                    "sample_id": sample.id,
                    "project_name": sample.project.name,
                    "project_id": sample.project.id,
                }
            )

    if list_contents:
        logging.info(pool_contents_msg)

    df = pd.DataFrame(rows)
    table_str = tabulate(df, headers=df.columns)
    indented_table_str = "\n".join(["\t" + line for line in table_str.split("\n")])
    if print_dataframe:
        logging.info(
            f"Sample-level information compiled for library '{ont_library.name}':\n{indented_table_str}"
        )

    return df


def get_kit_string(process: Process) -> str:
    """Combine prep kit and expansion kit UDFs (if any) into space-separated string"""
    prep_kit = process.udf.get("ONT prep kit")
    expansion_kit = process.udf.get("ONT expansion kit")

    if expansion_kit != "None":
        prep_kit += f" {expansion_kit.replace('.','-')}"

    return prep_kit


def sanitize_string(string: str) -> str:
    """Remove potentially problematic characters from string."""

    # Patterns
    disallowed_characters = re.compile("[^a-zA-Z0-9_-]")
    consecutive_underscores = re.compile("__+")

    # Replace any disallowed characters with underscores
    string = disallowed_characters.sub("_", string)
    # Remove any consecutive underscores
    string = consecutive_underscores.sub("_", string)
    # Remove heading/trailing underscores
    string = string.strip("_")

    return string


def write_minknow_csv(df: pd.DataFrame, file_path: str):
    columns = [
        "flow_cell_id",
        "position_id",
        "sample_id",
        "experiment_id",
        "flow_cell_product_code",
        "kit",
    ]

    if df.position_id[0] == "None":
        columns.remove("position_id")

    if "alias" in df.columns and "barcode" in df.columns:
        columns.append("alias")
        columns.append("barcode")

    df_csv = df.loc[:, columns]

    df_csv.to_csv(file_path, index=False)


def generate_MinKNOW_samplesheet(process: Process, args: Namespace):
    """=== Sample sheet columns ===

    flow_cell_id                E.g. 'PAM96489'
    position_id                 Only included for PromethION runs: '1A', '1B', ... '3G'
    sample_id                   LIMS sample/pool name, stripped of problematic characters, optionally prepended with "QC_"
    experiment_id               LIMS process ID, optionally prepended with "QC_"
    flow_cell_product_code      E.g. 'FLO-MIN106D'
    kit                         Product codes separated by spaces, e.g. 'SQK-LSK109 EXP-NBD196'
    alias                       Only included for barcoded pools, LIMS sample name stripped, of problematic characters
    barcode                     E.g. 'barcode01', 'barcode02', etc, fetched from LIMS

    === Constraints ===

    Must be the same across sheet:
    - kit
    - flow_cell_product_code
    - experiment_id

    Must be unique within the same flowcell
    - alias
    - barcode

    """
    errors = []

    valid_flowcell_type_strings = [
        "FLO-PRO002 (PromethION R9.4.1)",
        "FLO-MIN106D (MinION R9.4.1)",
        "FLO-FLG001 (Flongle R9.4.1)",
        "FLO-PRO114M (PromethION R10.4.1)",
        "FLO-MIN114 (MinION R10.4.1)",
        "FLO-FLG114 (Flongle R10.4.1)",
    ]

    ont_libraries = [art for art in process.all_outputs() if art.type == "Analyte"]
    ont_libraries.sort(key=lambda art: art.id)

    rows = []
    for ont_library in ont_libraries:
        # In case of errors, skip to next artifact
        try:
            # Perform sample-level demultiplexing on the library to see if it's a QC run
            # i.e. if the ONT pooling inputs consist of Illumina pools
            library_df = get_ont_library_contents(
                ont_library=ont_library,
                ont_pooling_step_name=args.pooling_step,
            )
            qc = True if "illumina_index" in library_df.columns else False
            logging.info(f"'{ont_library.name}' parsed as {'QC' if qc else 'user'} run")
            ont_barcodes = True if "ont_barcode" in library_df.columns else False
            logging.info(
                f"'{ont_library.name}' parsed as containing {'' if ont_barcodes else 'no '}ONT barcodes"
            )

            # Assert flowcell type is written in a valid format
            assert (
                process.udf["ONT flow cell type"] in valid_flowcell_type_strings
            ), f"Invalid flow cell type {process.udf['ONT flow cell type']}."

            # Parse flowcell product code
            flowcell_product_code = process.udf["ONT flow cell type"].split(" ", 1)[0]
            flow_cell_type = (
                process.udf["ONT flow cell type"].split(" ", 1)[1].strip("()")
            )

            # Start building the row in the samplesheet corresponding to the current artifact
            row = {
                "experiment_id": process.id if not qc else f"QC_{process.id}",
                "sample_id": sanitize_string(ont_library.name)
                if not qc
                else f"QC_{sanitize_string(ont_library.name)}",
                "flow_cell_product_code": flowcell_product_code,
                "flow_cell_type": flow_cell_type,
                "kit": get_kit_string(process),
                "flow_cell_id": ont_library.udf["ONT flow cell ID"],
                "position_id": ont_library.udf["ONT flow cell position"],
            }

            # Assert position makes sense with the flowcell type
            if "PromethION" in row["flow_cell_type"]:
                assert (
                    row["position_id"] != "None"
                ), "Positions must be specified for PromethION flow cells."
            else:
                assert (
                    row["position_id"] == "None"
                ), "Positions must be unassigned for non-PromethION flow cells."

            # If there are ONT barcodes, append rows for each barcode
            if process.udf.get("ONT expansion kit") == "None":
                assert not ont_barcodes, f"ONT expansion kit is set to 'None', but library '{ont_library.name}' appears to contain ONT barcodes."
                rows.append(row)

            else:
                assert ont_barcodes, f"ONT expansion kit is assigned, but no ONT barcodes were found within library {ont_library.name}"

                # Append rows for each barcode
                alias_column_name = "illumina_pool_name" if qc else "ont_pool_name"
                barcode_rows_data = (
                    library_df[[alias_column_name, "ont_barcode"]]
                    .drop_duplicates()
                    .sort_values(by=alias_column_name)
                    .to_dict(orient="records")
                )

                for barcode_row_data in barcode_rows_data:
                    row["alias"] = sanitize_string(barcode_row_data[alias_column_name])
                    barcode_label_match = re.match(
                        ONT_BARCODE_LABEL_PATTERN,
                        barcode_row_data["ont_barcode"],
                    )
                    barcode_id = barcode_label_match.group(2)
                    row["barcode"] = f"barcode{barcode_id}"

                    assert re.match(r"barcode\d{2}", row["barcode"])
                    assert "" not in row.values(), "All fields must be populated."

                    rows.append(row.copy())

        except AssertionError as e:
            logging.error(str(e), exc_info=True)
            logging.warning(f"Skipping {ont_library.name} due to error.")
            errors.append(ont_library.name)
            continue

    # Abort on errors processing samples, else compile samplesheet
    if errors:
        raise AssertionError(f"Errors occurred when parsing artifacts {errors}")

    df = pd.DataFrame(rows)

    # Samplesheet-wide assertions
    if len(ont_libraries) > 1:
        assert all(
            ["PromethION" in fc_type for fc_type in df.flow_cell_type.unique()]
        ), "Only PromethION flowcells can be grouped together in the same sample sheet."
        assert (
            len(ont_libraries) <= 24
        ), "Only up to 24 PromethION flowcells may be started at once."
    elif len(ont_libraries) == 1 and "MinION" in df.flow_cell_type[0]:
        assert (
            df.position_id[0] == "None"
        ), "MinION flow cells should not have a position assigned."
    assert (
        len(df.flow_cell_product_code.unique()) == len(df.kit.unique()) == 1
    ), "All rows must have the same flow cell type and kits"
    assert (
        len(df.position_id.unique())
        == len(df.flow_cell_id.unique())
        == len(ont_libraries)
    ), "All rows must have different flow cell positions and IDs"

    # Generate samplesheet
    file_name = f"MinKNOW_samplesheet_{process.id}_{TIMESTAMP}_{process.technician.name.replace(' ','')}.csv"
    write_minknow_csv(df, file_name)

    return file_name


def main():
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
    parser.add_argument(
        "--file",
        required=True,
        type=str,
        help="Samplesheet file slot",
    )
    parser.add_argument(
        "--pooling_step",
        required=True,
        type=str,
        help="Name of ONT pooling step",
    )
    args = parser.parse_args()

    # Set up LIMS
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    process = Process(lims, id=args.pid)

    # Set up logging
    log_filename: str = (
        "_".join(
            [
                SCRIPT_NAME,
                process.id,
                TIMESTAMP,
                process.technician.name.replace(" ", ""),
            ]
        )
        + ".log"
    )

    logging.basicConfig(
        filename=log_filename,
        filemode="w",
        format="%(levelname)s: %(message)s",
        level=logging.INFO,
    )

    # Start logging
    logging.info(f"Script '{SCRIPT_NAME}' started at {TIMESTAMP}.")
    logging.info(
        f"Launched in step '{process.type.name}' ({process.id}) by {process.technician.name}."
    )
    args_str = "\n\t".join([f"'{arg}': {getattr(args, arg)}" for arg in vars(args)])
    logging.info(f"Script called with arguments: \n\t{args_str}")

    try:
        file_name = generate_MinKNOW_samplesheet(process=process, args=args)
        logging.info("Uploading samplesheet to LIMS...")
        upload_file(
            file_name,
            args.file,
            process,
            lims,
        )

        logging.info("Moving samplesheet to ngi-nas-ns...")
        try:
            shutil.move(
                file_name,
                f"/srv/ngi-nas-ns/samplesheets/nanopore/{dt.now().year}/{file_name}",
            )
        except:
            logging.error("Failed to move samplesheet to ngi-nas-ns.", exc_info=True)
        else:
            logging.info("Samplesheet moved to ngi-nas-ns.")

    except Exception as e:
        # Post error to LIMS GUI
        logging.error(str(e), exc_info=True)
        logging.shutdown()
        upload_file(
            file_path=log_filename,
            file_slot=args.log,
            process=process,
            lims=lims,
        )
        os.remove(log_filename)
        sys.stderr.write(str(e))
        sys.exit(2)
    else:
        logging.info("")
        logging.info("Script completed successfully.")
        logging.shutdown()
        upload_file(
            file_path=log_filename,
            file_slot=args.log,
            process=process,
            lims=lims,
        )
        # Check log for errors and warnings
        log_content = open(log_filename).read()
        os.remove(log_filename)
        if "ERROR:" in log_content or "WARNING:" in log_content:
            sys.stderr.write(
                "Script finished successfully, but log contains errors or warnings, please have a look."
            )
            sys.exit(2)
        else:
            sys.exit(0)


if __name__ == "__main__":
    main()
