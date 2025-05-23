#!/usr/bin/env python

import logging
import os
import re
import shutil
from argparse import ArgumentParser
from datetime import datetime as dt

import pandas as pd
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Artifact, Process
from genologics.lims import Lims
from tabulate import tabulate

from data.ONT_barcodes import format_well, ont_label2dict
from scilifelab_epps.epp import (
    get_pool_sample_label_mapping,
    traceback_to_step,
    upload_file,
)
from scilifelab_epps.utils.udf_tools import fetch
from scilifelab_epps.wrapper import epp_decorator

DESC = """ Script to generate MinKNOW samplesheet for starting ONT runs.
"""

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")


def get_ont_library_contents(
    ont_library: Artifact,
    list_contents: bool = False,
    print_dataframe: bool = False,
) -> pd.DataFrame:
    """For an ONT sequencing library, compile a dataframe with sample-level information.

    Will backtrack the library to previous ONT pooling step (if any) to elucidate
    sample and index information and decide whether to demultiplex at the level of
    ONT barcodes, Illumina indices, both or neither.

    """

    # Link samples to reagent_labels via database queries, if applicable
    if len(ont_library.reagent_labels) > 0:
        sample2label = get_pool_sample_label_mapping(ont_library)

    logging.info(
        f"Compiling sample-level information for library '{ont_library.name}'..."
    )
    library_contents_msg = f"ONT sequencing library '{ont_library.name}' consists of:"

    # Instantiate list to collect dataframe rows
    rows = []

    # Remaining possibilities:
    # (1) ONT-barcodes and Illumina indexes
    # (2) ONT-barcodes only
    # (3) Illumina-indexes only
    # (4) No labels

    # See if library can be backtracked to an ONT pooling step
    ont_pooling_traceback = traceback_to_step(
        ont_library, re.compile(r"ONT.*Pooling"), allow_multiple_inputs=True
    )
    if ont_pooling_traceback is not None:
        # Remaining possibilities:
        # (1) ONT-barcodes and Illumina indexes
        # (2) ONT-barcodes only

        _ont_pooling_step, ont_pooling_inputs, ont_pooling_output = (
            ont_pooling_traceback
        )

        library_contents_msg += f"\n - '{ont_pooling_output.name}': ONT-barcoded pool"

        # Iterate across ONT pooling inputs
        for ont_pooling_input in ont_pooling_inputs:
            if len(ont_pooling_input.samples) > 1:
                # Remaining possibilities:
                # (1) ONT-barcodes and Illumina indexes

                # ONT barcodes on Illumina libraries are assigned via UDFs rather than reagent labels, since LIMS can't handle double demultiplexing
                udf_ont_barcode_well = fetch(
                    ont_pooling_input, "ONT Barcode Well", on_fail=None
                )
                assert udf_ont_barcode_well, (
                    f"Pooling input '{ont_pooling_input.name}' consists of multiple samples, but has not been assigned an ONT barcode."
                )

                # Link ONT barcode well to ONT barcode
                ont_barcode_well2label = {}
                for label, label_properties in ont_label2dict.items():
                    ont_barcode_well2label[label_properties["well"]] = label
                sanitized_well = format_well(udf_ont_barcode_well)
                ont_barcode = ont_barcode_well2label[sanitized_well]

                library_contents_msg += f"\n\t - '{ont_pooling_input.name}': Illumina indexed pool with ONT-barcode '{ont_barcode}'"

                for sample in ont_pooling_input.samples:
                    library_contents_msg += f"\n\t\t - '{sample.name}': Illumina sample with index '{sample2label[sample.name]}'."
                    rows.append(
                        {
                            "sample_name": sample.name,
                            "sample_id": sample.id,
                            "project_name": sample.project.name,
                            "project_id": sample.project.id,
                            "illumina_index": sample2label[sample.name],
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
                assert len(ont_pooling_input.reagent_labels) == 1, (
                    f"ONT-pooling input '{ont_pooling_input.name}' lacks any reagent labels. Mixing barcoded and non-barcoded samples is not allowed."
                )

                # ONT barcode-level demultiplexing
                for ont_sample in ont_pooling_input.samples:
                    library_contents_msg += f"\n\t - '{ont_pooling_input.name}': ONT sample with barcode '{sample2label[ont_sample.name]}'"
                    rows.append(
                        {
                            "sample_name": ont_sample.name,
                            "sample_id": ont_sample.id,
                            "project_name": ont_sample.project.name,
                            "project_id": ont_sample.project.id,
                            "ont_barcode": sample2label[ont_sample.name],
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

            for sample in ont_library.samples:
                library_contents_msg += f"\n - '{sample.name}': Illumina sample with index '{sample2label[sample.name]}'."
                rows.append(
                    {
                        "sample_name": sample.name,
                        "sample_id": sample.id,
                        "project_name": sample.project.name,
                        "project_id": sample.project.id,
                        "illumina_index": sample2label[sample.name],
                        "illumina_pool_name": ont_library.name,
                        "illumina_pool_id": ont_library.id,
                    }
                )

        else:
            # Remaining possibilities:
            # (4) No labels
            sample = ont_library.samples[0]

            library_contents_msg += f"\n - {sample.name}: Non-labeled sample"
            rows.append(
                {
                    "sample_name": sample.name,
                    "sample_id": sample.id,
                    "project_name": sample.project.name,
                    "project_id": sample.project.id,
                }
            )

    if list_contents:
        logging.info(library_contents_msg)

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
        prep_kit += f" {expansion_kit.replace('.', '-')}"

    return prep_kit


def sanitize_string(string: str) -> str:
    """Remove parenthesized content and potentially problematic characters from string."""

    # Patterns
    parenthesized_content = re.compile(r"\([^()]*\)")
    disallowed_characters = re.compile("[^a-zA-Z0-9_-]")
    consecutive_underscores = re.compile("__+")

    # Remove parenthesized content
    string = parenthesized_content.sub("", string)
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


def generate_MinKNOW_samplesheet(process):
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

    qc = True if "QC" in process.type.name else False
    logging.info(f"QC run: {qc}")

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
                list_contents=True,
                print_dataframe=True,
            )
            ont_barcodes = True if "ont_barcode" in library_df.columns else False
            logging.info(
                f"'{ont_library.name}' parsed as containing {'' if ont_barcodes else 'no '}ONT barcodes"
            )

            # Assert flowcell type is written in a valid format
            assert process.udf["ONT flow cell type"] in valid_flowcell_type_strings, (
                f"Invalid flow cell type {process.udf['ONT flow cell type']}."
            )

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
                assert row["position_id"] != "None", (
                    "Positions must be specified for PromethION flow cells."
                )
            else:
                assert row["position_id"] == "None", (
                    "Positions must be unassigned for non-PromethION flow cells."
                )

            # 1) Barcodes implied from kit selection
            if process.udf.get("ONT expansion kit") != "None" or process.udf.get(
                "ONT prep kit"
            ) in ["SQK-PCB114-24"]:
                # Assert barcodes are found within library
                assert ont_barcodes, (
                    f"ONT barcodes are implied from kit selection, but no ONT barcodes were found within library {ont_library.name}"
                )

                # Append rows for each barcode
                alias_column_name = "illumina_pool_name" if qc else "sample_name"
                barcode_rows_data = (
                    library_df[[alias_column_name, "ont_barcode"]]
                    .drop_duplicates()
                    .sort_values(by=alias_column_name)
                    .to_dict(orient="records")
                )

                for barcode_row_data in barcode_rows_data:
                    row["alias"] = sanitize_string(barcode_row_data[alias_column_name])
                    barcode_id = ont_label2dict[barcode_row_data["ont_barcode"]]["num"]
                    row["barcode"] = f"barcode{str(barcode_id).zfill(2)}"

                    assert re.match(r"barcode\d{2}", row["barcode"])
                    assert "" not in row.values(), "All fields must be populated."

                    rows.append(row.copy())

            # 2) No barcodes implied from kit selection
            else:
                # Assert barcodes are not found within library
                assert not ont_barcodes, (
                    f"Library '{ont_library.name}' appears to contain ONT barcodes, but no ONT barcodes are implied from the kit selection."
                )

                # Append single row
                rows.append(row)

        except AssertionError as e:
            logging.error(str(e), exc_info=True)
            logging.warning(f"Skipping '{ont_library.name}' due to error.")
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
        assert len(ont_libraries) <= 24, (
            "Only up to 24 PromethION flowcells may be started at once."
        )
    elif len(ont_libraries) == 1 and "MinION" in df.flow_cell_type[0]:
        assert df.position_id[0] == "None", (
            "MinION flow cells should not have a position assigned."
        )
    assert len(df.flow_cell_product_code.unique()) == len(df.kit.unique()) == 1, (
        "All rows must have the same flow cell type and kits"
    )
    assert (
        len(df.position_id.unique())
        == len(df.flow_cell_id.unique())
        == len(ont_libraries)
    ), "All rows must have different flow cell positions and IDs"

    # Generate samplesheet
    file_name = f"ONT_ss_{process.id}_{TIMESTAMP}_{process.technician.name.replace(' ', '')}.csv"
    write_minknow_csv(df, file_name)

    return file_name


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    file_name = generate_MinKNOW_samplesheet(process)

    logging.info("Uploading samplesheet to LIMS...")
    upload_file(
        file_name,
        args.file,
        process,
        lims,
    )

    logging.info("Moving samplesheet to ngi-nas-ns...")
    try:
        dst = f"/srv/ngi-nas-ns/samplesheets/nanopore/{dt.now().year}"
        if not os.path.exists(dst):
            logging.info(f"Happy new year! Creating {dst}")
            os.mkdir(dst)
        shutil.copyfile(
            file_name,
            f"{dst}/{file_name}",
        )
        os.remove(file_name)
    except:
        logging.error("Failed to move samplesheet to ngi-nas-ns.", exc_info=True)
    else:
        logging.info("Samplesheet moved to ngi-nas-ns.")


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
    parser.add_argument(
        "--file",
        required=True,
        type=str,
        help="Samplesheet file slot",
    )
    args = parser.parse_args()

    main(args)
