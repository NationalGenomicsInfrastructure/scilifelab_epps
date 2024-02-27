#!/usr/bin/env python

import logging
import os
import re
import shutil
import sys
from argparse import ArgumentParser
from datetime import datetime as dt

import pandas as pd
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from data.Chromium_10X_indexes import Chromium_10X_indexes as idxs_10x
from epp_utils import udf_tools
from epp_utils.formula import well_name2num_96plate as well2num

DESC = """ Script for EPP "Generate ONT Sample Sheet" and file slot(s) "ONT sample sheet" (and optionally "Anglerfish sample sheet").
Used to generate MinKNOW (and Anglerfish) samplesheets.
"""

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")
SCRIPT_NAME: str = os.path.basename(__file__).split(".")[0]


def generate_samplesheets(process, lims, args):
    """
    === Sample sheet columns ===

    flow_cell_id                e.g. PAM96489
    position_id                 [1-3A-G] for PromethION, else None
    sample_id                   lims-sample-name, stripped of problematic characters
    experiment_id               lims-step-id
    flow_cell_product_code      e.g. FLO-MIN106D
    kit                         Product codes separated by spaces, e.g. SQK-LSK109 EXP-NBD196
    alias                       Only included for barcoded pools, sample name e.g. P12345_101
    barcode                     barcode01, barcode02, etc, fetched from LIMS

    === Constraints ===

    Must be the same across sheet:
    - kit
    - flow_cell_product_code
    - experiment_id

    Must be unique within the same flowcell
    - alias
    - barcode

    === Flowcell product codes ===

    FLO-PRO002 (PromethION R9.4.1)
    FLO-MIN106D (MinION R9.4.1)
    FLO-FLG001 (Flongle R9.4.1)
    FLO-PRO114M (PromethION R10.4.1)
    FLO-MIN114 (MinION R10.4.1)
    FLO-FLG114 (Flongle R10.4.1)

    === Outputs ===

    ONT_samplesheet_lims-step_yymmdd_hhmmss.csv
    """

    minknow_samplesheet_file = (
        minknow_samplesheet_for_qc(process)
        if "ONT QC" in process.type.name
        else minknow_samplesheet_default(process)
    )
    upload_file(
        minknow_samplesheet_file,
        "ONT sample sheet",
        process,
        lims,
    )
    shutil.copyfile(
        minknow_samplesheet_file,
        f"/srv/ngi-nas-ns/samplesheets/nanopore/{dt.now().year}/{minknow_samplesheet_file}",
    )
    os.remove(minknow_samplesheet_file)

    if "ONT QC" in process.type.name:
        anglerfish_samplesheet_file = anglerfish_samplesheet(process)
        upload_file(
            anglerfish_samplesheet_file,
            "Anglerfish sample sheet",
            process,
            lims,
        )
        shutil.copyfile(
            anglerfish_samplesheet_file,
            f"/srv/ngi-nas-ns/samplesheets/anglerfish/{dt.now().year}/{anglerfish_samplesheet_file}",
        )
        os.remove(anglerfish_samplesheet_file)


def minknow_samplesheet_default(process):
    arts = [art for art in process.all_outputs() if art.type == "Analyte"]
    arts.sort(key=lambda art: art.id)

    rows = []
    for art in arts:
        row = {
            "flow_cell_id": art.udf.get("ONT flow cell ID"),
            "position_id": art.udf.get("ONT flow cell position"),
            "sample_id": strip_characters(art.name),
            "experiment_id": f"{process.id}",
            "flow_cell_product_code": process.udf["ONT flow cell type"].split(" ")[0],
            "flow_cell_type": process.udf["ONT flow cell type"]
            .split(" ")[1]
            .strip("()"),
            "kit": get_kit_string(process),
        }

        if "PromethION" in row["flow_cell_type"]:
            assert (
                row["position_id"] != "None"
            ), "Positions must be specified for PromethION flow cells."

        # Add extra columns for barcodes
        if len(art.reagent_labels) > 0:
            assert (
                process.udf.get("ONT expansion kit") != "None"
            ), f"Barcodes found in pool {art.name}, but no barcode kit specified."
        if process.udf.get("ONT expansion kit") != "None":
            assert (
                len(art.reagent_labels) > 0
            ), f"No barcodes found within pool {art.name}"
            label_tuples = [(e[0], e[1]) for e in zip(art.samples, art.reagent_labels)]
            label_tuples.sort(key=str)
            for sample, label in label_tuples:
                row["alias"] = strip_characters(sample.name)
                row["barcode"] = strip_characters(
                    "barcode" + label[0:2]
                )  # TODO double check extraction of barcode number
                rows.append(row.copy())
        else:
            rows.append(row)

        assert "" not in row.values(), "All fields must be populated."

    df = pd.DataFrame(rows)

    if len(arts) > 1:
        assert all(
            ["PromethION" in fc_type for fc_type in df.flow_cell_type.unique()]
        ), "Only PromethION flowcells can be grouped together in the same sample sheet."
        assert (
            len(arts) <= 24
        ), "Only up to 24 PromethION flowcells may be started at once."
    elif len(arts) == 1 and "MinION" in df.flow_cell_type[0]:
        assert (
            df.position_id[0] == "None"
        ), "MinION flow cells should not have a position assigned."

    assert (
        len(df.flow_cell_product_code.unique()) == len(df.kit.unique()) == 1
    ), "All rows must have the same flow cell type and kits"
    assert (
        len(df.position_id.unique()) == len(df.flow_cell_id.unique()) == len(arts)
    ), "All rows must have different flow cell positions and IDs"

    file_name = write_minknow_csv(
        df,
        f"MinKNOW_samplesheet_{process.id}_{TIMESTAMP}_{process.technician.name.replace(' ','')}.csv",
    )

    return file_name


def minknow_samplesheet_for_qc(process):
    measurements = []

    # Differentiate file outputs from measurements outputs by name, i.e. "P12345_101" vs "Scilifelab SampleSheet"
    sample_pattern = re.compile(r"P\d{5}_\d{3,4}")
    for art in process.all_outputs():
        if re.search(sample_pattern, art.name):
            measurements.append(art)

    # Build an input output map objects omitting the files
    art_tuples = []
    for art_tuple in process.input_output_maps:
        if art_tuple[1]["uri"].id in [m.id for m in measurements]:
            art_tuples.append(art_tuple)
        else:
            pass

    rows = []

    # Iterate through the input Illumina pools one by one
    for pool in process.all_inputs():
        # Find all outputs belonging to the current Illumina pool
        pool_samples = [
            art_tuple[1]["uri"]
            for art_tuple in art_tuples
            if art_tuple[0]["uri"].id == pool.id
        ]

        # Assert ONT barcode wells are correctly populated
        barcode_wells_in_pool = [
            udf_tools.fetch(art, "ONT Barcode Well", on_fail=None)
            for art in pool_samples
        ]

        assert (
            len(set(barcode_wells_in_pool)) == 1
        ), f"All ONT barcodes must be the same within a pool, not the case for pool {pool.name}"
        barcode_well = barcode_wells_in_pool[0]

        # Assert well looks like a well, e.g. "A:11", "G4", "C:1"
        barcode_well_pattern = re.compile(r"^[A-H]:?([1-9]$|(1[0-2])$)")

        if barcode_well:
            assert (
                process.udf.get("ONT expansion kit") != "None"
            ), "ONT Barcodes have been assigned, but no 'ONT expansion kit' is specified."

            assert re.match(
                barcode_well_pattern, barcode_well
            ), f"The 'ONT Barcode Well' entry '{barcode_well}' in pool {pool.name} doesn't look like a plate well."

            if barcode_well not in well2num:
                barcode_well = barcode_well[0] + ":" + barcode_well[1:]
            barcode_int = well2num[barcode_well]
        else:
            assert (
                process.udf.get("ONT expansion kit") == "None"
            ), "ONT Barcodes have not been assigned."

        row = {
            "position_id": "None",
            "flow_cell_id": process.udf["ONT flow cell ID"],
            "sample_id": f"QC_{art.name}",
            "experiment_id": f"{process.id}",
            "flow_cell_product_code": process.udf["ONT flow cell type"].split(" ")[0],
            "flow_cell_type": process.udf["ONT flow cell type"]
            .split(" ")[1]
            .strip("()"),
            "kit": get_kit_string(process),
        }

        if barcode_well:
            row["alias"] = strip_characters(pool.name)
            row["barcode"] = "barcode" + str(barcode_int).zfill(2)

        rows.append(row)

    df = pd.DataFrame(rows)

    if "barcode" in df.columns:
        assert all(
            df.barcode.notna()
        ), "Nanopore barcodes must be specified for either ALL samples, or NONE."

        assert len(df.barcode.unique()) == len(
            process.all_inputs()
        ), "Nanopore barcodes are shared between Illumina pools"

    file_name = write_minknow_csv(
        df,
        f"ONT_samplesheet_{df.experiment_id.unique()[0]}_{TIMESTAMP}.csv",
    )
    return file_name


def anglerfish_samplesheet(process):
    measurements = []

    # Differentiate file outputs from measurements outputs by name, i.e. "P12345_101" vs "Scilifelab SampleSheet"
    sample_pattern = re.compile(r"P\d{5}_\d{3,4}")
    for art in process.all_outputs():
        if re.search(sample_pattern, art.name):
            measurements.append(art)

    ont_barcode_bools = [
        udf_tools.fetch(art, "ONT Barcode Well", on_fail=None) is not None
        for art in measurements
    ]

    if all(ont_barcode_bools):
        ont_barcodes = True
    elif not any(ont_barcode_bools):
        ont_barcodes = False
    else:
        raise AssertionError(
            "ONT barcodes must be present either for all samples or for none."
        )

    rows = []

    # Iterate through the samples
    for sample in measurements:
        if ont_barcodes:
            barcode_well = udf_tools.fetch(sample, "ONT Barcode Well")

            if barcode_well not in well2num:
                barcode_well = barcode_well[0] + ":" + barcode_well[1:]
            barcode_int = well2num[barcode_well]

            fastq_path = f"./fastq_pass/barcode{str(barcode_int).zfill(2)}/*.fastq.gz"  # Assuming the Anglerfish working dir is the ONT run dir TODO

        elif not ont_barcodes:
            fastq_path = "./fastq_pass/*.fastq.gz"

        index_seq_list, adaptors_name = get_index_info(sample)

        # For multi-index samples, append multiple rows
        for index_seq in index_seq_list:
            row = {
                "sample_name": sample.name,
                "adaptors": adaptors_name,
                "index": index_seq,
                "fastq_path": fastq_path,
            }

            rows.append(row)

    df = pd.DataFrame(rows)
    df.sort_values(by="sample_name", inplace=True)

    file_name = f"Anglerfish_samplesheet_{process.id}_{TIMESTAMP}.csv"
    df.to_csv(
        file_name,
        header=False,
        index=False,
    )

    return file_name


def get_index_info(sample):
    """
    Input: LIMS API measurement object

    Output: tuple(
        List of indexes (either i7 or i7-i5),
        The name of the adaptors as defined in Anglerfish config
        )
    """

    index_seq = None

    assert (
        len(sample.reagent_labels) == 1
    ), f"Multiple reagent labels found for sample {sample.name}"

    label = sample.reagent_labels[0]

    index_pattern = re.compile("[ACTG]{4,}-?[ACTG]{4,}")

    ### Get the index sequence ####

    # 1) Look for idx sequence contained directly in .reagent_labels attribute
    index_search = re.search(index_pattern, label)

    if index_search:
        index_seq = index_search.group()

    else:
        # 2) Look for idx among 10X idxs
        if label in idxs_10x:
            idx_10x_list = idxs_10x[label]

            if len(idx_10x_list) == 2:
                # Return i7-i5
                index_seq = "-".join(idx_10x_list)
            elif len(idx_10x_list) == 4:
                # Return list of combination i7 idxs
                index_seq = idx_10x_list
            else:
                raise AssertionError("Unrecognized format of 10X index.")

    ### Get the name of the adaptors ###

    # For now, only support truseq and truseq_dual adaptors TODO
    if "-" in index_seq:
        adaptors_name = "truseq_dual"
    else:
        adaptors_name = "truseq"

    # Return

    if index_seq:
        if not isinstance(index_seq, list):
            index_seq = [index_seq]
        return index_seq, adaptors_name

    else:
        assert index_search, f"No index information found for sample {sample.name}"


def upload_file(file_name, file_slot, process, lims):
    for out in process.all_outputs():
        if out.name == file_slot:
            for f in out.files:
                lims.request_session.delete(f.uri)
            lims.upload_new_file(out, file_name)


def write_minknow_csv(df, file_name):
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

    df_csv.to_csv(file_name, index=False)

    return file_name


def strip_characters(input_string):
    """Remove potentially problematic characters from string."""

    allowed_characters = re.compile("[^a-zA-Z0-9_-]")
    # Replace any disallowed characters with underscores
    subbed_string = allowed_characters.sub("_", input_string)

    # Remove any consecutive underscores
    string_to_shorten = re.compile("__+")
    shortened_string = string_to_shorten.sub("_", subbed_string)

    return shortened_string


def get_kit_string(process):
    """Combine prep kit and expansion kit UDFs (if any) into space-separated string"""
    prep_kit = process.udf.get("ONT prep kit")
    expansion_kit = process.udf.get("ONT expansion kit")

    if expansion_kit != "None":
        prep_kit += f" {expansion_kit.replace('.','-')}"

    return prep_kit


def main():
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", type=str, help="Lims ID for current Process")
    parser.add_argument("--log", type=str, help="Which log file slot to use")
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
        generate_samplesheets(process, lims, args)
    except Exception as e:
        # Post error to LIMS GUI
        logging.error(str(e), exc_info=True)
        logging.shutdown()
        upload_file(
            file_name=log_filename,
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
            file_name=log_filename,
            file_slot=args.log,
            process=process,
            lims=lims,
        )
        # Check log for errors and warnings
        log_content = open(log_filename).read()
        os.remove(log_filename)
        if "ERROR:" in log_content or "WARNING:" in log_content:
            sys.stderr.write(
                "Script finished successfully, but log contains erros or warnings, please have a look."
            )
            sys.exit(2)
        else:
            sys.exit(0)


if __name__ == "__main__":
    main()