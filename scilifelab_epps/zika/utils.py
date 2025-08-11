#!/usr/bin/env python

DESC = """
This module contains functions used to generate worklists for the Mosquito X1
instrument Zika using sample data fetched from Illumina Clarity LIMS.

The functions are written with the intention of being re-useable for different
applications of the instrument.

Written by Alfred Kedhammar
"""

import sys
from datetime import datetime as dt

import numpy as np
import pandas as pd
from genologics.entities import Process

from scilifelab_epps.utils.udf_tools import fetch_last


def verify_step(currentStep, targets=None):
    """
    Given a LIMS step and a list of targets, check whether they match. Workflow information unfortunately needs to be excavated from the samples.

    The "targets" consist of a list of tuples, whose elements are partial string matches of a workflow and step, respectively.
    Empty strings will match any workflow or step.
    """

    if currentStep.instrument.name == "Zika":
        if not targets:
            # Instrument is correct and no workflows or steps are specified
            return True

        elif any(
            target_tuple[1] in currentStep.type.name and target_tuple[0] == ""
            for target_tuple in targets
        ):
            # Instrument and step are correct and no workflow is specified
            return True

        else:
            # Need to check all samples match at least one ongoing workflow / step combo of the targets
            sample_bools = []
            for art in [
                art for art in currentStep.all_inputs() if art.type == "Analyte"
            ]:
                active_stages = [
                    stage_tuple
                    for stage_tuple in art.workflow_stages_and_statuses
                    if stage_tuple[1] == "IN_PROGRESS"
                ]
                sample_bools.append(
                    # Sample has at least one ongoing target step in the target workflow
                    any(
                        workflow_string in active_stage[0].workflow.name
                        and step_string in active_stage[2]
                        for active_stage in active_stages
                        for workflow_string, step_string in targets
                    )
                )

            return all(sample_bools)

    else:
        return False


class CheckLog(Exception):
    def __init__(self, log, log_filename, lims, currentStep):
        write_log(log, log_filename)
        upload_log(currentStep, lims, log_filename)

        sys.stderr.write("ERROR: Check log for more info.")
        sys.exit(2)


def fetch_sample_data(currentStep: Process, to_fetch: dict) -> pd.DataFrame:
    """
    Given a LIMS step and a dictionary detailing which info to fetch, this function
    will go through all analyte input/output tuples of the step (or previous steps)
    and try to fetch the relevant information in a Pandas dataframe.

    In the dictionary "to_fetch":
    - Dict keys will be the column names in the returned df
    - Dict items are either...
       1) an expression to be evaulated to fetch the info
       2) the name of a UDF to fetch recursively

    Examples of dictionary contents:
    to_fetch = {
        "vol"   : "art_tuple[0]['uri'].udf['Final Volume (uL)']",       # Eval string
        "conc"  : "art_tuple[0]['uri'].udf['Final Concentration']",     # Eval string
        "size"  : 'Size (bp)'                                           # UDF name, to fetch recursively
    }
    """

    # Fetch all input/output sample tuples
    art_tuples = [
        art_tuple
        for art_tuple in currentStep.input_output_maps
        if art_tuple[0]["uri"].type == art_tuple[1]["uri"].type == "Analyte"
    ]

    # Fetch all target data
    rows = []
    for art_tuple in art_tuples:
        row = {}
        for col_name, udf_query in to_fetch.items():
            if "art_tuple" in udf_query:
                try:
                    row[col_name] = eval(udf_query)
                except KeyError:
                    row[col_name] = None
            else:
                row[col_name] = fetch_last(
                    target_art=art_tuple[0]["uri"],
                    target_udfs=udf_query,
                    include_current=True,
                )
        rows.append(row)

    # Transform to dataframe
    df = pd.DataFrame(rows)

    return df


def format_worklist(df, deck):
    """
    - Add columns in Mosquito-intepretable format

    - Split transfers exceeding max pipetting volume.
      Create splits of 5000 nl at a time until the remaining volume is >5000 and <= 10000,
      then split it in half.

    - Sort by buffer/sample, dst col, dst row
    """

    # Add columns for plate positions
    df["src_pos"] = df["src_name"].apply(lambda x: deck[x])
    df["dst_pos"] = df["dst_name"].apply(lambda x: deck[x])

    # Convert volumes to whole nl
    df["transfer_vol"] = round(df.transfer_vol * 1000, 0)
    df["transfer_vol"] = df["transfer_vol"].astype(int)

    # Convert well names to r/c coordinates
    df["src_row"], df["src_col"] = well2rowcol(df.src_well)
    df["dst_row"], df["dst_col"] = well2rowcol(df.dst_well)

    # Sort df
    try:
        # Normalization, buffer first, work column-wise dst
        df.sort_values(by=["src_type", "dst_col", "dst_row"], inplace=True)
    except KeyError:
        # Pooling, sort by column-wise dst (pool), then by descending transfer volume
        df.sort_values(
            by=["dst_col", "dst_row", "transfer_vol"],
            ascending=[True, True, False],
            inplace=True,
        )
    df.reset_index(inplace=True, drop=True)

    # Split >5000 nl transfers

    assert all(df.transfer_vol < 180000), "Some transfer volumes exceed 180 ul"
    max_vol = 5000

    # We need to split the row across multiple sub-transfers.
    # Make a list to house the sub-transfers as dicts.
    subtransfers = []

    # Iterate across rows
    for _idx, row in df.iterrows():
        # If transfer volume of current row exceeds max
        if row.transfer_vol > max_vol:
            # Create a row corresponding to the max permitted volume
            max_vol_transfer = row.copy().to_dict()
            max_vol_transfer["transfer_vol"] = max_vol

            # As long as the transfer volume of the current row exceeds twice the max
            while row.transfer_vol > 2 * max_vol:
                # Add a max-volume sub-transfer and deduct the same volume from the current row
                subtransfers.append(max_vol_transfer)
                row.transfer_vol -= max_vol

            # The remaining volume is higher than the max but lower than twice the max. Split this volume across two transfers.
            final_split = row.copy().to_dict()
            final_split["transfer_vol"] = round(row.transfer_vol / 2)
            # Append both
            for i in range(2):
                subtransfers.append(final_split)

        else:
            subtransfers.append(row.to_dict())

    # Format all resolved sub-transfers back into dataframe
    df_split = pd.DataFrame(subtransfers)

    return df_split


class VolumeOverflow(Exception):
    pass


def resolve_buffer_transfers(
    df=None,
    wl_comments=None,
    buffer_strategy="adaptive",
    well_dead_vol=5,
    well_max_vol=180,
    zika_max_vol=5,
):
    """
    Melt buffer and sample information onto separate rows to
    produce a "one row <-> one transfer" dataframe.
    """

    # Pivot buffer transfers
    df.rename(columns={"sample_vol": "sample", "buffer_vol": "buffer"}, inplace=True)
    to_pivot = ["sample", "buffer"]
    to_keep = ["src_name", "src_well", "dst_name", "dst_well"]
    df = df.melt(
        value_vars=to_pivot,
        var_name="src_type",
        value_name="transfer_vol",
        id_vars=to_keep,
    )

    # Sort df
    split_dst_well = df.dst_well.str.split(":", expand=True)
    df["dst_well_row"] = split_dst_well[0]
    df["dst_well_col"] = split_dst_well[1].apply(int)

    df.sort_values(by=["src_type", "dst_well_col", "dst_well_row"], inplace=True)

    # Remove zero-vol transfers
    df = df[df.transfer_vol > 0]

    # Re-set index
    df = df.reset_index(drop=True)

    df.loc[df.src_type == "buffer", "src_name"] = "buffer_plate"
    df.loc[df.src_type == "buffer", "src_well"] = np.nan

    # Assign buffer src wells
    if buffer_strategy == "first_column":
        # Keep rows, but only use column 1
        df.loc[df["src_type"] == "buffer", "src_well"] = df.loc[
            df["src_type"] == "buffer", "src_well"
        ].apply(lambda x: x[0:-1] + "1")

    elif buffer_strategy == "adaptive":
        # Make well iterator
        wells = []
        for col in range(1, 13):
            for row in list("ABCDEFGH"):
                wells.append(f"{row}:{col}")
        well_iter = iter(wells)

        # Start "filling up" buffer wells based on transfer list
        try:
            # Start at first well
            current_well = next(well_iter)
            current_well_vol = well_dead_vol

            for idx, row in df[df.src_type == "buffer"].iterrows():
                # How many subtransfers will be needed?
                n_transfers = (row.transfer_vol // zika_max_vol) + 1
                # Estimate 0.2 ul loss per transfer due to overaspiration
                vol_to_add = row.transfer_vol + 0.2 * n_transfers

                # TODO support switching buffer wells in the middle of subtransfer block
                if current_well_vol + vol_to_add > well_max_vol:
                    # Start on the next well
                    current_well = next(well_iter)
                    current_well_vol = well_dead_vol

                current_well_vol += vol_to_add
                df.loc[idx, "src_well"] = current_well

        except StopIteration:
            raise AssertionError("Total buffer volume exceeds plate capacity.")

        wl_comments.append(
            f"Fill up the buffer plate column-wise up to well {current_well} with {well_max_vol} uL buffer."
        )

    else:
        raise Exception("No buffer strategy defined")

    return df, wl_comments


def well2rowcol(well_iter):
    """
    Translates iterable of well names to list of row/column integer tuples to specify
    well location in Mosquito worklists.
    """

    # In an advanced worklist: startcol, endcol, row
    rows = []
    cols = []
    for well in well_iter:
        [row_letter, col_number] = str.split(well, sep=":")
        rowdict = {}
        for l, n in zip("ABCDEFGH", "12345678"):
            rowdict[l] = n
        rows.append(int(rowdict[row_letter]))
        cols.append(int(col_number))
    return rows, cols


def get_filenames(method_name, pid):
    timestamp = dt.now().strftime("%y%m%d_%H%M%S")

    wl_filename = "_".join(["zika_worklist", method_name, pid, timestamp]) + ".csv"
    log_filename = "_".join(["zika_log", method_name, pid, timestamp]) + ".log"

    return wl_filename, log_filename


def write_worklist(df, deck, wl_filename, comments=None, max_transfers_per_tip=10):
    """
    Write a Mosquito-interpretable advanced worklist.
    """

    # Replace all commas with semi-colons, so they can be printed without truncating the worklist
    for c, is_string in zip(df.columns, df.map(type).eq(str).all()):
        if is_string:
            df[c] = df[c].apply(lambda x: x.replace(",", ";"))

    # Format comments for printing into worklist
    if comments:
        comments = ["COMMENT, " + e for e in comments]

    # Default transfer type is simple copy
    df["transfer_type"] = "COPY"

    # PRECAUTION Keep tip change strategy variable definitions immutable
    tip_strats = {"always": "[VAR1]", "never": "[VAR2]"}

    # Initially, set all transfers to always change tips
    df["tip_strat"] = tip_strats["always"]

    # As default, keep tips between buffer transfers
    df.loc[df.src_name == "buffer_plate", "tip_strat"] = tip_strats["never"]
    # Add tip changes every x buffer transfers
    n_transfers = 0
    for i, r in df.iterrows():
        if (
            r.tip_strat == tip_strats["never"]
            and n_transfers < max_transfers_per_tip - 1
        ):
            n_transfers += 1
        elif (
            r.tip_strat == tip_strats["never"]
            and n_transfers >= max_transfers_per_tip - 1
        ):
            df.loc[i, "tip_strat"] = tip_strats["always"]
            n_transfers = 0
        elif r.tip_strat != tip_strats["never"]:
            n_transfers = 0
        else:
            raise AssertionError("Unpredicted case")

    df.sort_index(inplace=True)
    df.reset_index(inplace=True, drop=True)

    # Convert all data to strings
    for c in df:
        df.loc[:, c] = df[c].apply(str)

    # Write worklist
    with open(wl_filename, "w") as wl:
        wl.write("worklist,\n")

        # Define variables
        variable_definitions = []
        for tip_strat in [
            tip_strat
            for tip_strat in tip_strats.items()
            if tip_strat[1] in df.tip_strat.unique()
        ]:
            variable_definitions.append(f"{tip_strat[1]}TipChangeStrategy")
            variable_definitions.append(tip_strat[0])
        wl.write(",".join(variable_definitions) + "\n")

        # Add barcode settings
        wl.write("[PROTOCOLOPTION]SCAN_BARCODES,TRUE\n")
        wl.write("[PROTOCOLOPTION]ABORT_ON_BAD_BARCODE,TRUE\n")

        # Write header
        wl.write(f"COMMENT, This is the worklist {wl_filename}\n")
        if comments:
            for line in comments:
                wl.write(line + "\n")
        wl.write(get_deck_comment(deck))

        # Write transfers
        for i, r in df.iterrows():
            if r.transfer_type == "COPY":
                wl.write(
                    ",".join(
                        [
                            r.transfer_type,
                            r.src_pos,
                            r.src_col,
                            r.src_col,
                            r.src_row,
                            r.dst_pos,
                            r.dst_col,
                            r.dst_row,
                            r.transfer_vol,
                            r.tip_strat,
                        ]
                    )
                    + "\n"
                )
            elif r.transfer_type == "MULTI_ASPIRATE":
                wl.write(
                    ",".join(
                        [
                            r.transfer_type,
                            r.src_pos,
                            r.src_col,
                            r.src_row,
                            "1",
                            r.transfer_vol,
                        ]
                    )
                    + "\n"
                )
            elif r.transfer_type == "CHANGE_PIPETTES":
                wl.write(r.transfer_type + "\n")
            else:
                raise AssertionError("No transfer type defined")

        wl.write("COMMENT, Done")


def get_deck_comment(deck):
    """Convert the plate:position 'decktionary' into a worklist comment"""

    pos2plate = {pos: plate for plate, pos in deck.items()}

    l = [
        pos2plate[i].replace(",", "") if i in pos2plate else "[Empty]"
        for i in range(1, 6)
    ]

    deck_comment = "COMMENT, Set up layout:    " + "     ".join(l) + "\n"

    return deck_comment


def write_log(log, log_filename):
    with open(log_filename, "w") as logContext:
        logContext.write("\n".join(log))


def upload_log(currentStep, lims, log_filename):
    for out in currentStep.all_outputs():
        if out.name == "Mosquito Log":
            for f in out.files:
                lims.request_session.delete(f.uri)
            lims.upload_new_file(out, log_filename)


def upload_csv(currentStep, lims, wl_filename):
    for out in currentStep.all_outputs():
        if out.name == "Mosquito CSV File":
            for f in out.files:
                lims.request_session.delete(f.uri)
            lims.upload_new_file(out, wl_filename)
