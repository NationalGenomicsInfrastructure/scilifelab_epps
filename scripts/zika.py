#!/usr/bin/env python

DESC = """
This module contains functions used to generate worklists for the Mosquito X1
instrument Zika using sample data fetched from Illumina Clarity LIMS.

The functions are written with the intention of being re-useable for different
applications of the instrument.

Written by Alfred Kedhammar
"""

import pandas as pd
import numpy as np
from datetime import datetime as dt
from genologics.lims import Lims
from genologics.config import BASEURI, USERNAME, PASSWORD
from genologics.entities import Process


def verify_step(lims, currentStep, target_instrument, target_workflow, target_step):
    """Verify the instrument, workflow and step for a given process is correct"""

    workflows = [
        art.workflow_stages[0].workflow.name for art in currentStep.all_inputs()
    ]

    verify_bools = [
        currentStep.instrument.name == target_instrument,
        currentStep.type.name == target_step,
        all([wf == target_workflow for wf in workflows]),
    ]

    if all(verify_bools):
        return True
    else:
        return False


def fetch_sample_data(currentStep, to_fetch):
    """
    Within this function is the dictionary key2expr, its keys being the given name of a particular piece of information linked to a
    transfer input/output sample and it's values being the string that when evaluated will yield the desired info from whatever
    the variable "art_tuple" is currently pointing at.

    Given the positional arguments of a LIMS transfer process and a list of keys, it will return a dataframe containing the information
    fetched from each transfer based on the keys.
    """

    key2expr = {
        "sample_name": "art_tuple[0]['uri'].name",
        "source_fc": "art_tuple[0]['uri'].location[0].name",
        "source_well": "art_tuple[0]['uri'].location[1]",
        "conc_units": "art_tuple[0]['uri'].samples[0].artifact.udf['Conc. Units']",
        "conc": "art_tuple[0]['uri'].samples[0].artifact.udf['Concentration']",
        "vol": "art_tuple[0]['uri'].samples[0].artifact.udf['Volume (ul)']",
        "amt": "art_tuple[0]['uri'].samples[0].artifact.udf['Amount (ng)']",
        "dest_fc": "art_tuple[1]['uri'].location[0].id",
        "dest_well": "art_tuple[1]['uri'].location[1]",
        "dest_fc_name": "art_tuple[1]['uri'].location[0].name",
        "target_vol": "art_tuple[1]['uri'].udf['Total Volume (uL)']",
        "target_amt": "art_tuple[1]['uri'].udf['Amount taken (ng)']",
    }

    assert all(
        [k in key2expr.keys() for k in to_fetch]
    ), "fetch_sample_data() did not recognize key"

    l = []
    art_tuples = [
        art_tuple
        for art_tuple in currentStep.input_output_maps
        if art_tuple[0]["uri"].type == art_tuple[1]["uri"].type == "Analyte"
    ]

    for art_tuple in art_tuples:
        key2val = {}
        for k in to_fetch:
            key2val[k] = eval(key2expr[k])

        l.append(key2val)

    df = pd.DataFrame(l)

    return df


def format_worklist(df, deck, buffer_strategy):
    """
    Pivot buffer transfers and add columns in Mosquito-intepretable format.
    """

    # Pivot transfers to elucidate buffer transfers and remove zero-volume transfers
    df = pivot_buffer_transfers(df, buffer_strategy)

    # Add columns for plate positions
    df["src_pos"] = df["source_fc"].apply(lambda x: deck[x])
    df["dst_pos"] = df["dest_fc"].apply(lambda x: deck[x])

    # Convert volumes to whole nl
    df["transfer_vol"] = round(df.transfer_vol * 1000, 0)
    df["transfer_vol"] = df["transfer_vol"].astype(int)

    # Convert well names to r/c coordinates
    df["src_row"], df["src_col"] = well2rowcol(df.source_well)
    df["dst_row"], df["dst_col"] = well2rowcol(df.dest_well)

    return df


def pivot_buffer_transfers(df, buffer_strategy):
    """
    Melt buffer and sample information onto separate rows to
    produce a "one row <-> one transfer" dataframe.
    """

    # Pivot buffer transfers
    to_pivot = ["sample", "buffer"]
    to_keep = ["source_fc", "source_well", "dest_fc", "dest_well"]
    df = df.melt(
        value_vars=to_pivot,
        var_name="src_type",
        value_name="transfer_vol",
        id_vars=to_keep,
    ).sort_values(by=["dest_well", "src_type"])

    # Remove zero-vol transfers
    df = df[df.transfer_vol > 0]

    # Re-set index
    df = df.reset_index(drop=True)

    # Assign buffer transfers to buffer plate
    df.loc[df["src_type"] == "buffer", "source_fc"] = "buffer_plate"

    # Assign buffer source wells
    if buffer_strategy == "column":
        # Keep rows, but only use column 1
        df.loc[df["src_type"] == "buffer", "source_well"] = df.loc[
            df["src_type"] == "buffer", "source_well"
        ].apply(lambda x: x[0:-1] + "1")
    else:
        raise Exception("No buffer strategy defined")

    return df


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
        rows.append(rowdict[row_letter])
        cols.append(col_number)
    return rows, cols


def get_filenames(method_name, pid):

    timestamp = dt.now().strftime("%y%m%d_%H%M%S")

    wl_filename = "_".join(["zika_worklist", method_name, pid, timestamp]) + ".csv"
    log_filename = "_".join(["zika_log", method_name, pid, timestamp]) + ".csv"

    return wl_filename, log_filename


def write_worklist(df, deck, wl_filename, comments=None, strategy=None):
    """
    Write a Mosquito-interpretable advanced worklist.

    Strategies (optional):
    multi-aspirate -- If a buffer transfer is followed by a sample transfer
                      to the same well, and the sum of their volumes
                      is <= 5000 nl, use multi-aspiration.
    """

    # Format comments for printing into worklist
    if comments:
        comments = ["COMMENT, " + e for e in comments]

    # Default transfer type is simple copy
    df["transfer_type"] = "COPY"

    if strategy == "multi-aspirate":
        filter = np.all(
            [
                # End position of next transfer is the same
                df.dst_pos == df.shift(-1).dst_pos,
                # End well of the next transfer is the same
                df.dest_well == df.shift(-1).dest_well,
                # This transfer is buffer
                df.source_fc == "buffer_plate",
                # Next transfer is not buffer
                df.shift(-1).source_fc != "buffer_plate",
                # Sum of this and next transfer is >= 5 ul
                df.transfer_vol + df.shift(-1).transfer_vol <= 5000,
            ],
            axis=0,
        )
        df.loc[filter, "transfer_type"] = "MULTI_ASPIRATE"

    # PRECAUTION Keep tip change strategy in a single dict to avoid mix-ups
    tip_strats = {
        "always": ("[VAR1]", "TipChangeStrategy,always"),
        "never": ("[VAR2]", "TipChangeStrategy,never"),
    }

    # Convert all data to strings
    for c in df:
        df.loc[:, c] = df[c].apply(str)

    # Write worklist
    with open(wl_filename, "w") as wl:

        # Define worklist and variables
        wl.write("worklist,\n")
        for k in tip_strats:
            wl.write("".join(tip_strats[k]) + "\n")

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
                            tip_strats["always"][0],
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
            else:
                raise AssertionError("No transfer type defined")


def get_deck_comment(deck):
    """Convert the plate:position 'decktionary' into a worklist comment."""

    pos2plate = dict([(pos, plate) for plate, pos in deck.items()])

    l = [pos2plate[i] if i in pos2plate else "[Empty]" for i in range(1, 6)]

    deck_comment = "COMMENT, Set up layout:    " + "     ".join(l) + "\n"

    return deck_comment


def upload_log(currentStep, lims, log, log_filename):
    with open(log_filename, "w") as logContext:
        logContext.write("\n".join(log))

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