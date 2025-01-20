#!/usr/bin/env python


import sys
from argparse import ArgumentParser
from datetime import datetime as dt

import pandas as pd
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from ont_send_reloading_info_to_db import parse_run
from tabulate import tabulate

from scilifelab_epps.utils import udf_tools

DESC = """Script for the EPP "Log fields" and file slot "Field log".

Use this script to create / update an append-only log where the step UDF info is stored.
Step-specific input assertions may be added where indicated.
"""


def main(lims, args):
    try:
        currentStep = Process(lims, id=args.pid)
        udfs_to_log = args.udfs

        timestamp = dt.now().strftime("%y%m%d_%H%M%S")

        # Collect last log file, if there is one
        file_art = [op for op in currentStep.all_outputs() if op.name == "Field log"][0]
        if file_art.files:
            file_bytes = lims.get_file_contents(uri=file_art.files[0].uri).read()
            file_str = file_bytes.decode("utf-8")
        else:
            file_str = None

        # Parse outputs and their UDFs
        if udf_tools.get_art_tuples(currentStep) == []:
            arts = [art for art in currentStep.all_inputs() if art.type == "Analyte"]
        else:
            arts = [art for art in currentStep.all_outputs() if art.type == "Analyte"]

        # Step-specific assertions
        if "ONT Sequencing and Reloading" in currentStep.type.name:
            for art in currentStep.all_inputs():
                parse_run(art)

        # Start building log dataframe
        rows = []
        for art in arts:
            row = {}
            row["Sample"] = art.name
            for udf in udfs_to_log:
                try:
                    row[udf] = art.udf[udf]
                except KeyError:
                    row[udf] = None

            rows.append(row)

        df = pd.DataFrame(rows).set_index("Sample")
        df.sort_index(inplace=True)

        new_log_name = f"UDF_log_{currentStep.id}_{timestamp}.txt"
        with open(new_log_name, "w") as f:
            # Write previous file contents
            if file_str:
                f.write(file_str)

            # Write new file contents
            f.write(timestamp + "\n")
            f.write(tabulate(df, headers="keys"))
            f.write("\n\n")

        for out in currentStep.all_outputs():
            if out.name == "Field log":
                for f in out.files:
                    lims.request_session.delete(f.uri)
                lims.upload_new_file(out, new_log_name)

    except AssertionError as e:
        sys.stderr.write(str(e))
        sys.exit(2)


if __name__ == "__main__":
    r"""Example script call:

        python \
        log_udfs.py \
        --pid {processLuid} \
        --udfs \
        'ONT run ID' \
        'ONT reload run time (hh:mm)'  \
        'ONT reload amount (fmol)'  \
        'ONT reload wash kit'
    """
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    parser.add_argument(
        "--udfs", metavar="U", type=str, nargs="+", help="UDFs to log, as strings"
    )
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
