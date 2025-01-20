#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from datetime import datetime as dt

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from scilifelab_epps.utils import formula, udf_tools

DESC = """ EPP "ONT Update Amounts".

Calculate and populate the ng and fmol amount of the output UDFs of the current step.

The calculation is based on the "Concentration", "Conc. Units" and "Volume (ul)" output UDFs as well as the last recorded "Size (bp)" UDF of the sample.

Alfred Kedhammar, NGI SciLifeLab
"""


def main(lims, args):
    try:
        currentStep = Process(lims, id=args.pid)

        log = []
        art_tuples = udf_tools.get_art_tuples(currentStep)

        for art_tuple in art_tuples:
            art_in = art_tuple[0]["uri"]
            art_out = art_tuple[1]["uri"]

            log.append(f"Input {art_in.name} --> Output {art_out.name}")

            # Get size
            if udf_tools.is_filled(art_out, "Size (bp)"):
                size_bp = udf_tools.fetch(art_out, "Size (bp)")
                log.append(f"'Size (bp)': {size_bp}")
            else:
                # Fetch recursively, if appropriate
                if (
                    "ONT End-Prep" in currentStep.type.name
                    or "ONT Barcoding" in currentStep.type.name
                ):
                    size_bp, size_bp_history = udf_tools.fetch_last(
                        target_art=art_out,
                        target_udfs="Size (bp)",
                        log_traceback=True,
                        return_traceback=True,
                        on_fail=None,
                    )
                    log.append(f"'Size (bp)': {size_bp}\n{size_bp_history}")
                else:
                    raise AssertionError(f"Size is not provided for {art_out.name}")

            # Get current metrics
            vol = udf_tools.fetch(art_out, "Volume (ul)")
            log.append(f"'Volume (ul)': {vol}")
            conc_units = udf_tools.fetch(art_out, "Conc. Units")
            assert conc_units in [
                "ng/ul",
                "nM",
            ], f'Unsupported conc. units "{conc_units}" for art {art_out.name}'
            log.append(f"'Conc. Units': {conc_units}")

            # Fetch or calculate conc in ng/ul
            if conc_units == "nM" and size_bp:
                conc_nM = udf_tools.fetch(art_out, "Concentration")
                log.append(f"'Concentration': {conc_nM}")
                conc_ng_ul = formula.nM_to_ng_ul(nM=conc_nM, bp=size_bp)
                log.append(f"--> Concentration (ng/ul): {conc_ng_ul}")
            elif conc_units == "ng/ul":
                conc_ng_ul = udf_tools.fetch(art_out, "Concentration")
                log.append(f"'Concentration': {conc_ng_ul}")
            else:
                raise AssertionError(f"Cannot parse concentration of {art_out.name}")

            # Calculate and put ng amount
            amount_ng = round(conc_ng_ul * vol, 2)
            udf_tools.put(art_out, "Amount (ng)", round(amount_ng, 2), on_fail=None)
            log.append(f"--> 'Amount (ng)': {amount_ng}")

            # Calculate and put fmol amount
            amount_fmol = round(
                formula.ng_to_fmol(amount_ng, size_bp),
                2,
            )
            log.append(f"--> 'Amount (fmol)': {amount_fmol}")
            udf_tools.put(
                art_out,
                "Amount (fmol)",
                round(
                    formula.ng_to_fmol(amount_ng, size_bp),
                    2,
                ),
                on_fail=None,
            )
            log.append("\n")

        # Write log
        timestamp = dt.now().strftime("%y%m%d_%H%M%S")
        log_filename = (
            "_".join(["ont_update_amounts_log", currentStep.id, timestamp]) + ".txt"
        )
        with open(log_filename, "w") as logContext:
            logContext.write("\n".join(log))

        # Upload log
        for out in currentStep.all_outputs():
            if out.name == "ONT Update Amounts log":
                for f in out.files:
                    lims.request_session.delete(f.uri)
                lims.upload_new_file(out, log_filename)

    except Exception as e:
        sys.stderr.write(str(e))
        sys.exit(2)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
