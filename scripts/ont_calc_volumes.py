from __future__ import division
from argparse import ArgumentParser
from genologics.lims import Lims
from genologics.config import BASEURI, USERNAME, PASSWORD
from genologics.entities import Process
from utils import udf_tools, formula
import sys
from datetime import datetime as dt

DESC = """
EPP "ONT calculate volumes"

Given any:

...output UDF(s)
- Amount (fmol)
- Amount (ng)
- Volume to take (uL)

...input UDF(s)
and last known UDFs
- Final Volume (uL) / Volume (uL)
- Concentration

...and last known UDF(s)
- Size (bp)

Will use ONE of the output UDFs (prioritized in the listed order) to calculate all three output UDFs.
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

            # Get last known length
            size_bp, size_bp_history = udf_tools.fetch_last(
                currentStep, art_tuple, "Size (bp)", on_fail=None, print_history=True
            )
            log.append(f"'Size (bp): {size_bp}\n{size_bp_history}")

            # Get current stats
            vol = udf_tools.fetch(art_in, "Volume (ul)")
            log.append(f"'Volume (ul)': {vol}")
            conc = udf_tools.fetch(art_in, "Concentration")
            log.append(f"'Concentration': {conc}")
            conc_units = udf_tools.fetch(art_in, "Conc. Units")
            log.append(f"'Conc. Units': {conc_units}")
            assert conc_units in [
                "ng/ul",
                "nM",
            ], f'Unsupported conc. units "{conc_units}" for art {art_in.name}'

            # Calculate volume to take, based on supplied info
            if udf_tools.is_filled(art_out, "Amount (fmol)"):
                log.append(
                    f"Basing calculations on 'Amount (fmol): {udf_tools.fetch(art_out, 'Amount (fmol)')}'"
                )
                if conc_units == "nM":
                    vol_to_take = min(
                        udf_tools.fetch(art_out, "Amount (fmol)") / conc, vol
                    )
                elif conc_units == "ng/ul":
                    vol_to_take = min(
                        formula.fmol_to_ng(
                            udf_tools.fetch(art_out, "Amount (fmol)"), size_bp
                        )
                        / conc,
                        vol,
                    )
            elif udf_tools.is_filled(art_out, "Amount (ng)"):
                log.append(
                    f"Basing calculations on 'Amount (ng)': {udf_tools.fetch(art_out, 'Amount (ng)')}"
                )
                if conc_units == "ng/ul":
                    vol_to_take = min(
                        udf_tools.fetch(art_out, "Amount (ng)") / conc, vol
                    )
                elif conc_units == "nM":
                    vol_to_take = min(
                        formula.ng_to_fmol(
                            udf_tools.fetch(art_out, "Amount (ng)"), size_bp
                        )
                        / conc,
                        vol,
                    )
            elif udf_tools.is_filled(art_out, "Volume to take (uL)"):
                log.append(
                    f"Basing calculations on 'Volume to take (uL)': {udf_tools.fetch(art_out, 'Volume to take (uL)')}"
                )
                vol_to_take = min(udf_tools.fetch(art_out, "Volume to take (uL)"), vol)
            else:
                raise AssertionError(f"No target metrics specified for {art_out.name}")

            # Based on volume to take, calculate corresponding amounts
            if conc_units == "nM":
                amt_taken_fmol = conc * vol_to_take
                amt_taken_ng = formula.fmol_to_ng(amt_taken_fmol, size_bp)
            elif conc_units == "ng/ul":
                amt_taken_ng = conc * vol_to_take
                amt_taken_fmol = formula.ng_to_fmol(amt_taken_ng, size_bp)

            log.append(f"--> 'Volume to take (uL)': {vol_to_take}")
            log.append(f"--> 'Amount (ng)': {amt_taken_ng}")
            log.append(f"--> 'Amount (fmol)': {amt_taken_fmol}")

            # Populate fields
            udf_tools.put(art_out, "Amount (fmol)", round(amt_taken_fmol, 1))
            udf_tools.put(art_out, "Amount (ng)", round(amt_taken_ng, 1))
            udf_tools.put(art_out, "Volume to take (uL)", round(vol_to_take, 1))

            log.append("\n")

        # Write log
        timestamp = dt.now().strftime("%y%m%d_%H%M%S")
        log_filename = (
            "_".join(["ont_volume_calc_log", currentStep.id, timestamp]) + ".txt"
        )
        with open(log_filename, "w") as logContext:
            logContext.write("\n".join(log))

        # Upload log
        for out in currentStep.all_outputs():
            if out.name == "Volume Calculation Log":
                for f in out.files:
                    lims.request_session.delete(f.uri)
                lims.upload_new_file(out, log_filename)

    except BaseException as e:
        sys.stderr.write(str(e))
        sys.exit(2)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
