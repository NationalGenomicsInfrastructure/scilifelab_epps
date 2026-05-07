#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from datetime import datetime as dt

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from numpy import minimum
from tabulate import tabulate

from scilifelab_epps.utils import formula
from scilifelab_epps.zika.utils import fetch_sample_data

DESC = """
EPP "ONT pooling", file slot "ONT pooling log".

For each pool, given either a target amount (UDF 'Amount (fmol)') or a target volume (UDF 'Final Volume (uL)'),
will calculate the other UDF (if both are specified, amount will overwrite volume) as well as corresponding
volumes, amounts and molar percentages of samples. These values are tabulated and returned as a log.
"""


def main(lims, args):
    try:
        currentStep = Process(lims, id=args.pid)
        log = []

        pools = [art for art in currentStep.all_outputs() if art.type == "Analyte"]
        pools.sort(key=lambda pool: pool.name)

        to_fetch = {
            # Search within step
            "sample_name": "art_tuple[0]['uri'].name",
            "dst_name": "art_tuple[1]['uri'].name",
            "vol_ul": "art_tuple[0]['uri'].udf['Volume (ul)']",
            "conc": "art_tuple[0]['uri'].udf['Concentration']",
            "conc_units": "art_tuple[0]['uri'].udf['Conc. Units']",
            # Seach recursively
            "size_bp": "Size (bp)",
        }

        df = fetch_sample_data(currentStep, to_fetch)

        assert all([i in ["ng/ul", "nM"] for i in df.conc_units]), (
            "Some of the pool inputs have invalid concentration units."
        )

        df["conc_nM"] = df.apply(
            lambda x: (
                x
                if x["conc_units"] == "nM"
                else formula.ng_ul_to_nM(x["conc"], x["size_bp"])
            ),
            axis=1,
        )

        for pool in pools:
            log.append(f"{pool.name}")

            # Subset data
            df_pool = df[df.dst_name == pool.name].copy()

            # Get molar proportions between samples within pool
            df_pool["prop_nM"] = df_pool.conc_nM / sum(df_pool.conc_nM)
            df_pool["prop_nM_inv"] = (1 / df_pool["prop_nM"]) / sum(
                1 / df_pool["prop_nM"]
            )

            # If amount is specified, use for calculations and ignore target vol
            try:
                pool_target_amt_fmol = pool.udf["Amount (fmol)"]
                pool_target_vol = None

                log.append(f"Target amt: {round(pool_target_amt_fmol, 1)} fmol")

                target_amt_fmol = pool_target_amt_fmol / len(df_pool)

                # Apply molar proportions to target amount to get transfer volumes
                df_pool["transfer_vol_ul"] = minimum(
                    target_amt_fmol / df_pool.conc_nM, df_pool.vol_ul
                )
                pool_transfer_vol = sum(df_pool.transfer_vol_ul)
                df_pool["transfer_amt_fmol"] = df_pool.transfer_vol_ul * df_pool.conc_nM
                pool_transfer_amt = sum(df_pool.transfer_amt_fmol)

            # If amount is omitted but not target volume, use target volume for calculations
            except KeyError:
                pool_target_amt_fmol = None
                pool_target_vol = pool.udf["Final Volume (uL)"]

                log.append(f"Target vol: {round(pool_target_vol, 1)} uL")

                # Apply molar proportions to target volume to get transfer amounts
                df_pool["transfer_vol_ul"] = minimum(
                    pool_target_vol * df_pool.prop_nM_inv, df_pool.vol_ul
                )
                pool_transfer_vol = sum(df_pool.transfer_vol_ul)
                df_pool["transfer_amt_fmol"] = df_pool.transfer_vol_ul * df_pool.conc_nM
                pool_transfer_amt = sum(df_pool.transfer_amt_fmol)

            # Evaluate target fraction
            df_pool["molar_percentage"] = (
                df_pool.transfer_amt_fmol / sum(df_pool.transfer_amt_fmol) * 100
            )

            # Update UDFs
            pool.udf["Amount (fmol)"] = round(pool_transfer_amt, 1)
            pool.udf["Final Volume (uL)"] = round(pool_transfer_vol, 1)
            pool.put()

            # Create table
            df_to_print = df_pool.loc[
                :,
                [
                    "sample_name",
                    "vol_ul",
                    "conc",
                    "conc_units",
                    "size_bp",
                    "conc_nM",
                    "transfer_vol_ul",
                    "transfer_amt_fmol",
                    "molar_percentage",
                ],
            ]
            df_to_print.reset_index(inplace=True, drop=True)
            df_to_print = df_to_print.round(1)
            log.append(tabulate(df_to_print, headers="keys"))
            log.append(f"\nFinal amt: {round(pool_transfer_amt, 1)} fmol")
            log.append(f"Final vol: {round(pool_transfer_vol, 1)} uL")
            log.append("\n")

        # Write log
        timestamp = dt.now().strftime("%y%m%d_%H%M%S")
        log_filename = "_".join(["ont_pooling_log", currentStep.id, timestamp]) + ".txt"
        with open(log_filename, "w") as logContext:
            logContext.write("\n".join(log))

        # Upload log
        for out in currentStep.all_outputs():
            if out.name == "ONT pooling log":
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
