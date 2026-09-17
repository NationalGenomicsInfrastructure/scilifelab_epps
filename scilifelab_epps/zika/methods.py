#!/usr/bin/env python

DESC = """
This module contains methods for normalization and pooling on the Mosquito Zika instrument.

Written by Alfred Kedhammar
"""

import sys

import numpy as np
import pandas as pd

from scilifelab_epps import zika
from scilifelab_epps.utils.udf_tools import is_filled


def pool_fixed_vol(
    currentStep=None,
    lims=None,
    # Volume constraints
    zika_min_vol=0.5,  # 0.5 lowest validated, 0.1 lowest possible
    well_max_vol=180,  # TwinTec96
):
    # Write log header
    log = []
    for e in [
        f"LIMS process {currentStep.id}\n\n=== Volume constraints ===",
        f"Minimum pipetting volume: {zika_min_vol} ul",
        f"Maximum allowed dst well volume: {well_max_vol} ul",
    ]:
        log.append(e)

    # Get transfer vol
    fixed_vol_step_udf = "Transfer Volume for Pooling (uL)"
    fixed_vol = currentStep.udf[fixed_vol_step_udf]
    assert type(fixed_vol) in [int, float], f"'{fixed_vol_step_udf}' must be a number."
    assert zika_min_vol <= fixed_vol <= well_max_vol, (
        f"'{fixed_vol_step_udf}' must be between {zika_min_vol} and {well_max_vol} ul."
    )
    log.append(f"Fixed transfer volume: {fixed_vol} ul")

    # Get pools
    pools = [art for art in currentStep.all_outputs() if art.type == "Analyte"]
    pools.sort(key=lambda pool: pool.name)

    # Get sample dataframe
    to_fetch = {
        # Input sample
        "sample_name": "art_tuple[0]['uri'].name",
        "src_name": "art_tuple[0]['uri'].location[0].name",
        "src_id": "art_tuple[0]['uri'].location[0].id",
        "src_well": "art_tuple[0]['uri'].location[1]",
        # Output pool
        "target_name": "art_tuple[1]['uri'].name",
        "dst_name": "art_tuple[1]['uri'].location[0].name",
        "dst_id": "art_tuple[1]['uri'].location[0].id",
        "dst_well": "art_tuple[1]['uri'].location[1]",
    }
    df_all = zika.utils.fetch_sample_data(currentStep, to_fetch)

    # Define deck, a dictionary mapping plate names to deck positions
    assert len(df_all.src_id.unique()) <= 4, "Only one to four input plates allowed"
    assert len(df_all.dst_id.unique()) == 1, "Only one output plate allowed"
    deck = {}
    deck[df_all.dst_name.unique()[0]] = 3
    available = [2, 4, 1, 5][0 : len(df_all.src_name.unique())]
    # TODO assign deck positions to minimize travel distance
    for plate, pos in zip(df_all.src_name.unique(), available):
        deck[plate] = pos

    # Populate worklist
    df_wl = pd.DataFrame()
    for pool in pools:
        df_pool = df_all[df_all.target_name == pool.name].copy()
        df_pool["transfer_vol"] = fixed_vol
        df_wl = pd.concat([df_wl, df_pool], axis=0)

    # Format worklist
    df_formatted = zika.utils.format_worklist(df_wl.copy(), deck)
    wl_filename, log_filename = zika.utils.get_filenames(
        method_name="pool", pid=currentStep.id
    )

    # Write the output files
    zika.utils.write_worklist(
        df=df_formatted.copy(),
        deck=deck,
        wl_filename=wl_filename,
    )
    zika.utils.write_log(log, log_filename)

    # Upload files
    zika.utils.upload_csv(currentStep, lims, wl_filename)
    zika.utils.upload_log(currentStep, lims, log_filename)

    # Issue warnings, if any
    if any("WARNING" in entry for entry in log):
        sys.stderr.write(
            "CSV-file generated with warnings, please check the Log file\n"
        )
        sys.exit(2)


def pool(
    currentStep=None,
    lims=None,
    # Volume constraints
    zika_min_vol=0.5,  # 0.5 lowest validated, 0.1 lowest possible
    well_dead_vol=5,  # 5 ul generous estimate of dead volume in TwinTec96
    well_max_vol=180,  # TwinTec96
    # Input and output metrics
    udfs={
        # Different steps may use different UDFs in different contexts
        # Here, ambiguity is eliminated within the script
        "target_amt": None,  # Per sample
        "target_vol": None,  # Pool
        "target_conc": None,  # Pool
        "final_amt": None,  # Per sample
        "final_vol": None,  # Pool
        "final_conc": None,  # Pool
    },
):
    """
    Pool samples.

    Input UDFs:
    - The target amount taken per sample (ng) OR target pool concentration (nM)
    - The target final pool volume (ul)

    Calculations:
    1) The inputs are translated to the desired pool concentration and volume.
    2) The desired pool concentration and volume is translated to the target transfer amount for each sample.
    3) The minimum and maximum transferrable amounts are calculated for each sample based on the Zika minimum transfer volume and the accessible sample volume

    Cases:

    - If a sample has very low or negligible concentration
        --> Set concentration to 0.01
        --> # TODO add option to use everything

    - If the highest minimum transferrable amount is lower than the lowest maximum transfer amount, an even pool can be created
        A) An even pool can be created:
            --> Set concentration and volume as close as possible to the specified targets
            --> If necessary, expand the amount taken of all samples and the total volume to maintain the initially specified pool concentration
        B) An even pool cannot be created:
            --> Use the highest minimum transferrable amount as the target
                If this causes volume overflow
                    --> The worklist will not be generated
            - If a sample does not have enough accessible volume to reach the target representation in the pool
                --> Let it be under-represented

    """

    try:
        # Write log header
        log = []
        for e in [
            f"LIMS process {currentStep.id}\n\n=== Volume constraints ===",
            f"Minimum pipetting volume: {zika_min_vol} ul",
            f"Applied dead volume: {well_dead_vol} ul",
            f"Maximum allowed dst well volume: {well_max_vol} ul",
        ]:
            log.append(e)

        pools = [art for art in currentStep.all_outputs() if art.type == "Analyte"]
        pools.sort(key=lambda pool: pool.name)

        # Supplement df with additional info
        to_fetch = {
            # Input sample
            "sample_name": "art_tuple[0]['uri'].name",
            "vol": "art_tuple[0]['uri'].udf['Volume (ul)']",
            "conc": "art_tuple[0]['uri'].udf['Concentration']",
            "conc_units": "art_tuple[0]['uri'].udf['Conc. Units']",
            "src_name": "art_tuple[0]['uri'].location[0].name",
            "src_id": "art_tuple[0]['uri'].location[0].id",
            "src_well": "art_tuple[0]['uri'].location[1]",
            # Output pool
            "target_name": "art_tuple[1]['uri'].name",
            "dst_name": "art_tuple[1]['uri'].location[0].name",
            "dst_id": "art_tuple[1]['uri'].location[0].id",
            "dst_well": "art_tuple[1]['uri'].location[1]",
        }

        for k, v in udfs.items():
            if v:
                to_fetch[k] = f"art_tuple[1]['uri'].udf['{v}']"

        df_all = zika.utils.fetch_sample_data(currentStep, to_fetch)

        # All samples should have accessible volume
        assert all(df_all.vol > well_dead_vol), (
            f"The minimum required source volume is {well_dead_vol} ul"
        )

        assert all(df_all.target_vol <= well_max_vol), (
            f"All target volumes must be at or below {well_max_vol} uL"
        )

        # Adjust for dead volume
        df_all["full_vol"] = df_all.vol.copy()
        df_all.loc[:, "vol"] = df_all.vol - well_dead_vol

        # Define deck, a dictionary mapping plate names to deck positions
        assert len(df_all.src_id.unique()) <= 4, "Only one to four input plates allowed"
        assert len(df_all.dst_id.unique()) == 1, "Only one output plate allowed"
        deck = {}
        deck[df_all.dst_name.unique()[0]] = 3
        available = [2, 4, 1, 5][0 : len(df_all.src_name.unique())]
        # TODO assign deck positions to minimize travel distance
        for plate, pos in zip(df_all.src_name.unique(), available):
            deck[plate] = pos

        # Work through the pools one at a time
        df_wl = pd.DataFrame()
        buffer_vols = {}
        errors = False
        for pool in pools:
            try:
                # === PREPARE CALCULATION INPUTS ===

                # Subset data to current pool
                df_pool = df_all[df_all.target_name == pool.name].copy()

                # Find target parameters, amount and conentration will be either in ng and ng/ul or fmol and nM
                target_pool_vol = df_pool.target_vol.unique()[0]
                if udfs["target_conc"] == "Pool Conc. (nM)":
                    target_pool_conc = df_pool.target_conc.values[0]
                    target_amt_taken = target_pool_conc * target_pool_vol / len(df_pool)
                    amt_unit = "fmol"
                    conc_unit = "nM"
                elif udfs["target_amt"] == "Amount for prep (ng)":
                    target_amt_taken = df_pool.target_amt.unique()[0]
                    target_pool_conc = target_amt_taken * len(df_pool) / target_pool_vol
                    amt_unit = "ng"
                    conc_unit = "ng/ul"
                else:
                    raise AssertionError("Could not make sense of input UDFs")
                assert all(df_all.conc_units == conc_unit), (
                    "Samples and pools have different conc units"
                )

                # Append target parameters to log
                log.append(f"\n\nPooling {len(df_pool)} samples into {pool.name}...")
                log.append("Target parameters:")
                log.append(
                    f" - Amount per sample: {round(target_amt_taken, 2)} {amt_unit}"
                )
                log.append(f" - Pool volume: {round(target_pool_vol, 1)} ul")
                log.append(
                    f" - Pool concentration: {round(target_pool_conc, 2)} {conc_unit}"
                )

                # Set any negative or negligible concentrations to 0.01 and flag in log
                conc_floor = 0.01
                if not df_pool.loc[df_pool.conc < conc_floor, "conc"].empty:
                    neg_conc_sample_names = df_pool.loc[
                        df_pool.conc < conc_floor, "sample_name"
                    ].sort_values()
                    df_pool.loc[df_pool.conc < conc_floor, "conc"] = conc_floor
                    log.append(
                        f"\nWARNING: The following {len(neg_conc_sample_names)} sample(s) fell short of, and will be treated as, "
                        + f"{conc_floor} {conc_unit}: {', '.join(neg_conc_sample_names)}"
                    )
                    log.append(
                        "Low concentration samples will warrant high transfer volumes and may cause pool overflow."
                    )

                # === CALCULATE SAMPLE RANGES ===

                # Calculate the range of transferrable amount for each sample
                df_pool["min_amount"] = zika_min_vol * df_pool.conc
                df_pool["max_amount"] = df_pool.vol * df_pool.conc

                # === CALCULATE POSSIBLE OUTCOMES AND MAKE ADJUSTMENTS ===

                # Isolate highest concentrated sample
                highest_conc_sample = df_pool.sort_values(
                    by="conc", ascending=False
                ).iloc[0]

                # Given the input samples, can an even pool be produced? I.e. is there an overlap in the transfer amount ranges of all samples?
                even_pool_is_possible = max(df_pool.min_amount) < min(
                    df_pool.max_amount
                )

                if even_pool_is_possible:
                    lowest_common_amount = max(df_pool.min_amount)
                    highest_common_amount = min(df_pool.max_amount)

                    # Calculate pool limits given samples
                    pool_min_amt = lowest_common_amount * len(df_pool)
                    pool_min_sample_vol = sum(lowest_common_amount / df_pool.conc)
                    pool_max_sample_vol = sum(highest_common_amount / df_pool.conc)
                    if pool_max_sample_vol < well_max_vol:
                        pool_max_sample_amt = highest_common_amount * len(df_pool)
                    else:
                        # If the max amount corresponds to a volume higher than max, scale it down accordingly
                        pool_max_sample_amt = (
                            highest_common_amount
                            * len(df_pool)
                            * well_max_vol
                            / pool_max_sample_vol
                        )
                    pool_min_conc = pool_min_amt / well_max_vol
                    pool_max_conc = (
                        pool_min_amt / pool_min_sample_vol
                    )  # also equals pool_max_amt / pool_max_sample_vol because amt / vol proportions are the same between samples

                    # Ensure that pool will not overflow
                    if pool_min_sample_vol > well_max_vol:
                        log.append(
                            f"\nERROR: Overflow in {pool.name}. Decrease number of samples or dilute highly concentrated outliers"
                        )
                        log.append(
                            f"Highest concentrated sample: {highest_conc_sample.sample_name} at {round(highest_conc_sample.conc, 2)} {conc_unit}"
                        )
                        log.append(
                            f"Pooling cannot be normalized to less than {round(pool_min_sample_vol, 1)} ul"
                        )

                        errors = True
                        raise zika.utils.VolumeOverflow

                    log.append(
                        "\nAn even pool can be created within the following parameter ranges:"
                    )
                    log.append(
                        f" - Amount per sample {round(lowest_common_amount, 2)} - {round(pool_max_sample_amt / len(df_pool), 2)} {amt_unit}"
                    )
                    log.append(
                        f" - Pool volume {round(pool_min_sample_vol, 1)} - {round(well_max_vol, 1)} ul"
                    )
                    log.append(
                        f" - Pool concentration {round(pool_min_conc, 2)} - {round(pool_max_conc, 2)} {conc_unit}"
                    )

                    # Nudge conc, if necessary
                    if target_pool_conc > pool_max_conc:
                        # Pool conc. has to be decreased from target to the maximum possible, given samples
                        pool_conc = pool_max_conc
                    elif target_pool_conc < pool_min_conc:
                        # Pool conc. has to be increased from target to the minimum possible, given samples
                        pool_conc = pool_min_conc
                    else:
                        # Pool conc. can be set to target
                        pool_conc = target_pool_conc

                    # Nudge vol, if necessary
                    pool_min_vol_given_conc = min(
                        pool_min_amt / pool_conc, well_max_vol
                    )
                    pool_max_vol_given_conc = min(
                        highest_common_amount * len(df_pool) / pool_conc, well_max_vol
                    )
                    if target_pool_vol < pool_min_vol_given_conc:
                        # Pool vol has to be increased from target to minimum possible, given samples
                        pool_vol = pool_min_vol_given_conc
                    elif target_pool_vol > pool_max_vol_given_conc:
                        # Pool vol has to be decreased from target to maximum possible, given samples
                        pool_vol = pool_max_vol_given_conc
                    else:
                        # Pool vol can be set to target
                        pool_vol = target_pool_vol

                    target_transfer_amt = pool_vol * pool_conc / len(df_pool)

                else:
                    # There is no common transfer amount, and sample volumes can NOT be expanded without worsening the even-ness of the pool

                    # Use the minimum transfer amount of the most concentrated sample as the common transfer amount
                    target_transfer_amt = max(df_pool.min_amount)

                    df_low = df_pool[df_pool.max_amount < target_transfer_amt]

                    log.append("\nWARNING: The samples cannot be evenly pooled!")
                    log.append(
                        f"The minimum transfer amount of the highest concentrated sample {highest_conc_sample.sample_name} ({round(highest_conc_sample.conc, 2)} {highest_conc_sample.conc_units}) exceeds the maximum transfer amount of the following samples:"
                    )
                    for i, r in df_low.iterrows():
                        log.append(
                            f"{r.sample_name} ({round(r.conc, 2)} {r.conc_units}, {round(r.vol, 2)} uL accessible volume)"
                        )
                    log.append(
                        "The above samples will be depleted and under-represented in the final pool."
                    )

                    # Calculate pool limits taking into account sample depletion
                    pool_real_min_amt = sum(
                        np.minimum(target_transfer_amt, df_pool.max_amount)
                    )
                    pool_real_min_sample_vol = sum(
                        np.minimum(target_transfer_amt / df_pool.conc, df_pool.vol)
                    )
                    pool_real_max_conc = pool_real_min_amt / pool_real_min_sample_vol
                    pool_real_min_conc = pool_real_min_amt / well_max_vol

                    # Ensure that pool will not overflow
                    if pool_real_min_sample_vol > well_max_vol:
                        log.append(
                            f"\nERROR: Overflow in {pool.name}. Decrease number of samples or dilute highly concentrated outliers"
                        )
                        log.append(
                            f"Highest concentrated sample: {highest_conc_sample.sample_name} at {round(highest_conc_sample.conc, 2)} {conc_unit}"
                        )
                        log.append(
                            f"Pooling cannot be normalized to less than {round(pool_real_min_sample_vol, 1)} ul"
                        )

                        errors = True
                        raise zika.utils.VolumeOverflow

                    log.append(
                        "\nWill try to create a pool that is as even as possible. Accounting for sample depletion, a pool can be created with the following parameter ranges: "
                    )
                    log.append(
                        f" - Target amount per sample {round(target_transfer_amt, 2)}"
                    )
                    log.append(
                        f" - Pool volume {round(pool_real_min_sample_vol, 1)}-{round(well_max_vol, 1)} ul"
                    )
                    log.append(
                        f" - Pool concentration {round(pool_real_min_conc, 2)}-{round(pool_real_max_conc, 2)} {conc_unit}"
                    )

                    # Nudge conc, if necessary
                    # Use the flawed target parameters for comparison and ignore sample depletion
                    if target_pool_conc > pool_real_max_conc:
                        pool_conc = pool_real_max_conc
                    elif target_pool_conc < pool_real_min_conc:
                        pool_conc = pool_real_min_conc
                    else:
                        pool_conc = target_pool_conc

                    # No volume expansion is allowed, so pool volume is set to the minimum, given the conc
                    pool_vol = pool_real_min_sample_vol

            except zika.utils.VolumeOverflow:
                continue

            # === STORE FINAL CALCULATION RESULTS ===

            # Append transfer volumes and corresponding fraction of target conc. for each sample
            df_pool["transfer_vol"] = np.minimum(
                target_transfer_amt / df_pool.conc, df_pool.vol
            )
            df_pool["transfer_amt"] = df_pool.transfer_vol * df_pool.conc
            df_pool["final_amt_fraction"] = round(
                (df_pool.transfer_vol * df_pool.conc / pool_vol)
                / (target_transfer_amt / pool_vol),
                2,
            )

            # Report adjustments in log
            log.append("\nAdjustments:")
            if round(target_pool_conc, 2) != round(pool_conc, 2):
                log.append(
                    f" - WARNING: Target pool concentration is adjusted from {round(target_pool_conc, 2)} --> {round(pool_conc, 2)} {conc_unit}"
                )
            if round(target_pool_vol, 1) != round(pool_vol, 1):
                log.append(
                    f" - WARNING: Target pool volume is adjusted from {round(target_pool_vol, 1)} --> {round(pool_vol, 1)} ul"
                )
            if round(target_pool_conc, 2) == round(pool_conc, 2) and round(
                target_pool_vol, 1
            ) == round(pool_vol, 1):
                log.append("Pooling OK")
            if round(target_transfer_amt, 2) != round(target_amt_taken, 2):
                log.append(
                    f" - INFO: Amount taken per sample is adjusted from {round(target_amt_taken, 2)} --> {round(target_transfer_amt, 2)} {amt_unit}"
                )

            # Calculate and store pool buffer volume
            total_sample_vol = sum(df_pool["transfer_vol"])
            buffer_vol = (
                pool_vol - total_sample_vol
                if pool_vol - total_sample_vol > zika_min_vol
                else 0
            )
            buffer_vols[pool.name] = buffer_vol
            log.append(
                f"\nThe final pool volume is {round(pool_vol, 1)} ul ({round(total_sample_vol, 1)} ul sample + {round(buffer_vol, 1)} ul buffer)"
            )

            # === REPORT DEVIATING SAMPLES ===

            # Report deviating conc samples
            outlier_samples = df_pool[
                np.logical_or(
                    df_pool.final_amt_fraction < 0.995,
                    df_pool.final_amt_fraction > 1.005,
                )
            ][["sample_name", "final_amt_fraction"]].sort_values("sample_name")
            if not outlier_samples.empty:
                log.append(
                    "\nThe following samples deviate from the target representation within the pool:"
                )
                log.append("Sample\tFraction")
                for name, frac in outlier_samples.values:
                    log.append(f" - {name}\t{round(frac, 2)}")

            df_wl = pd.concat([df_wl, df_pool], axis=0)

            # Update UDFs
            pool.udf["Final Volume (uL)"] = float(round(pool_vol, 1))
            if amt_unit == "fmol":
                pool.udf["Pool Conc. (nM)"] = float(round(pool_conc, 2))
            elif amt_unit == "ng":
                if even_pool_is_possible:
                    pool.udf["Amount for prep (ng)"] = float(
                        round(df_pool["transfer_amt"].unique()[0], 2)
                    )
                else:
                    pool.udf["Amount for prep (ng)"] = float(
                        round(target_transfer_amt, 2)
                    )
            pool.put()

        # Get filenames and upload log if errors
        wl_filename, log_filename = zika.utils.get_filenames(
            method_name="pool", pid=currentStep.id
        )
        if errors:
            raise zika.utils.CheckLog(log, log_filename, lims, currentStep)

        # Format worklist
        df_formatted = zika.utils.format_worklist(df_wl.copy(), deck)

        # Comments to attach to the worklist header
        comments = [
            f"This worklist will enact pooling of {len(df_all)} samples",
            "For detailed parameters see the worklist log",
        ]
        for pool in pools:
            if buffer_vols[pool.name] > 0:
                comments.append(
                    f"Add {round(buffer_vols[pool.name], 1)} ul buffer to pool {pool.name} (well {pool.location[1]})"
                )

        # Write the output files
        zika.utils.write_worklist(
            df=df_formatted.copy(),
            deck=deck,
            wl_filename=wl_filename,
            comments=comments,
        )
        zika.utils.write_log(log, log_filename)

        # Upload files
        zika.utils.upload_csv(currentStep, lims, wl_filename)
        zika.utils.upload_log(currentStep, lims, log_filename)

        # Issue warnings, if any
        if any("WARNING" in entry for entry in log):
            sys.stderr.write(
                "CSV-file generated with warnings, please check the Log file\n"
            )
            sys.exit(2)

    except AssertionError as e:
        sys.stderr.write(str(e))
        sys.exit(2)


# =====================================================================================================


def norm(
    # LIMS info
    currentStep=None,
    lims=None,
    # Dilution strategy
    volume_expansion=True,  # For samples that are too concentrated, increase target volume to obtain correct conc
    # Volume constraints
    zika_min_vol=0.5,  # 0.5 lowest validated, 0.1 lowest possible
    well_dead_vol=5,  # 5 ul generous estimate of dead volume in TwinTec96
    well_max_vol=180,  # TwinTec96
    # Input and output metrics
    use_customer_metrics=False,
    allow_multi_plate=False,  # Allow up to 3 source plates; by default only one is allowed
    udfs={
        # Different steps may use different UDFs in different contexts
        # Here, ambiguity is eliminated within the script
        "target_amt": None,
        "target_vol": None,
        "target_conc": None,
        "final_amt": None,
        "final_vol": None,
        "final_conc": None,
    },
):
    """
    Normalize to target amount and volume.

    Cases:
    1) Not enough sample       --> Decrease amount, flag
    2) Enough sample           --> OK
    3) Sample too concentrated --> if volume_expansion:
                                    Increase volume to obtain target concentration
                                   else:
                                    Maintain target volume and allow sample to be above target concentration
    """

    try:
        # Write log header

        log = []
        for e in [
            f"LIMS process {currentStep.id}\n\n=== Dilution strategy ===",
            f"Expand volume to obtain target conc: {volume_expansion}",
            f"Base calculations on user measurements: {use_customer_metrics}",
            "\n=== Volume constraints ===",
            f"Minimum pipetting volume: {zika_min_vol} ul",
            f"Applied dead volume: {well_dead_vol} ul",
            f"Maximum allowed dst well volume: {well_max_vol} ul",
        ]:
            log.append(e)

        # Import output artifacts
        outputs = {
            art.name: art for art in currentStep.all_outputs() if art.type == "Analyte"
        }

        # Assert required UDFs are populated in step
        for output_name, output in outputs.items():
            assert is_filled(output, udfs["target_amt"]), (
                f"UDF '{udfs['target_amt']}' missing for {output.name}"
            )
            assert is_filled(output, udfs["target_vol"]), (
                f"UDF '{udfs['target_vol']}' missing for {output.name}"
            )

        # Fetch sample data

        to_fetch = {
            # Input sample
            "sample_name": "art_tuple[0]['uri'].name",
            "src_name": "art_tuple[0]['uri'].location[0].name",
            "src_id": "art_tuple[0]['uri'].location[0].id",
            "src_well": "art_tuple[0]['uri'].location[1]",
            # Output sample
            "dst_name": "art_tuple[1]['uri'].location[0].name",
            "dst_id": "art_tuple[1]['uri'].location[0].id",
            "dst_well": "art_tuple[1]['uri'].location[1]",
        }

        if use_customer_metrics:
            to_fetch["conc"] = "art_tuple[0]['uri'].samples[0].udf['Customer Conc']"
            to_fetch["vol"] = "art_tuple[0]['uri'].samples[0].udf['Customer Volume']"
        else:
            to_fetch["conc_units"] = "art_tuple[0]['uri'].udf['Conc. Units']"
            to_fetch["conc"] = "art_tuple[0]['uri'].udf['Concentration']"
            to_fetch["vol"] = "art_tuple[0]['uri'].udf['Volume (ul)']"

        for k, v in udfs.items():
            if v:
                to_fetch[k] = f"art_tuple[1]['uri'].udf['{v}']"

        df = zika.utils.fetch_sample_data(currentStep, to_fetch)

        # Cast numeric fields and skip inputs that cannot be normalized.
        # This keeps hard validation for valid samples while allowing controls
        # without source metrics to be ignored instead of crashing the step.
        for numeric_col in ["conc", "vol", "target_amt", "target_vol"]:
            df[numeric_col] = pd.to_numeric(df[numeric_col], errors="coerce")

        invalid_rows = (
            df.conc.isna()
            | df.vol.isna()
            | df.target_amt.isna()
            | df.target_vol.isna()
            | (df.conc <= 0)
            | (df.target_amt <= 0)
            | (df.target_vol <= 0)
            | (df.vol <= well_dead_vol)
        )
        if invalid_rows.any():
            skipped = df.loc[
                invalid_rows, ["sample_name", "conc", "vol", "target_amt", "target_vol"]
            ].copy()
            for _, skipped_row in skipped.iterrows():
                log.append(
                    f"WARNING-SKIPPED: Sample {skipped_row.sample_name} skipped due to missing/insufficient source metrics "
                    f"(conc={skipped_row.conc}, vol={skipped_row.vol} uL, target_amt={skipped_row.target_amt}, target_vol={skipped_row.target_vol} uL; minimum required source volume is {well_dead_vol} uL)"
                )
            df = df.loc[~invalid_rows].copy()
            df = df.reset_index(drop=True)

        assert not df.empty, (
            "No valid samples left after filtering missing/insufficient source metrics"
        )

        conc_unit = "ng/ul" if use_customer_metrics else df.conc_units.iloc[0]
        amt_unit = "ng" if conc_unit == "ng/ul" else "fmol"

        # Assertions
        assert all(df.target_vol <= well_max_vol), (
            f"All target volumes must be at or below {well_max_vol} uL"
        )
        df["full_vol"] = df.vol.copy()
        df.loc[:, "vol"] = df.vol - well_dead_vol

        # Define deck
        assert len(df.dst_id.unique()) == 1, "Only one output plate allowed"
        if allow_multi_plate:
            assert len(df.src_id.unique()) <= 3, (
                "Only one to three input plates allowed"
            )
            deck = {}
            deck[df.dst_name.unique()[0]] = 3
            available = [2, 4, 1, 5][: len(df.src_name.unique())]
            for plate, pos in zip(df.src_name.unique(), available):
                deck[plate] = pos
            deck["buffer_plate"] = next(
                p for p in [4, 1, 5, 6] if p not in deck.values()
            )
        else:
            assert len(df.src_id.unique()) == 1, "Only one input plate allowed"
            deck = {
                df.src_name.unique()[0]: 2,
                df.dst_name.unique()[0]: 3,
                "buffer_plate": 4,
            }

        # Make calculations
        df["target_conc"] = df.target_amt / df.target_vol
        df["min_transfer_amt"] = np.minimum(df.vol, zika_min_vol) * df.conc
        df["max_transfer_amt"] = np.minimum(df.vol, df.target_vol) * df.conc

        # Iterate across samples
        d = {
            "sample_vol": [],
            "buffer_vol": [],
            "tot_vol": [],
            "sample_amt": [],
            "final_conc": [],
        }
        for i, r in df.iterrows():
            log.append(
                f"\n{r.sample_name} (conc {round(r.conc, 2)} {conc_unit}, vol {round(r.vol, 1)} ul)"
            )

            # Cases

            # 1) Not enough sample --> Conc below target
            if round(r.max_transfer_amt, 2) < round(r.target_amt, 2):
                if r.conc < r.target_conc:
                    log.append(
                        "WARNING: Sample concentration is less than target concentration"
                    )
                if r.vol < r.target_vol:
                    log.append("WARNING: Sample is depleted")
                sample_vol = min(r.vol, r.target_vol)
                log.append(f"INFO: Using maximum sample live volume {sample_vol} ul")
                tot_vol = r.target_vol
                buffer_vol = tot_vol - sample_vol

            # 2) Ideal case
            elif (
                round(r.min_transfer_amt, 2)
                <= round(r.target_amt, 2)
                <= round(r.max_transfer_amt, 2)
            ):
                sample_vol = r.target_amt / r.conc
                buffer_vol = r.target_vol - sample_vol
                tot_vol = sample_vol + buffer_vol

            # 3) Sample too concentrated -> Increase final volume if possible
            elif round(r.min_transfer_amt, 2) > round(r.target_amt, 2):
                log.append("WARNING: Sample is too concentrated")

                if volume_expansion:
                    if r.min_transfer_amt / r.target_conc <= well_max_vol:
                        tot_vol = r.min_transfer_amt / r.target_conc
                    else:
                        tot_vol = well_max_vol
                    log.append(
                        f"INFO: Expanding total volume to {round(tot_vol, 1)} ul"
                    )
                    sample_vol = zika_min_vol
                    buffer_vol = tot_vol - sample_vol

                else:
                    log.append(f"INFO: Using minimum sample volume {zika_min_vol} ul")
                    sample_vol = zika_min_vol
                    tot_vol = r.target_vol
                    buffer_vol = tot_vol - sample_vol

            # Adress cases where buffer volume is lower than the minimum transfer amount
            if 0 < buffer_vol < zika_min_vol:
                log.append(
                    f"WARNING: Required buffer volume ({round(buffer_vol, 1)} ul) is less than minimum transfer volume {zika_min_vol} ul"
                )
                log.append("INFO: Omitting buffer")
                tot_vol -= buffer_vol
                buffer_vol = 0

            # Finalize calculations
            final_amt = sample_vol * r.conc
            final_conc = final_amt / tot_vol
            final_conc_frac = final_conc / r.target_conc
            if round(final_conc_frac, 2) > 1:
                log.append("WARNING: Final concentration is above target")
            elif round(final_conc_frac, 2) < 1:
                log.append("WARNING: Final concentration is below target")
            log.append(
                f"--> Diluting {round(sample_vol, 1)} ul ({round(final_amt, 2)} {amt_unit}) to {round(tot_vol, 1)} ul ({round(final_conc, 2)} {conc_unit}, {round(final_conc_frac * 100, 1)}% of target)"
            )

            # Append calculation results to dict

            d["sample_amt"].append(final_amt)
            d["sample_vol"].append(sample_vol)
            d["buffer_vol"].append(buffer_vol)
            d["tot_vol"].append(tot_vol)
            d["final_conc"].append(final_conc)

            # Update UDFs
            op = outputs[r.sample_name]
            op.udf[udfs["final_amt"]] = float(round(final_amt, 2))
            op.udf[udfs["final_vol"]] = float(round(tot_vol, 1))
            if round(final_amt, 2) < round(r.target_amt, 2):
                op.udf[udfs["target_amt"]] = float(round(final_amt, 2))
            op.put()

        log.append("\nDone.\n")

        # Join dict to dataframe
        df = df.join(pd.DataFrame(d))

        # Comments to attach to the worklist header
        wl_comments = []

        # Resolve buffer transfers
        df_buffer, wl_comments = zika.utils.resolve_buffer_transfers(
            df=df.copy(), wl_comments=wl_comments
        )

        # Format worklist
        df_formatted = zika.utils.format_worklist(df_buffer.copy(), deck=deck)
        wl_comments.append(
            f"This worklist will enact normalization of {len(df)} samples. For detailed parameters see the worklist log"
        )

        # Write files

        wl_filename, log_filename = zika.utils.get_filenames(
            method_name="norm", pid=currentStep.id
        )

        zika.utils.write_worklist(
            df=df_formatted.copy(),
            deck=deck,
            wl_filename=wl_filename,
            comments=wl_comments,
        )

        zika.utils.write_log(log, log_filename)

        # Upload files
        zika.utils.upload_csv(currentStep, lims, wl_filename)
        zika.utils.upload_log(currentStep, lims, log_filename)

        # Issue warnings, if any
        if any("WARNING" in entry for entry in log):
            sys.stderr.write(
                "CSV-file generated with warnings, please check the Log file\n"
            )
            sys.exit(2)

    except AssertionError as e:
        sys.stderr.write(str(e))
        sys.exit(2)
