#!/usr/bin/env python

import json
import logging
import os
import re
import shutil
from argparse import ArgumentParser, Namespace
from datetime import datetime as dt
from zipfile import ZipFile

import pandas as pd
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from Levenshtein import distance

from data.Chromium_10X_indexes import Chromium_10X_indexes
from scilifelab_epps.epp import get_pool_sample_label_mapping, upload_file
from scilifelab_epps.wrapper import epp_decorator

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")

# Pre-compile regexes in global scope:
IDX_PAT = re.compile("([ATCG]{4,}N*)-?([ATCG]*)")
TENX_SINGLE_PAT = re.compile("SI-(?:GA|NA)-[A-H][1-9][0-2]?")
TENX_DUAL_PAT = re.compile("SI-(?:TT|NT|NN|TN|TS)-[A-H][1-9][0-2]?")
SMARTSEQ_PAT = re.compile("SMARTSEQ[1-9]?-[1-9][0-9]?[A-P]")

# Set up Element PhiX control sets, keys are options in LIMS dropdown UDF
PHIX_SETS = {
    "PhiX Control Library, Adept": {
        "nickname": "PhiX_Adept",
        "indices": [
            ("ATGTCGCTAG", "CTAGCTCGTA"),
            ("CACAGATCGT", "ACGAGAGTCT"),
            ("GCACATAGTC", "GACTACTAGC"),
            ("TGTGTCGACA", "TGTCTGACAG"),
        ],
    },
    "Cloudbreak PhiX Control Library, Elevate": {
        "nickname": "PhiX_Elevate",
        "indices": [
            ("ACGTGTAGC", "GCTAGTGCA"),
            ("CACATGCTG", "AGACACTGT"),
            ("GTACACGAT", "CTCGTACAG"),
            ("TGTGCATCA", "TAGTCGATC"),
        ],
    },
    "Cloudbreak Freestyle PhiX Control, Third Party": {
        "nickname": "PhiX_Third",
        "indices": [
            ("ATGTCGCTAG", "CTAGCTCGTA"),
            ("CACAGATCGT", "ACGAGAGTCT"),
            ("GCACATAGTC", "GACTACTAGC"),
            ("TGTGTCGACA", "TGTCTGACAG"),
        ],
    },
}

# Load SS3 indexes
SMARTSEQ3_indexes_json = (
    "/opt/gls/clarity/users/glsai/repos/scilifelab_epps/data/SMARTSEQ3_indexes.json"
)
with open(SMARTSEQ3_indexes_json) as file:
    SMARTSEQ3_INDEXES = json.loads(file.read())


def revcomp(seq: str) -> str:
    """Reverse-complement a DNA string."""
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def idxs_from_label(label: str) -> list[str | tuple[str, str]]:
    """From a LIMS reagent label, return list whose elements are
    single indices or tuples of dual index pairs.
    """

    # Initialize result
    idxs: list[str | tuple[str, str]] = []

    # Expand 10X single indexes
    if TENX_SINGLE_PAT.findall(label):
        match = TENX_SINGLE_PAT.findall(label)[0]
        for tenXidx in Chromium_10X_indexes[match]:
            idxs.append(tenXidx)
    # Case of 10X dual indexes
    elif TENX_DUAL_PAT.findall(label):
        match = TENX_DUAL_PAT.findall(label)[0]
        i7_idx = Chromium_10X_indexes[match][0]
        i5_idx = Chromium_10X_indexes[match][1]
        idxs.append((i7_idx, revcomp(i5_idx)))
    # Case of SS3 indexes
    elif SMARTSEQ_PAT.findall(label):
        match = SMARTSEQ_PAT.findall(label)[0]
        for i7_idx in SMARTSEQ3_INDEXES[match][0]:
            for i5_idx in SMARTSEQ3_INDEXES[match][1]:
                idxs.append((i7_idx, revcomp(i5_idx)))
    # NoIndex cases
    elif label.replace(",", "").upper() == "NOINDEX" or (
        label.replace(",", "").upper() == ""
    ):
        raise AssertionError("NoIndex cases not allowed.")
    # Ordinary indexes
    elif IDX_PAT.findall(label):
        match = IDX_PAT.findall(label)[0]
        if "-" in match:
            idx1, idx2 = match.split("-")
            idxs.append((idx1, revcomp(idx2)))
        else:
            idx1 = match
            idxs.append(idx1)
    else:
        raise AssertionError(f"Could not parse index from '{label}'.")
    return idxs


def get_flowcell_id(process: Process) -> str:
    """Get the Element flowcell ID from the process."""
    flowcell_ids = [
        op.container.name for op in process.all_outputs() if op.type == "Analyte"
    ]

    assert len(set(flowcell_ids)) == 1, "Expected one flowcell ID."
    flowcell_id = flowcell_ids[0]

    if "-" in flowcell_id:
        logging.warning(
            f"Container name {flowcell_id} contains a dash, did you forget to set the name of the LIMS container to the flowcell ID?"
        )

    return flowcell_id


def dict_to_manifest_col(d: dict) -> str:
    """Turn a list of key-value pairs into a string fitting into a manifest column."""
    for k, v in d.items():
        for char in [",", ":", " "]:
            assert char not in k, f"Character '{char}' not allowed in manifest columns."
            assert char not in v, f"Character '{char}' not allowed in manifest columns."

    s = " ".join([f"{k}:{v}" for k, v in d.items()])

    return s


def get_manifests(
    process: Process, manifest_root_name: str
) -> list[tuple[str, str | None]]:
    """Generate multiple manifests, grouping samples by index multiplicity and length,
    adding PhiX controls of appropriate lengths as needed.
    """

    # Assert output analytes loaded on flowcell
    arts_out = [op for op in process.all_outputs() if op.type == "Analyte"]
    assert len(arts_out) == 1 or len(arts_out) == 2, (
        "Expected one or two output analytes."
    )

    # Assert lanes
    lanes = [art_out.location[1].split(":")[0] for art_out in arts_out]
    assert set(lanes) == {"1"} or set(lanes) == {
        "1",
        "2",
    }, "Expected a single-lane or dual-lane flowcell."

    # Iterate over pool / lane
    sample_rows = []
    for pool, lane in sorted(zip(arts_out, lanes), key=lambda x: x[1]):
        # Get sample-label linkage via database
        sample2label: dict[str, str] = get_pool_sample_label_mapping(pool)
        assert len(set(pool.reagent_labels)) == len(pool.reagent_labels), (
            f"Detected non-unique reagent labels in lane {lane}"
        )

        # Record PhiX UDFs for each output artifact
        phix_loaded: float = pool.udf.get("% phiX", 0)
        phix_set_name: str = pool.udf.get("Element PhiX Set", "")
        if phix_loaded != 0:
            assert phix_set_name != "", "PhiX controls loaded but no kit specified."
        else:
            assert phix_set_name == "", "PhiX controls specified but not loaded."

        # Collect rows for each sample
        for sample in pool.samples:
            # Include project name and sequencing setup
            if sample.project:
                project = sample.project.name.replace(".", "__").replace(",", "")
                seq_setup = sample.project.udf.get("Sequencing setup", "0-0")
                user_library = (
                    True
                    if sample.project.udf["Library construction method"]
                    == "Finished library (by user)"
                    else False
                )
            else:
                project = "Control"
                seq_setup = "0-0"

            # Add row(s), depending on index type
            lims_label = sample2label[sample.name]
            for idx in idxs_from_label(lims_label):
                row = {}
                row["SampleName"] = sample.name
                if isinstance(idx, tuple):
                    row["Index1"], row["Index2"] = idx
                    # Special cases to reverse-complement index2
                    if not user_library or (
                        user_library
                        and (
                            TENX_DUAL_PAT.findall(lims_label)
                            or SMARTSEQ_PAT.findall(lims_label)
                        )
                    ):
                        logging.info(f"Reverse-complementing index2 of {sample.name}.")
                        row["Index2"] = revcomp(row["Index2"])
                else:
                    row["Index1"] = idx
                    # Assume long idx2 from recipe + no idx2 from label means idx2 is UMI
                    if int(process.udf.get("Index Read 2", 0)) > 12:
                        row["Index2"] = "N" * int(process.udf["Index Read 2"])
                    else:
                        row["Index2"] = ""
                row["Lane"] = lane
                row["Project"] = project
                row["Recipe"] = seq_setup
                row["phix_loaded"] = phix_loaded
                row["phix_set_name"] = phix_set_name
                row["lims_label"] = lims_label

                # Add special case settings
                row_settings = {}
                if TENX_SINGLE_PAT.findall(lims_label):
                    # For 10X 8-mer single indexes (e.g. SI-NA-A1) it is usually required that
                    #  index 1 sequences shall be written as a separate FastQ file (I1).
                    # In this case we need the additional option I1Fastq,TRUE.
                    row_settings["I1Fastq"] = "True"
                row["settings"] = dict_to_manifest_col(row_settings)

                sample_rows.append(row)

    # Compile sample dataframe
    df_samples = pd.DataFrame(sample_rows)

    # Add PhiX controls
    df_samples_and_controls = df_samples.copy()
    for lane, group in df_samples.groupby(["Lane"]):
        if group["phix_loaded"].any():
            phix_set_name = group["phix_set_name"].iloc[0]
            phix_set = PHIX_SETS[phix_set_name]

            # Add row for each PhiX index pair
            for phix_idx_pair in phix_set["indices"]:
                row = {}
                row["SampleName"] = phix_set["nickname"]
                row["Index1"] = phix_idx_pair[0]
                row["Index2"] = phix_idx_pair[1]
                row["Lane"] = group["Lane"].iloc[0]
                row["Project"] = "Control"
                row["Recipe"] = "0-0"

                df_samples_and_controls = pd.concat(
                    [df_samples_and_controls, pd.DataFrame([row])], ignore_index=True
                )

    df_samples_and_controls.sort_values(by=["Lane", "SampleName"], inplace=True)
    df_samples_and_controls.reset_index(drop=True, inplace=True)

    # Check for index collision per lane, across samples and PhiX
    for lane, group in df_samples_and_controls.groupby("Lane"):
        rows_to_check = group.to_dict(orient="records")
        check_distances(rows_to_check)

    # Start building manifests
    manifests: list[tuple[str, str | None]] = []
    for manifest_type in ["untrimmed", "trimmed", "phix", "empty"]:
        manifest_name, manifest_contents = make_manifest(
            df_samples_and_controls,
            process,
            manifest_root_name,
            manifest_type,
        )
        manifests.append((manifest_name, manifest_contents))

    return manifests


def make_manifest(
    df_samples_and_controls: pd.DataFrame,
    process: Process,
    manifest_root_name: str,
    manifest_type: str,
) -> tuple[str, str | None]:
    logging.info(f"Building {manifest_type} manifest...")

    # Get the index cycles from the step fields
    idx1_cycles = int(process.udf.get("Index Read 1"))
    idx2_cycles = int(process.udf.get("Index Read 2"))

    # Make copy of input df and subset columns to include in manifest
    df = df_samples_and_controls[
        [
            "SampleName",
            "Index1",
            "Index2",
            "Lane",
            "Project",
            "Recipe",
            "lims_label",
            "settings",
        ]
    ].copy()

    file_name = f"{manifest_root_name}_{manifest_type}.csv"

    runValues_section = "\n".join(
        [
            "[RUNVALUES]",
            "KeyName, Value",
            f"lims_step_name, {process.type.name}",
            f"lims_step_id, {process.id}",
            f"manifest_file, {file_name}",
        ]
    )

    # Build the [SAMPLES] section of the manifest, depending on the manifest type.
    if manifest_type == "untrimmed":
        samples_section = f"[SAMPLES]\n{df.to_csv(index=None, header=True)}"

    elif manifest_type in ["trimmed", "phix"]:
        if manifest_type == "phix":
            # Subset to PhiX controls
            df = df[df["Project"] == "Control"]

        # Make index lengths conform to number of cycles
        for idx, cycles in zip(["Index1", "Index2"], [idx1_cycles, idx2_cycles]):
            # If there are any indexes shorter than the number of cycles
            if not df[df[idx].apply(len) < cycles].empty:
                for row in df[df[idx].apply(len) < cycles].to_dict(orient="records"):
                    logging.error(
                        f"'{row['SampleName']}' has {idx} '{row[idx]}' of length {len(row[idx])} shorter than {cycles} cycles."
                    )
                logging.error(
                    f"Could not generate {manifest_type} manifest because indexes are shorter than the number of index cycles. Skipping."
                )
                return (file_name, None)
            # If there are any indexes longer than the number of cycles
            if not df[df[idx].apply(len) > cycles].empty:
                # For each one, log how it's trimmed
                for row in df[df[idx].apply(len) > cycles].to_dict(orient="records"):
                    logging.info(
                        f"Trimming '{row['SampleName']}' {idx} '{row[idx]}' of length {len(row[idx])} to {cycles} cycles."
                    )
            df[idx] = df[idx].apply(lambda x: x[:cycles])

        samples_section = f"[SAMPLES]\n{df.to_csv(index=None, header=True)}"

    elif manifest_type == "empty":
        samples_section = ""

    else:
        raise AssertionError("Invalid manifest type.")

    settings_section = "\n".join(
        [
            "[SETTINGS]",
            "SettingName, Value",
        ]
    )

    # Customize mismatch thresholds, if necessary
    if manifest_type not in ["untrimmed", "empty"]:
        try:
            logging.info(
                f"Getting custom mismatch thresholds for {manifest_type} manifest..."
            )
            i1_mismatch, i2_mismatch = get_custom_mistmatch_thresholds(df)
        except AssertionError as e:
            logging.error(e)
            logging.error(
                f"Could not generate {manifest_type} manifest without index collisions. Skipping."
            )
            return (file_name, None)

        settings_section += "\n" + "\n".join(
            [
                f"I1MismatchThreshold, {i1_mismatch}",
                f"I2MismatchThreshold, {i2_mismatch}",
            ]
        )

    # Write manifest
    manifest_contents = "\n\n".join(
        [runValues_section, settings_section, samples_section]
    )

    return (file_name, manifest_contents)


def get_custom_mistmatch_thresholds(df: pd.DataFrame) -> tuple[int, int]:
    # Defaults, according to Element documentation
    i1MismatchThreshold = 1
    i2MismatchThreshold = 1

    # Collect distances
    idx1_dists = []
    idx2_dists = []
    total_dists = []
    # Iterate across all sample pairings per lane
    for lane in df["Lane"].unique():
        df_lane = df[df["Lane"] == lane]
        df_lane.reset_index(drop=True, inplace=True)
        for i in range(0, len(df_lane)):
            for j in range(i + 1, len(df_lane)):
                # TODO skip NNN
                idx1_dist = distance(df_lane["Index1"][i], df_lane["Index1"][j])
                idx2_dist = distance(df_lane["Index2"][i], df_lane["Index2"][j])

                # Collect distances between all sample pairings on index and index-pair level
                idx1_dists.append(idx1_dist)
                idx2_dists.append(idx2_dist)
                total_dists.append(idx1_dist + idx2_dist)

    if min(total_dists) == 0:
        raise AssertionError("Total index distance of 0 detected.")
    if min(idx1_dists) <= 2:
        logging.warning(
            "Minimum distance between Index1 sequences is at or below 2. Reducing allowed mismatches from 1 to 0."
        )
        i1MismatchThreshold = 0
    if min(idx2_dists) <= 2:
        logging.warning(
            "Minimum distance between Index2 sequences is at or below 2. Reducing allowed mismatches from 1 to 0."
        )
        i2MismatchThreshold = 0

    return (i1MismatchThreshold, i2MismatchThreshold)


def check_distances(rows: list[dict], threshold=2) -> None:
    """Iterator function to check index distances between all pairs of samples."""
    for i in range(len(rows)):
        row = rows[i]

        for row_comp in rows[i + 1 :]:
            check_pair_distance(row, row_comp, threshold=threshold)


def check_pair_distance(row, row_comp, check_flips: bool = False, threshold: int = 3):
    """Distance check between two index pairs.

    row                     dict   manifest row of sample A
    row_comp                dict   manifest row of sample B
    check_flips             bool   check all reverse-complement combinations
    threshold               int    trigger warning for distances at or below this value

    """

    if check_flips:
        flips: list[tuple[int, str, str]] = []
        for s1i1, s1i1_name in zip(
            [row["Index1"], revcomp(row["Index1"])],
            ["Index1", "Index1_rc"],
        ):
            for s1i2, s1i2_name in zip(
                [row["Index2"], revcomp(row["Index2"])],
                ["Index2", "Index2_rc"],
            ):
                for s2i1, s2i1_name in zip(
                    [row_comp["Index1"], revcomp(row_comp["Index1"])],
                    ["Index1", "Index1_rc"],
                ):
                    for s2i2, s2i2_name in zip(
                        [row_comp["Index2"], revcomp(row_comp["Index2"])],
                        ["Index2", "Index2_rc"],
                    ):
                        flips.append(
                            (
                                distance(s1i1, s2i1) + distance(s1i2, s2i2),
                                f"{s1i1}-{s1i2} {s2i1}-{s2i2}",
                                f"{s1i1_name}-{s1i2_name} {s2i1_name}-{s2i2_name}",
                            )
                        )
        dist, compared_seqs, flip_conf = min(flips, key=lambda x: x[0])

    else:
        dist = distance(
            row["Index1"] + row["Index2"], row_comp["Index1"] + row_comp["Index2"]
        )
        compared_seqs = (
            f"{row['Index1']}-{row['Index2']} {row_comp['Index1']}-{row_comp['Index2']}"
        )

    if dist <= threshold:
        # Build a warning message for the pair
        warning_lines = [
            f"Hamming distance {dist} between {row['SampleName']} and {row_comp['SampleName']}"
        ]
        # If the distance is derived from a flip, show the original and the flipped conformation
        if check_flips:
            warning_lines.append(
                f"Given: {row['Index1']}-{row['Index2']} <-> {row_comp['Index1']}-{row_comp['Index2']}"
            )
            warning_lines.append(f"Distance: {dist} when flipped to {flip_conf}")
        # If the index lengths are equal, add a simple visual representation
        if len(row["Index1"]) + len(row["Index2"]) == len(row_comp["Index1"]) + len(
            row_comp["Index2"]
        ):
            warning_lines.append(show_match(*compared_seqs.split()))

        warning = "\n".join(warning_lines)
        logging.warning(warning)

        # For identical collisions, kill the process
        if dist == 0:
            raise AssertionError("Identical indices detected.")


def show_match(seq1: str, seq2: str) -> str:
    """Visualize base-by-base match between sequences of equal length."""

    assert len(seq1) == len(seq2)

    m = ""
    for seq1_base, seq2_base in zip(seq1, seq2):
        if seq1_base == seq2_base:
            m += "|"
        else:
            m += "X"

    lines = "\n".join([seq1, m, seq2])
    return lines


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args: Namespace):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    # Create manifest root name
    flowcell_id = get_flowcell_id(process)
    manifest_root_name = f"AVITI_run_manifest_{flowcell_id}_{process.id}_{TIMESTAMP}_{process.technician.name.replace(' ', '')}"

    # Create manifest(s)
    manifests: list[tuple[str, str | None]] = get_manifests(process, manifest_root_name)

    # Write and zip manifest(s)
    zip_file = f"{manifest_root_name}.zip"
    with ZipFile(zip_file, "w") as zip_stream:
        for file, content in manifests:
            if content:
                open(file, "w").write(content)
                zip_stream.write(file)
                os.remove(file)
            else:
                logging.warning(f"Not writing {file} due to missing contents.")

    # Upload manifest(s)
    logging.info("Uploading run manifest to LIMS...")
    upload_file(
        zip_file,
        args.file,
        process,
        lims,
    )

    # Move manifest(s)
    logging.info("Moving run manifest to ngi-nas-ns...")
    try:
        dst = f"/srv/ngi-nas-ns/samplesheets/Aviti/{dt.now().year}"
        if not os.path.exists(dst):
            logging.info(f"Happy new year! Creating {dst}")
            os.mkdir(dst)
        shutil.copyfile(
            zip_file,
            f"{dst}/{zip_file}",
        )
        os.remove(zip_file)
    except:
        logging.error("Failed to move run manifest to ngi-nas-ns.", exc_info=True)
    else:
        logging.info("Run manifest moved to ngi-nas-ns.")


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser()
    parser.add_argument(
        "--pid",
        required=True,
        type=str,
        help="Lims ID for current Process.",
    )
    parser.add_argument(
        "--log",
        required=True,
        type=str,
        help="Which file slot to use for the script log.",
    )
    parser.add_argument(
        "--file",
        required=True,
        type=str,
        help="Which file slot to use for the run manifest.",
    )
    args = parser.parse_args()

    main(args)
