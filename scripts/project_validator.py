#!/usr/bin/env python

import re
import smtplib
import sys
from argparse import ArgumentParser
from email.message import Message
from email.mime.text import MIMEText
from typing import TypedDict

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Project
from genologics.lims import Lims

from data.Chromium_10X_indexes import Chromium_10X_indexes
from scilifelab_epps.utils.get_epp_user import get_epp_user

DESC = """EPP used to validate a project including checking sample name format, index format and index distance in library pool.
"""

# Pre-compile regexes in global scope:
NGISAMPLE_PAT = re.compile("P[0-9]+_[0-9]+")
INDEX_PAT = re.compile(
    r"^([ATGC]{6,24}(-[ATGC]{6,12})?|SI-[A-Z0-9]{2,4}-[A-Z]\d{1,2}|NB\d{2}|BC\d{2}|NoIndex)$"
)
IDX_PAT = re.compile("([ATCG]{4,}N*)-?([ATCG]*)")
VALIDBASES_PAT = re.compile(r"^[ATCG\-]+$")
TENX_SINGLE_PAT = re.compile("SI-(?:GA|NA)-[A-H][1-9][0-2]?")
TENX_DUAL_PAT = re.compile("SI-(?:TT|NT|NN|TN|TS)-[A-H][1-9][0-2]?")
SMARTSEQ_PAT = re.compile("SMARTSEQ[1-9]?-[1-9][0-9]?[A-P]")
compl = {"A": "T", "C": "G", "G": "C", "T": "A"}


class IndexPair(TypedDict):
    idx1: str
    idx2: str


class WellData(TypedDict):
    count: int
    labels: dict[str, IndexPair]  # Maps sample_id to resolved IndexPair


def verify_samplename(sample_name: str, proj_id: str) -> list[str]:
    message: list[str] = []
    if not NGISAMPLE_PAT.findall(sample_name):
        message.append(f"SAMPLE NAME WARNING: Bad sample name format {sample_name}")
    else:
        if sample_name.split("_")[0] != proj_id:
            message.append(
                f"SAMPLE NAME WARNING: Sample name {sample_name} does not match project ID {proj_id}"
            )
    return message


def get_index_format_error(index: str) -> str | None:
    if TENX_SINGLE_PAT.fullmatch(index):
        return None
    if TENX_DUAL_PAT.fullmatch(index):
        return None
    if SMARTSEQ_PAT.fullmatch(index):
        return None

    if not IDX_PAT.fullmatch(index):
        return "does not match known index patterns"
    parts = index.split("-")
    if len(parts) > 2:
        return "too many parts (expected single or dual index)"
    if any(part == "" for part in parts):
        return "empty part around '-'"

    if not VALIDBASES_PAT.fullmatch(index):
        return "contains invalid characters (allowed: A, T, C, G, -)"

    return None


def my_distance(idx_a: str, idx_b: str) -> int:
    diffs = 0
    short = min((idx_a, idx_b), key=len)
    lon = idx_a if short == idx_b else idx_b
    for i, c in enumerate(short):
        if c != lon[i]:
            diffs += 1
    return diffs


def validate_reagent_label(
    reagent_label: str, sample_id: str, pool: str, data: dict[str, WellData]
) -> list[str]:
    """Validate a single reagent label and check index distance. Returns messages."""
    message: list[str] = []
    curr_idx: IndexPair | None = None

    if reagent_label == "NoIndex":
        if data[pool]["count"] > 1:
            message.append(
                f"NOINDEX ERROR: Pool {pool} has NoIndex but contains more than one sample"
            )
    else:
        reason = get_index_format_error(reagent_label)
        if reason:
            message.append(
                f"INDEX FORMAT ERROR: Sample {sample_id} with index '{reagent_label}' has a bad format: {reason}"
            )
        else:
            is_tenx_index = TENX_SINGLE_PAT.findall(
                reagent_label
            ) or TENX_DUAL_PAT.findall(reagent_label)
            is_smartseq_index = SMARTSEQ_PAT.findall(reagent_label)
            if is_tenx_index:
                if TENX_SINGLE_PAT.findall(reagent_label):
                    message.append(
                        f"INDEX FORMAT WARNING: Sample {sample_id} with index '{reagent_label}' is a TENX single index, skipping detailed checks"
                    )
                    return message
                else:
                    idx_1 = Chromium_10X_indexes[reagent_label][0].replace(",", "")
                    idx_2 = "".join(
                        reversed(
                            [
                                compl.get(b, b)
                                for b in Chromium_10X_indexes[reagent_label][1]
                                .replace(",", "")
                                .upper()
                            ]
                        )
                    )
            # skipping checks for SMARTSEQ indexes for now
            elif is_smartseq_index:
                message.append(
                    f"INDEX FORMAT WARNING: Sample {sample_id} with index '{reagent_label}' is a SMARTSEQ index, skipping detailed checks"
                )
                return message
            else:
                idxs_matches = IDX_PAT.findall(reagent_label)
                if not idxs_matches:
                    message.append(
                        f"INDEX FORMAT ERROR: Sample {sample_id} with index '{reagent_label}' could not be parsed"
                    )
                    return message

                idxs = idxs_matches[0]
                idx_1 = idxs[0]
                idx_2 = idxs[1] if len(idxs) > 1 else ""

            if curr_idx is None:
                curr_idx = {
                    "idx1": idx_1,
                    "idx2": idx_2,
                }
                data[pool]["labels"][sample_id] = curr_idx

                # Check index distance from previous samples in pool
                for prev_sample_id, prev_idx in data[pool]["labels"].items():
                    if prev_sample_id == sample_id:
                        continue  # Skip self-comparison
                    dist = my_distance(prev_idx["idx1"], curr_idx["idx1"])
                    if prev_idx.get("idx2", "") and curr_idx.get("idx2", ""):
                        dist += my_distance(prev_idx["idx2"], curr_idx["idx2"])

                    if dist < 2:
                        idx_a = f"{prev_idx.get('idx1', '')}-{prev_idx.get('idx2', '')}"
                        idx_b = f"{curr_idx.get('idx1', '')}-{curr_idx.get('idx2', '')}"
                        if dist == 0:
                            message.append(
                                f"INDEX COLLISION ERROR: Index {idx_a} and Index {idx_b} (for sample {sample_id}) in pool {pool}"
                            )
                        else:
                            message.append(
                                f"SIMILAR INDEX WARNING: Index {idx_a} and Index {idx_b} (for sample {sample_id}) in pool {pool}"
                            )

    return message


def validate_plate_sequences(plates: dict[str, list[int]]) -> list[str]:
    """Validate sample numbering within each plate."""
    message = []

    for plate_digit in sorted(plates.keys()):
        plate_suffixes = sorted(plates[plate_digit])
        first_suffix_num = plate_suffixes[0]

        # Compute expected first suffix based on the number of digits
        num_digits = len(str(first_suffix_num))
        expected_first = int(f"{plate_digit}{'0' * (num_digits - 2)}1")
        if first_suffix_num != expected_first:
            message.append(
                f"SAMPLE SEQUENCE WARNING: Plate {plate_digit} first "
                f"sample should be _{expected_first}, but got _{first_suffix_num}"
            )

        # Check gaps within plate
        for curr, next_val in zip(plate_suffixes, plate_suffixes[1:]):
            if next_val - curr != 1:
                message.append(
                    f"SAMPLE SEQUENCE WARNING: Gap detected in plate "
                    f"{plate_digit} numbering between {curr} and "
                    f"{next_val}. Verify missing samples in the "
                    f"uploaded CSV file."
                )

    return message


def validate_pool_label_lengths(pool: str, pool_data: WellData) -> list[str]:
    """Validate that all resolved indices in a pool have the same combined length."""
    message: list[str] = []

    if not pool_data["labels"]:
        return message

    # Calculate lengths for all samples
    length_map: dict[int, list[str]] = {}  # Maps length -> [sample_ids]
    for sample_id, index_pair in pool_data["labels"].items():
        resolved_length = len(index_pair["idx1"]) + len(index_pair["idx2"])
        if resolved_length not in length_map:
            length_map[resolved_length] = []
        length_map[resolved_length].append(sample_id)

    # If only one length, all is good
    if len(length_map) == 1:
        return message

    # Find the most common length
    most_common_length = max(length_map, key=lambda k: len(length_map[k]))

    # Report any indices with different lengths
    for resolved_length, samples in length_map.items():
        if resolved_length != most_common_length:
            for sample_id in samples:
                index_pair = pool_data["labels"][sample_id]
                idx_str = f"{index_pair.get('idx1', '')}-{index_pair.get('idx2', '')}"
                message.append(
                    f"LABEL LENGTH WARNING: Pool {pool}, Sample {sample_id} has resolved index length {resolved_length} ({idx_str}), "
                    f"but majority of samples have length {most_common_length}"
                )

    return message


# Verify sample IDs
def verify_samples(lims: Lims, project: Project) -> list[str]:
    """Validate project sample IDs for format, count, and sequence."""
    message: list[str] = []
    data: dict[str, WellData] = {}
    samples = lims.get_samples(projectname=project.name)
    if not samples:
        message.append("SAMPLE COUNT WARNING: No samples found for the project")
        return message

    # Group samples by plate digit (first digit of suffix: _1XXX, _2XXX, etc.)
    plates: dict[str, list[int]] = {}

    # Validate sample name format and collect data
    for sample in sorted(samples, key=lambda s: s.name):
        sample_id = sample.name
        customer_name = sample.udf.get("Customer Name")
        pool = sample.udf.get("Pooling", "")
        if not customer_name:
            message.append(
                f"SAMPLE NAME WARNING: Sample {sample_id} has no customer name"
            )

        # Validate format and prefix match
        message.extend(verify_samplename(sample_id, project.id))
        sample_suffix = sample_id.split("_")[1] if "_" in sample_id else None
        if sample_suffix:
            try:
                plate_digit = sample_suffix[0]
                suffix_num = int(sample_suffix)
            except ValueError:
                message.append(
                    f"SAMPLE SEQUENCE WARNING: Sample {sample_id} has a non-numeric suffix {sample_suffix}"
                )
                continue

            plates.setdefault(plate_digit, []).append(suffix_num)

        if not sample.artifact.reagent_labels:
            message.append(f"INDEX WARNING: Sample {sample.name} has no label")
        else:
            data.setdefault(pool, {"count": 0, "labels": {}})
            data[pool]["count"] += 1

            reagent_label = sample.artifact.reagent_labels[0].strip()
            if not reagent_label:
                message.append(f"INDEX WARNING: Sample {sample.name} has no label")
            else:
                # Validate the reagent label and get the resolved index
                validation_msgs = validate_reagent_label(
                    reagent_label, sample_id, pool, data
                )
                message.extend(validation_msgs)

    # Validate label lengths in each pool
    for pool, pool_data in data.items():
        message.extend(validate_pool_label_lengths(pool, pool_data))

    # Validate sample numbering within plates
    message.extend(validate_plate_sequences(plates))

    return message


def email_responsible(
    message: str,
    resp_email: str,
    subject: str,
) -> None:
    msg: Message
    body = "Samplesheet validation Errors: \n" + message
    body += "\n\n--\nThis is an automatically generated error notification"
    msg = MIMEText(body)
    msg["Subject"] = subject

    msg["From"] = "Lims_monitor"
    msg["To"] = resp_email

    with smtplib.SMTP("localhost") as s:
        s.sendmail("genologics-lims@scilifelab.se", msg["To"], msg.as_string())


def main(lims: Lims, pid: str, auto: bool) -> None:
    messages = []
    project = Project(lims, id=pid)
    # Get all samples in the project
    # Validate sample IDs
    if project.udf.get("Library construction method") == "Finished library (by user)":
        messages.append(
            f"Running checks for Project {pid}: Project is marked as 'Finished library (by user)'"
        )
        messages = verify_samples(lims, project)
    else:
        print(f"Project {pid}: Project is not a user library, skipping sample checks")

    if messages:
        resp_email = get_epp_user(lims, project_id=pid).email
        if auto:
            if not resp_email:
                print(
                    "**Errors exist in the samplesheet: **\n"
                    "Email with the error could not be sent as no email address was found for the EPP user.\n"
                    + "\n".join(messages),
                    file=sys.stderr,
                )
            else:
                print(
                    "**Errors exist in the samplesheet: **\n"
                    f"Email with the error has been sent to {resp_email}.\n"
                )
                email_responsible(
                    message="\n".join(messages),
                    resp_email=resp_email,
                    subject=f"[Error] Project {pid} failed sample sheet validation",
                )
        else:
            sys.stderr.write("; ".join(messages))
    else:
        print("No issue detected with indexes or placement")

    with open(args.log, "w") as logContext:
        logContext.write("\n".join(messages))
    # Throw red warning message when it is not automatically run
    if not auto and not messages:
        sys.exit(0)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Project ID for current Project")
    parser.add_argument(
        "--log",
        dest="log",
        default="index_checker.log",
        help=("File name for standard log file, for runtime information and problems."),
    )
    parser.add_argument(
        "--auto",
        action="store_true",
        help=("Used when the script is running automatically in LIMS."),
    )
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args.pid, args.auto)
