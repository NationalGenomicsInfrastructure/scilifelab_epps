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


class IndexPair(TypedDict):
    idx1: str
    idx2: str


class WellData(TypedDict):
    count: int
    labels: list[IndexPair]
    label_length: set[int]


def verify_samplename(sample_name: str, proj_id: str) -> list[str]:
    message = []
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
    """Validate a single reagent label and check index distance."""
    message = []

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
            idxs = (
                TENX_SINGLE_PAT.findall(reagent_label)
                or TENX_DUAL_PAT.findall(reagent_label)
                or SMARTSEQ_PAT.findall(reagent_label)
            )
            # We'll skip TENX and SMARTSEQ indexes for now
            if not idxs:
                idxs_matches = IDX_PAT.findall(reagent_label)
                if not idxs_matches:
                    message.append(
                        f"INDEX FORMAT ERROR: Sample {sample_id} with index '{reagent_label}' could not be parsed"
                    )
                    return message

                idxs = idxs_matches[0]
                curr_idx: IndexPair = {
                    "idx1": idxs[0],
                    "idx2": idxs[1] if len(idxs) > 1 else "",
                }
                data[pool]["labels"].append(curr_idx)

                # Check index distance from previous samples in pool
                for prev_idx in data[pool]["labels"][:-1]:
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


# Verify sample IDs
def verify_samples(lims: Lims, project: Project) -> list[str]:
    """Validate project sample IDs for format, count, and sequence."""
    message = []
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
            data.setdefault(pool, {"count": 0, "labels": [], "label_length": set()})
            data[pool]["count"] += 1

            reagent_label = sample.artifact.reagent_labels[0].strip()
            if not reagent_label:
                message.append(f"INDEX WARNING: Sample {sample.name} has no label")
            else:
                # Track label length and warn if mismatched
                if (
                    data[pool]["label_length"]
                    and len(reagent_label) not in data[pool]["label_length"]
                ):
                    common_index = next(iter(data[pool]["label_length"]))
                    message.append(
                        f"LABEL LENGTH WARNING: Multiple label lengths noticed in pool {pool} for Sample {sample_id}, length {len(reagent_label)} is different from {common_index}"
                    )
                else:
                    data[pool]["label_length"].add(len(reagent_label))

                # Validate the reagent label
                message.extend(
                    validate_reagent_label(reagent_label, sample_id, pool, data)
                )

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
    messages = verify_samples(lims, project)

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

    with open("index_checker.log", "w") as logContext:
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
