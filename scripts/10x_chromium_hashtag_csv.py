#!/usr/bin/env python

DESC = """This script is intended to run as a project automation,
taking a CSV whose columns are sample names and semi-colon separated
antibody product numbers.

Generated at https://genomics-status.scilifelab.se/10X_chromium_hashtag_csv

Example CSV "P31907_20250403_1334.csv"

    Sample_ID,Antibodies
    P31907_1001,394663;682205;682207
    P31907_1002,394667;394669
    P31907_1003,394669
    P31907_1004,

The script will scan the project files for such a csv, and upon finding it will
write the 2nd column to the matching sample's UDF "Antibodies".
"""

import re
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Project
from genologics.lims import Lims


def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    proj = Project(lims, id=args.proj_id)

    antibody_csv_pattern = re.compile(proj.id + r"_\d{8}_\d{4}\.csv")
    antibodies_udf_name = "Antibodies"

    # Search project files
    matching_files = []
    for f in proj.files:
        f_name = f.original_location
        if re.match(antibody_csv_pattern, f_name):
            matching_files.append(f)

    if len(matching_files) == 0:
        print(
            f"No files matching pattern '{antibody_csv_pattern.pattern}' found in project."
        )
        sys.exit(2)
    elif len(matching_files) > 1:
        print(
            f"Multiple files matching pattern '{antibody_csv_pattern.pattern}' found in project."
        )
        sys.exit(2)

    # Parse and assert csv file
    csv = matching_files[0]
    csv_name = csv.original_location
    csv_contents = lims.get_file_contents(csv.id)
    lines = csv_contents.splitlines()
    header = lines[0]
    assert header == "Sample_ID,Antibodies", (
        f"Unrecognized CSV header, expected 'Sample_ID,Antibodies', got '{header}'"
    )
    rows = lines[1:]

    # Map sample names to ABs string
    sample2abs = {}
    for row in rows:
        sample_name, abs = row.split(",")
        sample2abs[sample_name] = abs

    samples = lims.get_samples(projectlimsid=proj.id)
    assert set(sample2abs.keys()) == set([s.name for s in samples]), (
        "Sample names in project do not match sample names in file."
    )

    for sample in samples:
        if sample2abs[sample.name]:
            sample.udf[antibodies_udf_name] = sample2abs[sample.name]
            sample.put()

    print(
        f"Used file '{csv_name}' to update sample UDF '{antibodies_udf_name}' of {len(samples)} samples."
    )


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument(
        "--proj_id",
        type=str,
        help="Lims ID for Project",
    )
    args = parser.parse_args()

    main(args)
