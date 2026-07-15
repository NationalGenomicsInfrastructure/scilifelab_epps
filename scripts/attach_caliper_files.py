#!/usr/bin/env python
DESC = """EPP script to fetch and upload Caliper image files for Clarity LIMS.
Searches the directory given by the path argument for filenames matching
a specific pattern ending with:
${INPUT.CONTAINER.PLACEMENT}_${INPUT.NAME}_${INPUT.CONTAINER.LIMSID}_${INPUT.LIMSID}.
This is done for each artifact of type ResultFile, that is of type PerInput. Any
file found matching is copied to the current working directory with a name suffixed with the
output artifact it is connected to. When executed as an EPP, this will cause the
Clarity LIMS EPP wrapper to associate the file with this artifact.

Written by Johannes Alneberg, Science for Life Laboratory, Stockholm, Sweden.
"""

import logging
import os
import re
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Artifact, Process
from genologics.lims import Lims

from scilifelab_epps.epp import EppLogger, attach_file


def normalize_sample_name(sample_name):
    """Normalize sample names so '-' and '_' variants map to the same value."""
    return re.sub(r"[^A-Za-z0-9]", "", sample_name or "").lower()


def build_sample_suffix_regex(sample_name):
    """Build a regex that matches sample name at end of filename stem."""
    tokens = [
        re.escape(token)
        for token in re.split(r"[^A-Za-z0-9]+", sample_name or "")
        if token
    ]
    if not tokens:
        return None

    core = r"[^A-Za-z0-9]*".join(tokens)
    pattern = rf"(?<![A-Za-z0-9]){core}$"
    return re.compile(pattern, re.IGNORECASE)


def normalize_well(well):
    """Normalize well names so A:1, A1 and A01 are considered equal."""
    compact_well = (well or "").replace(":", "")
    match = re.match(r"^([A-Za-z])(\d+)$", compact_well)
    if not match:
        return compact_well.upper()

    row, column = match.groups()
    return f"{row.upper()}{int(column)}"


def parse_fragment_analyzer_filename(filename):
    """Parse FA filename broadly by finding final sample token and preceding well."""
    stem, _ = os.path.splitext(filename)
    # Keep prefix matching broad; sample is expected at end of stem.
    sample_suffix_regex = getattr(
        parse_fragment_analyzer_filename, "sample_suffix_regex", None
    )
    if sample_suffix_regex is None:
        return None

    sample_match = sample_suffix_regex.search(stem)
    if not sample_match:
        return None

    sample_from_file = sample_match.group(0)
    prefix = stem[: sample_match.start()].rstrip()
    # Well is expected immediately before sample token.
    well_match = re.search(r"([A-Za-z]:?\d{1,2})\s*$", prefix)
    if not well_match:
        return None

    well_from_file = well_match.group(1)
    return {
        "well": normalize_well(well_from_file),
        "sample": normalize_sample_name(sample_from_file),
    }


def main(lims, args, epp_logger):
    p = Process(lims, id=args.pid)

    if not args.path:
        args.path = os.getcwd()

    file_list = os.listdir(args.path)

    # Find all per input result files
    io = p.input_output_maps
    io_filtered = [x for x in io if x[1]["output-generation-type"] == "PerInput"]
    io_filtered = [x for x in io_filtered if x[1]["output-type"] == "ResultFile"]

    artifact_missing_file = []
    artifact_multiple_file = []
    found_files = []

    for input, output in io_filtered:
        i_a = Artifact(lims, id=input["limsid"])
        o_a = Artifact(lims, id=output["limsid"])

        # Input Well, Input Container
        i_w, i_c = i_a.location[1], i_a.location[0]

        # Well is typed without colon in filename:
        i_w = "".join(i_w.split(":"))

        # Use a reguluar expression to find the file name given
        # the container and sample. This is all assuming the driver template name ends with:
        # ${INPUT.CONTAINER.PLACEMENT}_${INPUT.NAME}_${INPUT.CONTAINER.LIMSID}_${INPUT.LIMSID}
        # For Fragment Analyzer, filenames are matched using sample-at-end and preceding well.
        if args.instrument == "fragment_analyzer":
            expected_well = normalize_well(o_a.location[1])
            # Use the output artifact name (e.g. "Fragment Analyzer P40753P1-A1") to
            # derive the sample token written in filenames by the FA instrument.
            # This is stable even if the underlying LIMS sample name is renamed,
            # because frag_an_driver_gen.py uses output.name.split("Fragment Analyzer ")[1].
            fa_prefix = "Fragment Analyzer "
            if o_a.name.startswith(fa_prefix):
                sample_token = o_a.name[len(fa_prefix):]
            else:
                sample_token = o_a.name
            expected_sample = normalize_sample_name(sample_token)
            parse_fragment_analyzer_filename.sample_suffix_regex = (
                build_sample_suffix_regex(sample_token)
            )

            logging.info(
                "Looking for FA file for artifact id: %s, well: %s, token: %s",
                i_a.id,
                expected_well,
                sample_token,
            )

            fns = []
            for fn in file_list:
                parsed = parse_fragment_analyzer_filename(fn)
                if not parsed:
                    continue

                has_name_match = expected_sample and parsed["sample"] == expected_sample
                has_well_match = parsed["well"] == expected_well

                if has_name_match and has_well_match:
                    fns.append(fn)
        else:
            info = {"well": i_w, "container_id": i_c.id, "input_artifact_id": i_a.id}
            re_str = ".*{well}_.*_.*{container_id}_.*{input_artifact_id}".format(**info)
            logging.info(
                (
                    "Looking for file for artifact id: {input_artifact_id} "
                    "from container with id: {container_id}."
                ).format(**info)
            )
            im_file_r = re.compile(re_str)
            fns = list(filter(im_file_r.match, file_list))

        if len(fns) == 0:
            logging.warning(f"No image file found for artifact with id {i_a.id}")
            artifact_missing_file.append(i_a)
        elif len(fns) > 1:
            logging.warning(
                f"Multiple image files found for artifact with id {i_a.id}, "
                "please attach files manually"
            )
            artifact_multiple_file.append(i_a)
        else:
            fn = fns[0]
            found_files.append(fn)
            logging.info(f"Found image file {fn} for artifact with id {i_a.id}")
            fp = os.path.join(args.path, fn)

            # Attach file to the LIMS
            location = attach_file(fp, o_a)
            logging.debug(f"Moving {fp} to {location}")

    warning = ""
    if len(artifact_missing_file):
        warning = (
            f"Did not find any file for {len(artifact_missing_file)} artifact(s). "
        )

    if len(artifact_multiple_file):
        warning += f"Found multiple files for {len(artifact_multiple_file)} artifact(s), none of these were uploaded."

    if warning:
        warning = "Warning: " + warning

    abstract = f"Uploaded {len(found_files)} file(s). {warning}"
    print(abstract, file=sys.stderr)  # stderr will be logged and printed in GUI


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    parser.add_argument("--log", help="Log file for runtime info and errors")
    parser.add_argument("--path", help="Path where image files are located")
    parser.add_argument(
        "--instrument",
        default="caliper",
        help="instrument deciding the file regex format",
    )
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()

    with EppLogger(args.log, lims=lims, prepend=True) as epp_logger:
        main(lims, args, epp_logger)
