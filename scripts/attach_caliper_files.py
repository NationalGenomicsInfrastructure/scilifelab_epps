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
        # However, names are excluded to improve robustness.
        if args.instrument == "fragment_analyzer":
            info = {
                "well": o_a.location[1].replace(":", ""),
                "output_artifact_name": o_a.samples[0].name,
            }
            re_str = ".*{well}.*{output_artifact_name}".format(**info)
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
            logging.warning(
                f"No image file found for artifact with id {i_a.id}"
            )
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
            logging.info(
                f"Found image file {fn} for artifact with id {i_a.id}"
            )
            fp = os.path.join(args.path, fn)

            # Attach file to the LIMS
            location = attach_file(fp, o_a)
            logging.debug(f"Moving {fp} to {location}")

    warning = ""
    if len(artifact_missing_file):
        warning = "Did not find any file for {} artifact(s). ".format(
            len(artifact_missing_file)
        )

    if len(artifact_multiple_file):
        warning += "Found multiple files for {} artifact(s), none of these were uploaded.".format(
            len(artifact_multiple_file)
        )

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
