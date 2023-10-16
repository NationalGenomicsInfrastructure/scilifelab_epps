#!/usr/bin/env python
DESC="""EPP script to fetch and upload Caliper image files for Clarity LIMS.
Searches the directory given by the path argument for filenames matching
a specific pattern ending with:
${INPUT.CONTAINER.PLACEMENT}_${INPUT.NAME}_${INPUT.CONTAINER.LIMSID}_${INPUT.LIMSID}.
This is done for each artifact of type ResultFile, that is of type PerInput. Any
file found matching is copied to the current working directory with a name suffixed with the
output artifact it is connected to. When executed as an EPP, this will cause the
Clarity LIMS EPP wrapper to associate the file with this artifact.

Written by Johannes Alneberg, Science for Life Laboratory, Stockholm, Sweden.
"""

from genologics.lims import Lims
from genologics.entities import Artifact, Process,Container, Sample
from genologics.config import BASEURI,USERNAME,PASSWORD

from epp_utils.nextcloud_webdav import NextCloudClient

from argparse import ArgumentParser
import os
import re
import sys
import logging

from scilifelab_epps.epp import attach_file,EppLogger, unique_check, EmptyError

def main(lims, args, epp_logger):
    p = Process(lims, id=args.pid)

    if args.nextcloud:
        client = NextCloudClient()
        file_list = client.listdir(args.path)
    else:
        if not args.path:
            args.path = os.getcwd()

        file_list = os.listdir(args.path)

    # Find all per input result files
    io = p.input_output_maps
    io_filtered = [x for x in io if x[1]['output-generation-type']=='PerInput']
    io_filtered = [x for x in io_filtered if x[1]['output-type']=='ResultFile']

    artifact_missing_file = []
    artifact_multiple_file = []
    found_files = []



    for input, output in io_filtered:
        i_a = Artifact(lims,id=input['limsid'])
        o_a = Artifact(lims,id=output['limsid'])

        # Input Well, Input Container
        i_w, i_c = i_a.location[1], i_a.location[0]

        # Well is typed without colon in filename:
        i_w = ''.join(i_w.split(':'))


        # Use a regular expression to find the file name given
        # the container and sample. This is all assuming the driver template name ends with:
        # ${INPUT.CONTAINER.PLACEMENT}_${INPUT.NAME}_${INPUT.CONTAINER.LIMSID}_${INPUT.LIMSID}
        # However, names are excluded to improve robustness.
        if args.instrument == "fragment_analyzer":
            info = {'well':o_a.location[1].replace(":", ""),
                'output_artifact_name':o_a.samples[0].name}
            re_str = '.*{well}.*{output_artifact_name}'\
                                       .format(**info)
        else:
            info = {'well':i_w,
                'container_id':i_c.id,
                'input_artifact_id':i_a.id}
            re_str = '.*{well}_.*_.*{container_id}_.*{input_artifact_id}'\
                                       .format(**info)
            logging.info(("Looking for file for artifact id: {input_artifact_id} "
                      "from container with id: {container_id}.").format(**info))

        im_file_r = re.compile(re_str)
        fns = list(filter(im_file_r.match, file_list))

        if len(fns) == 0:
            logging.warning("No image file found for artifact with id {0}".format(i_a.id))
            artifact_missing_file.append(i_a)
        elif len(fns) > 1:
            logging.warning(("Multiple image files found for artifact with id {0}, "
                            "please attach files manually").format(i_a.id))
            artifact_multiple_file.append(i_a)
        else:
            file_path = fns[0]
            file_name = file_path.split("/")[-1]

            found_files.append(file_name)
            logging.info("Found image file {0} for artifact with id {1}".format(file_name, i_a.id))
            fp = os.path.join(args.path, file_name)


            if args.nextcloud:
                client.download_file(file_path, file_name)
                location = attach_file(file_name, o_a)
            else:
                # Attach file to the LIMS
                location = attach_file(fp, o_a)

            logging.debug("Moving {0} to {1}".format(fp,location))

    warning = ""
    if len(artifact_missing_file):
        warning = "Did not find any file for {0} artifact(s). ".format(len(artifact_missing_file))

    if len(artifact_multiple_file):
        warning += "Found multiple files for {0} artifact(s), none of these were uploaded.".format(len(artifact_multiple_file))

    if warning:
       warning = "Warning: " + warning

    abstract = "Uploaded {0} file(s). {1}".format(len(found_files), warning)
    print(abstract, file=sys.stderr) # stderr will be logged and printed in GUI


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument('--pid',
                        help='Lims id for current Process')
    parser.add_argument('--log',
                        help='Log file for runtime info and errors')
    parser.add_argument('--path',
                        help='Path where image files are located')
    parser.add_argument('--instrument',default="caliper",
                        help='instrument deciding the file regex format')
    parser.add_argument('--nextcloud', action='store_true', help="Files are fetched from nextcloud instead of local file system")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()

    with EppLogger(args.log, lims=lims, prepend=True) as epp_logger:
        main(lims, args, epp_logger)
