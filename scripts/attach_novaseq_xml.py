#!/usr/bin/env python

import re
import os
import sys
import glob

from argparse import ArgumentParser
from datetime import datetime
from genologics.lims import Lims
from genologics.entities import Process
from genologics.config import BASEURI, USERNAME, PASSWORD

DESC = """EPP for attaching RunInfo.xml and RunParameters.xml from NovaSeq run dir, and copying run parameters from the previous step
Author: Chuan Wang, Science for Life Laboratory, Stockholm, Sweden
"""


def main(lims, args):
    log = []
    content = None
    process = Process(lims, id=args.pid)

    # Copy Read and index parameter from the step "Load to Flowcell (NovaSeq 6000 v2.0)"
    UDF_to_copy = ["Read 1 Cycles", "Read 2 Cycles", "Index Read 1", "Index Read 2"]
    for i in UDF_to_copy:
        if process.parent_processes()[0].udf.get(i):
            process.udf[i] = process.parent_processes()[0].udf[i]
    process.put()

    # Fetch Flowcell ID
    FCID = process.parent_processes()[0].output_containers()[0].name

    for outart in process.all_outputs():
        if outart.type == "ResultFile" and outart.name == "Run Info":
            try:
                lims.upload_new_file(
                    outart,
                    max(
                        glob.glob("/srv/mfs/NovaSeq_data/*{}/RunInfo.xml".format(FCID)),
                        key=os.path.getctime,
                    ),
                )
            except:
                raise RuntimeError("No RunInfo.xml Found!")
        elif outart.type == "ResultFile" and outart.name == "Run Parameters":
            try:
                lims.upload_new_file(
                    outart,
                    max(
                        glob.glob(
                            "/srv/mfs/NovaSeq_data/*{}/RunParameters.xml".format(FCID)
                        ),
                        key=os.path.getctime,
                    ),
                )
            except:
                raise RuntimeError("No RunParameters.xml Found!")


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
