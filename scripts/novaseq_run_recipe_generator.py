#!/usr/bin/env python

import os
import sys
import json

from argparse import ArgumentParser
from datetime import datetime
from genologics.lims import Lims
from genologics.entities import Process
from genologics.config import BASEURI, USERNAME, PASSWORD


DESC = """EPP used to create run recipe for NovaSeq sequencing
Author: Chuan Wang, Science for Life Laboratory, Stockholm, Sweden
"""


def main(lims, args):
    log = []
    thisyear = datetime.now().year
    process = Process(lims, id=args.pid)

    # Fetch FC ID
    for out in process.all_outputs():
        try:
            if out.type == "Analyte":
                fc_name = out.location[0].name
        except Exception as e:
            log.append(str(e))

    # Fetch required run step UDFs
    run_mode = process.udf.get("Run Mode", "")
    sample_loading_type = (
        "NovaSeqXp"
        if process.udf.get("Loading Workflow Type") == "NovaSeq Xp"
        else "NovaSeqStandard"
    )
    workflow_type = process.udf.get("Workflow Type", "").replace(" ", "")
    librarytube_ID = process.udf.get("Library Tube Barcode", "")
    paired_end = True if process.udf.get("Paired End") == "True" else False
    read1 = process.udf.get("Read 1 Cycles", 0)
    read2 = process.udf.get("Read 2 Cycles", 0)
    index_read1 = process.udf.get("Index Read 1", 0)
    index_read2 = process.udf.get("Index Read 2", 0)
    output_folder = "\\\\172.16.1.6\\novaseqdata\\Runs\\"
    attachment = "\\\\172.16.1.6\\samplesheets\\novaseq\\{}\\\\{}.csv".format(
        thisyear, fc_name
    )
    basespace_mode = process.udf.get("BaseSpace Sequence Hub Configuration")
    if basespace_mode == "Not Used":
        use_basespace = False
    else:
        use_basespace = True
    use_custom_read1_primer = process.udf.get("Use Custom Read 1 Primer")
    use_custom_read2_primer = process.udf.get("Use Custom Read 2 Primer")
    use_custom_index_read1_primer = process.udf.get("Use Custom Index Read 1 Primer")

    # Prepare json file
    output = {
        "run_name": fc_name,
        "run_mode": workflow_type,
        "workflow_type": workflow_type,
        "sample_loading_type": sample_loading_type,
        "librarytube_ID": librarytube_ID,
        "flowcell_ID": fc_name,
        "rehyb": False,
        "paired_end": paired_end,
        "read1": read1,
        "read2": read2,
        "index_read1": index_read1,
        "index_read2": index_read2,
        "output_folder": output_folder,
        "attachment": attachment,
        "use_basespace": use_basespace,
        "basespace_mode": basespace_mode,
        "use_custom_read1_primer": use_custom_read1_primer,
        "use_custom_read2_primer": use_custom_read2_primer,
        "use_custom_index_read1_primer": use_custom_index_read1_primer,
    }

    # Write json file
    if os.path.exists("/srv/mfs/NovaSeq_data/gls_recipe_novaseq/"):
        try:
            with open(
                "/srv/mfs/NovaSeq_data/gls_recipe_novaseq/{}.json".format(fc_name), "w"
            ) as sf:
                json.dump(output, sf, separators=(",", ":"))
        except Exception as e:
            log.append(str(e))

    for out in process.all_outputs():
        if out.name == "Run Recipe":
            ss_art = out
        elif out.name == "Run Recipe Log":
            log_id = out.id

    with open("{}.json".format(fc_name), "w", 0o664) as sf:
        json.dump(output, sf, separators=(",", ":"))
    os.chmod("{}.json".format(fc_name), 0o664)
    for f in ss_art.files:
        lims.request_session.delete(f.uri)
    lims.upload_new_file(ss_art, "{}.json".format(fc_name))

    # Write log
    if log:
        with open("{}_{}_Error.log".format(log_id, fc_name), "w") as f:
            f.write("\n".join(log))
        sys.stderr.write("Errors were met, check the log.")
        sys.exit(1)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
