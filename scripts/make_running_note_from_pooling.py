#!/usr/bin/env python
DESC = """EPP used to create running notes from the pooling step"""

import os
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process, Project
from genologics.lims import Lims
from write_notes_to_couchdb import write_note_to_couch

from scilifelab_epps.epp import attach_file


def main(lims, args):
    p = Process(lims, id=args.pid)
    log = []
    datamap = {}
    wsname = None
    user_email = p.technician.email
    for art in p.all_inputs():
        if len(art.samples) != 1:
            log.append(f"Warning : artifact {art.id} has more than one sample")
        for sample in art.samples:
            # take care of lamda DNA
            if sample.project:
                if sample.project.id not in datamap:
                    datamap[sample.project.id] = [sample.name]
                else:
                    datamap[sample.project.id].append(sample.name)

    for art in p.all_outputs():
        try:
            wsname = art.location[0].name
            break
        except:
            pass

    for pid in datamap:
        pj = Project(lims, id=pid)
        if len(datamap[pid]) > 1:
            rnt = f"{len(datamap[pid])} samples planned for {wsname}"
        else:
            rnt = f"{len(datamap[pid])} sample planned for {wsname}"

        running_note = {}
        running_note["note"] = rnt
        running_note["email"] = user_email
        running_note["categories"] = ["Workset"]
        running_note["note_type"] = "project"
        write_note_to_couch(pid, running_note, lims.get_uri())
        log.append(
            f"Updated project {pid} : {pj.name}, {len(datamap[pid])} samples in this workset"
        )

    with open("EPP_Notes.log", "w") as flog:
        flog.write("\n".join(log))
    for out in p.all_outputs():
        # attach the log file
        if out.name == "RNotes Log":
            attach_file(os.path.join(os.getcwd(), "EPP_Notes.log"), out)

    sys.stderr.write(f"Updated {len(list(datamap.keys()))} projects successfully")


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
