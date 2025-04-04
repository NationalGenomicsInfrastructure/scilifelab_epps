#!/usr/bin/env python
DESC = """EPP script to summarize sequencing start results to the projects running notes
"""
import datetime
import logging
import re
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from write_notes_to_couchdb import write_note_to_couch

from scilifelab_epps.utils.get_epp_user import get_epp_user
from scilifelab_epps.wrapper import epp_decorator

regex_projectid_line = re.compile(r"^(P\d+)[\s]*:")
# regex_projectname_line = re.compile(r"^([A-Za-z]+\.[A-Za-z]+_\d{2}_\d{2})[\s]*:")
TIMESTAMP: str = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
"""
Example running note

Comment from Load to Flowcell (NovaSeqXPlus) v1.0 (LIMS) :
Sequencing started 2025-03-31
Pool 'POOLNAME' in lane 1, xpM, y% PhiX,
Pool 'POOLNAME' in lane 2, xpM, y% PhiX,
1.5B FC=22XDFCRT, on Inst A, 1_3_4_2

/Technician Name
"""


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    pro = Process(lims, id=args.pid)
    pro_udfs = {}
    pools = {}
    general_comments = []
    project_specific_comments = {}
    projects = {}

    for name, value in pro.udf.items():
        pro_udfs[name] = value
    # Use .get() for non-mandatory UDFs
    seq_setup = f"{pro_udfs['Read 1 Cycles']}_{pro_udfs['Index Read 1']}_{pro_udfs.get('Index Read 2', 'x')}_{pro_udfs.get('Read 2 Cycles', 'x')}"
    inst = f"{pro_udfs['Instrument']} {pro_udfs['Side']}"

    date_started = datetime.datetime.fromisoformat(pro.step.date_started).date()

    pool_artifacts = pro.analytes()[0]
    for pool_artifact in pool_artifacts:
        pool_obj = {}
        for name, value in pool_artifact.udf.items():
            pool_obj[name] = value

        pool_obj["pool_name"] = pool_artifact.name
        pools[pool_artifact.id] = pool_obj
        for sample in pool_artifact.samples:
            if sample.project.id not in projects:
                projects[sample.project.id] = {pool_artifact.id}
            else:
                projects[sample.project.id].add(pool_artifact.id)

    # The default value for the "Comments" step field in LIMS contains lines that start with '//' which should be filtered out if they are not removed
    # They could also be used to save info in the step that does not have to be in Genstat
    for line in pro_udfs.get("Comments").splitlines():
        if not line.startswith("//"):
            if regex_projectid_line.match(line):
                result = regex_projectid_line.search(line)
                if result.group(1) not in projects:
                    logging.warning(
                        f"Project {result.group(1)} not in list of projects for this run. Skipping comment."
                    )
                if result.group(1) not in project_specific_comments:
                    project_specific_comments[result.group(1)] = []
                project_specific_comments[result.group(1)].append(
                    line[result.end() :].strip()
                )
            else:
                general_comments.append(line.strip())

    an_analyte_container = pro.output_containers()[0]
    container_name = an_analyte_container.name
    container_type = an_analyte_container.type.name
    for well, pool_artifact in an_analyte_container.placements.items():
        if pool_artifact.id in pools:
            pools[pool_artifact.id]["lane"] = well.split(":")[0]

    note_creation_date = datetime.datetime.now()
    epp_initiator = get_epp_user(lims, pro.id)
    for project in projects:
        project_comments = "\n".join(project_specific_comments.get(project, []))
        pool_text = ""
        for pool_id in projects[project]:
            pool = pools[pool_id]
            pool_text += f"Pool '{pool['pool_name']}' in lane {pool['lane']}, {pool['Loading Conc. (pM)']}pM, {pool['% phiX']}% PhiX, \n"
        note = (
            f"Comment from {pro.type.name} ([LIMS]({BASEURI}/clarity/work-details/{pro.id.split('-')[1]})) : \n"
            f"**Sequencing started {date_started} ** by {pro.technician.name}\n"
            f"{pool_text}"
            f"{container_type} FC={container_name}, on {inst}, {seq_setup} \n"
            f"{project_comments} \n"
            f"{'/n'.join(general_comments)} \n"
            f"/{epp_initiator.name}"
        )
        note_obj = {
            "_id": f"{project}:{datetime.datetime.timestamp(note_creation_date)}",
            "categories": ["Lab"],
            "note_type": "project",
            "parent": project,
            "created_at_utc": note_creation_date.isoformat(),
            "updated_at_utc": note_creation_date.isoformat(),
            "projects": [project],
            "user": epp_initiator.name,
            "email": epp_initiator.email,
            "note": note,
        }
        write_note_to_couch(project, note_creation_date, note_obj, BASEURI)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", type=str, help="Lims id for current Process")
    parser.add_argument("--log", type=str, help="Which log file slot to use")
    args = parser.parse_args()

    main(args)
