#!/usr/bin/env python
DESC = """EPP script to summarize sequencing start results to the projects running notes
"""
import datetime
import re
import xml.etree.ElementTree as ET
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from write_notes_to_couchdb import write_note_to_couch

from scilifelab_epps.wrapper import epp_decorator

regex_project_line = re.compile(r"^\[(P\d+)\][\s]*:")
TIMESTAMP: str = datetime.datetime.now().strftime("%y%m%d_%H%M%S")


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    pro = Process(lims, id=args.pid)
    pro_udfs = {}
    samples = {}
    general_comments = []
    project_specific_comments = {}

    for name, value in pro.udf.items():
        pro_udfs[name] = value
    seq_setup = f"{pro_udfs.get('Read 1 Cycles')}_{pro_udfs.get('Index Read 1')}_{pro_udfs.get('Index Read 2')}_{pro_udfs.get('Read 2 Cycles')}"
    inst = f"{pro_udfs.get('Instrument')} {pro_udfs.get('Side')}"

    for line in pro_udfs.get("Comment").splitlines():
        if not line.startswith("//"):
            if regex_project_line.match(line):
                result = regex_project_line.search(line)
                if result.group(1) not in project_specific_comments:
                    project_specific_comments[result.group(1)] = []
                project_specific_comments[result.group(1)].append(line)
            else:
                general_comments.append(line)

    step_xml = ET.fromstring(pro.step.xml())
    date_started = datetime.datetime.fromisoformat(
        step_xml.find("date-started").text
    ).date()

    sample_artifacts = pro.all_outputs(unique=True)
    for art in sample_artifacts:
        if art.type == "Analyte":
            sample_item = {}
            sample_item["projects"] = set()
            for sample in art.samples:
                sample_item["projects"].add(sample.project.id)

            for name, value in art.udf.items():
                sample_item[name] = value

            sample_item["pool"] = art.name
            samples[art.id] = sample_item

    container_name = sample_artifacts[0].container.name
    container_type = sample_artifacts[0].container.type.name
    for well, art in sample_artifacts[0].container.placements.items():
        if art.id in samples:
            samples[art.id]["lane"] = well.split(":")[0]

    note_creation_date = datetime.datetime.now()
    kit_size = ""  # TODO: get kit size
    for sample in samples.values():
        for project in sample["projects"]:
            note_obj = {
                "_id": f"{project}:{datetime.datetime.timestamp(note_creation_date)}",
                "categories": ["Lab"],
                "note_type": "project",
                "parent": project,
                "created_at_utc": note_creation_date.isoformat(),
                "updated_at_utc": note_creation_date.isoformat(),
                "projects": [project],
                "user": pro.technician.name,
                "email": pro.technician.email,
            }
            project_comments = project_specific_comments.get(project, []).join("\n")
            note_obj["note"] = (
                f"Comment from {pro.type.name} ([LIMS]({BASEURI}/clarity/work-details/{pro.id.split('-')[1]}) : \n \
                                **Sequencing started {date_started} ** \n \
                                Pool {sample['pool']} in lane {sample['lane']}, {sample['Loading Conc. (pM)']}pM, {sample['% phiX']}% PhiX, \n \
                                {container_type}-{kit_size} FC = {container_name}, on {inst}, {seq_setup} \n \
                                {project_comments} \n \
                                {general_comments.join('/n')} \n \
                                /{pro.technician.name}"
            )
            write_note_to_couch(project, note_creation_date, note_obj, BASEURI)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    main(args)
