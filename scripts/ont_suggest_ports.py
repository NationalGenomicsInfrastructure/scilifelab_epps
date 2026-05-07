#!/usr/bin/env python


import re
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from ont_send_reloading_info_to_db import get_ONT_db

DESC = """ Script for EPP "Suggest PromethION ports".
Use StatusDB to find the least used ports and populate the positions UDFs with them.
"""


def main(lims, args):
    try:
        currentStep = Process(lims, id=args.pid)
        outputs = [op for op in currentStep.all_outputs() if op.type == "Analyte"]

        # Get database
        client, db_name = get_ONT_db()

        # Instantiate dict for counting port usage
        ports = {}
        for c in list("123"):
            for r in list("ABCDEFGH"):
                ports[c + r] = 0

        pattern = re.compile(
            r"/\d{8}_\d{4}_([1-8][A-H])_"
        )  # Matches start of run name, capturing position as a group

        view = client.post_view(
            db=db_name,
            ddoc="info",
            view="all_stats",
            include_docs=False,
        ).get_result()
        # Count port usage
        for row in view["rows"]:
            if re.search(pattern, row["value"]["TACA_run_path"]):
                position = re.search(pattern, row["value"]["TACA_run_path"]).groups()[0]
                ports[position] += 1

        # Sort ports (a sort of port sort, if you will)
        ports_list = list(ports.items())
        ports_list.sort(key=lambda x: x[1])

        # Collect which ports are already specified in UDFs
        ports_used = []
        for output in outputs:
            try:
                if output.udf["ONT flow cell position"] != "None":
                    assert output.udf["ONT flow cell position"] in ports.keys(), (
                        f"{output.udf['ONT flow cell position']} is not a valid position"
                    )
                    ports_used.append(output.udf["ONT flow cell position"])
            except KeyError:
                continue

        # Populate non-specified port UDFs with least used ports
        for output in outputs:
            # If port is already specified, skip to next output
            try:
                if output.udf["ONT flow cell position"] != "None":
                    continue
                else:
                    pass
            except KeyError:
                pass

            # Find the next non-used port
            for port_tuple in ports_list:
                if port_tuple[0] in ports_used:
                    continue
                else:
                    break

            output.udf["ONT flow cell position"] = port_tuple[0]
            output.put()
            ports_used.append(port_tuple[0])

        # Print ports to stdout, starting with the least used
        message = f"Listing ports, from least to most used, down-prioritizing column 3: {', '.join([port[0] for port in ports_list])}"
        sys.stdout.write(message)

    except AssertionError as e:
        sys.stderr.write(str(e))
        sys.exit(2)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
