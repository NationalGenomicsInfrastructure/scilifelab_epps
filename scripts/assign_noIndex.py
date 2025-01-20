#!/usr/bin/env python
import logging
from argparse import ArgumentParser
from datetime import datetime as dt
from xml.etree import ElementTree as ET

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from scilifelab_epps.wrapper import epp_decorator

DESC = """Script to automatically assign the 'noIndex' reagent label
to all unlabeled input or output samples of a step.

The development of this script was made necessary by the need to pass
unlabeled ONT samples through a LIMS demultiplexing step.

Running the script at the start of a demultiplexing step will unfortunately
not work as intended, since LIMS performs the demultiplexing magic
inherent to the step prior to running any EPPs. Running the script on the
step prior to the demux step will work as intended.
"""

TIMESTAMP: str = dt.now().strftime("%y%m%d_%H%M%S")


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    # Set up LIMS
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    # Determine whether to process input or output analytes
    if args.label == "input":
        arts_list = process.all_inputs()
    elif args.label == "output":
        arts_list = process.all_outputs()
    else:
        raise ValueError(f"Invalid value '{args.label}' for argument 'label'")

    # Filter for unlabeled analytes
    unlabeled_analytes = []
    for analyte in arts_list:
        if not analyte.type == "Analyte":
            continue
        if analyte.reagent_labels:
            logging.info(
                f"{args.label.capitalize()} '{analyte.name}' is already labeled ({analyte.reagent_labels}), skipping."
            )
            continue
        unlabeled_analytes.append(analyte)
    logging.info(f"Found {len(unlabeled_analytes)} unlabeled {args.label} analytes.")

    xml_element_noIndex = ET.Element("reagent-label", name="NoIndex")
    failed_analytes = []
    for analyte in unlabeled_analytes:
        try:
            analyte.root.append(xml_element_noIndex)
            analyte.put()
            logging.info(
                f"Assigned 'noIndex' reagent label to {args.label} '{analyte.name}'"
            )
        except Exception as e:
            logging.error(e)
            failed_analytes.append(analyte)

    if failed_analytes:
        logging.error(
            f"Failed to assign 'noIndex' reagent label to the following {args.label}s: {', '.join(art.name for art in failed_analytes)}"
        )


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", type=str, help="Lims ID for current Process")
    parser.add_argument("--log", type=str, help="Which log file slot to use")
    parser.add_argument("--label", type=str, help="Either 'input' or 'output'")

    args = parser.parse_args()

    main(args)
