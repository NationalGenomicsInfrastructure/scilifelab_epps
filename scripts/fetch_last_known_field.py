#!/usr/bin/env python
import logging
from argparse import ArgumentParser
from datetime import datetime as dt

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Artifact, Process
from genologics.lims import Lims

from scilifelab_epps.utils import udf_tools
from scilifelab_epps.wrapper import epp_decorator

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    """This script will get the name of an artifact UDF from a master step field,
    and for every sample artifact in the current step:

    - Use API calls to recursively back-trace the sample history using
        input-output links until it finds an artifact with the specified UDF
    - Copy the value of the specified UDF from the found artifact to the
        artifact of the current step

    Example use-case:
    - For Nanopore libraries in the Aggregate QC step of the Library Validation protocol,
        fetch the last recorded artifact UDF "Size (bp)" from the library prep for all samples.
    """
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    # Get the name of the target UDF from the step field
    target_udf = process.udf.get(args.step_udf, None)
    assert target_udf is not None and target_udf != "None", (
        f"No target UDF supplied from step field '{args.step_udf}'"
    )

    # Check whether process has output artifacts, not the case for e.g. QC steps
    no_outputs: bool = True if udf_tools.get_art_tuples(process) == [] else False

    # Load input artifacts
    arts_in: list[Artifact] = [
        art for art in process.all_inputs() if art.type == "Analyte"
    ]

    # Find target output artifacts, if any
    if no_outputs:
        logging.info("Step has no output artifacts. Assigning to input artifacts.")
    else:
        art_tuples: list[tuple[dict]] = process.input_output_maps
        art_in2out: dict[str:Artifact] = {
            i["uri"].id: o["uri"]
            for i, o in art_tuples
            if i["uri"].type == "Analyte" and o["uri"].type == "Analyte"
        }

    steps_used = []
    for art_in in arts_in:
        if no_outputs:
            target_artifact = art_in
        else:
            target_artifact = art_in2out[art_in.id]
        logging.info(
            f"Looking for last recorded UDF '{target_udf}' of {'input' if no_outputs else 'output'} artifact '{target_artifact.name}'..."
        )
        udf_value, traceback = udf_tools.fetch_last(
            target_art=target_artifact,
            target_udfs=target_udf,
            log_traceback=True,
            return_traceback=True,
            on_fail=None,
        )

        steps_used.append(
            f"'{traceback[-1]['Artifact']['Parent Step']['Name']}'"
            + f" ({traceback[-1]['Artifact']['Parent Step']['ID']})"
        )

        if udf_value is not None:
            target_artifact.udf[target_udf] = udf_value
            target_artifact.put()
            logging.info(
                f"Updated UDF '{target_udf}' for {'input' if no_outputs else 'output'} artifact '{target_artifact.name}' to '{udf_value}'"
            )
        else:
            logging.warning(
                f"Could not traceback UDF '{target_udf}' for {'input' if no_outputs else 'output'} artifact '{target_artifact.name}'"
            )
            logging.info(f"Traceback:\n{traceback}")

        # Write to stdout for the green banner
        print(
            f"UDF '{target_udf}' pulled from steps: {' ,'.join(set(steps_used))}. Please double check the values."
        )


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser()
    parser.add_argument(
        "--pid",
        required=True,
        type=str,
        help="Lims ID for current Process.",
    )
    parser.add_argument(
        "--log",
        required=True,
        type=str,
        help="Which file slot to use for the script log.",
    )
    parser.add_argument(
        "--step_udf",
        required=True,
        type=str,
        help="The name of the step UDF listing the target artifact UDF.",
    )
    args = parser.parse_args()

    main(args)
