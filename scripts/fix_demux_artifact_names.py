#!/usr/bin/env python
DESC = """
This script addresses a specific LIMS bug in which the demux artifacts of a
demultiplexing step are associated to multiple samples, instead of a single sample.

The bug is described in deviation #312 on AM.

I.e. a list of demux artifacts for a pool which should be e.g.
    P34304_1004 (FASTQ reads)
    P34304_1005 (FASTQ reads)
    P34304_1006 (FASTQ reads)
instead becomes e.g.
    P34304_1004 + P34304_1005 + P34304_1006 (FASTQ reads)
    P34304_1005 + P34304_1006 (FASTQ reads)
    P34304_1005 + P34304_1006 (FASTQ reads)

The identity of each demux artifact is dubious, but they do appear to maintain
the correct reagent label.

This script
 - Iterates through each pool in the process
    - Creates an SQL-derived sample-label mapping for the pool contents
    - Iterates through each demux artifact associated with the pool
        - Uses the sample-label mapping to find the "correct" sample name
        - Renames the demux artifact to the "correct" sample name

The demux artifacts will still have multiple samples associated with them, but at least
the name will be corrected to facilitate downstream processing.

"""
import logging
from argparse import ArgumentParser, Namespace
from datetime import datetime as dt

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from scilifelab_epps.epp import get_pool_sample_label_mapping
from scilifelab_epps.wrapper import epp_decorator

TIMESTAMP = dt.now().strftime("%y%m%d_%H%M%S")


@epp_decorator(__file__, TIMESTAMP)
def main(args: Namespace):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    process = Process(lims, id=args.pid)

    pools = process.all_inputs()
    for pool in pools:
        # Collect demux artifacts associated with pool
        demux_arts = []
        for art_tuple in process.input_output_maps:
            if (
                art_tuple[0]["uri"].id == pool.id
                and art_tuple[1]["output-generation-type"] == "PerReagentLabel"
            ):
                demux_arts.append(art_tuple[1]["uri"])

        try:
            logging.info(f"Checking pool '{pool.name}' ({pool.id}).")
            for demux_art in demux_arts:
                assert len(demux_art.samples) == 1
            logging.info(f"Pool '{pool.name}' ({pool.id}) looks OK.")
        except AssertionError:
            logging.info(f"Pool '{pool.name}' ({pool.id}) needs fixing.")
            logging.info(
                f"Getting sample-label mapping for pool '{pool.name}' ({pool.id})"
            )
            sample2label = get_pool_sample_label_mapping(pool)

            # Iterate across demux arts
            for demux_art in demux_arts:
                # Get reagent label
                assert len(demux_art.reagent_labels) == 1, (
                    f"Demux artifact '{demux_art.name}' ({demux_art.id}) "
                    + " has multiple labels. That should not happen."
                )
                label = demux_art.reagent_labels[0]

                # Use sample-label mapping to find the "correct" sample name
                correct_sample_name = None
                for sample in demux_art.samples:
                    if (
                        sample.name in sample2label
                        and sample2label[sample.name] == label
                    ):
                        correct_sample_name = sample.name
                        break
                if correct_sample_name is None:
                    msg = (
                        "Could not find correct sample name"
                        + f" for demux artifact '{demux_art.name}' ({demux_art.id}) "
                        + f" in pool {pool.id}. Skipping."
                    )
                    logging.error(msg)
                    continue

                correct_demux_art_name = f"{correct_sample_name} (FASTQ reads)"
                if correct_demux_art_name == demux_art.name:
                    logging.info(
                        f"Demux artifact '{demux_art.name}' ({demux_art.id}) "
                        + "already has the correct name. Skipping."
                    )
                    continue
                else:
                    logging.info(
                        f"Renaming '{demux_art.name}' ({demux_art.id}) -> '{correct_demux_art_name}'"
                    )
                    demux_art.name = correct_demux_art_name
                    demux_art.put()


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
    args = parser.parse_args()

    main(args)
