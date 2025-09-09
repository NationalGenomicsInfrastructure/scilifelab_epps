#!/usr/bin/env python
import datetime
import logging
import subprocess
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from scilifelab_epps.epp import upload_file
from scilifelab_epps.wrapper import epp_decorator

TIMESTAMP = datetime.datetime.now().strftime("%y%m%d_%H%M%S")


def build_zpl_format(format_lines, data_lines, copies=1):
    """Helper to build ZPL format and data for a label."""
    lines = []
    lines.append("^XA")  # Start format definition
    lines.append("^DFFORMAT^FS")  # Delete previous format named FORMAT
    lines.append("^LH0,0")  # Set label home position
    lines.extend(format_lines)  # Add format (layout) commands
    lines.append("^XZ")  # End format definition
    for _ in range(copies):
        lines.append("^XA")  # Start label instance
        lines.append("^XFFORMAT^FS")  # Recall the format defined above
        lines.extend(data_lines)  # Add data fields for this label
        lines.append("^XZ")  # End label instance
    return lines


def make_container_label(plateid, copies=1):
    """Construct label with container id as human readable and barcode"""
    format_lines = [
        "^FO360,30^AFN 78,39^FN1^FS",  # Field origin, font, field number 1 (human readable)
        "^FO70,10^BCN,70,N,N^FN2^FS",  # Field origin, barcode, field number 2 (barcode)
    ]
    data_lines = [
        f"^FN1^FD{plateid}^FS",  # Assign plateid to field 1 (human readable)
        f"^FN2^FD{plateid}^FS",  # Assign plateid to field 2 (barcode)
    ]
    return build_zpl_format(format_lines, data_lines, copies)


def makeContainerNameBarcode(plate_name, copies=1):
    """Construct label with container name as human readable"""
    format_lines = []
    # Adjust font size and position based on name length
    if len(plate_name) > 21:
        format_lines.append("^FO20,40^AFN 54,30^FN1^FS")  # Smaller font for long names
    else:
        format_lines.append("^FO20,30^AFN 78,39^FN1^FS")  # Larger font for short names
    data_lines = [
        f"^FN1^FD{plate_name}^FS"  # Assign plate_name to field 1 (human readable)
    ]
    return build_zpl_format(format_lines, data_lines, copies)


def makeOperatorAndDateBarcode(operator, date, copies=1):
    """Construct label with operator name and date in human readable format"""
    format_lines = [
        "^FO420,35^ADN,36,20^FN1^FS",  # Field for date (right side)
        "^FO20,35^ADN,36,20^FN2^FS",  # Field for operator (left side)
    ]
    if len(operator) > 19:
        operator = operator[:19]  # Truncate operator name if too long
    data_lines = [
        f"^FN1^FD{date}^FS",  # Assign date to field 1
        f"^FN2^FD{operator}^FS",  # Assign operator to field 2
    ]
    return build_zpl_format(format_lines, data_lines, copies)


def makeProcessNameBarcode(process_name, copies=1):
    """Construct label with process name as human readable"""
    format_lines = []
    # Adjust font size and position based on process name length
    if len(process_name) > 21:
        format_lines.append("^FO20,40^ADN 54,30^FN1^FS")  # Smaller font for long names
    else:
        format_lines.append("^FO20,30^AFN 78,39^FN1^FS")  # Larger font for short names
    data_lines = [
        f"^FN1^FD{process_name}^FS"  # Assign process_name to field 1 (human readable)
    ]
    return build_zpl_format(format_lines, data_lines, copies)


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    process = Process(lims, id=args.pid)

    # Build a list of ZPL (=Zebra Programming Language) lines, corresponding to 4 labels per output container
    zpl_code = []
    for container in process.output_containers():
        logging.info(
            f"Making label for container ID with barcode: <barcode> {container.id}"
        )
        zpl_code += make_container_label(container.id)

        logging.info(f"Making label for container ID: {container.id}")
        zpl_code += makeContainerNameBarcode(container.name)

        logging.info(
            f"Making label for operator and date: {process.technician.name} {str(datetime.date.today())}"
        )
        zpl_code += makeOperatorAndDateBarcode(
            process.technician.name, str(datetime.date.today())
        )

        logging.info(f"Making label for step name: {process.type.name}")
        zpl_code += makeProcessNameBarcode(process.type.name)

    # Build args list to label printer command
    lp_args = ["lp"]
    lp_args += ["-h", "homer2.scilifelab.se:631"]
    lp_args += ["-d", "zebrabarcode"]
    lp_args.append("-")  # make lp command read from stdin
    logging.info(f"Using command: '{' '.join(lp_args)}'")

    # Call label printer command
    logging.info("Calling command...")
    if not args.test:
        lp_process = subprocess.Popen(
            lp_args,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="utf8",
        )
        logging.info("Piping ZPL contents...")
        lp_process.stdin.write(str("\n".join(zpl_code)))
        stdout, stderr = lp_process.communicate()  # Will wait for subprocess to finish
        logging.info(f"lp stdout: {stdout}")
        logging.info(f"lp stderr: {stderr}")
        logging.info("Command finished, closing subprocess.")
        lp_process.stdin.close()
    else:
        logging.info("Just kidding. This is a test run.")

    # Upload file with ZPL contents, will persist after finishing step, useful for re-prints and doing LIMS from home
    filename = f"barcodes_{process.id}_{TIMESTAMP}_{process.technician.name.replace(' ', '')}.txt"
    logging.info(f"Uploading ZPL contents as {filename}")

    with open(filename, "w") as f:
        f.write(str("\n".join(zpl_code)))

    upload_file(
        filename, args.file, process, lims, remove=True, fail_on_missing_file_slot=False
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--pid", help="The process LIMS id.")
    parser.add_argument("--file", help="LIMS file slot name to use for barcode file.")
    parser.add_argument("--log", help="LIMS file slot name to use for log file.")
    parser.add_argument(
        "--test",
        action="store_true",
        default=False,
        help="Test run, suppress actual label printing.",
    )
    args = parser.parse_args()

    main(args)
