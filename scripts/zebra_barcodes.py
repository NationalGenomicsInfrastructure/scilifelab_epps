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
DEFAULT_CONFIG = {
    "printer": {
        "command": "lp",
        "host": "ipp.sys.kth.se:631",
        "destination": "zebrabarcode",
        "stdin_arg": "-",
    },
    "labels": {
        "container_id": {
            "copies": 1,
            "format_lines": [
                "^FO360,30^AFN 78,39^FN1^FS",
                "^FO70,10^BCN,70,N,N^FN2^FS",
            ],
        },
        "container_name": {
            "copies": 1,
            "max_length": 21,
            "short_format_line": "^FO20,30^AFN 78,39^FN1^FS",
            "long_format_line": "^FO20,40^AFN 54,30^FN1^FS",
        },
        "operator_date": {
            "copies": 1,
            "operator_max_length": 19,
            "format_lines": [
                "^FO420,35^ADN,36,20^FN1^FS",
                "^FO20,35^ADN,36,20^FN2^FS",
            ],
        },
        "process_name": {
            "copies": 1,
            "max_length": 21,
            "short_format_line": "^FO20,30^AFN 78,39^FN1^FS",
            "long_format_line": "^FO20,40^ADN 54,30^FN1^FS",
        },
    },
}


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


def make_container_label(plateid):
    """Construct label with container id as human readable and barcode"""
    label_config = DEFAULT_CONFIG["labels"]["container_id"]
    format_lines = label_config["format_lines"]
    data_lines = [
        f"^FN1^FD{plateid}^FS",  # Assign plateid to field 1 (human readable)
        f"^FN2^FD{plateid}^FS",  # Assign plateid to field 2 (barcode)
    ]
    return build_zpl_format(format_lines, data_lines, label_config["copies"])


def makeNameBarcode(plate_name, type):
    """Construct label with container name as human readable"""
    label_config = DEFAULT_CONFIG["labels"][type]
    format_lines = []
    # Adjust font size and position based on name length
    if len(plate_name) > label_config["max_length"]:
        format_lines.append(
            label_config["long_format_line"]
        )  # Smaller font for long names
    else:
        format_lines.append(
            label_config["short_format_line"]
        )  # Larger font for short names
    data_lines = [
        f"^FN1^FD{plate_name}^FS"  # Assign plate_name to field 1 (human readable)
    ]
    return build_zpl_format(format_lines, data_lines, label_config["copies"])


def makeOperatorAndDateBarcode(operator, date):
    """Construct label with operator name and date in human readable format"""
    label_config = DEFAULT_CONFIG["labels"]["operator_date"]
    format_lines = label_config["format_lines"]
    if len(operator) > label_config["operator_max_length"]:
        operator = operator[:19]  # Truncate operator name if too long
    data_lines = [
        f"^FN1^FD{date}^FS",  # Assign date to field 1
        f"^FN2^FD{operator}^FS",  # Assign operator to field 2
    ]
    return build_zpl_format(format_lines, data_lines, label_config["copies"])


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
        zpl_code += makeNameBarcode(container.name, "container_name")

        logging.info(
            "Making label for operator and date: "
            f"{process.technician.name} {str(datetime.date.today())}"
        )
        zpl_code += makeOperatorAndDateBarcode(
            process.technician.name, str(datetime.date.today())
        )

        logging.info(f"Making label for step name: {process.type.name}")
        zpl_code += makeNameBarcode(process.type.name, "process_name")

    # Build args list to label printer command
    printer_config = DEFAULT_CONFIG["printer"]
    lp_args = [printer_config["command"]]
    lp_args += ["-h", printer_config["host"]]
    lp_args += ["-d", printer_config["destination"]]
    lp_args.append(printer_config["stdin_arg"])  # make lp command read from stdin
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
