#!/usr/bin/env python

import logging
import os
import sys

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process, Researcher
from genologics.lims import Lims

from scilifelab_epps.epp import upload_file
from scilifelab_epps.utils.get_epp_user import get_epp_user


class TrackingRootLogger(logging.RootLogger):
    """A root logger that tracks whether any errors or warnings have been emitted."""

    def __init__(self, level):
        """Initialize with the given level and errors_or_warnings set to False."""
        super().__init__(level)
        self.errors_or_warnings = False

    def handle(self, record):
        """Handle a log record and track if it's an error or warning."""
        # Check if the record is a warning or error before handling it
        if record.levelno >= logging.WARNING:
            self.errors_or_warnings = True

        # Let the parent class handle the record normally
        return super().handle(record)


def epp_decorator(script_path: str, timestamp: str):
    """This top-level decorator is meant to be used on EPP scripts' main functions.

    It receives the script path (__file__) and timestamp (yymmdd_hhmmss) as arguments to
    pass on to it's children which wrap the main function to handle logging and graceful failure.
    """
    script_name: str = os.path.basename(script_path).split(".")[0]

    def _epp_decorator(script_main):
        def epp_wrapper(args):
            """General wrapper for EPP scripts."""

            # Set up LIMS
            lims = Lims(BASEURI, USERNAME, PASSWORD)
            lims.check_version()
            process = Process(lims, id=args.pid)

            # Get EPP user
            try:
                epp_user: Researcher = get_epp_user(lims, args.pid)
            except ValueError:
                epp_user = None

            # Name log file
            log_filename: str = (
                "_".join(
                    [
                        script_name,
                        process.id,
                        timestamp,
                        (epp_user or process.technician).name.replace(" ", ""),
                    ]
                )
                + ".log"
            )

            # Set up logging

            # Set custom subclass as root logger
            logger = TrackingRootLogger(level=logging.INFO)
            logging.root = logger

            # Clear any existing handlers (to avoid duplicates)
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)

            # Create file handler
            file_handler = logging.FileHandler(log_filename, mode="w")
            file_handler.setLevel(logging.INFO)
            file_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))

            # Create stdout handler
            stdout_handler = logging.StreamHandler(sys.stdout)
            stdout_handler.setLevel(logging.INFO)
            stdout_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))

            # Add both handlers to logger
            logger.addHandler(file_handler)
            logger.addHandler(stdout_handler)

            # Start logging
            if not epp_user:
                logging.warning(f"No EPP user found for process ID {args.pid}.")
            logging.info(
                f"Script '{script_name}' started at {timestamp} by {(epp_user.name if epp_user else 'unknown')}."
            )
            logging.info(
                f"Launched in step '{process.type.name}' ({process.id}) opened by {process.technician.name}."
            )
            args_str = "\n\t".join(
                [f"'{arg}': {getattr(args, arg)}" for arg in vars(args)]
            )
            logging.info(f"Script called with arguments: \n\t{args_str}")

            # Run EPP, check for errors
            script_exception = None
            try:
                script_main(args)
            except Exception as e:
                # Post error to LIMS GUI
                script_exception = e
                logger.error(e)

            # Check for error / warning messages in log
            errors_or_warnings_in_log = logger.errors_or_warnings
            logging.shutdown()

            # Check for log upload success
            try:
                upload_file(
                    file_path=log_filename,
                    file_slot=args.log,
                    process=process,
                    lims=lims,
                    remove=True,
                )
                log_upload_failed = False
            except:
                log_upload_failed = True

            # Build final message
            final_message_lines = []
            if script_exception:
                final_message_lines.append("Script failed.")
            else:
                final_message_lines.append("Script finished successfully.")
            if errors_or_warnings_in_log:
                final_message_lines.append(
                    "Errors or warnings occurred, please check log."
                )
            if log_upload_failed:
                final_message_lines.append("Log upload failed.")

            # Post final message and exit
            sys.stderr.write(" ".join(final_message_lines))
            if errors_or_warnings_in_log:
                sys.exit(1)
            else:
                sys.exit(0)

        return epp_wrapper

    return _epp_decorator
