"""Contains useful and reusable code for EPP scripts.

Classes, methods and exceptions.

Johannes Alneberg, Science for Life Laboratory, Stockholm, Sweden.
Copyright (C) 2013 Johannes Alneberg
"""

import csv
import logging
import os
import re
import sys
from importlib.metadata import PackageNotFoundError, version
from logging.handlers import RotatingFileHandler
from shutil import copy
from time import localtime, strftime

import psycopg2
import yaml
from genologics.config import MAIN_LOG
from genologics.entities import Artifact, Process
from genologics.lims import Lims
from requests import HTTPError


def attach_file(src, resource):
    """Attach file at src to given resource

    Copies the file to the current directory, EPP node will upload this file
    automatically if the process output is properly set up"""
    original_name = os.path.basename(src)
    new_name = resource.id + "_" + original_name
    dir = os.getcwd()
    location = os.path.join(dir, new_name)
    copy(src, location)
    return location


class EmptyError(ValueError):
    "Raised if an iterator is unexpectedly empty."

    pass


class NotUniqueError(ValueError):
    "Raised if there are unexpectedly more than 1 item in an iterator"

    pass


def unique_check(to_check, msg):
    "Check that l is of length 1, otherwise raise error, with msg appended"
    if len(to_check) == 0:
        raise EmptyError(f"No item found for {msg}")
    elif len(to_check) != 1:
        raise NotUniqueError(f"Multiple items found for {msg}")


def set_field(element):
    try:
        element.put()
    except (TypeError, HTTPError) as e:
        logging.warning(f"Error while updating element: {e}")


class EppLogger:
    """Context manager for logging module useful for EPP script execution.

    This context manager (CM) automatically logs what script that is executed,
    with what parameters it was executed and what version (including) commit
    hash of the genologics package used. Since EPP scripts are often ran
    automatically by the genologics LIMS client, the stdout and stderr is
    captured and logged within this CM. Stderr is duplicated so that the
    last line can be shown in the GUI. In order to track multiple runs
    of the same process from the genologics LIMS GUI, the previous log
    files can be prepended. Also a main log file can be used that is
    supposed to be common for all scripts executed on the server.

    """

    PACKAGE = "genologics"

    def __enter__(self):
        logging.info(f"Executing file: {sys.argv[0]}")
        logging.info(f"with parameters: {sys.argv[1:]}")
        try:
            logging.info(f"Version of {self.PACKAGE}: " + version(self.PACKAGE))
        except PackageNotFoundError as e:
            logging.error(e)
            logging.error(f"Make sure you have the {self.PACKAGE} package installed")
            sys.exit(-1)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # If no exception has occurred in block, turn off logging.
        if not exc_type:
            logging.shutdown()
            sys.stderr = self.saved_stderr
            sys.stdout = self.saved_stdout
        # Do not repress possible exception
        return False

    def __init__(self, log_file=None, level=logging.INFO, lims=None, prepend=False):
        """Initialize the logger with custom settings.

        Arguments:
        log_file  -- file to write individual log to

        Keyword Arguments:
        level   -- Logging level, default logging.INFO
        lims    -- Lims instance, needed for prepend to work
        prepend -- If True, prepend old log file to new, requires lims
        """
        self.lims = lims
        self.log_file = log_file
        self.level = level
        self.prepend = prepend

        if prepend and self.log_file:
            self.prepend_old_log()

        # Loggers that will capture stdout and stderr respectively
        stdout_logger = logging.getLogger("STDOUT")
        self.slo = self.StreamToLogger(stdout_logger, logging.INFO)
        self.saved_stdout = sys.stdout
        sys.stdout = self.slo

        stderr_logger = logging.getLogger("STDERR")
        self.saved_stderr = sys.stderr
        # Duplicate stderr stream to log
        self.sle = self.StreamToLogger(stderr_logger, logging.INFO, self.saved_stderr)
        sys.stderr = self.sle

        # Root logger with filehandler(s)
        self.logger = logging.getLogger()
        self.logger.setLevel(self.level)
        formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(name)s:%(message)s")
        if self.log_file:
            individual_fh = logging.FileHandler(self.log_file, mode="a")
            individual_fh.setFormatter(formatter)
            self.logger.addHandler(individual_fh)

        if MAIN_LOG:
            # Rotating file handler, that will create up to 10 backup logs,
            # each no bigger than 100MB.
            main_fh = RotatingFileHandler(
                MAIN_LOG, mode="a", maxBytes=1e8, backupCount=10
            )
            main_fh.setFormatter(formatter)
            self.logger.addHandler(main_fh)
        else:
            self.logger.warning("No main log file found.")

    def prepend_old_log(self, external_log_file=None):
        """Prepend the old log to the new log.

        The location of the old log file is retrieved through the REST api.
        In order to work, the script should be executed on the LIMS server
        since the location on the disk is parsed out from the sftp string
        and then used for local copy of file.

        This method does not use logging since that could mess up the
        logging settings, instead warnings are printed to stderr."""
        if external_log_file:
            log_file_name = external_log_file
        else:
            log_file_name = self.log_file

        local_log_path = os.path.join(os.getcwd(), log_file_name)
        if not os.path.isfile(local_log_path):
            try:
                log_artifact = Artifact(self.lims, id=log_file_name)
                log_artifact.get()
                if log_artifact.files:
                    log_path = log_artifact.files[0].content_location.split(
                        self.lims.baseuri.split(":")[1]
                    )[1]
                    if not log_path.startswith("/"):
                        log_path = f"/{log_path}"
                    copy(log_path, local_log_path)
                    with open(local_log_path, "a") as f:
                        f.write("=" * 80 + "\n")
            except HTTPError:  # Probably no artifact found, skip prepending
                print(
                    (f"No log file artifact found for id: {log_file_name}"),
                    file=sys.stderr,
                )
            except OSError as e:  # Probably some path was wrong in copy
                print(
                    (
                        "Log could not be prepended, "
                        f"make sure {log_path} and {log_file_name} are "
                        "proper paths."
                    ),
                    file=sys.stderr,
                )
                raise e

    class StreamToLogger:
        """Fake file-like stream object that redirects writes to a logger
        instance.

        source:
        http://www.electricmonk.nl/log/2011/08/14/
        redirect-stdout-and-stderr-to-a-logger-in-python/
        """

        def __init__(self, logger, log_level=logging.INFO, stream=None):
            self.logger = logger
            self.log_level = log_level
            self.linebuf = ""
            self.stream = stream

        def write(self, buf):
            if self.stream:
                self.stream.write(buf)
            for line in buf.rstrip().splitlines():
                self.logger.log(self.log_level, line.rstrip())


class ReadResultFiles:
    """Class to read pars different kinds of result files from a process.
    The class stores the parsed content of all shared result files in a
    dictionary 'shared_files'. The data is parsed as lists of lists."""

    def __init__(self, process):
        self.process = process
        self.shared_files = self._pars_file("SharedResultFile")
        self.perinput_files = self._pars_file("ResultFile")

    def get_file_path(self, artifact):
        if len(artifact.files) > 0:
            file = artifact.files[0]
            file_path = file.content_location.split("scilifelab.se")[1]
            if len(file_path.split(".")) > 1:
                return file_path
        return None

    def _pars_file(self, output_type):
        """Reads a csv or txt into a list of lists, where sub lists are lines
        of the csv."""
        outs = self.process.all_outputs()
        outarts = [a for a in outs if a.output_type == output_type]
        parsed_files = {}
        for outart in outarts:
            file_path = self.get_file_path(outart)
            if file_path:
                of = open(file_path)
                file_ext = file_path.split(".")[-1]
                if file_ext == "csv":
                    pf = [row for row in csv.reader(of.read().splitlines())]
                    parsed_files[outart.name] = pf
                elif file_ext == "txt":
                    pf = [row.strip().strip("\\").split("\t") for row in of.readlines()]
                    parsed_files[outart.name] = pf
                of.close()
        return parsed_files

    def format_file(
        self,
        parsed_file,
        name="",
        first_header=None,
        header_row=None,
        root_key_col=0,
        find_keys=[],
    ):
        """Function to format a parsed csv or txt file.

        Arguments and Output:
            parsed_file     A list of lists where sublists are rows of the csv.
            name            Name of parsed file.
            first_header    First column of the heather section in the file.
                            default value is 'None'
            root_key_col    If you want the root keys to be given by some other
                            column than the first one, set root_key_col to the
                            column number.
            header_row      Instead of specifying first_header you can choose
                            from what line to reed by setting header_row to the
                            row number where you want to start reading.
            find_keys       List of row names to look for. Will exclude all
                            others.
            file_info       Dict of dicts. Keys of root dict are the first
                            column in the csv starting from the line after the
                            heather line. Keys of sub dicts are the columns of
                            the heather line."""
        file_info = {}
        keys = []
        error_message = ""
        duplicated_lines = []
        exeptions = ["Sample", "Fail", ""]
        if not isinstance(first_header, list):
            if first_header:
                first_header = [first_header]
            else:
                first_header = []
        for row, line in enumerate(parsed_file):
            if keys and len(line) == len(keys):
                root_key = line[root_key_col]
                cond1 = find_keys == [] and root_key not in exeptions
                cond2 = root_key in find_keys
                if root_key in file_info:
                    duplicated_lines.append(root_key)
                elif cond1 or cond2:
                    file_info[root_key] = {}
                    if not duplicated_lines:
                        for col in range(len(keys)):
                            if keys[col] != "":
                                file_info[root_key][keys[col]] = line[col]
                            elif keys[col - 1] != "":
                                tupl = (file_info[root_key][keys[col - 1]], line[col])
                                file_info[root_key][keys[col - 1]] = tupl

            head = line[root_key_col] if len(line) > root_key_col else None
            if first_header and head in first_header:
                keys = line
            elif header_row and row == header_row:
                keys = line
        if duplicated_lines:
            error_message = (
                "Row names {} occurs more than once in file {}. "
                "Fix the file to continue. "
            ).format(",".join(duplicated_lines), name)
        if not file_info:
            error_message = error_message + f"Could not format parsed file {name}."
        if error_message:
            print(error_message, file=sys.stderr)
            sys.exit(-1)
        return file_info


class CopyField:
    """Class to copy any filed (or udf) from any lims element to any
    udf on any other lims element

    arguments:

    s_elt           source element - instance of a type
    d_elt           destination element - instance of a type
    s_field_name    name of source field (or udf) to be copied
    d_udf_name      name of destination udf name. If not specified
                    s_field_name will be used.

    The copy_udf() function takes a log file as optional argument.
    If this is given the changes will be logged there.

    Written by Maya Brandi and Johannes Alnberg
    """

    def __init__(self, s_elt, d_elt, s_field_name, d_udf_name=None):
        if not d_udf_name:
            d_udf_name = s_field_name
        self.s_elt = s_elt
        self.s_field_name = s_field_name
        self.s_field = self._get_field(s_elt, s_field_name)
        self.d_elt = d_elt
        self.d_type = d_elt._URI
        self.d_udf_name = d_udf_name
        self.old_dest_udf = self._get_field(d_elt, d_udf_name)

    def _current_time(self):
        return strftime("%Y-%m-%d %H:%M:%S", localtime())

    def _get_field(self, elt, field):
        if field in elt.udf:
            return elt.udf[field]
        else:
            try:
                return elt.field
            except:
                return None

    def _set_udf(self, elt, udf_name, val):
        try:
            elt.udf[udf_name] = val
            elt.put()
            return True
        except (TypeError, HTTPError) as e:
            print(f"Error while updating element: {e}", file=sys.stderr)
            sys.exit(-1)
            return False

    def _log_before_change(self, changelog_f=None):
        if changelog_f:
            d = {
                "ct": self._current_time(),
                "s_udf": self.s_field_name,
                "sn": self.d_elt.name,
                "si": self.d_elt.id,
                "su": self.old_dest_udf,
                "nv": self.s_field,
                "d_elt_type": self.d_type,
            }

            changelog_f.write(
                (
                    "{ct}: udf: '{s_udf}' on {d_elt_type}: '{sn}' ("
                    "id: {si}) is changed from '{su}' to '{nv}'.\n"
                ).format(**d)
            )

        logging.info(
            f"Copying from element with id: {self.s_elt.id} to element with "
            f" id: {self.d_elt.id}"
        )

    def _log_after_change(self):
        d = {
            "s_udf": self.s_field_name,
            "d_udf": self.d_udf_name,
            "su": self.old_dest_udf,
            "nv": self.s_field,
            "d_elt_type": self.d_type,
        }

        logging.info(
            "Updated {d_elt_type} udf: {d_udf}, from {su} to {nv}.".format(**d)
        )

    def copy_udf(self, changelog_f=None):
        if self.s_field != self.old_dest_udf:
            self._log_before_change(changelog_f)
            log = self._set_udf(self.d_elt, self.d_udf_name, self.s_field)
            self._log_after_change()
            return log
        else:
            return False


def get_well_number(art: Artifact, count_per: str) -> int:
    """Convert well names (e.g. 'A:1', 'H:12') to well numbers.

    Choose between column-wise or row-wise counting.
    """

    assert count_per in ["row", "col"], "Invalid function argument"

    # Ensure container well names match classical convention
    assert (
        art.container.type.y_dimension["is_alpha"]
        and art.container.type.y_dimension["offset"] == 0
        and not art.container.type.x_dimension["is_alpha"]
        and art.container.type.x_dimension["offset"] == 1
    ), "Can't convert well name --> well number for invalid container"

    # Get dimensions of artifact container
    n_cols = art.container.type.x_dimension["size"]
    n_rows = art.container.type.y_dimension["size"]

    # Get simple dict translating letters to numbers
    letter2num = {}
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i, letter in zip(range(1, len(alphabet) + 1), alphabet):
        letter2num[letter] = i

    # Collect well data from LIMS
    well_name = art.location[1]
    row_letter, col_num = well_name.split(":")
    row_num = letter2num[row_letter]

    if count_per == "row":
        well_num = (int(row_num) - 1) * n_cols + int(col_num)
    elif count_per == "col":
        well_num = (int(col_num) - 1) * n_rows + int(row_num)
    else:
        raise AssertionError

    return well_num


def get_matching_inputs(
    process: Process, output_artifact: Artifact
) -> list[Artifact] | None:
    """Get the input artifacts coupled to an output artifact."""
    input_arts = []
    for io_tuple in process.input_output_maps:
        input_art = io_tuple[0]["uri"]
        output_art = io_tuple[1]["uri"]
        if input_art.type == "Analyte" and output_art.id == output_artifact.id:
            input_arts.append(io_tuple[0]["uri"])

    if input_arts:
        return input_arts
    else:
        return None


def traceback_to_step(
    art: Artifact, step_name_pattern: re.Pattern, allow_multiple_inputs: bool = False
) -> tuple[Process, list[Artifact], Artifact] | None:
    """Try to backtrack an artifact to a target step based on a supplied step name pattern.

    Returns:
    - The target step
    - A list of it's matching input Artifacts
    - It's matching output Artifact

    If a step is reached where the artifact has multiple linked inputs, linear traceback is not possible.
    This will either return None or raise an error, depending on the value of allow_multiple_inputs.

    Example:

        To backtrack an ONT sequencing library to the step in which it was pooled,
        either "ONT Pooling" or "ONT QC Pooling".

            traceback_to_step(
                art=ont_library_artifact,
                step_name_pattern=re.compile(r"ONT.*Pooling"),
                allow_multiple_inputs=True,
            )

        - If the ONT sequencing library did not pass through a matching pooling step,
            the function will return None.
        - If we set allow_multiple_inputs=False and the library backtracks to a non-matching pooling step,
            it will throw an error instead.

    """

    current_art = art
    logging.info(
        f"Attempting to backtrack artifact '{current_art.name}' to step matching pattern {step_name_pattern.pattern}, from step '{current_art.parent_process.type.name}'"
    )

    # Loop until return, or as long as there is a parent process
    while current_art.parent_process is not None:
        current_pp = current_art.parent_process
        logging.info(f"Backtracking to parent process '{current_pp.type.name}'")

        input_arts = get_matching_inputs(current_pp, current_art)
        assert input_arts is not None, "No matching input artifacts found."

        match = re.match(step_name_pattern, current_pp.type.name)

        if match:
            logging.info(f"Found matching step '{current_pp.type.name}'. Returning.")
            return (current_pp, input_arts, current_art)
        elif len(input_arts) > 1:
            msg = f"Output artifact {current_art.name} in step '{current_pp.type.name}' has multiple inputs. Can't traceback further."
            logging.info(msg)
            if allow_multiple_inputs:
                logging.info("Target step not found, returning None.")
                return None
            else:
                raise AssertionError(msg)
        else:
            # Continue backtracking
            current_art = input_arts[0]

    logging.info(
        f"Traceback reached the beginning of the process tree ('{current_pp.type.name}'), returning None."
    )
    return None


def upload_file(
    file_path: str,
    file_slot: str,
    process: Process,
    lims: Lims,
    remove: bool = True,
    fail_on_missing_file_slot: bool = True,
):
    matching_file_slots = [
        output_artifact
        for output_artifact in process.all_outputs()
        if output_artifact.name == file_slot
    ]

    if not matching_file_slots:
        msg = f"Found no matching file slots '{file_slot}' in process '{process.type.name}' ({process.id})."
        logging.warning(msg)

        if remove:
            os.remove(file_path)
            logging.info(f"Removed local instance of '{file_path}'")

        if fail_on_missing_file_slot:
            raise AssertionError(msg)

    else:
        correct_file_slot: Artifact = matching_file_slots[0]
        for f in correct_file_slot.files:
            lims.request_session.delete(f.uri)
        lims.upload_new_file(correct_file_slot, file_path)

        logging.info(
            f"'{file_path}' uploaded to LIMS file slot '{correct_file_slot.name}' ({correct_file_slot.id})."
        )
        if remove:
            os.remove(file_path)
            logging.info(f"Removed local instance of '{file_path}'")


def get_pool_sample_label_mapping(pool: Artifact) -> dict[str, str]:
    """Given a pool artifact containing labeled samples, use database queries to
    build a dictionary mapping each sample name to its reagent label.
    """
    with open("/opt/gls/clarity/users/glsai/config/genosqlrc.yaml") as f:
        config = yaml.safe_load(f)

    # Setup DB connection
    connection = psycopg2.connect(
        user=config["username"],
        host=config["url"],
        database=config["db"],
        password=config["password"],
    )
    cursor = connection.cursor()

    """Supply a pool artifact ID and a sample name:
    1. Find the ancestor artifacts of the pool artifact.
    2. Filter for derived sample artifacts
    3. Filter for artifacts sharing a name with the target sample
    4. Filter for artifacts with reagent labels
    """
    query = """--sql
        SELECT
            DISTINCT rl.name
        FROM
            -- Table mapping artifact IDs to ancestor artifact IDs
            artifact_ancestor_map aam
            -- Join artifact information on ancestor artifact IDs
            JOIN artifact parent ON aam.ancestorartifactid = parent.artifactid
            -- Join reagent label information on ancestor artifact IDs
            LEFT JOIN artifact_label_map alm ON parent.artifactid = alm.artifactid
            LEFT JOIN reagentlabel rl ON rl.labelid = alm.labelid
        WHERE
            -- The pool artifact ID is used to find its ancestors
            aam.artifactid = {}
            -- Filter for derived sample artifacts
            AND parent.artifacttypeid = 2
            -- Filter for artifacts sharing a name with the target sample
            AND parent.name = '{}'
            -- Filter for artifacts with reagent labels
            AND rl.name IS NOT NULL;
    """

    errors = False
    sample2label = {}
    pool_db_id = int(pool.id.split("-")[1])
    for sample in pool.samples:
        try:
            cursor.execute(query.format(pool_db_id, sample.name))
            query_results = cursor.fetchall()

            assert len(query_results) != 0, (
                f"No reagent labels found for sample '{sample.name}'."
            )
            assert len(query_results) == 1, (
                f"Multiple reagent labels found for sample '{sample.name}'."
            )

            label = query_results[0][0]
            sample2label[sample.name] = label
        except AssertionError as e:
            logging.error(str(e), exc_info=True)
            logging.warning(f"Skipping sample '{sample.name}' due to error.")
            errors = True
            continue

    if errors:
        raise AssertionError(
            "Errors occurred when linking samples and indices. Please report this error."
        )
    else:
        return sample2label
