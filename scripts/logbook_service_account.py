#!/usr/bin/env python

DESC = """EPP script for automatically logging the record of instrument use
into the electronic logbooks in Google Doc.
Author: Chuan Wang, Science for Life Laboratory, Stockholm, Sweden
"""

import datetime
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError

from data.logbook_data import GDoc_logbook, lims_process_record
from scilifelab_epps.epp import EppLogger

# If modifying these scopes, delete your previously saved credentials
# at ~/.credentials/*.json
SCOPES = ["https://www.googleapis.com/auth/spreadsheets"]
SERVICE_ACCOUNT_FILE = "/opt/gls/clarity/users/glsai/.credentials/ngi-lims-epp-001.json"

# In coordinates, the values correspond to LineID of the header, columns for date, operator, instrument name and details. The sheet name should be "Logbook"
LOGBOOK_COORDINATES = [22, "B", "C", "D", "E"]
SAMPLE_LEVEL_FALLBACK_UDFS = {"PCR Cycler"}


def get_udf_value(pro, udf_name):
    """Return process UDF value, or sample-level fallback for selected UDFs."""
    value = pro.udf.get(udf_name)
    if value:
        return value

    if udf_name not in SAMPLE_LEVEL_FALLBACK_UDFS:
        return None

    values = set()
    for art in pro.all_outputs():
        if art.type != "Analyte":
            continue
        sample_value = art.udf.get(udf_name)
        if sample_value:
            values.add(sample_value)

    if values:
        return ", ".join(sorted(str(v) for v in values))
    return None


def get_credentials():
    """Gets valid user credentials from storage.

    If nothing has been stored, or if the stored credentials are invalid,
    the OAuth2 flow is completed to obtain the new credentials.

    Returns:
        Credentials, the obtained credential.
    """
    try:
        credentials = service_account.Credentials.from_service_account_file(
            SERVICE_ACCOUNT_FILE, scopes=SCOPES
        )
        return credentials
    except FileNotFoundError:
        sys.exit("Missing credentials for service account")


def write_record(content, dest_file):
    credentials = get_credentials()
    service = build("sheets", "v4", credentials=credentials)

    spreadsheetId = GDoc_logbook[dest_file]["File"]

    # Insert empty line
    batch_update_spreadsheet_request_body = {
        "requests": [
            {
                "insertDimension": {
                    "range": {
                        "sheetId": 1,
                        "dimension": "ROWS",
                        "startIndex": LOGBOOK_COORDINATES[0],
                        "endIndex": LOGBOOK_COORDINATES[0] + 1,
                    },
                    "inheritFromBefore": False,
                }
            },
        ]
    }
    request = service.spreadsheets().batchUpdate(
        spreadsheetId=spreadsheetId, body=batch_update_spreadsheet_request_body
    )
    try:
        request.execute()
    except HttpError:
        sys.exit("Service account has no editing access to the logbook")

    # Fill in values
    rangeName = "".join(
        [
            "Logbook",
            "!",
            str(LOGBOOK_COORDINATES[1]),
            str(LOGBOOK_COORDINATES[0] + 1),
            ":",
            str(LOGBOOK_COORDINATES[4]),
            str(LOGBOOK_COORDINATES[0] + 1),
        ]
    )
    values = [content]
    data = [{"range": rangeName, "majorDimension": "ROWS", "values": values}]
    body = {"valueInputOption": "USER_ENTERED", "data": data}
    try:
        (
            service.spreadsheets()
            .values()
            .batchUpdate(spreadsheetId=spreadsheetId, body=body)
            .execute()
        )
    except HttpError:
        sys.exit("Service account has no editing access to the logbook")


# Fetch UDF details
def get_details(record_instrument, pro):
    udf_detail = []
    for item in record_instrument.get("details", []):
        if item == "Processname":
            udf_detail.append(pro.type.name)
        elif pro.udf.get(item):
            udf_detail.append(item + ":" + pro.udf.get(item))
    return ",".join(udf_detail) if udf_detail else "-"


# All logics about logging
def main(lims, pid, epp_logger):
    pro = Process(lims, id=pid)
    log = []
    time = datetime.datetime.now().strftime("%Y-%m-%d")
    log.append(time)
    user = f"{pro.technician.first_name} {pro.technician.last_name}"
    log.append(user)
    log_tmp = log[:]
    record = lims_process_record[pro.type.name]

    for instrument in record:
        log = log_tmp[:]
        if instrument == "default":
            log.append("-")
        elif instrument.startswith("lims_"):
            if pro.instrument.type == "Manual":
                continue
            log.append(pro.instrument.name)
        elif instrument.startswith("udf_"):
            udf_value = get_udf_value(pro, instrument[4:])
            if udf_value:
                log.append(udf_value)
            else:
                continue
        details = get_details(record[instrument], pro)
        log.append(details)
        if isinstance(record[instrument]["dest_file"], list):
            write_record(log, pro.instrument.type)
        else:
            write_record(log, record[instrument]["dest_file"])


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", dest="pid", help="Lims id for current Process")
    parser.add_argument(
        "--log",
        dest="log",
        help=("File name for standard log file, for runtime information and problems."),
    )
    args = parser.parse_args()
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    with EppLogger(log_file=args.log, lims=lims, prepend=True) as epp_logger:
        main(lims, args.pid, epp_logger)
