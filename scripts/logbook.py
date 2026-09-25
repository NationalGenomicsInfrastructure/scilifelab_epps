#!/usr/bin/env python

DESC = """EPP script for automatically logging the record of instrument use
into the electronic logbooks in Google Doc.
Author: Chuan Wang, Science for Life Laboratory, Stockholm, Sweden
"""

import argparse
import datetime
import os

import httplib2
from apiclient import discovery
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
from oauth2client import client, tools
from oauth2client.file import Storage

from data.logbook_data import GDoc_logbook, lims_process_record
from scilifelab_epps.epp import EppLogger

# If modifying these scopes, delete your previously saved credentials
# at ~/.credentials/*.json
SCOPES = ["https://www.googleapis.com/auth/spreadsheets"]
CLIENT_SECRET_FILE = "~/.credentials/client_secret.json"
APPLICATION_NAME = "GDOC API FOR LIMS EPP"
flags = argparse.ArgumentParser(parents=[tools.argparser]).parse_args(
    "--auth_host_name localhost --logging_level INFO".split()
)


def get_credentials():
    """Gets valid user credentials from storage.

    If nothing has been stored, or if the stored credentials are invalid,
    the OAuth2 flow is completed to obtain the new credentials.

    Returns:
        Credentials, the obtained credential.
    """
    home_dir = os.path.expanduser("~")
    credential_dir = os.path.join(home_dir, ".credentials")
    if not os.path.exists(credential_dir):
        os.makedirs(credential_dir)
    credential_path = os.path.join(
        credential_dir, "sheets.googleapis.com-python-limsepp.json"
    )
    store = Storage(credential_path)
    credentials = store.get()

    if not credentials or credentials.invalid:
        flow = client.flow_from_clientsecrets(CLIENT_SECRET_FILE, SCOPES)
        flow.user_agent = APPLICATION_NAME
        if flags:
            credentials = tools.run_flow(flow, store, flags)
        else:  # Needed only for compatibility with Python 2.6
            credentials = tools.run(flow, store)
        print("Storing credentials to " + credential_path)
    return credentials


# In coordinates, the values correspond to LineID of the header, columns for date, operator, instrument name and details. The sheet name shoud be "Logbook"
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


def write_record(content, dest_file):
    credentials = get_credentials()
    http = credentials.authorize(httplib2.Http())
    discoveryUrl = "https://sheets.googleapis.com/$discovery/rest?version=v4"
    service = discovery.build(
        "sheets", "v4", http=http, discoveryServiceUrl=discoveryUrl
    )
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
    request.execute()

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

    (
        service.spreadsheets()
        .values()
        .batchUpdate(spreadsheetId=spreadsheetId, body=body)
        .execute()
    )


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
    parser = argparse.ArgumentParser(description=DESC)
    parser.add_argument(
        "--pid", default="24-38458", dest="pid", help="Lims id for current Process"
    )
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
